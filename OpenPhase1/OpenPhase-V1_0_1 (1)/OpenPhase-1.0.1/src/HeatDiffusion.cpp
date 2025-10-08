/*
 *   This file is part of the OpenPhase (R) software library.
 *
 *   Copyright (c) 2009-2022 Ruhr-Universitaet Bochum,
 *                 Universitaetsstrasse 150, D-44801 Bochum, Germany
 *             AND 2018-2022 OpenPhase Solutions GmbH,
 *                 Universitaetsstrasse 136, D-44799 Bochum, Germany.
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   File created :   2014
 *   Main contributors :   Oleg Shchyglo; Philipp Engels; Raphael Schiedung;
 *                         Marvin Tegeler; Helge Schaar
 *
 */

#include "HeatDiffusion.h"
#include "Info.h"
#include "Settings.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"
#include "Composition.h"
#include "FluidDynamics/FlowSolverLBM.h"

namespace openphase
{
using namespace std;

void HeatDiffusion::Initialize(Settings& locSettings)
{
    thisclassname = "HeatDiffusion";

    TotalNx = locSettings.TotalNx;
    OffsetX = locSettings.OffsetX;
    TotalNy = locSettings.TotalNy;
    OffsetY = locSettings.OffsetY;
    TotalNz = locSettings.TotalNz;
    OffsetZ = locSettings.OffsetZ;

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;

    dNx = locSettings.dNx;
    dNy = locSettings.dNy;
    dNz = locSettings.dNz;

    Nphases = locSettings.Nphases;
    dx = locSettings.dx;

    Tolerance = 1.0e-6;
    MaxIterations = 10000;

    PhaseThermalConductivity.Allocate(Nphases);
    PhaseVolumetricHeatCapacity.Allocate(Nphases);
    PhaseSpecificHeatCapacity.Allocate(Nphases);

    LatentHeat.Allocate(Nphases,Nphases);

    size_t Bcells = locSettings.Bcells;
    EffectiveThermalConductivity.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
    EffectiveHeatCapacity.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);

    TxOld.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
    dTx.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
    Qdot.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void HeatDiffusion::ReadInput(const std::string InputFileName)
{
    Info::WriteLineInsert("HeatDiffusion input");
    Info::WriteStandard("Source", InputFileName.c_str());

    fstream inpF(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inpF)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };
    stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();
    ReadInput(inp);
}

void HeatDiffusion::ReadInput(std::stringstream& inp)
{
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    for(size_t n = 0; n < Nphases; n++)
    {
        stringstream converter;
        converter << n;
        string counter = converter.str();

        Info::WriteBlankLine();

        PhaseThermalConductivity[n]    = UserInterface::ReadParameterD(inp, moduleLocation, string("ThermalConductivity_") + counter);
        PhaseVolumetricHeatCapacity[n] = UserInterface::ReadParameterD(inp, moduleLocation, string("VolumetricHeatCapacity_") + counter, false, 0.0);
        PhaseSpecificHeatCapacity[n] = UserInterface::ReadParameterD(inp, moduleLocation, string("SpecificHeatCapacity_") + counter, false, 0.0);
    }

    for (size_t n = 0; n < Nphases; n++)
    for (size_t m = n; m < Nphases; m++)
    {
        if(n != m)
        {
            stringstream idx;

            idx << "_" << n << "_" << m;
            Info::WriteBlankLine();

            LatentHeat(n,m) = UserInterface::ReadParameterD(inp, moduleLocation, "LatentHeat" + idx.str(), false, 0.0);
            LatentHeat(m,n) = -LatentHeat(n,m);
        }
    }

    Tolerance     = UserInterface::ReadParameterD(inp, moduleLocation, string("Tolerance"), false, Tolerance);
    MaxIterations = UserInterface::ReadParameterI(inp, moduleLocation, string("MaxIterations"), false, MaxIterations);

    MaxThermalConductivity = 0.0;
    for(size_t n = 0; n < Nphases; n++)
    {
        MaxThermalConductivity = max(MaxThermalConductivity, PhaseThermalConductivity[n]);
    }

    Info::WriteLine();
    Info::WriteBlankLine();
}

void HeatDiffusion::Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)
{
    EffectiveThermalConductivity.Reallocate(newNx, newNy, newNz);
    EffectiveHeatCapacity.Reallocate(newNx, newNy, newNz);
    TxOld.Reallocate(newNx, newNy, newNz);
    dTx.Reallocate(newNx, newNy, newNz);
    Qdot.Reallocate(newNx, newNy, newNz);

    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    Info::WriteStandard(thisclassname, "Remeshed");
}

void HeatDiffusion::SetLatentHeat(const PhaseField& Phase)
{
    /** This function accounts for the latent heat release due to
        phase-transformations. */

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Qdot,0,)
    if(Phase.Interface(i,j,k))
    {
        for(auto it  = Phase.FieldsDot(i,j,k).cbegin();
                 it != Phase.FieldsDot(i,j,k).cend(); ++it)
        if(fabs(it->value1) > DBL_EPSILON)
        {
            size_t PIdxA = Phase.FieldsStatistics[it->indexA].Phase;
            size_t PIdxB = Phase.FieldsStatistics[it->indexB].Phase;

            Qdot(i,j,k) += LatentHeat(PIdxA,PIdxB)*it->value1;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void HeatDiffusion::SetEffectiveProperties(const PhaseField& Phase, const Temperature& Tx)
{
    SetLatentHeat(Phase);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EffectiveThermalConductivity,EffectiveThermalConductivity.Bcells(),)
    {
        EffectiveHeatCapacity(i,j,k) = 0.0;
        EffectiveThermalConductivity(i,j,k) = 0.0;

        for(size_t alpha = 0; alpha < Phase.Nphases; ++alpha)
        if(Phase.Fractions(i,j,k)({alpha}) > 0.0)
        {
            EffectiveHeatCapacity(i,j,k)
            += PhaseVolumetricHeatCapacity[alpha]*Phase.Fractions(i,j,k)({alpha});

            EffectiveThermalConductivity(i,j,k)
            += PhaseThermalConductivity[alpha]*Phase.Fractions(i,j,k)({alpha});
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void HeatDiffusion::SetEffectiveProperties(const PhaseField& Phase, const Temperature& Tx, const FlowSolverLBM& FL)
{
    SetLatentHeat(Phase);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EffectiveThermalConductivity,EffectiveThermalConductivity.Bcells(),)
    {
        EffectiveHeatCapacity(i,j,k) = 0.0;
        EffectiveThermalConductivity(i,j,k) = 0.0;

        for(size_t alpha = 0; alpha < Phase.Nphases; ++alpha)
        if(Phase.Fractions(i,j,k)({alpha}) > 0.0)
        {
        	if (Phase.PhaseAggregateStates[alpha] == AggregateStates::Liquid or
        			Phase.PhaseAggregateStates[alpha] == AggregateStates::Gas)
        	{
        	    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
        	    {
        	    	PhaseVolumetricHeatCapacity[alpha] = PhaseSpecificHeatCapacity[alpha] * FL.DensityWetting(i,j,k)({n});
        	    	PhaseThermalConductivity[alpha]    = FL.nut(i,j,k)({n}) * PhaseVolumetricHeatCapacity[alpha] /Tx.Pr;

        	    }
        	}
        	EffectiveHeatCapacity(i,j,k)
            		+= PhaseVolumetricHeatCapacity[alpha]*Phase.Fractions(i,j,k)({alpha});

        	EffectiveThermalConductivity(i,j,k)
        		+= PhaseThermalConductivity[alpha]*Phase.Fractions(i,j,k)({alpha});

        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

int HeatDiffusion::SolveImplicit(const PhaseField& Phase,
                                  const BoundaryConditions& BC, Temperature& Temp,
                                  const double dt)
{
    /** This function calculates the heat conductivity in the simulation domain
    using an implicit solver method. */

    /* Populate "old" temperature for the iterations. */

    if(Temp.ExtensionsActive)
    {
        if(Temp.ExtensionX0.isActive()) {Temp.ExtensionX0.store_temporary();}
        if(Temp.ExtensionXN.isActive()) {Temp.ExtensionXN.store_temporary();}
        if(Temp.ExtensionY0.isActive()) {Temp.ExtensionY0.store_temporary();}
        if(Temp.ExtensionYN.isActive()) {Temp.ExtensionYN.store_temporary();}
        if(Temp.ExtensionZ0.isActive()) {Temp.ExtensionZ0.store_temporary();}
        if(Temp.ExtensionZN.isActive()) {Temp.ExtensionZN.store_temporary();}
    }
    Temp.SetMinMaxAvg();
    TxOld = Temp.Tx;

    int iteration = 0;                                                          // Iterations counter
    double residual = 0.0;                                                      // Residual for calculating the accuracy of the converged solution
    const double dt_dx2 = dt/(dx*dx);
    const double dimension = 2.0*double(dNx + dNy + dNz);                       // Laplacian stencil dimension parameter

    do
    {
        iteration++;
        residual = 0.0;

        if(Temp.ExtensionsActive)
        {
            /* Calculation of heat diffusion in 1D extension using Gauss–Seidel implicit method. */
#ifdef MPI_PARALLEL
            if(MPI_RANK == 0)
#endif
            if (Temp.ExtensionX0.isActive())
            {
                Temp.ExtensionX0.PerformImplicitIteration(Temp, *this, residual, dt);
                Temp.ExtensionX0.setBC();
            }
#ifdef MPI_PARALLEL
            if(MPI_RANK == MPI_SIZE - 1)
#endif
            if (Temp.ExtensionXN.isActive())
            {
                Temp.ExtensionXN.PerformImplicitIteration(Temp, *this, residual, dt);
                Temp.ExtensionXN.setBC();
            }
            if (Temp.ExtensionY0.isActive())
            {
                Temp.ExtensionY0.PerformImplicitIteration(Temp, *this, residual, dt);
                Temp.ExtensionY0.setBC();
            }
            if (Temp.ExtensionYN.isActive())
            {
                Temp.ExtensionYN.PerformImplicitIteration(Temp, *this, residual, dt);
                Temp.ExtensionYN.setBC();
            }
            if (Temp.ExtensionZ0.isActive())
            {
                Temp.ExtensionZ0.PerformImplicitIteration(Temp, *this, residual, dt);
                Temp.ExtensionZ0.setBC();
            }
            if (Temp.ExtensionZN.isActive())
            {
                Temp.ExtensionZN.PerformImplicitIteration(Temp, *this, residual, dt);
                Temp.ExtensionZN.setBC();
            }
        }

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Temp.Tx,0,reduction(max:residual))
        {
            /* Calculation of heat diffusion using Jacobi implicit method. */

            double RhoCp      = EffectiveHeatCapacity(i,j,k);                   // Volumetric heat capacity [J/(m^3 K)]
            double Lambda     = EffectiveThermalConductivity(i,j,k);            // Thermal conductivity [J/(m s K)]
            double locQdot    = Qdot(i,j,k);                                    // Heat source [J/(m^3 s)]
            double locStencil = 0.0;                                            // Temperature stencil from the updated field [K]

            if (dNx > 0) locStencil += (Temp(i+1,j,k) + Temp(i-1,j,k));
            if (dNy > 0) locStencil += (Temp(i,j+1,k) + Temp(i,j-1,k));
            if (dNz > 0) locStencil += (Temp(i,j,k+1) + Temp(i,j,k-1));

            dTx(i,j,k) = (RhoCp*TxOld(i,j,k) + Lambda*locStencil*dt_dx2 + locQdot*dt)
                        /(RhoCp + dimension*Lambda*dt_dx2) - Temp(i,j,k);

            residual = max(residual,dTx(i,j,k)*dTx(i,j,k));
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        /* Updating the temperature with the calculated increment.*/

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Temp.Tx,0,)
        {
            Temp(i,j,k) += dTx(i,j,k);
            dTx(i,j,k) = 0.0;
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        /* Assigning values for the boundary cells. */
        Temp.SetBoundaryConditions(BC);

        /* Warning output, if number of iterations exceeds MaxIterations. */

        residual = sqrt(residual)/(Temp.Tavg + Tolerance); // Added Tolerance in denominator in case if Tavg == 0

#ifdef MPI_PARALLEL
        double rresidual = residual;
        MPI_Allreduce(&rresidual, &residual, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

        if(iteration == MaxIterations)
        {
            stringstream value;
            value.precision(16);
            value << residual;

            string tmp = "Maximum number of " + to_string(MaxIterations) + " iterations reached.\n"
                       + "Solver failed to converge to the requested tolerance of " + to_string(Tolerance) + ".\n"
                       + "Current convergence residual is " + value.str() + ".\n"
                       + "Simulation will continue!\n";

            Info::WriteWarning(tmp,thisclassname,"SolveImplicit()");
            break;
        }
    }
    while (residual >= Tolerance);

    /* Evaluate temperature statistics. */
    Temp.SetMinMaxAvg();

    if(Temp.ExtensionsActive)
    {
        if(Temp.ExtensionX0.isActive()) {Temp.ExtensionX0.Qdot = 0.0;}
        if(Temp.ExtensionXN.isActive()) {Temp.ExtensionXN.Qdot = 0.0;}
        if(Temp.ExtensionY0.isActive()) {Temp.ExtensionY0.Qdot = 0.0;}
        if(Temp.ExtensionYN.isActive()) {Temp.ExtensionYN.Qdot = 0.0;}
        if(Temp.ExtensionZ0.isActive()) {Temp.ExtensionZ0.Qdot = 0.0;}
        if(Temp.ExtensionZN.isActive()) {Temp.ExtensionZN.Qdot = 0.0;}
    }
    Qdot.Clear();
    return iteration;
}
}// namespace openphase
