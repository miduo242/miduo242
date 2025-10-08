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
 *   Main contributors :   Oleg Shchyglo; Dmitri Medvedev; Amol Subhedar;
 *                         Marvin Tegeler; Raphael Schiedung
 *
 */

#include "Base/CommonFunctions.h"
#include "Base/UserInterface.h"
#include "BoundaryConditions.h"
#include "Composition.h"
#include "FluidDynamics/BenziGas.h"
#include "FluidDynamics/D3Q27.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "FluidDynamics/VanDerWaalsGas.h"
#include "Info.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Temperature.h"
#include "VTK.h"
#include "Velocities.h"

namespace openphase
{
using namespace std;

FlowSolverLBM::FlowSolverLBM(const Settings& locSettings, double in_dt, const std::string InputFileName)
{
    Initialize(locSettings,in_dt);
    ReadInput(InputFileName);
}

void FlowSolverLBM::Initialize(const Settings& locSettings, double in_dt)
{
    thisclassname = "FlowSolverLBM";

    Nx      = locSettings.Nx;
    Ny      = locSettings.Ny;
    Nz      = locSettings.Nz;

    dNx     = locSettings.dNx;
    dNy     = locSettings.dNy;
    dNz     = locSettings.dNz;

    dx      = locSettings.dx;
    dt      = in_dt;
    Nphases = locSettings.Nphases;
    Ncomp   = (locSettings.Ncomp > 0) ? locSettings.Ncomp - 1 : 0;

    Pth    = 1.0e05;
    PthOld = 1.0e05;
    Pth0   = 1.0e05;
    Poutlet = 1e-6;

    cs2 = lbcs2*dx*dx/dt/dt;

    ObstaclesChanged = false;

    ReadInput(DefaultInputFileName); // NOTE N_Fluid_Comp needs to be read

    size_t Bcells = locSettings.Bcells;
    Obstacle.Allocate               (Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
    ObstacleAppeared.Allocate       (Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
    ObstacleChangedDensity.Allocate (Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
    ObstacleVanished.Allocate       (Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
    DensityWetting.Allocate         (Nx, Ny, Nz, dNx, dNy, dNz, {N_Fluid_Comp}, Bcells);
    ForceDensity.Allocate           (Nx, Ny, Nz, dNx, dNy, dNz, {N_Fluid_Comp}, Bcells);
    MomentumDensity.Allocate        (Nx, Ny, Nz, dNx, dNy, dNz, {N_Fluid_Comp}, Bcells);
    lbPopulations.Allocate          (Nx, Ny, Nz, dNx, dNy, dNz, {N_Fluid_Comp}, Bcells);
    lbPopulationsTMP.Allocate       (Nx, Ny, Nz, dNx, dNy, dNz, {N_Fluid_Comp}, Bcells);
    nut.Allocate                    (Nx, Ny, Nz, dNx, dNy, dNz, {N_Fluid_Comp}, Bcells);
    HydroDynPressure.Allocate       (Nx, Ny, Nz, dNx, dNy, dNz, {N_Fluid_Comp}, Bcells);
    DivergenceVel.Allocate          (Nx, Ny, Nz, dNx, dNy, dNz, {N_Fluid_Comp}, Bcells);

    switch(locSettings.ActiveDimensions())
    {
        case 1:
        {
            for(int ii = 0; ii < 3; ii++)
            for(int jj = 0; jj < 3; jj++)
            for(int kk = 0; kk < 3; kk++)
            {
                lbWeights[ii][jj][kk] = LBStencil1D[ii][jj][kk];
            }
            break;
        }
        case 2:
        {
            for(int ii = 0; ii < 3; ii++)
            for(int jj = 0; jj < 3; jj++)
            for(int kk = 0; kk < 3; kk++)
            {
                lbWeights[ii][jj][kk] = LBStencil2D[ii][jj][kk];
            }
            break;
        }
        case 3:
        {
            for(int ii = 0; ii < 3; ii++)
            for(int jj = 0; jj < 3; jj++)
            for(int kk = 0; kk < 3; kk++)
            {
                lbWeights[ii][jj][kk] = LBStencil3D[ii][jj][kk];
            }
            break;
        }
        default:
        {

        }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    {
        Obstacle               (i,j,k) = false;
        ObstacleAppeared       (i,j,k) = false;
        ObstacleChangedDensity (i,j,k) = false;
        ObstacleVanished       (i,j,k) = false;

        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            DensityWetting      (i,j,k)({n}) = 0.0;
            HydroDynPressure    (i,j,k)({n}) = 0.0;
            DivergenceVel       (i,j,k)({n}) = 0.0;
            ForceDensity        (i,j,k)({n}).set_to_zero();
            MomentumDensity     (i,j,k)({n}).set_to_zero();
            lbPopulations       (i,j,k)({n}).set_to_zero();
            lbPopulationsTMP    (i,j,k)({n}).set_to_zero();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    VTKDir = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;

    initialized = true;
    Info::Write(thisclassname, "Initialized");
}

void FlowSolverLBM::ReadInput(const std::string InputFileName)
{
    Info::WriteLineInsert("FlowSolverLBM input");
    Info::WriteStandard("Source", InputFileName);

    fstream inpF(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inpF)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };
    stringstream inp;
    inp << inpF.rdbuf();

    ReadInput(inp);

    inpF.close();
}

void FlowSolverLBM::ReadInput(std::stringstream& inp)
{
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    N_Fluid_Comp = UserInterface::ReadParameterD(inp, moduleLocation, std::string("N_FLUID_COMP"));
    nu.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        std::stringstream NUstr;
        NUstr << "NU[" << n << "]";
        nu[n] = UserInterface::ReadParameterD(inp, moduleLocation, NUstr.str());
    }

    Do_BounceBack        = UserInterface::ReadParameterB(inp, moduleLocation, std::string("BOUNCEBACK"));
    Do_BounceBackElastic = UserInterface::ReadParameterB(inp, moduleLocation, std::string("BBELASTIC"), false, false);
    Do_Buoyancy          = UserInterface::ReadParameterB(inp, moduleLocation, std::string("BUOYANCY"));
    Do_Drag              = UserInterface::ReadParameterB(inp, moduleLocation, std::string("DRAG"));
    Do_Gravitation       = UserInterface::ReadParameterB(inp, moduleLocation, std::string("GRAVITY"));
    Do_ThermalComp       = UserInterface::ReadParameterB(inp, moduleLocation, std::string("THERMALCOMP"));
    Do_SolidSolid        = UserInterface::ReadParameterB(inp, moduleLocation, std::string("SOLIDSOLID"));
    Do_StickySolids      = UserInterface::ReadParameterB(inp, moduleLocation, std::string("STICKYSOLIDS"), Do_SolidSolid, false);
    Do_Kupershtokh       = UserInterface::ReadParameterB(inp, moduleLocation, std::string("KUPERSHTOKH"), false, false);
    Do_Benzi             = UserInterface::ReadParameterB(inp, moduleLocation, std::string("BENZI"),       false, false);
    Do_TwoPhase          = UserInterface::ReadParameterB(inp, moduleLocation, std::string("TWOPHASE"),  !Do_Benzi and !Do_Kupershtokh, Do_Benzi or Do_Kupershtokh);
    Do_Benzi             = Do_Benzi or (Do_TwoPhase and !Do_Kupershtokh);
    Do_GuoForcing        = UserInterface::ReadParameterB(inp, moduleLocation, std::string("GUO_FORCING"), false,  Do_TwoPhase);
    Do_EDForcing         = UserInterface::ReadParameterB(inp, moduleLocation, std::string("ED_FORCING"),  false, !Do_TwoPhase);
    Do_FixPopulations    = UserInterface::ReadParameterB(inp, moduleLocation, std::string("bFixPop"),  false, true);
    Do_FluidRedistribution_Apearing  = UserInterface::ReadParameterB(inp, moduleLocation, std::string("bFluidRedistributionApearing"),   false, true);
    Do_FluidRedistribution_Vanishing = UserInterface::ReadParameterB(inp, moduleLocation, std::string("bFluidRedistributionVanishing"),  false, true);
    FluidRedistributionRange         = UserInterface::ReadParameterI(inp, moduleLocation, std::string("FluidRedistributionRange"),       false, Do_TwoPhase ? 20 : 0);

    h_star    = UserInterface::ReadParameterD(inp, moduleLocation, std::string("H_STAR"), Do_Drag, 0.0);
    ParaKuper = UserInterface::ReadParameterD(inp, moduleLocation, std::string("ParaKuper"), false, -0.0152);

    // Determine density discretization
    // NOTE that dRho is the equilibrium fluid density [Kg/m^2] if no liquid-vapor phase separation is considered!
    dRho = UserInterface::ReadParameterD(inp, moduleLocation, std::string("dRho"));
    dM   = dRho*dx*dx*dx;
    dP   = dM/dx/dt/dt;
    dm   = dM/dx/dx/dt;
    df   = dM/dx/dx/dt/dt;
    dnu  = dx*dx/dt;

    if (Ncomp > 0)
    {
        drhodc.resize(Ncomp);
        for (size_t n = 0; n < Ncomp; ++n)
        {
            std::stringstream DRHODCstr;
            DRHODCstr << "DRHODC[" << n << "]";
            drhodc[n] = UserInterface::ReadParameterD(inp, moduleLocation, DRHODCstr.str(), Do_Buoyancy,0.0);
        }
    }

    for (size_t i = 0; i < 3; i++)
    {
        std::stringstream GAstr;
        GAstr << "G0[" << i << "]";
        GA[i] = UserInterface::ReadParameterD(inp, moduleLocation, GAstr.str(), Do_Gravitation, 0.0);
    }

    // Read Benzi interaction parameter
    Gb.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        Gb[n].resize(N_Fluid_Comp);
        for (size_t m = 0; m < N_Fluid_Comp; ++m)
        {
            std::stringstream Gbstr;
            Gbstr << "GB[" << n << "][" << m << "]";
            Gb[n][m] = UserInterface::ReadParameterD(inp, moduleLocation, Gbstr.str(), Do_Benzi, (n==m) ? 1 : 0);
        }
    }

    // Read Kupershtokh interaction parameter
    lbGK.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        lbGK[n].resize(N_Fluid_Comp);
        for (size_t m = 0; m < N_Fluid_Comp; ++m)
        {
            std::stringstream GKstr;
            GKstr << "lbGK[" << n << "][" << m << "]";
            lbGK[n][m] = UserInterface::ReadParameterD(inp, moduleLocation, GKstr.str(), Do_Kupershtokh, (n==m) ? 1 : 0);
        }
    }

    rho_0.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        std::stringstream RHO0str;
        RHO0str << "RHO0[" << n << "]";
        rho_0[n] = UserInterface::ReadParameterD(inp, moduleLocation, RHO0str.str(), Do_Benzi, dRho);
    }

    CriticalDensity.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        std::stringstream RHOCstr;
        RHOCstr << "RHOC[" << n << "]";
        CriticalDensity[n] = UserInterface::ReadParameterD(inp, moduleLocation, RHOCstr.str(), Do_Kupershtokh, dRho);
    }

    lbTemperature.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        std::stringstream TEMPstr;
        TEMPstr << "TEMP[" << n << "]";
        lbTemperature[n] = UserInterface::ReadParameterD(inp, moduleLocation, TEMPstr.str(), Do_Kupershtokh, 0.0);
    }

    lbCriticalTemperature.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
       std::stringstream TEMPCstr;
       TEMPCstr << "TEMPC[" << n << "]";
       lbCriticalTemperature[n] = UserInterface::ReadParameterD(inp, moduleLocation, TEMPCstr.str(), Do_Kupershtokh, 0.0);
    }

    GasParameter.resize(N_Fluid_Comp);
    lbCriticalPressure.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        std::stringstream TEMPCstr;
        TEMPCstr << "CRITP[" << n << "]";
        lbCriticalPressure[n] = UserInterface::ReadParameterD(inp, moduleLocation, TEMPCstr.str(), Do_Kupershtokh, 0.0);
        GasParameter[n]       = lbCriticalPressure[n]/CriticalDensity[n]*dRho*(dt*dt)/(dx*dx);
        Info::Write("Gas Parameter["+std::to_string(n)+"]", GasParameter[n]);
    }

    LiquidDensity .resize(N_Fluid_Comp);
    VaporDensity  .resize(N_Fluid_Comp);
    SurfaceTension.resize(N_Fluid_Comp);
    InterfaceWidth.resize(N_Fluid_Comp);
    if (Do_Benzi)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        const BenziGas::EquilibriumValues_t Equilibrium =
            BenziGas::EquilibriumValues(Gb[n][n]/dRho);

        LiquidDensity    [n] = Equilibrium.lbLiquidDensity*dRho;
        VaporDensity     [n] = Equilibrium.lbVaporDensity*dRho;
        SurfaceTension   [n] = Equilibrium.lbSurfaceTension*dM/dt/dt;
        InterfaceWidth   [n] = 2.0; //TODO calculate

        Info::Write("Vapor density   ["+std::to_string(n)+"]", VaporDensity   [n]);
        Info::Write("Liquid density  ["+std::to_string(n)+"]", LiquidDensity  [n]);
        Info::Write("Surface tension ["+std::to_string(n)+"]", SurfaceTension [n]);
        Info::Write("Interface width ["+std::to_string(n)+"]", InterfaceWidth [n]);
    }
    else if (Do_Kupershtokh)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        const VanDerWaalsGas::EquilibriumValues_t Equilibrium =
            VanDerWaalsGas::EquilibriumValues(lbTemperature[n], lbCriticalTemperature[n]);

        LiquidDensity    [n] = Equilibrium.LiquidDensity*dRho;
        VaporDensity     [n] = Equilibrium.VaporDensity*dRho;
        SurfaceTension   [n] = 0.0*dM/dt/dt; //TODO
        InterfaceWidth   [n] = 0.402819*std::pow(GasParameter[n],-0.470509);

        Info::Write("Vapor density   ["+std::to_string(n)+"]", VaporDensity   [n]);
        Info::Write("Liquid density  ["+std::to_string(n)+"]", LiquidDensity  [n]);
        Info::Write("Surface tension ["+std::to_string(n)+"]", SurfaceTension [n]);
        Info::Write("Interface width ["+std::to_string(n)+"]", InterfaceWidth [n]);
    }
    else
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        LiquidDensity  [n] = dRho;
        VaporDensity   [n] = dRho;
        InterfaceWidth [n] = 1.0;
    }

    U0X  = UserInterface::ReadParameterD(inp, moduleLocation, std::string("U0X"), false, 0.0);
    U0Y  = UserInterface::ReadParameterD(inp, moduleLocation, std::string("U0Y"), false, 0.0);
    U0Z  = UserInterface::ReadParameterD(inp, moduleLocation, std::string("U0Z"), false, 0.0);
    Pth0 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Pth0"), false, 1.0e5);
    Poutlet = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Poutlet"), false, 1.0e5);
    Pth    = Pth0;
    PthOld = Pth0;
    Poutlet = Poutlet - 1.0e5 + 1.0e-6;

    Wetting.resize(N_Fluid_Comp);
    if (Do_TwoPhase)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        Wetting[n].resize(Nphases);
        for (size_t m = 0; m < Nphases; ++m)
        {
            std::stringstream TEMPCstr;
            TEMPCstr << "Wetting[" << n << "][" << m << "]";
            Wetting[n][m] = UserInterface::ReadParameterD(inp, moduleLocation, TEMPCstr.str());
        }
    }

    FluidMass.resize(N_Fluid_Comp, 0.0);

    // NOTE: the minimum and maximum kinematic viscosity values are based
    // on experience and may be adjusted.
    const double nu_min = 0.05*dx*dx/dt;
    const double nu_max = 1.00*dx*dx/dt;
    const double nu_opt = 0.5/3.0*dx*dx/dt; // optimal value !

    // Set relaxation parameter tau
    lbtau.resize(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        if ((nu[n] < nu_min) or (nu[n] > nu_max))
        {
            std::stringstream message;
            message << "Bad kinematic viscosity"
                    <<            "nu_lb = "<< std::scientific << nu[n]
                    << "\n\t  (Min nu_lb = "<< std::scientific << nu_min
                    <<      "; Opt nu_lb = "<< std::scientific << nu_opt
                    <<      "; Max nu_lb = "<< std::scientific << nu_max
                    <<")!\n\t  Simulation may be unstable! Adjust dx or dt.";
            Info::WriteWarning(message.str(), thisclassname, "Initialize");
        }

        // Set relaxation parameter tau
        lbtau [n] = 3*nu[n]/dx/dx*dt + 0.5;
        Info::WriteStandard("lbtau["+std::to_string(n)+"]", lbtau[n]);
    }

    Info::WriteLine();
    Info::WriteBlankLine();
}

void FlowSolverLBM::Remesh(int newNx, int newNy, int newNz,
                           const BoundaryConditions& BC)
{
    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    DensityWetting  .Reallocate(Nx, Ny, Nz);
    Obstacle        .Reallocate(Nx, Ny, Nz);
    MomentumDensity .Reallocate(Nx, Ny, Nz);
    ForceDensity    .Reallocate(Nx, Ny, Nz);
    lbPopulations   .Reallocate(Nx, Ny, Nz);
    lbPopulationsTMP.Reallocate(Nx, Ny, Nz);
}

void FlowSolverLBM::Read(const BoundaryConditions& BC, const int tStep)
{

#ifdef MPI_PARALLEL
    std::string FileName = UserInterface::MakeFileName(RawDataDir, thisclassname+'_'+std::to_string(MPI_RANK)+"_", tStep, ".dat");
#else
    std::string FileName = UserInterface::MakeFileName(RawDataDir, thisclassname+'_', tStep, ".dat");
#endif

    std::ifstream inp(FileName.c_str(), std::ios::in | std::ios::binary);

    if (!inp)
    {
        std::stringstream message;
        message << "File \"" << FileName << "\" could not be opened";
        Info::WriteExit(message.str(), thisclassname, "Read()");
#ifdef DEBUG
        throw std::runtime_error(message.str());
#else
        std::exit(EXIT_FAILURE);
#endif
    };

    int locNx = Nx;
    int locNy = Ny;
    int locNz = Nz;
    inp.read(reinterpret_cast<char*>(&Nx), sizeof(int));
    inp.read(reinterpret_cast<char*>(&Ny), sizeof(int));
    inp.read(reinterpret_cast<char*>(&Nz), sizeof(int));

    if(locNx != Nx or locNy != Ny or locNz != Nz)
    {
        std::stringstream message;
        message << "Inconsistent system dimensions!\n"
                << "Input data dimensions: (" << locNx
                << ", " << locNy << ", " << locNz << ") grid points.\n"
                << "Required data dimensions: (" << Nx
                << ", " << Ny << ", " << Nz << ") grid points.\n";
#ifdef DEBUG
        throw std::runtime_error(message.str());
#else
        std::exit(EXIT_FAILURE);
#endif
    }

    size_t locN_Fluid_Comp = N_Fluid_Comp;
    inp.read(reinterpret_cast<char*>(&N_Fluid_Comp), sizeof(size_t));

    if(locN_Fluid_Comp != N_Fluid_Comp)
    {
        std::stringstream message;
        message << "Inconsistent number of fluid components\n"
                << "Input number of fluid components: "
                << locN_Fluid_Comp << "\n"
                << "Required number of fluid components: "
                << N_Fluid_Comp << "\n";
#ifdef DEBUG
        throw std::runtime_error(message.str());
#else
        std::exit(EXIT_FAILURE);
#endif
    }

    STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,0)
    {
        for(size_t n = 0; n < N_Fluid_Comp; ++n)
        for(int ii = -dNx; ii <= dNx; ++ii)
        for(int jj = -dNy; jj <= dNy; ++jj)
        for(int kk = -dNz; kk <= dNz; ++kk)
        {
            inp.read(reinterpret_cast<char*>(&lbPopulations(i,j,k)({n})(ii,jj,kk)), sizeof(double));
        }
    }
    STORAGE_LOOP_END

    STORAGE_LOOP_BEGIN(i,j,k,ForceDensity,0)
    {
        for(size_t n = 0; n < N_Fluid_Comp; ++n)
        for(int m = 0; m < 3; ++m)
        {
            inp.read(reinterpret_cast<char*>(&ForceDensity(i,j,k)({n})[m]), sizeof(double));
        }
    }
    STORAGE_LOOP_END

    STORAGE_LOOP_BEGIN(i,j,k,Obstacle,0)
    {
        inp.read(reinterpret_cast<char*>(&Obstacle(i,j,k)), sizeof(bool));
    }
    STORAGE_LOOP_END

    CalculateDensityAndMomentum();
    SetBoundaryConditions(BC);

    Info::WriteStandard(thisclassname, "Binary input loaded");
}

void FlowSolverLBM::SetUniformVelocity(const BoundaryConditions& BC,
        const dVector3 U0)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        DensityWetting  (i,j,k)({n}) = dRho;
        MomentumDensity (i,j,k)({n}) = U0*DensityWetting(i,j,k)({n});
        lbPopulations   (i,j,k)({n}) = EquilibriumDistribution(1.0, lbWeights, MomentumDensity(i,j,k)({n})/dm);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    SetBoundaryConditions(BC);
}

void FlowSolverLBM::Write(const int tStep) const
{
#ifdef MPI_PARALLEL
    std::string FileName = UserInterface::MakeFileName(RawDataDir, thisclassname+'_'+std::to_string(MPI_RANK)+"_", tStep, ".dat");
#else
    std::string FileName = UserInterface::MakeFileName(RawDataDir, thisclassname+'_', tStep, ".dat");
#endif

    std::ofstream out(FileName.c_str(), std::ios::out | std::ios::binary);

    if (!out)
    {
        std::stringstream message;
        message << "File \"" << FileName
            << "\" could not be created! Terminating!!!" << std::endl;
        throw std::runtime_error(message.str());
    };

    out.write(reinterpret_cast<const char*>(&Nx),           sizeof(int));
    out.write(reinterpret_cast<const char*>(&Ny),           sizeof(int));
    out.write(reinterpret_cast<const char*>(&Nz),           sizeof(int));
    out.write(reinterpret_cast<const char*>(&N_Fluid_Comp), sizeof(size_t));

    STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,0)
    {
        for(size_t n = 0; n < N_Fluid_Comp; ++n)
        for(int ii = -dNx; ii <= dNx; ++ii)
        for(int jj = -dNy; jj <= dNy; ++jj)
        for(int kk = -dNz; kk <= dNz; ++kk)
        {
            const double value = lbPopulations(i,j,k)({n})(ii,jj,kk);
            out.write(reinterpret_cast<const char*>(&value), sizeof(double));
        }
    }
    STORAGE_LOOP_END

    STORAGE_LOOP_BEGIN(i,j,k,ForceDensity,0)
    {
        for(size_t n = 0; n < N_Fluid_Comp; ++n)
        for(int m = 0; m < 3; ++m)
        {
            const double value = ForceDensity(i,j,k)({n})[m];
            out.write(reinterpret_cast<const char*>(&value), sizeof(double));
        }
    }
    STORAGE_LOOP_END

    STORAGE_LOOP_BEGIN(i,j,k,Obstacle,0)
    {
        const bool value = Obstacle(i,j,k);
        out.write(reinterpret_cast<const char*>(&value), sizeof(bool));
    }
    STORAGE_LOOP_END

    out.close();
}

void FlowSolverLBM::WriteVTK(const int tStep, const PhaseField& Phase, const Settings& locSettings, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        ListOfFields.push_back((VTK::Field_t) {"DensityWetting_" +std::to_string(n), [n, this, &Phase](int i,int j,int k){return DensityWetting(i,j,k)({n})/dRho;}});
        ListOfFields.push_back((VTK::Field_t) {"Density_"        +std::to_string(n), [n, this, &Phase](int i,int j,int k){return Density(Phase,i,j,k,n);}});
        ListOfFields.push_back((VTK::Field_t) {"Velocity_"       +std::to_string(n), [n, this]        (int i,int j,int k){return Velocity(i,j,k,n);}});
    }
    ListOfFields.push_back((VTK::Field_t) {"Pressure", [this](int i,int j,int k){return Pressure(i,j,k);}});
    ListOfFields.push_back((VTK::Field_t) {"Obstacle", [this](int i,int j,int k){return int(Obstacle(i,j,k));}});

    std::string Filename = UserInterface::MakeFileName(VTKDir, thisclassname+'_', tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void FlowSolverLBM::lbWriteVTK(const int tStep, const PhaseField& Phase, const Settings& locSettings, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        ListOfFields.push_back((VTK::Field_t) {"DensityWetting_" +std::to_string(n), [n, this, &Phase](int i,int j,int k){return DensityWetting(i,j,k)({n})/dRho;}});
        ListOfFields.push_back((VTK::Field_t) {"Density_"        +std::to_string(n), [n, this, &Phase](int i,int j,int k){return Density(Phase,i,j,k,n)/dRho;}});
        ListOfFields.push_back((VTK::Field_t) {"Velocity_"       +std::to_string(n), [n, this]        (int i,int j,int k){return Velocity(i,j,k,n)*dt/dx;}});
    }
    ListOfFields.push_back((VTK::Field_t) {"Pressure", [this](int i,int j,int k){return Pressure(i,j,k)/dP;}});
    ListOfFields.push_back((VTK::Field_t) {"Obstacle", [this](int i,int j,int k){return int(Obstacle(i,j,k));}});

    std::string Filename = UserInterface::MakeFileName(VTKDir, thisclassname+"_lattice_units_", tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void FlowSolverLBM::SetBoundaryConditions(const BoundaryConditions& BC)
{
    if(dNx) BC.SetXVector(lbPopulations);
    if(dNy) BC.SetYVector(lbPopulations);
    if(dNz) BC.SetZVector(lbPopulations);

    if(dNx) BC.SetX(DensityWetting);
    if(dNy) BC.SetY(DensityWetting);
    if(dNz) BC.SetZ(DensityWetting);

    if(dNx) BC.SetXVector(MomentumDensity);
    if(dNy) BC.SetYVector(MomentumDensity);
    if(dNz) BC.SetZVector(MomentumDensity);
}

void FlowSolverLBM::CalculateForceTwoPhase(const int i, const int j,
        const int k, PhaseField& Phase)
{
    for(int ii = -dNx; ii <= dNx; ++ii)
    for(int jj = -dNy; jj <= dNy; ++jj)
    for(int kk = -dNz; kk <= dNz; ++kk)
    {
        const dVector3 vel {ii*dx/dt,jj*dx/dt,kk*dx/dt};

        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        for (size_t m = 0; m < N_Fluid_Comp; ++m)
        {
            double tmp = 0.0;
            if (Do_Benzi)
            {
                // Benzi et. al. (2006) - Mesoscopic modeling of a two-phase flow
                // in the presence of boundaries
                tmp = - Gb[n][m]/dRho * psi(DensityWetting(i,j,k)({n}), rho_0[n]) * lbWeights[ii+1][jj+1][kk+1] * psi(DensityWetting(i+ii,j+jj,k+kk)({m}), rho_0[m]);
            }
            else if (Do_Kupershtokh)
            {
                // Calculate phase sepration force according to:
                // Kupershtokh, Medvedev and Karpov, 2009, On equations of state in a lattice Boltzmann method
                // Kupershtokh, Medvedev and Gribanov, 2017, Thermal lattice Boltzmann method for multiphase flows

                // Evaluate equation of state (Van der Waals)
                const double p_000 = VanDerWaalsGas::ReducedPressure(DensityWetting(i   ,j   ,k   )({n})/dRho, lbTemperature[n], CriticalDensity[n]/dRho, lbCriticalTemperature[n]);
                const double p_ijk = VanDerWaalsGas::ReducedPressure(DensityWetting(i+ii,j+jj,k+kk)({m})/dRho, lbTemperature[m], CriticalDensity[m]/dRho, lbCriticalTemperature[m]);

                // Calculate interaction potential
                const double Phi_000 = Phi(DensityWetting(i,   j,   k   )({n})/dRho, p_000, GasParameter[n]);
                const double Phi_ijk = Phi(DensityWetting(i+ii,j+jj,k+kk)({m})/dRho, p_ijk, GasParameter[m]);

                // Interpolate between both way
                tmp = 6.0*lbWeights[ii+1][jj+1][kk+1]*lbGK[n][m]*(ParaKuper*(Phi_ijk*Phi_ijk) + (1.0-2.0*ParaKuper)*(Phi_000*Phi_ijk));
            }

            const dVector3 BenziForceDensity = {tmp*ii*df, tmp*jj*df, tmp*kk*df};

            ForceDensity(i,j,k)({n}) += BenziForceDensity;

            if (!Obstacle(i,j,k) && Obstacle(i+ii,j+jj,k+kk))
            for(auto& it : Phase.Fields(i+ii,j+jj,k+kk))
            {
                Grain& grain = Phase.FieldsStatistics[it.index];
                if(grain.State == AggregateStates::Solid)
                {
                    const dVector3 pos = {double(i+ii), double(j+jj), double(k+kk)};
                    dVector3 distanceCM;
                    CommonFunctions::CalculateDistancePeriodic(pos, grain.Rcm, distanceCM, Nx, Ny, Nz);
                    const dVector3 R = distanceCM * dx;

                    #ifdef _OPENMP
                    #pragma omp critical
                    #endif
                    {
                        grain.Force  -= BenziForceDensity *  dx*dx*dx;
                        grain.Torque -= R.cross(BenziForceDensity) * dx*dx*dx;
                    }
                }
            }
        }
    }
}

void FlowSolverLBM::CalculateForceGravitation(PhaseField& Phase)
{
    // Calculate gravity on solid bodies
    for (size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    if (Phase.FieldsStatistics[idx].State == AggregateStates::Solid and
        Phase.FieldsStatistics[idx].Mobile)
    {
        Phase.FieldsStatistics[idx].Acm += GA;
    }

    // Calculate gravity on fluid nodes
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    if (not Obstacle(i,j,k))
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        ForceDensity(i,j,k)({n}) += GA * DensityWetting(i,j,k)({n});
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::CalculateDensityTC(Temperature& Tx, Composition& Cx, const BoundaryConditions& BC)
{
    double Ts=273.0;
    double mus=1.68e-05; //it's for air
    double S=110.5;
    double AW_air = Cx.AtomicWeightMixture;
    double R= PhysicalConstants::R;     ///<  Universal gas constant  j/mol k
    double Rm=R/AW_air;                 ///<  gas constant for air j/kg k

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        DensityWetting(i,j,k)({n}) = Pth/(Rm * Tx.Tx(i,j,k));
        double mu=mus*(Tx.Tx(i,j,k)/Ts)*sqrt(Tx.Tx(i,j,k)/Ts)*(Ts+S)/(Tx.Tx(i,j,k)+S); ///< Dynamic viscosity is calculated based on the Sutherland's Law
        nut(i,j,k)({n}) = mu/DensityWetting(i,j,k)({n});
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::SetInitialPopulationsTC(BoundaryConditions& BC)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    if (!Obstacle(i,j,k))
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        HydroDynPressure(i,j,k)({n}) = Poutlet;
        MomentumDensity(i,j,k)({n}).set_to_zero();

        dVector3 vel = MomentumDensity(i,j,k)({n})/dm/(DensityWetting(i,j,k)({n})/dRho*lbcs2) + ForceDensity(i,j,k)({n})*df * 0.5/(DensityWetting(i,j,k)({n})/dRho);

        double u2 = (vel*vel);
        double Feq;
        double geq;

        for(int ii = -dNx; ii <= dNx; ++ii)
        for(int jj = -dNy; jj <= dNy; ++jj)
        for(int kk = -dNz; kk <= dNz; ++kk)
        {
            double cu = ii*vel[0] + jj*vel[1] + kk*vel[2];
            Feq = DensityWetting(i,j,k)({n})/dRho * lbWeights[ii+1][jj+1][kk+1]*(1.0 - 1.5*u2 + cu*(3.0 + 4.5*cu));
            geq = Feq * lbcs2 + (HydroDynPressure(i,j,k)({n}) / dP - DensityWetting(i,j,k)({n})/dRho * lbcs2) * lbWeights[ii+1][jj+1][kk+1] ;

            lbPopulations(i,j,k)({n})(ii,jj,kk) = geq;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    if(dNx) BC.SetXVector(lbPopulations);
    if(dNy) BC.SetYVector(lbPopulations);
    if(dNz) BC.SetZVector(lbPopulations);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,HydroDynPressure,HydroDynPressure.Bcells(),)
    if (!Obstacle(i,j,k))
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        //for (int ii = -dNx; ii <= dNx; ++ii)
        //for (int jj = -dNy; jj <= dNy; ++jj)
        //for (int kk = -dNz; kk <= dNz; ++kk)
        //{
            HydroDynPressure (i,j,k)({n}) = Poutlet;
        //}
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::CalculateHydrodynamicPressureAndMomentum( Temperature& Tx,  Composition& Cx, Velocities& Vel)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells()-1,)
    if(!Obstacle(i,j,k))
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        HydroDynPressure(i,j,k)({n}) = 0.0;
        MomentumDensity(i,j,k)({n}).set_to_zero();

        for (int ii = -dNx; ii <= dNx; ++ii)
        for (int jj = -dNy; jj <= dNy; ++jj)
        for (int kk = -dNz; kk <= dNz; ++kk)
        {
            double rr = lbPopulations(i,j,k)({n})(ii,jj,kk);
            HydroDynPressure (i,j,k)({n})  += rr*dP;
            MomentumDensity(i,j,k)({n})[0] += rr*ii*dm;
            MomentumDensity(i,j,k)({n})[1] += rr*jj*dm;
            MomentumDensity(i,j,k)({n})[2] += rr*kk*dm;
        }

        dVector3 vel = (MomentumDensity(i,j,k)({n}) * 3.0 + ForceDensity(i,j,k)({n})*dt/2.0)/DensityWetting(i,j,k)({n});
        cs2 = lbcs2 * dx/dt * dx/dt;
        
        HydroDynPressure (i,j,k)({n})   =   cs2*dt/2.0 * (
                                            (vel[0]* (DensityWetting(i+1,j,k)({n})-DensityWetting(i-1,j,k)({n}))/2.0/dx )
                                        +   (vel[1]* (DensityWetting(i,j+1,k)({n})-DensityWetting(i,j-1,k)({n}))/2.0/dx )
                                        +   (vel[2]* (DensityWetting(i,j,k+1)({n})-DensityWetting(i,j,k-1)({n}))/2.0/dx )
                                        +             DensityWetting(i,j,k)({n}) * DivergenceVel(i,j,k)({n})        
                                                      ) + HydroDynPressure (i,j,k)({n});
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}


void FlowSolverLBM::CollisionTC(Temperature& Tx,  Composition& Cx, Velocities& Vel, PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    if (!Obstacle(i,j,k))
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        dVector3 lbvel = ((MomentumDensity(i,j,k)({n}) * 3.0 + ForceDensity(i,j,k)({n})*dt/2.0)/DensityWetting(i,j,k)({n})) *dt/dx;
        double u2 = lbvel[0]*lbvel[0]+lbvel[1]*lbvel[1]+lbvel[2]*lbvel[2];
        lbtau [n] = nut(i,j,k)({n})/dnu/lbcs2 + 0.5;
        for(int ii = -dNx; ii <= dNx; ii++)
        for(int jj = -dNy; jj <= dNy; jj++)
        for(int kk = -dNz; kk <= dNz; kk++)
        {
            double cu = ii*lbvel[0] + jj*lbvel[1] + kk*lbvel[2];
            double Feq = DensityWetting(i,j,k)({n})/dRho * lbWeights[ii+1][jj+1][kk+1]*(1.0 + cu/lbcs2- u2/(2.0*lbcs2) + cu*cu/(2.0*lbcs2*lbcs2) );
            double geq = Feq*lbcs2 + ( HydroDynPressure(i,j,k)({n})/dP - DensityWetting(i,j,k)({n})/dRho * lbcs2 ) * lbWeights[ii+1][jj+1][kk+1];

            lbPopulations(i,j,k)({n})(ii,jj,kk) = lbPopulations(i,j,k)({n})(ii,jj,kk) + 1.0/lbtau[n] * ( geq  - lbPopulations(i,j,k)({n})(ii,jj,kk) )
            + ((
           (ii-lbvel[0]) * ( ((DensityWetting(i+1,j,k)({n})-DensityWetting(i-1,j,k)({n}))/dRho/2.0)
                   * lbcs2 *(Feq/DensityWetting(i,j,k)({n})*dRho - lbWeights[ii+1][jj+1][kk+1]) + GA[0]*dt*dt/dx*Feq )
       +   (jj-lbvel[1]) * ( ((DensityWetting(i,j+1,k)({n})-DensityWetting(i,j-1,k)({n}))/dRho/2.0)
                   * lbcs2 *(Feq/DensityWetting(i,j,k)({n})*dRho - lbWeights[ii+1][jj+1][kk+1]) + GA[1]*dt*dt/dx*Feq )
       +   (kk-lbvel[2]) * ( ((DensityWetting(i,j,k+1)({n})-DensityWetting(i,j,k-1)({n}))/dRho/2.0)
                   * lbcs2 *(Feq/DensityWetting(i,j,k)({n})*dRho - lbWeights[ii+1][jj+1][kk+1]) + GA[2]*dt*dt/dx*Feq )
                                     )
                                    + DensityWetting(i,j,k)({n})/dRho * lbcs2 * DivergenceVel(i,j,k)({n})*dt *lbWeights[ii+1][jj+1][kk+1])* (1.0-1.0/(2.0*lbtau[n]));
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::CalculateForceGravitation(PhaseField& Phase,
        const Composition& Cx, const Temperature& Tx)
{
    double R= PhysicalConstants::R;     ///<  Universal gas constant  j/mol k
    double Aw_air = Cx.AtomicWeightMixture;
    double Rm=R/Aw_air;             ///< gas constant for air j/kg k
    double Tref=(Tx.TBC0X+Tx.TSphere)/2.0;

    double Rhoref= Pth0 /(Rm *Tref);
    if(Tx.TBC0X==Tx.TSphere)
    {
        Rhoref=0.0;
    }

    // Calculate gravity on fluid nodes bodies
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    {
        if (not Obstacle(i,j,k))
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            if (Do_ThermalComp)
            {
                ForceDensity(i,j,k)({n}) = GA * (DensityWetting(i,j,k)({n})- Rhoref);
            }
            else if (!Do_ThermalComp)
            {
                ForceDensity(i,j,k)({n}) = GA * DensityWetting(i,j,k)({n});
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::ApplyForces(PhaseField& Phase, const Velocities& Vel,
        const Composition& Cx, Temperature& Tx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ForceDensity,0,)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        ForceDensity(i,j,k)({n}).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    if (!Obstacle(i,j,k))
    {
        if (Do_TwoPhase) CalculateForceTwoPhase(i,j,k, Phase);
        if (Do_Drag and Phase.Interface(i,j,k))
        {
            CalculateForceDrag(i,j,k, Phase, Vel);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    if (Do_Gravitation) CalculateForceGravitation(Phase, Cx, Tx);
}

void FlowSolverLBM::CalculateForceBuoyancy(const int i, const int j,
        const int k, const PhaseField& Phase, const Composition& Cx)
{
    for(auto it = Phase.Fields(i,j,k).cbegin();
             it != Phase.Fields(i,j,k).cend(); it++)
    {
        const Grain& grain = Phase.FieldsStatistics[it->index];
        if(grain.State == AggregateStates::Liquid)
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        for (size_t Comp = 0; Comp < Ncomp; ++Comp)
        {
            ForceDensity(i,j,k)({n}) += GA * drhodc[Comp]*(Cx.MoleFractions(i,j,k)({grain.Phase, Comp}) - Cx.Initial({grain.Phase, Comp}));
        }
    }
}

void FlowSolverLBM::CalculateForceDrag(const int i,const int j,const int k,
        PhaseField& Phase, const Velocities& Vel)
{
    double dx3 = dx*dx*dx;

    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    if(DensityWetting(i,j,k)({n}) != 0.0)
    {
        double locSolidFraction = SolidFraction(i,j,k,Phase);
        const double mu = nu[n]*DensityWetting(i,j,k)({n});  //Dynamic viscosity [kg/m/s]
        for(auto& it : Phase.Fields(i,j,k))
        {
            Grain& grain = Phase.FieldsStatistics[it.index];
            if(grain.State == AggregateStates::Solid)
            {
                const double factor =
                    (it.value * it.value * (1.0 - locSolidFraction))/
                    (Phase.iWidth * Phase.iWidth);
                const dVector3 DragForceDensity = (Velocity(i,j,k,n) - Vel.Phase(i,j,k)({grain.Phase})) * mu * factor * h_star;
                const dVector3 pos = {double(i), double(j), double(k)};
                dVector3 distanceCM;
                CommonFunctions::CalculateDistancePeriodic(pos, grain.Rcm, distanceCM, Nx, Ny, Nz);
                const dVector3 locR = distanceCM * dx;
                #ifdef _OPENMP
                #pragma omp critical
                #endif
                {
                    grain.Force  += DragForceDensity * dx3;
                    grain.Torque += locR.cross(DragForceDensity) * dx3;
                }
            }
        }
    }
}

double FlowSolverLBM::BounceBack(const int i, const int j, const int k,
        const int ii, const int jj, const int kk, const size_t n,
        PhaseField& Phase, const BoundaryConditions& BC, double& lbDensityChange)
{
    double dx3 = dx*dx*dx;

    double NewPopulation = lbPopulations(i,j,k)({n})(-ii,-jj,-kk);
    if (Do_BounceBack)
    {
        for(auto it : Phase.Fields(i-ii,j-jj,k-kk))
        {
            Grain& grain = Phase.FieldsStatistics[it.index];
            if(grain.State == AggregateStates::Solid)
            {
                dMatrix3x3 W;
                W(0,0) = 0.0;
                W(1,1) = 0.0;
                W(2,2) = 0.0;
                W(0,1) = -grain.aVel[2];
                W(0,2) =  grain.aVel[1];
                W(1,2) = -grain.aVel[0];
                W(1,0) =  grain.aVel[2];
                W(2,0) = -grain.aVel[1];
                W(2,1) =  grain.aVel[0];

                const dVector3 pos = {double(i-0.5*ii), double(j-0.5*jj), double(k-0.5*kk)};
                dVector3 distanceCM;
                CommonFunctions::CalculateDistancePeriodic(pos, grain.Rcm, distanceCM, Nx, Ny, Nz);
                const dVector3 locR = distanceCM * dx;
                const dVector3 Vel = grain.Vcm + W * locR;

                const double tmp = lbWeights[ii+1][jj+1][kk+1] *
                    DensityWetting(i,j,k)({n})/dRho * (Vel[0]*ii + Vel[1]*jj + Vel[2]*kk) * dt/dx;

                const double lbBBDensity =
                    2.0 * (lbPopulations(i,j,k)({n})(-ii,-jj,-kk) + 3.0 * tmp);

                NewPopulation += it.value * 6.0 * tmp;

                // NOTE: Bounce Back Momentum Density is now in physical units
                dVector3 lbBounceBackForceDensity;
                lbBounceBackForceDensity[0] = - ii * lbBBDensity * it.value;
                lbBounceBackForceDensity[1] = - jj * lbBBDensity * it.value;
                lbBounceBackForceDensity[2] = - kk * lbBBDensity * it.value;

                const dVector3 BounceBackForceDensity = lbBounceBackForceDensity * df;

                #ifdef _OPENMP
                #pragma omp critical
                #endif
                {
                    grain.Force  += BounceBackForceDensity * dx3;
                    grain.Torque += locR.cross(BounceBackForceDensity) * dx3;
                }
            }
        }
    }
    else if (Do_BounceBackElastic)
    {
        // Author: raphael.schiedung@rub.de
        const double rho1    = lbPopulations(i,j,k)({n})(-ii,-jj,-kk);
        const double ci      = std::sqrt(ii*ii+jj*jj+kk*kk);
        const double ici     = 1./ci;
        const double p1      = rho1*ci;
        const dVector3 nn    = {-ii*ici, -jj*ici, -kk*ici};
        const double Elastic = 1.00;  // 0.0 models fully inelastic impact

        NewPopulation = 0.0;
        for(auto it : Phase.Fields(i-ii,j-jj,k-kk))
        {
            Grain& grain = Phase.FieldsStatistics[it.index];
            if(grain.State == AggregateStates::Solid)
            {
                dMatrix3x3 W;
                W(0,0) = 0.0;
                W(1,1) = 0.0;
                W(2,2) = 0.0;
                W(0,1) = -grain.aVel[2];
                W(0,2) =  grain.aVel[1];
                W(1,2) = -grain.aVel[0];
                W(1,0) =  grain.aVel[2];
                W(2,0) = -grain.aVel[1];
                W(2,1) =  grain.aVel[0];

                const dVector3 pos = {double(i-0.5*ii), double(j-0.5*jj), double(k-0.5*kk)};
                dVector3 distanceCM;
                CommonFunctions::CalculateDistancePeriodic(pos, grain.Rcm,distanceCM, Nx, Ny, Nz);
                const dVector3 locR   = distanceCM * dx;
                const dVector3 Vel = grain.Vcm + W * locR;

                const double v2   = Vel*nn*dt/dx;
                const double rho2 = grain.Density/dRho;
                const double p2   = rho2*v2;
                const double p1p  = (p1*(rho1-rho2*Elastic) + p2*(rho1+rho1*Elastic))/(rho1+rho2);
                const double p2p  = (p1*(rho2+rho2*Elastic) + p2*(rho2-rho1*Elastic))/(rho1+rho2);

                if (p1p < 0)
                {
                    NewPopulation -= it.value*p1p*ici;
                }
                else
                {
                    // The solid bounce back from the fluid (locally).
                    // Check grain Densities !
                    std::cerr << "Warning fluid momentum is too high!"
                              << "Check parameter (eg. Densities).\n";
                }

                dVector3 BounceBackForceDensity = {0,0,0};
                BounceBackForceDensity += nn*it.value*(p2p-p2)*df;
                #ifdef _OPENMP
                #pragma omp critical
                #endif
                {
                    grain.Force  += BounceBackForceDensity * dx3;
                    grain.Torque += locR.cross(BounceBackForceDensity) * dx3;
                }
            }
        }
    }
    lbDensityChange += NewPopulation - lbPopulations(i,j,k)({n})(-ii,-jj,-kk);
    return NewPopulation;
}

void FlowSolverLBM::Propagation(PhaseField& Phase, const BoundaryConditions& BC)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulationsTMP,0,)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        lbPopulationsTMP(i,j,k)({n}).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,0,)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        double lbDensityChange = 0.0;
        for(int ii = -dNx; ii <= dNx; ++ii)
        for(int jj = -dNy; jj <= dNy; ++jj)
        for(int kk = -dNz; kk <= dNz; ++kk)
        {
            if (Obstacle(i, j, k))
            {
                lbPopulationsTMP(i,j,k)({n})(ii,jj,kk) =
                    lbPopulations(i,j,k)({n})(ii,jj,kk);
            }
            else if (Obstacle(i-ii, j-jj, k-kk))
            {
                lbPopulationsTMP(i,j,k)({n})(ii,jj,kk) =
                    BounceBack(i, j, k, ii, jj, kk, n, Phase, BC, lbDensityChange);
            }
            else
            {
                lbPopulationsTMP(i,j,k)({n})(ii,jj,kk) =
                    lbPopulations(i-ii, j-jj, k-kk)({n})(ii,jj,kk);
            }
        }
        lbPopulationsTMP(i,j,k)({n})(0,0,0) -= lbDensityChange; //NoSlip
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,0,)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        lbPopulations(i,j,k)({n}) = lbPopulationsTMP(i,j,k)({n});
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    if(dNx) BC.SetXVector(lbPopulations);
    if(dNy) BC.SetYVector(lbPopulations);
    if(dNz) BC.SetZVector(lbPopulations);
}

void FlowSolverLBM::ApplyForces(PhaseField& Phase, const Velocities& Vel, const Composition& Cx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ForceDensity,0,)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        ForceDensity(i,j,k)({n}).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    if (!Obstacle(i,j,k))
    {
        if (Do_TwoPhase) CalculateForceTwoPhase(i,j,k, Phase);
        if (Do_Buoyancy) CalculateForceBuoyancy(i,j,k, Phase, Cx);
        if (Do_Drag and Phase.Interface(i,j,k))
        {
            CalculateForceDrag(i,j,k, Phase, Vel);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    if (Do_Gravitation) CalculateForceGravitation(Phase);
}

void FlowSolverLBM::ApplyForces(PhaseField& Phase, const Velocities& Vel)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ForceDensity,0,)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        ForceDensity(i,j,k)({n}).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    if (!Obstacle(i,j,k))
    {
        if (Do_TwoPhase) CalculateForceTwoPhase(i,j,k, Phase);
        if (Do_Drag and Phase.Interface(i,j,k))
        {
            CalculateForceDrag(i,j,k, Phase, Vel);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    if (Do_Gravitation) CalculateForceGravitation(Phase);
}

void FlowSolverLBM::Collision()
{
    // Apply body force
    if (Do_GuoForcing) // Guo, Zheng, and, Shi, Phy. Rev. E (2002)
    {
        // Apply body force with force projection method
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
        if (!Obstacle(i,j,k))
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            const dVector3 EqMomentum = MomentumDensity(i,j,k)({n}) + ForceDensity(i,j,k)({n})*dt/2;
            lbPopulations(i,j,k)({n}) += (EquilibriumDistribution(DensityWetting(i,j,k)({n})/dRho, lbWeights, EqMomentum/dm) - lbPopulations(i,j,k)({n}))/lbtau[n];
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        // Apply body force with force projection method
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
        if (!Obstacle(i,j,k))
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            assert(DensityWetting(i,j,k)({n})/dRho > DBL_EPSILON);
            assert(DensityWetting(i,j,k)({n}) < DBL_MAX);

            const dVector3 lbVel = Velocity(i,j,k,n)*dt/dx;

            D3Q27 lblocForcing;
            lblocForcing.set_to_zero();
            for(int ii = -dNx; ii <= dNx; ++ii)
            for(int jj = -dNy; jj <= dNy; ++jj)
            for(int kk = -dNz; kk <= dNz; ++kk)
            {
                lblocForcing(ii,jj,kk) =
                    (1.0-0.5/lbtau[n]) * lbWeights[ii+1][jj+1][kk+1]*
                    (3.0*((ii-lbVel[0])*ForceDensity(i,j,k)({n})[0]/df +
                          (jj-lbVel[1])*ForceDensity(i,j,k)({n})[1]/df +
                          (kk-lbVel[2])*ForceDensity(i,j,k)({n})[2]/df) +
                     9.0*(lbVel[0]*ii + lbVel[1]*jj + lbVel[2]*kk)*
                     (ii*ForceDensity(i,j,k)({n})[0]/df +
                      jj*ForceDensity(i,j,k)({n})[1]/df +
                      kk*ForceDensity(i,j,k)({n})[2]/df));
            }

            lbPopulations(i,j,k)({n}) += lblocForcing;
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    else if (Do_EDForcing)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
        if (!Obstacle(i,j,k))
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            const double q = 1.0 - 1.0/lbtau[n];
            lbPopulations(i,j,k)({n}) = (lbPopulations(i,j,k)({n}) -
                                         EquilibriumDistribution(DensityWetting(i,j,k)({n})/dRho, lbWeights,  MomentumDensity(i,j,k)({n})/dm))*q +
                                         EquilibriumDistribution(DensityWetting(i,j,k)({n})/dRho, lbWeights, (MomentumDensity(i,j,k)({n})/dm + ForceDensity(i,j,k)({n})/df*dt));
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

size_t FlowSolverLBM::CountObstacleNodes(void) const
{
    size_t locObstacleNodes = 0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0, reduction(+:locObstacleNodes))
    {
        if (Obstacle(i,j,k))
        {
            locObstacleNodes++;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    #ifdef MPI_PARALLEL
    size_t tmp = locObstacleNodes;
    MPI_Allreduce(&tmp, &locObstacleNodes, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    #endif

    return locObstacleNodes;
}

size_t FlowSolverLBM::CountFluidNodes(void) const
{
    size_t locFluidNodes = 0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0, reduction(+:locFluidNodes))
    {
        if (not Obstacle(i,j,k)) locFluidNodes++;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    #ifdef MPI_PARALLEL
    size_t tmp = locFluidNodes;
    MPI_Allreduce(&tmp, &locFluidNodes, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    #endif

    return locFluidNodes;
}

void FlowSolverLBM::EnforceMassConservation(void)
{
    // Calculate Fluid mass
    std::vector<double> DeltaDensity = CalculateFluidMass();
    size_t FluidNodes = CountFluidNodes();

    if (FluidNodes > 0)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    if (FluidMass[n] != 0)
    {
        DeltaDensity[n] -= FluidMass[n];
        DeltaDensity[n] /= FluidNodes;

        // Fix density
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
        {
            if (not Obstacle(i,j,k))
            {
                lbPopulations(i,j,k)({n}) = lbPopulations(i,j,k)({n}) -
                    EquilibriumDistribution(DensityWetting(i,j,k)({n})/dRho + DeltaDensity[n], lbWeights) +
                    EquilibriumDistribution(DensityWetting(i,j,k)({n})/dRho, lbWeights);
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    FluidMass = CalculateFluidMass();
}

void FlowSolverLBM::EnforceSolidMomentum(PhaseField& Phase, const dVector3 value)
{
    // Fix Solid momentum
    const double dx = Phase.dx;
    const double dV = dx*dx*dx;

    dVector3 SolidMomentum = {0.0,0.0,0.0};
    double SolidMass = 0.0;
    for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    if (Phase.FieldsStatistics[idx].State == AggregateStates::Solid and
        Phase.FieldsStatistics[idx].Mobile)
    {
        const double GrainMass = Phase.FieldsStatistics[idx].Volume*dV*
                                 Phase.FieldsStatistics[idx].Density;
        SolidMass     += GrainMass;
        SolidMomentum += Phase.FieldsStatistics[idx].Vcm * GrainMass;
    }

    const dVector3 VcmFix = (SolidMomentum - value)/SolidMass;

    for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    if (Phase.FieldsStatistics[idx].State == AggregateStates::Solid and
        Phase.FieldsStatistics[idx].Mobile)
    {
        Phase.FieldsStatistics[idx].Vcm -= VcmFix;
    }
}

bool FlowSolverLBM::SingleSolid(const int i, const int j, const int k, const PhaseField& Phase) const
{
    for(auto it = Phase.Fields(i,j,k).cbegin();
            it != Phase.Fields(i,j,k).cend(); ++it)
    if(Phase.FieldsStatistics[it->index].State == AggregateStates::Solid)
    {
        if(it->value >= 0.95) return true;
    }
    return false;
}

double FlowSolverLBM::SolidFraction(const int i, const int j, const int k, const PhaseField& Phase) const
{
    double SolidFraction = 0.0;
    for(auto it = Phase.Fields(i,j,k).cbegin();
            it != Phase.Fields(i,j,k).cend(); ++it)
    if(Phase.FieldsStatistics[it->index].State == AggregateStates::Solid)
    {
        SolidFraction += it->value;
    }
    return SolidFraction;
}

bool FlowSolverLBM::LocalObstacle(const int i, const int j, const int k, const PhaseField& Phase) const
{
    if (Do_StickySolids) return (SolidFraction(i,j,k,Phase) >= 0.95);
    else return SingleSolid(i,j,k,Phase);
}

void FlowSolverLBM::DetectObstacles(const PhaseField& Phase)
{
    ObstaclesChanged = false;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    {
        if (LocalObstacle(i,j,k,Phase))
        {
            if(not Obstacle(i,j,k)) ObstaclesChanged = true;
            Obstacle(i,j,k) = true;
        }
        else
        {
            if(Obstacle(i,j,k)) ObstaclesChanged = true;
            Obstacle(i,j,k) = false;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    #ifdef MPI_PARALLEL
    int tmp = ObstaclesChanged;
    MPI_Allreduce(&tmp, &(ObstaclesChanged), 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    #endif
}

void FlowSolverLBM::DetectObstaclesAdvection(const PhaseField& Phase,
        const Velocities &Vel, const BoundaryConditions& BC)
{
    #ifdef DEBUG
    dVector3 FluidMomentumOld = CalculateFluidMomentum(Phase)[0];
    double FluidMassOld = CalculateFluidMass()[0];
    #endif

    std::size_t ObstaclesAppeared = 0;
    std::size_t ObstaclesVanished = 0;

    bool locObstaclesChanged = false;
    // TODO Problem with periodic boundaries!
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(), reduction (+:ObstaclesAppeared,ObstaclesVanished))
    {
        ObstacleVanished(i,j,k) = false;
        ObstacleAppeared(i,j,k) = false;
        if (LocalObstacle(i,j,k,Phase))
        {
            if(not Obstacle(i,j,k))
            {
                ObstacleAppeared(i,j,k) = true;
                locObstaclesChanged = true;
                ObstaclesAppeared++;
            }
            Obstacle(i,j,k) = true;
        }
        else
        {
            if(Obstacle(i,j,k))
            {
                ObstacleVanished(i,j,k) = true;
                locObstaclesChanged = true;
                ObstaclesVanished++;
            }
            Obstacle(i,j,k) = false;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    ObstaclesChanged = locObstaclesChanged;
    #ifdef MPI_PARALLEL
    MPI_Allreduce(&locObstaclesChanged, &ObstaclesChanged, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
    #endif

    if (ObstaclesAppeared + ObstaclesVanished > 0)
    {
        SetFluidNodesNearObstacle();
        SetObstacleNodes(Phase,Vel);
        CalculateDensityAndMomentum();
        SetBoundaryConditions(BC);

        #ifdef DEBUG
        dVector3 FluidMomentumNew = CalculateFluidMomentum(Phase)[0];
        dVector3 FluidMomentumDelta = FluidMomentumNew-FluidMomentumOld;
        double FluidMomentumDeltaMag = FluidMomentumDelta.abs();
        double FluidMassNew = CalculateFluidMass()[0];
        double FluidMassDelta = FluidMassNew - FluidMassOld;
        if (std::abs(FluidMassDelta)/dRho > 1.0e-10)//DBL_EPSILON)
        {
            Info::WriteExit("Mass conservation violated! Method needs to be debugged.", thisclassname, "DetectObstaclesAdvection");
            Info::Write("Obstacles Appeared", ObstaclesAppeared);
            Info::Write("Obstacles Vanished", ObstaclesVanished);
            Info::Write("Flud Mass Change [Kg]", FluidMassDelta);
            Info::Write("Flud Mass Change  [1]", FluidMassDelta/dRho);
            std::exit(EXIT_FAILURE);
        }
        if (FluidMomentumDeltaMag/dm/dx/dx/dx > 1.0e-10)//DBL_EPSILON)
        {
            Info::WriteExit("Momentum conservation violated! Method needs to be debugged.", thisclassname, "DetectObstaclesAdvection");
            Info::Write("Obstacles Appeared",ObstaclesAppeared);
            Info::Write("Obstacles Vanished",ObstaclesVanished);
            Info::Write("Fluid Momentum Change [Kg m/s]",FluidMomentumDeltaMag);
            Info::Write("Fluid Momentum Change      [1]",FluidMomentumDeltaMag/dm/dx/dx/dx);
            std::exit(EXIT_FAILURE);
        }
        #endif
    }
}

std::vector<dVector3> FlowSolverLBM::CalculateFluidMomentum(const PhaseField& Phase) const
{
    std::vector<dVector3> Momentum(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        std::vector<double> tmpMomentumX;
        std::vector<double> tmpMomentumY;
        std::vector<double> tmpMomentumZ;

        // Calculate fluid momentum
        #pragma omp parallel
        {
            std::vector<double> tmplocMomentumX;
            std::vector<double> tmplocMomentumY;
            std::vector<double> tmplocMomentumZ;
            #pragma omp for collapse(OMP_COLLAPSE_LOOPS)\
            schedule(dynamic,OMP_DYNAMIC_CHUNKSIZE)
            for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++)
            if (not Obstacle(i,j,k))
            for (int ii = -dNx; ii <= dNx; ii++)
            for (int jj = -dNy; jj <= dNy; jj++)
            for (int kk = -dNz; kk <= dNy; kk++)
            {
                double rr = lbPopulations(i,j,k)({n})(ii,jj,kk);
                tmplocMomentumX.push_back(rr*ii);
                tmplocMomentumY.push_back(rr*jj);
                tmplocMomentumZ.push_back(rr*kk);
            }

            std::sort(tmplocMomentumX.begin(),tmplocMomentumX.end(), [] (double a, double b) { return a < b;});
            std::sort(tmplocMomentumY.begin(),tmplocMomentumY.end(), [] (double a, double b) { return a < b;});
            std::sort(tmplocMomentumZ.begin(),tmplocMomentumZ.end(), [] (double a, double b) { return a < b;});

            double locMomentumX = 0.0;
            double locMomentumY = 0.0;
            double locMomentumZ = 0.0;

            for (auto value : tmplocMomentumX) locMomentumX += value;
            for (auto value : tmplocMomentumY) locMomentumY += value;
            for (auto value : tmplocMomentumZ) locMomentumZ += value;

            #pragma omp critical
            {
                tmpMomentumX.push_back(locMomentumX);
                tmpMomentumY.push_back(locMomentumY);
                tmpMomentumZ.push_back(locMomentumZ);
            }
        }

        std::sort(tmpMomentumX.begin(),tmpMomentumX.end(), [] (double a, double b) { return a < b;});
        std::sort(tmpMomentumY.begin(),tmpMomentumY.end(), [] (double a, double b) { return a < b;});
        std::sort(tmpMomentumZ.begin(),tmpMomentumZ.end(), [] (double a, double b) { return a < b;});

        for (auto value : tmpMomentumX) Momentum[n][0] += value;
        for (auto value : tmpMomentumY) Momentum[n][1] += value;
        for (auto value : tmpMomentumZ) Momentum[n][2] += value;

        #ifdef MPI_PARALLEL
        double tmpX = Momentum[n][0];
        double tmpY = Momentum[n][1];
        double tmpZ = Momentum[n][2];
        MPI_Allreduce(&tmpX, &(Momentum[n][0]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&tmpY, &(Momentum[n][1]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&tmpZ, &(Momentum[n][2]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        #endif

        Momentum[n] *= dx /dt *dRho*dx*dx*dx;
    }
    return Momentum;
}

void FlowSolverLBM::SetFluidNodesNearObstacle()
{
    if (lbPopulations.Bcells() < 2*FluidRedistributionRange+1)
    {
        Info::WriteExit("Too few boundary cells, at least "+std::to_string(2*FluidRedistributionRange+1)+" boundary cells are needed!", thisclassname, "DetectObstaclesAdvection");
        std::exit(EXIT_FAILURE);
    }
    //TODO use Initialisations::loop_sphere so that less boundary cells are needed

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,lbPopulations.Bcells()-FluidRedistributionRange,)
    {
        ObstacleChangedDensity(i,j,k) = false;
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            lbPopulationsTMP(i,j,k)({n}).set_to_zero();
            if (ObstacleVanished(i,j,k)) lbPopulations(i,j,k)({n}).set_to_zero();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    if (Do_FluidRedistribution_Vanishing)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,lbPopulations.Bcells()-FluidRedistributionRange,)
        if (ObstacleVanished(i,j,k))
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            // If an obstacle vanishes the present cell has to be filled with the
            // surrounding fluid from the neighbour cells

            // Calculate local average population
            double   SumWeights = 0;
            double   avDensity  = 0;
            dVector3 avMomentum = {0.0,0.0,0.0};
            int FluidNeighbours = 0;
            for (int ii = -FluidRedistributionRange*dNx; ii <= FluidRedistributionRange*dNx; ++ii)
            for (int jj = -FluidRedistributionRange*dNy; jj <= FluidRedistributionRange*dNy; ++jj)
            for (int kk = -FluidRedistributionRange*dNz; kk <= FluidRedistributionRange*dNz; ++kk)
            if (not Obstacle(i+ii,j+jj,k+kk)
                    and not ObstacleVanished(i+ii,j+jj,k+kk)
                    and ii*ii + jj*jj + kk*kk <= FluidRedistributionRange*FluidRedistributionRange
                    and ii + jj + kk != 0)
            {
                const double weight = DensityWetting(i+ii,j+jj,k+kk)({n});
                SumWeights += weight;
                avDensity  += DensityWetting (i+ii,j+jj,k+kk)({n})*weight;
                avMomentum += MomentumDensity(i+ii,j+jj,k+kk)({n})*weight;
                FluidNeighbours++;
            }

            if (FluidNeighbours)
            {
                avDensity  /= SumWeights;
                avMomentum /= SumWeights;

                DensityWetting  (i,j,k)({n}) = avDensity;
                MomentumDensity (i,j,k)({n}) = avMomentum;
                const D3Q27 DeltaPop = EquilibriumDistribution(avDensity/dRho, lbWeights, avMomentum/dm);
                lbPopulationsTMP(i,j,k)({n}) += DeltaPop;
                ObstacleChangedDensity(i,j,k) = true;

                // Subtract added fluid from neighbouring cells
                for (int ii = -FluidRedistributionRange*dNx; ii <= FluidRedistributionRange*dNx; ++ii)
                for (int jj = -FluidRedistributionRange*dNy; jj <= FluidRedistributionRange*dNy; ++jj)
                for (int kk = -FluidRedistributionRange*dNz; kk <= FluidRedistributionRange*dNz; ++kk)
                if (not Obstacle(i+ii,j+jj,k+kk)
                        and not ObstacleVanished(i+ii,j+jj,k+kk)
                        and ii*ii + jj*jj + kk*kk <= FluidRedistributionRange*FluidRedistributionRange
                        and ii + jj + kk != 0)
                {
                    #pragma omp critical
                    {
                        const double weight = DensityWetting(i+ii,j+jj,k+kk)({n});
                        lbPopulationsTMP(i+ii,j+jj,k+kk)({n}) -= DeltaPop*weight/SumWeights;
                        ObstacleChangedDensity(i+ii,j+jj,k+kk) = true;
                    }
                }
            }
            else std::cerr << "Warning no fluid Neighbours\n";
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    else
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,lbPopulations.Bcells()-FluidRedistributionRange,)
        if (ObstacleVanished(i,j,k))
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            ObstacleChangedDensity(i,j,k) = true;
            lbPopulations(i,j,k)({n}).set_to_zero();
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    if (Do_FluidRedistribution_Apearing)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,lbPopulations.Bcells()-FluidRedistributionRange,)
        if (ObstacleAppeared(i,j,k))
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            // If an obstacle appears the present fluid has to be moved to the
            // neighbouring cell

            // Count number of fluid neighbours
            double   SumWeights = 0;
            int FluidNeighbours = 0;
            for (int ii = -FluidRedistributionRange*dNx; ii <= FluidRedistributionRange*dNx; ++ii)
            for (int jj = -FluidRedistributionRange*dNy; jj <= FluidRedistributionRange*dNy; ++jj)
            for (int kk = -FluidRedistributionRange*dNz; kk <= FluidRedistributionRange*dNz; ++kk)
            if (not Obstacle(i+ii,j+jj,k+kk)
                    and not ObstacleAppeared(i+ii,j+jj,k+kk)
                    and ii*ii + jj*jj + kk*kk <= FluidRedistributionRange*FluidRedistributionRange
                    and ii + jj + kk != 0)
            {
                const double weight = DensityWetting(i+ii,j+jj,k+kk)({n});
                SumWeights += weight;
                FluidNeighbours++;
            }

            if (FluidNeighbours)
            {
                // Add fluid from this cell to the neighbouring cells
                const D3Q27 DeltaPop = lbPopulations(i,j,k)({n});
                for (int ii = -FluidRedistributionRange*dNx; ii <= FluidRedistributionRange*dNx; ++ii)
                for (int jj = -FluidRedistributionRange*dNy; jj <= FluidRedistributionRange*dNy; ++jj)
                for (int kk = -FluidRedistributionRange*dNz; kk <= FluidRedistributionRange*dNz; ++kk)
                if (not Obstacle(i+ii,j+jj,k+kk)
                       and not ObstacleAppeared(i+ii,j+jj,k+kk)
                       and ii*ii + jj*jj + kk*kk <= FluidRedistributionRange*FluidRedistributionRange
                       and ii + jj + kk != 0)
               {
                   #pragma omp critical
                   {
                       const double weight = DensityWetting(i+ii,j+jj,k+kk)({n});
                       lbPopulationsTMP(i+ii,j+jj,k+kk)({n}) += DeltaPop*weight/SumWeights;
                       ObstacleChangedDensity(i+ii,j+jj,k+kk) = true;
                   }
               }
           }
           else std::cerr << "Warning no fluid Neighbours\n";
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    // Add change in fluid nodes
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,0,)
    {
        if (ObstacleChangedDensity(i,j,k))
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            lbPopulations(i,j,k)({n}) += lbPopulationsTMP(i,j,k)({n});
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::FixPopulations(void)
{
    // Add change in fluid nodes
    if (Do_FixPopulations)
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,lbPopulations.Bcells(),)
    {
        if (not Obstacle(i,j,k))
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            bool FixPopulation = false;
            // Check for negative populations and fix them
            for (int ii = -dNx; ii <= dNx; ++ii)
            for (int jj = -dNy; jj <= dNy; ++jj)
            for (int kk = -dNz; kk <= dNz; ++kk)
            if (lbPopulations(i,j,k)({n})(ii,jj,kk) <= 0.0)
            {
                FixPopulation = true;
            }

            if (DensityWetting(i,j,k)({n}) <= 0.0)
            {
                FixPopulation = true;
                DensityWetting (i,j,k)({n}) = 0.0;
                MomentumDensity(i,j,k)({n}) = {0,0,0};
            }

            if (FixPopulation)
            {
                lbPopulations(i,j,k)({n}) =
                    EquilibriumDistribution(DensityWetting(i,j,k)({n})/dRho, lbWeights, MomentumDensity(i,j,k)({n})/dm);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::SetObstacleNodes(const PhaseField& Phase, const Velocities& Vel)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    if (Obstacle(i,j,k))
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        MomentumDensity(i,j,k)({n}).set_to_zero();
        if (Do_TwoPhase)
        {
            DensityWetting (i,j,k)({n}) = 0.0;
            for(auto it = Phase.Fields(i,j,k).cbegin();
                    it != Phase.Fields(i,j,k).cend(); it++)
            if (Phase.FieldsStatistics[it->index].State == AggregateStates::Solid)
            {
                size_t pIndex = Phase.FieldsStatistics[it->index].Phase;

                DensityWetting (i,j,k)({n}) += it->value * Wetting[n][pIndex];

                if (Phase.FieldsStatistics[it->index].Mobile)
                {
                    MomentumDensity(i,j,k)({n}) += Vel.Phase(i,j,k)({pIndex}) *
                        DensityWetting(i,j,k)({n}) * it->value;
                }
            }
        }
        else if (!Do_ThermalComp)
        {
            DensityWetting (i,j,k)({n}) = dRho;
        }
        lbPopulations(i,j,k)({n}) =
            EquilibriumDistribution(DensityWetting(i,j,k)({n})/dRho, lbWeights, MomentumDensity(i,j,k)({n})/dm);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::SetDivVelZeroNearObst()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    if (Obstacle(i,j,k))
    {
        for (int ii = -dNx; ii <= dNx; ii++)
        for (int jj = -dNy; jj <= dNy; jj++)
        for (int kk = -dNz; kk <= dNz; kk++)
        {
            DivergenceVel(i+ii,j+jj,k+kk)({0}) = 0.0;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::Solve(PhaseField& Phase, const Composition& Cx,
        Velocities& Vel, const BoundaryConditions& BC)
{
    for (size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    {
        Phase.FieldsStatistics[idx].Force  = {0,0,0};
        Phase.FieldsStatistics[idx].Torque = {0,0,0};
        Phase.FieldsStatistics[idx].Acm    = {0,0,0};
        Phase.FieldsStatistics[idx].aAcc   = {0,0,0};
    }

    DetectObstacles(Phase);
    SetObstacleNodes(Phase, Vel);
    if(ObstaclesChanged) CalculateDensityAndMomentum();

    Propagation(Phase, BC);
    CalculateDensityAndMomentum();

    // BEGIN Dmitry
    if(dNx) BC.SetX(DensityWetting);
    if(dNy) BC.SetY(DensityWetting);
    if(dNz) BC.SetZ(DensityWetting);
    // END Dmitry

    ApplyForces(Phase, Vel, Cx);
    Collision();
    CalculateDensityAndMomentum(); //Commented by Dmitry // Uncommented by Raphael
    SetBoundaryConditions(BC);

    CalculateFluidVelocities(Vel, Phase, BC);

    #ifdef MPI_PARALLEL
    for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    if(Phase.FieldsStatistics[idx].Exist and Phase.FieldsStatistics[idx].State == AggregateStates::Solid)
    {
        auto locForce = Phase.FieldsStatistics[idx].Force;
        MPI_Allreduce(&locForce[0], &(Phase.FieldsStatistics[idx].Force[0]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locForce[1], &(Phase.FieldsStatistics[idx].Force[1]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locForce[2], &(Phase.FieldsStatistics[idx].Force[2]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        auto locTorque = Phase.FieldsStatistics[idx].Torque;
        MPI_Allreduce(&locTorque[0], &(Phase.FieldsStatistics[idx].Torque[0]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locTorque[1], &(Phase.FieldsStatistics[idx].Torque[1]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locTorque[2], &(Phase.FieldsStatistics[idx].Torque[2]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    #endif
}

void FlowSolverLBM::Solve(PhaseField& Phase, Velocities& Vel,
        const BoundaryConditions& BC)
{
    for (size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    {
        Phase.FieldsStatistics[idx].Force  = {0,0,0};
        Phase.FieldsStatistics[idx].Torque = {0,0,0};
        Phase.FieldsStatistics[idx].Acm    = {0,0,0};
        Phase.FieldsStatistics[idx].aAcc   = {0,0,0};
    }

    DetectObstacles(Phase);
    SetObstacleNodes(Phase, Vel);
    if(ObstaclesChanged) CalculateDensityAndMomentum();

    Propagation(Phase,BC);
    CalculateDensityAndMomentum();

    // BEGIN Dmitry
    if(dNx) BC.SetX(DensityWetting);
    if(dNy) BC.SetY(DensityWetting);
    if(dNz) BC.SetZ(DensityWetting);
    // END Dmitry

    ApplyForces(Phase, Vel);
    Collision();
    CalculateDensityAndMomentum(); //Commented by Dmitry // Uncommented by Raphael
    SetBoundaryConditions(BC);

    CalculateFluidVelocities(Vel, Phase, BC);

    #ifdef MPI_PARALLEL
    for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    if(Phase.FieldsStatistics[idx].Exist and Phase.FieldsStatistics[idx].State == AggregateStates::Solid)
    {
        auto locForce = Phase.FieldsStatistics[idx].Force;
        MPI_Allreduce(&locForce[0], &(Phase.FieldsStatistics[idx].Force[0]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locForce[1], &(Phase.FieldsStatistics[idx].Force[1]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locForce[2], &(Phase.FieldsStatistics[idx].Force[2]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        auto locTorque = Phase.FieldsStatistics[idx].Torque;
        MPI_Allreduce(&locTorque[0], &(Phase.FieldsStatistics[idx].Torque[0]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locTorque[1], &(Phase.FieldsStatistics[idx].Torque[1]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locTorque[2], &(Phase.FieldsStatistics[idx].Torque[2]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    #endif
}

std::vector<double> FlowSolverLBM::CalculateFluidMass(void) const
{
    std::vector<double> Mass(N_Fluid_Comp);
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        std::vector<double> tmpMass;
        // Calculate fluid momentum
        #pragma omp parallel
        {
            std::vector<double> tmplocMass;
            #pragma omp for collapse(OMP_COLLAPSE_LOOPS) schedule(dynamic,OMP_DYNAMIC_CHUNKSIZE)
            for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++)
            if (not Obstacle(i,j,k))
            for (int ii = -dNx; ii <= dNx; ii++)
            for (int jj = -dNy; jj <= dNy; jj++)
            for (int kk = -dNz; kk <= dNz; kk++)
            {
                double rr = lbPopulations(i,j,k)({n})(ii,jj,kk);
                tmplocMass.push_back(rr);
            }

            std::sort(tmplocMass.begin(),tmplocMass.end(), [] (double a, double b) { return a < b;});

            double locMass = 0.0;
            for (auto value : tmplocMass) locMass += value;

            #pragma omp critical
            {
                tmpMass.push_back(locMass);
            }
        }

        std::sort(tmpMass.begin(),tmpMass.end(), [] (double a, double b) { return a < b;});

        for (auto value : tmpMass) Mass[n] += value;

        #ifdef MPI_PARALLEL
        auto tmp = Mass[n];
        MPI_Allreduce(&tmp, &(Mass[n]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        #endif

        Mass[n] *= dRho*dx*dx*dx;
    }
    return Mass;
}

void FlowSolverLBM::CalculateDensityAndMomentum(void)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    for (size_t n = 0; n < N_Fluid_Comp; ++n)
    {
        DensityWetting (i,j,k)({n}) = 0.0;
        MomentumDensity(i,j,k)({n}).set_to_zero();
        for (int ii = -dNx; ii <= dNx; ii++)
        for (int jj = -dNy; jj <= dNy; jj++)
        for (int kk = -dNz; kk <= dNz; kk++)
        {
            const dVector3 vel {ii*dx/dt,jj*dx/dt,kk*dx/dt};
            DensityWetting (i,j,k)({n}) +=     lbPopulations(i,j,k)({n})(ii,jj,kk)*dRho;
            MomentumDensity(i,j,k)({n}) += vel*lbPopulations(i,j,k)({n})(ii,jj,kk)*dRho;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    FixPopulations();
}

void FlowSolverLBM::CalculateFluidVelocities(Velocities& Vel,
        const PhaseField& Phase, const BoundaryConditions& BC) const
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    if(!Obstacle(i,j,k))
    {
        for(auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); it++)
        {
            const Grain& grain = Phase.FieldsStatistics[it->index];
            if(grain.State != AggregateStates::Solid)
            {
                Vel.Phase(i,j,k)({grain.Phase}) = Velocity(i,j,k,0);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Vel.CalculateAverage(Phase);
}

double FlowSolverLBM::Pressure(const int i, const int j, const int k) const
{
    double Pressure = 0.0;
    if(!Obstacle(i,j,k))
    {
        if (Do_Kupershtokh)
        {
            for (size_t n = 0; n < N_Fluid_Comp; ++n)
            if (DensityWetting(i,j,k)({n}) > DBL_EPSILON)
            {
                Pressure += lbCriticalPressure[n]*dP*VanDerWaalsGas::ReducedPressure(
                        DensityWetting(i,j,k)({n})/dRho, lbTemperature[n], CriticalDensity[n]/dRho,lbCriticalTemperature[n]);
            }
        }
        else if (Do_Benzi)
        {
            for (size_t n = 0; n < N_Fluid_Comp; ++n)
            if (DensityWetting(i,j,k)({n}) > DBL_EPSILON)
            {
                const double& rho   = DensityWetting(i,j,k)({n});
                const double locPsi = psi(rho,rho_0[n]);

                Pressure += cs2*rho + 0.5*dt*cs2*Gb[n][n]*locPsi*locPsi;
            }
        }
        else if (Do_ThermalComp)
        {
            for (size_t n = 0; n < N_Fluid_Comp; ++n)
            if (DensityWetting(i,j,k)({n}) > DBL_EPSILON)
            {
                Pressure += HydroDynPressure(i,j,k)({n});
            }
        }
        else
        {
            for (size_t n = 0; n < N_Fluid_Comp; ++n)
            if (DensityWetting(i,j,k)({n}) > DBL_EPSILON)
            {
                const double& rho = DensityWetting(i,j,k)({n});

                Pressure += cs2*rho;
            }
        }
    }
    return Pressure;
}

dMatrix3x3 FlowSolverLBM::PressureTensor(const int i, const int j, const int k) const
{
    dMatrix3x3 Pressure({0,0,0,0,0,0,0,0,0});
    if (Do_Kupershtokh)
    {
        dMatrix3x3 lbPressure({0,0,0,0,0,0,0,0,0});
        //Info::WriteWarning("Pressure tensor calculation is inaccurate", thisclassname, "PressureTensor");
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            double LaplacianPhi    = 0.0;
            double GradPhi[3]      = {0,0,0};
            const double weight[3] = {-0.5,0.0,0.5};
            for (int ii = -1; ii <= 1; ii++)
            {
                const double pi = VanDerWaalsGas::ReducedPressure(DensityWetting(i+ii,j,k)({n})/dRho, lbTemperature[n], CriticalDensity[n]/dRho, lbCriticalTemperature[n]);
                const double pj = VanDerWaalsGas::ReducedPressure(DensityWetting(i,j+ii,k)({n})/dRho, lbTemperature[n], CriticalDensity[n]/dRho, lbCriticalTemperature[n]);
                const double pk = VanDerWaalsGas::ReducedPressure(DensityWetting(i,j,k+ii)({n})/dRho, lbTemperature[n], CriticalDensity[n]/dRho, lbCriticalTemperature[n]);

                GradPhi[0] += weight[ii+1]*Phi(DensityWetting(i+ii,j,k)({n})/dRho, pi, GasParameter[n]);
                GradPhi[1] += weight[ii+1]*Phi(DensityWetting(i,j+ii,k)({n})/dRho, pj, GasParameter[n]);
                GradPhi[2] += weight[ii+1]*Phi(DensityWetting(i,j,k+ii)({n})/dRho, pk, GasParameter[n]);
            }

            for (int ii = -dNx; ii <= dNx; ii++)
            for (int jj = -dNy; jj <= dNy; jj++)
            for (int kk = -dNz; kk <= dNz; kk++)
            {
                const double pijk = VanDerWaalsGas::ReducedPressure(DensityWetting(i+ii,j+jj,k+kk)({n})/dRho, lbTemperature[n], CriticalDensity[n]/dRho, lbCriticalTemperature[n]);
                LaplacianPhi += LaplacianStencil3D_27a[ii+1][jj+1][kk+1]*
                    Phi(DensityWetting(i+ii,j+jj,k+kk)({n})/dRho, pijk , GasParameter[n]);
            }

            const double locPhi  = Phi(DensityWetting(i,j,k)({n})/dRho, VanDerWaalsGas::ReducedPressure(DensityWetting(i,j,k)({n})/dRho, lbTemperature[n], CriticalDensity[n]/dRho, lbCriticalTemperature[n]), GasParameter[n]);

            const double GradPhi2 =
                GradPhi[0]*GradPhi[0] +
                GradPhi[1]*GradPhi[1] +
                GradPhi[2]*GradPhi[2];

            const double lbp0 = lbcs2*DensityWetting(i,j,k)({n})/dRho
                - 6.0*Gb[n][n]/dRho*0.50*lbcs2*locPhi*locPhi
                - 6.0*Gb[n][n]/dRho*0.50*lbcs2*lbcs2*locPhi*LaplacianPhi
                - 6.0*Gb[n][n]/dRho*0.25*lbcs2*lbcs2*GradPhi2;

            lbPressure(0,0) += lbp0;
            lbPressure(1,1) += lbp0;
            lbPressure(2,2) += lbp0;

            for (int ii = 0; ii < 3; ii++)
            for (int jj = 0; jj < 3; jj++)
            {
                lbPressure(ii,jj) += 6.0*Gb[n][n]/dRho*0.5*lbcs2*lbcs2*GradPhi[ii]*GradPhi[jj];
            }
        }
        Pressure = lbPressure*dP;
    }
    else if (Do_Benzi)
    {
        //READ: Shan - 2008 - Pressure tensor calculation in a class of nonideal...
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            const double& rho = DensityWetting(i,j,k)({n});
            const double  p0  = cs2*rho;
            Pressure(0,0) += p0;
            Pressure(1,1) += p0;
            Pressure(2,2) += p0;

            for (size_t m = 0; m < N_Fluid_Comp; ++m)
            if (std::abs(Gb[n][m]) > DBL_EPSILON)
            for (int ii = -dNx; ii <= dNx; ii++)
            for (int jj = -dNy; jj <= dNy; jj++)
            for (int kk = -dNz; kk <= dNz; kk++)
            {
                const dVector3 vel {ii*dx/dt,jj*dx/dt,kk*dx/dt};

                Pressure += vel.dyadic(vel) * 0.5 * Gb[n][m] * psi(rho, rho_0[n]) * lbWeights[ii+1][jj+1][kk+1] * psi(DensityWetting(i+ii,j+jj,k)({m}), rho_0[m]);
            }
        }
    }
    else
    {
        for (size_t n = 0; n < N_Fluid_Comp; ++n)
        {
            const double p0 = cs2*DensityWetting(i,j,k)({n});
            Pressure(0,0) += p0;
            Pressure(1,1) += p0;
            Pressure(2,2) += p0;
        }
    }
    return Pressure;
}

double  FlowSolverLBM::OptimalParaKuper(const double ReducedTemperature,
        const double GasPrameter)
{
    const std::vector<std::vector<double>> OptimalValues =
        {{2.0000000000000018e-01,-1.4092641421153448e-01},
         {3.0000000000000016e-01,-1.4776128416548723e-01},
         {4.0000000000000013e-01,-1.4973625604844859e-01},
         {5.0000000000000011e-01,-1.5026704234216137e-01},
         {6.0000000000000009e-01,-1.5023912552780752e-01},
         {7.0000000000000007e-01,-1.4788430943930206e-01},
         {8.0000000000000004e-01,-1.4180410847718788e-01},
         {9.0000000000000002e-01,-7.8369694891999728e-02}};

    // Return smallest temperature if temperature is too small
    if (OptimalValues.front()[0] > ReducedTemperature) return OptimalValues.front()[1];

    // Search for optimal value
    std::vector<double> PreviousValue = OptimalValues[0];
    std::vector<double> Interpolated  = {0.0,0.0};
    for (auto& Value: OptimalValues)
    {
        if (Value[0] == ReducedTemperature)
        {
            return Value[1];
        }
        else if (ReducedTemperature < Value[0] and
                 ReducedTemperature > PreviousValue[0])
        {
            Interpolated = PreviousValue;

            const double dT = (Value[0]-PreviousValue[0]);
            const double DT = (ReducedTemperature-PreviousValue[0]);
            const double A  = DT/dT;

            Interpolated[0] = ReducedTemperature;
            Interpolated[1] += (Value[1] - PreviousValue[1])*A;

            return Interpolated[1];
        }
        PreviousValue = Value;
    }
    return PreviousValue[1];
}
D3Q27 FlowSolverLBM::EquilibriumDistribution(double lbDensity, const double weights[3][3][3],
                     dVector3 lbMomentum)
{
   dVector3 macroVel;
   if(lbDensity != 0)
   {
       macroVel = lbMomentum/lbDensity;
   }
   else
   {
       macroVel.set_to_zero();
   }
   double u2 = (macroVel*macroVel);
   D3Q27 locPopulations;
   for(int x = -1; x <= 1; x++)
   for(int y = -1; y <= 1; y++)
   for(int z = -1; z <= 1; z++)
   {
       double cu = x*macroVel[0] + y*macroVel[1] + z*macroVel[2];
       locPopulations(x,y,z) = lbDensity*weights[x+1][y+1][z+1]*(1.0 - 1.5*u2 + cu*(3.0 + 4.5*cu));
   }
   return locPopulations;
};

double FlowSolverLBM::DensityProfile(const double x, const size_t n) const
{
    if (Do_TwoPhase)
    {
        return 0.5 * (LiquidDensity[n]/dRho+VaporDensity[n]/dRho) -
               0.5 * (LiquidDensity[n]/dRho-VaporDensity[n]/dRho) *
               std::tanh(x/InterfaceWidth[n]);
    }
    else
    {
        return 1.0;
    }
}

dVector3 FlowSolverLBM::Velocity(const int i, const int j, const int k, const size_t comp) const
{
    if ((not Obstacle(i,j,k)) and (DensityWetting(i, j, k)({comp}) > DBL_EPSILON))
    {
        if (Do_ThermalComp)
        {
            return (MomentumDensity(i,j,k)({comp}) * 3.0 + ForceDensity(i,j,k)({comp})*dt/2)/DensityWetting(i,j,k)({comp});
        }
        else
        {
            return (MomentumDensity(i,j,k)({comp}) + ForceDensity(i,j,k)({comp})*dt/2)/DensityWetting(i,j,k)({comp});
        }
    }
    else return dVector3{0.0,0.0,0.0};
}

double FlowSolverLBM::Density(const PhaseField& Phase, const int i, const int j, const int k, const size_t comp) const
{
    double Density = 0.0;
    for (auto it = Phase.Fields(i,j,k).cbegin();
             it != Phase.Fields(i,j,k).cend(); it++)
    if (Phase.FieldsStatistics[it->index].State == AggregateStates::Solid)
    {
        Density += it->value * Phase.FieldsStatistics[it->index].Density;
    }

    if (not Obstacle(i,j,k))
    {
        for (auto it = Phase.Fields(i,j,k).cbegin();
                 it != Phase.Fields(i,j,k).cend(); it++)
        if (Phase.FieldsStatistics[it->index].State == AggregateStates::Liquid or
            Phase.FieldsStatistics[it->index].State == AggregateStates::Gas)
        {
           Density += it->value*DensityWetting(i, j, k)({comp});
        }
    }

    return Density;
}
}
