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
 */

#include "AdvectionHR/AdvectionHR.h"
#include "BoundaryConditions.h"
#include "DoubleObstacle.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "FluidDynamics/InteractionSolidFluid.h"
#include "FluidDynamics/InteractionSolidSolid.h"
#include "GrainInfo.h"
#include "Info.h"
#include "Initializations.h"
#include "InterfaceDiffusion.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "VTK.h"
#include "Velocities.h"

using namespace openphase;

class LocalLBM : public FlowSolverLBM                                           ///< Modification of FLowSolverLBM which are only used in this example
{
 public:
    LocalLBM(const Settings& locSettings, double in_dt): FlowSolverLBM(locSettings, in_dt){};

    void InitializeSingle();                                                    ///<  Initializes a single density
    void InitializeSphere(const double Radius, const int i0, const int j0,
            const int k0, const double* rho = nullptr);  ///<  Initializes a single sphere
};

int main(int argc, char *argv[])
{

    // Detects floating point errors on Linux machines
    //feenableexcept(FE_INEXACT);  // NOTE: THIS DETECTS ROUNDING ERRORS
    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID);    // domain error occurred
    feenableexcept(FE_OVERFLOW);   // the result was too large
    feenableexcept(FE_UNDERFLOW);  // the result was too small

    // Read benchmark specific input parameters
    Info::WriteBlankLine();
    Info::WriteLineInsert("LBGravity");
    Info::WriteStandard("Source", DefaultInputFileName);
    std::fstream inpF(DefaultInputFileName, std::ios::in);
    std::stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();
    int moduleLocation = UserInterface::FindModuleLocation(inp, "LBGravity");
    const bool   UseDoubleObstacle = UserInterface::ReadParameterB(inp, moduleLocation, "bUseDoubleObstacle");
    const bool   UseVolumeFix      = UserInterface::ReadParameterB(inp, moduleLocation, "bUseVolumeFix");
    const double SolidDensity      = UserInterface::ReadParameterD(inp, moduleLocation, "dSolidDensity");     //Lattice units
    const double VolumeFixAccuracy = UserInterface::ReadParameterD(inp, moduleLocation, "dVolumeFixAccuracy");
    const int    RadiusLiquid      = UserInterface::ReadParameterI(inp, moduleLocation, "iRadiusLiquid");     // Initial radius of solids
    const int    RadiusSolid       = UserInterface::ReadParameterI(inp, moduleLocation, "iRadiusSolid");      // Initial radius of solids

    Settings OPSettings;
    OPSettings.ReadInput();

    RunTimeControl       RTC   (OPSettings);
    AdvectionHR          ADHR  (OPSettings);
    BoundaryConditions   BC    (OPSettings);
    DoubleObstacle       DO    (OPSettings);
    InterfaceProperties  IP    (OPSettings);
    LocalLBM             LB    (OPSettings, RTC.dt);
    PhaseField           Phase (OPSettings);
    Velocities           Vel   (OPSettings);

    // Set initial conditions
    if(RTC.Restart)
    {
        Phase.Read(BC, RTC.tStart);
        LB   .Read(BC,RTC.tStart);
        Phase.SetBoundaryConditions(BC);
        LB   .SetBoundaryConditions(BC);
    }
    else
    {
        const int Nx = OPSettings.Nx;
        const int Ny = OPSettings.Ny;
        const int Nz = OPSettings.Nz;

        // Initialize Fluid phase (lattice Boltzmann representation)
        LB.InitializeSingle();

        // Initialize wall phase
        const size_t wall = Initializations::Single(Phase, 1, BC, OPSettings);
        Phase.FieldsStatistics[wall].Density  = SolidDensity;
        Phase.FieldsStatistics[wall].Mobile = false;

        // Initialize Fluid phase (Phase-Field representation)
        Initializations::Rectangular(Phase, 0, Nx-16, Ny-16, Nz, Nx/2, Ny/2, Nz/2, BC, OPSettings);

        // Initialize liquid sphere
        LB.InitializeSphere(RadiusLiquid, Nx/4, Ny/2, Nz/2, &LB.LiquidDensity[0]);

        // Initialize solid sphere
        const size_t sphere = Initializations::Sphere(Phase, 1, RadiusSolid, 3*Nx/4, Ny/2, Nz/2, BC, OPSettings);
        Phase.FieldsStatistics[sphere].Density  = SolidDensity;
        Phase.FieldsStatistics[sphere].Mobile = true;
        // Wet solid sphere
        LB.InitializeSphere(RadiusSolid, 3*Nx/4, Ny/2, Nz/2, &LB.Wetting[0][1]);

        // Finalize initialization
        Phase.SetBoundaryConditions(BC);
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,LB.DensityWetting,LB.DensityWetting.Bcells(),)
        {
            LB.lbPopulations(i,j,k)({0}) = FlowSolverLBM::EquilibriumDistribution(LB.DensityWetting(i,j,k)({0})/LB.dRho, LB.lbWeights);
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        LB.DetectObstacles(Phase);
        LB.SetObstacleNodes(Phase,Vel);
        LB.EnforceMassConservation();
        LB.SetBoundaryConditions(BC);
    }

    // Calculate initial center of mass
    InteractionSolidFluid::CollectGrainsStatistics(Phase, BC, OPSettings);

    // Start time loop
    std::vector<double> Forces;
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        if (RTC.WriteVTK())
        {
            Phase.WriteVTK(RTC.tStep, OPSettings);
            LB   .WriteVTK(RTC.tStep, Phase, OPSettings);
            LB .lbWriteVTK(RTC.tStep, Phase, OPSettings);
        }

        if (RTC.WriteRawData())
        {
            Phase.Write(RTC.tStep);
            LB   .Write(RTC.tStep);
        }

        if(RTC.WriteToScreen())
        {
            // Print out statistics to screen
            Info::WriteTimeStep(RTC.tStep, RTC.nSteps);
            double VolumeSolids = 0;
            int    grainNr      = 0;
            dVector3 MomentumSolids = {0.0,0.0,0.0};
            for (size_t idx = 0; idx < Phase.FieldsStatistics.size();idx++)
            if (Phase.FieldsStatistics[idx].State == AggregateStates::Solid and
                Phase.FieldsStatistics[idx].Mobile)
            {
                VolumeSolids += Phase.FieldsStatistics[idx].Volume;

                const double   dx       = OPSettings.dx;
                const double   Mass     = Phase.FieldsStatistics[idx].Volume*
                                          Phase.FieldsStatistics[idx].Density*dx*dx*dx;
                const dVector3 Momentum = Phase.FieldsStatistics[idx].Vcm*Mass;
                MomentumSolids += Momentum;

                Info::WriteSimple("Grain " + std::to_string(grainNr));
                Info::WriteLine("-");
                Info::Write("Position",     Phase.FieldsStatistics[idx].Rcm);
                Info::Write("Acceleration", Phase.FieldsStatistics[idx].Acm);
                Info::Write("Velocity",     Phase.FieldsStatistics[idx].Vcm);
                Info::Write("Force",        Phase.FieldsStatistics[idx].Force);
                Info::Write("Torque",       Phase.FieldsStatistics[idx].Torque);
                Info::WriteBlankLine();
            }

            const dVector3 FluidMomentum = LB.CalculateFluidMomentum(Phase)[0];
            const double   FluidMass     = LB.CalculateFluidMass()[0];

            Info::Write("Solid volume ",  VolumeSolids);
            Info::Write("Solid momentum", MomentumSolids);
            Info::WriteBlankLine();
            Info::Write("Fluid mass",     FluidMass);
            Info::Write("Fluid momentum", FluidMomentum);
            Info::WriteBlankLine();
        }

        LB.Solve(Phase, Vel, BC);

        const int    Order    = 4;     // solid-solid interaction potential order
        const int    Range    = 3;     // solid-solid interaction range
        const double Strength = 1.e+2; // solid-solid interaction strength
        const double Elastic  = 1.e-1; // 1.0 means fully elastic impact
        InteractionSolidSolid::Calculate(Phase, LB, RTC.dt, Order, Range, Strength, Elastic);
        for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
        {
            Phase.FieldsStatistics[idx].Torque = {0,0,0}; // NOTE Torque destroys simulation
            //Phase.FieldsStatistics[idx].Force[0] = 0.0;
            //Phase.FieldsStatistics[idx].Force[2] = 0.0;
        }
        InteractionSolidFluid::CalculateSolidVelocities(Phase, Vel, BC, OPSettings, RTC.dt);
        ADHR.AdvectPhaseField(Phase, Vel, BC, OPSettings.dx, RTC.dt, RTC.tStep, false);
        LB.DetectObstaclesAdvection(Phase, Vel, BC);
    }
    return 0;
}

void LocalLBM::InitializeSingle()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    {
        Obstacle       (i,j,k) = 0;
        DensityWetting (i,j,k)({0}) = VaporDensity[0];
        MomentumDensity(i,j,k)({0}) = {0.,0.,0.};
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LocalLBM::InitializeSphere(const double Radius, const int i0, const int j0,
        const int k0, const double* rho)
{
    const double loc_rho = (rho == nullptr) ? LiquidDensity[0] : *rho;
    const double D = 3;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    {
        const int ii = i - i0;
        const int jj = j - j0;
        const int kk = k - k0;
        const double rr = sqrt(ii*ii + jj*jj + kk*kk);

        DensityWetting(i,j,k)({0}) = std::max((0.5*(loc_rho+VaporDensity[0]/dRho) -0.5*(loc_rho-VaporDensity[0]/dRho)*tanh(2.*(rr-Radius)/(D)))*dRho, DensityWetting(i,j,k)({0}));
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
