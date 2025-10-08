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

#include "Base/Macros.h"
#include "AdvectionHR/AdvectionHR.h"
#include "BoundaryConditions.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "FluidDynamics/InteractionSolidFluid.h"
#include "Info.h"
#include "Initializations.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "VTK.h"
#include "Velocities.h"

#include "LocalLBM.h"

#include <complex>

using namespace openphase;

int main(int argc, char **argv)
{
    std::string InputFileName = DefaultInputFileName;

#ifdef MPI_PARALLEL
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);
    {
#endif

    // Detects floating point errors on Linux machines
    //feenableexcept(FE_INEXACT);  // NOTE: THIS DETECTS ROUNDING ERRORS
    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID);    // domain error occurred
    //feenableexcept(FE_OVERFLOW);  // the result was too large
    //feenableexcept(FE_UNDERFLOW); // the result was too small

    Settings           OPSettings;
    OPSettings.ReadInput(InputFileName);

    AdvectionHR        ADHR  (OPSettings);
    BoundaryConditions BC    (OPSettings);
    PhaseField         Phase (OPSettings);
    RunTimeControl     RTC   (OPSettings);
    Velocities         Vel   (OPSettings);
    LocalLBM           LB    (OPSettings,RTC.dt);

    Info::WriteBlankLine();
    Info::WriteLineInsert("Lattice Boltzmann Capillary Bridge Benchmark");
    Info::WriteStandard("Source", DefaultInputFileName);

    std::fstream inpF(InputFileName, std::ios::in);
    if (!inpF)
    {
        std::stringstream message;
        message << "File \"" << InputFileName << "\" could not be opened";
        throw std::runtime_error(message.str());
    };
    std::stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();

    int moduleLocation = UserInterface::FindModuleLocation(inp, "LBCapillaryBridge");

          double SolidDensity     = UserInterface::ReadParameterD(inp, moduleLocation, std::string("dSolidDensity"));
          int    tAdvect          = UserInterface::ReadParameterI(inp, moduleLocation, std::string("iTAdvect"));
    const bool   AdvectSolids     = UserInterface::ReadParameterB(inp, moduleLocation, std::string("bAdvect"));
    const bool   EnforceVolume    = UserInterface::ReadParameterB(inp, moduleLocation, std::string("bEnforceV"));
    const bool   EnforceDiameter  = UserInterface::ReadParameterB(inp, moduleLocation, std::string("bEnforceD"));
    const bool   EnforceCurvature = UserInterface::ReadParameterB(inp, moduleLocation, std::string("bEnforceC"));


    // Define short hand for system size
    const size_t TotalNx = OPSettings.TotalNx;
    const size_t Nx = OPSettings.Nx;
    const size_t Ny = OPSettings.Ny;
    const size_t Nz = OPSettings.Nz;

    // Determine the number of time-steps the system needs to reach the equilibrium
    const double dx = OPSettings.dx;
    const double dt = RTC.dt;

    // Determine the liquid radius in lattice units
    const double lbR0 = TotalNx/3;
    const double R0   = TotalNx/3*dx;

    // Correct extensive quantities in 2D
    const bool ThreeDim = Nx > 1 and Ny > 1 and Nz > 1;

    // Define phase indices
    size_t FluidPhaseIdx = 0;
    size_t WallPhaseIdx  = 1;

    // Define phase-field indices
    size_t FluidIdx = 0;
    size_t WallIdx  = 1;

    if(RTC.Restart)
    {
        // Read restart input
        LB.Read    (BC, RTC.tStart);
        Phase.Read (BC, RTC.tStart);
        Vel.Read   (BC, RTC.tStart);
    }
    else
    {
        // Set start time step
        RTC.tStart = 0;

        // Determine dimension of plate
        const double Lx = 4.00*lbR0 + 1.0*OPSettings.iWidth;
        const double Ly = 1.00*lbR0 + 1.0*OPSettings.iWidth;
        const double Lz = 4.00*lbR0 + 1.0*OPSettings.iWidth; // Thickness of the walls
        const double j0 = 1.50*lbR0;

        // Initialize phase fields
        FluidIdx = Initializations::Single(Phase, FluidPhaseIdx, BC, OPSettings);
        WallIdx  = Initializations::Rectangular(Phase, WallPhaseIdx, Lx , Ly, Lz, 0, j0, 0, BC, OPSettings);

        // Set phase states
        Phase.FieldsStatistics[FluidIdx].State  = AggregateStates::Liquid;
        Phase.FieldsStatistics[ WallIdx].Mobile = true;
        Phase.FieldsStatistics[ WallIdx].State  = AggregateStates::Solid;

        LB.SetInitialDF(Phase, OPSettings, Vel, BC, lbR0, 0, j0);
        Phase.CalculateVolumes();
        Phase.SetBoundaryConditions(BC);

        // Calculate initial center of mass
        InteractionSolidFluid::CollectGrainsStatistics(Phase, BC, OPSettings);
    }

    // Open log files
    std::string FileName = OPSettings.TextDir + "TimeLog.csv";
    std::fstream log(FileName, std::ios::trunc | std::ios::out);
    log << std::scientific << std::setprecision(16);

    double volumeW = Phase.FieldsStatistics[WallIdx].Volume*dx*dx*dx;
    if (ThreeDim) volumeW *= 4; else volumeW *= 2*lbR0;
    const double SolidMass = SolidDensity*volumeW;
    double tc = std::sqrt(SolidMass/LB.SurfaceTension[0]);
    double ErrorGap   = DBL_MAX;
    double ErrorForce = DBL_MAX;

    // Start time loop
    std::vector<double> WallForces;
    Info::WriteSimple("Entering the Time Loop!!!");
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        // NOTE: In the planar case
        volumeW = Phase.FieldsStatistics[WallIdx].Volume*dx*dx*dx;
        if (ThreeDim) volumeW *= 4; else volumeW *= 2*lbR0;

        // Set density of solid wall to ensure mass conservation
        SolidDensity = (volumeW!=0)? SolidMass/volumeW : 0.0;
        Phase.FieldsStatistics[WallIdx].Density = SolidDensity;

        if (RTC.WriteRawData())
        {
            LB   .Write(RTC.tStep);
            Phase.Write(RTC.tStep);
            Vel  .Write(RTC.tStep);
        }

        if(RTC.WriteVTK())
        {
            LB.WriteVTK(RTC.tStep, Phase, OPSettings);
        }

        // Calculate mid radius of the capillary bridge
        const Grain Wall  = Phase.FieldsStatistics[WallIdx];
        WallForces.push_back(Wall.Force[1]);

        //  Output to screen and do some diagnostics
        if(RTC.WriteToScreen())
        {
            const double t            = RTC.tStep*dt;
            const double t0           = tAdvect*dt;
            const double Time         = t-t0;
            const double TimeRescaled = Time/tc;

            // Calculate interface tensions (via density profile)
            const double sigma_lv = LB.CalculateSigmaX(0, 0);
            const double sigma_sv = LB.CalculateSigmaY(Wall.Rcm[1], Ny-1, 0);
            const double sigma_sl = LB.CalculateSigmaY(0, Wall.Rcm[1], 0);
            const double kappa    = LB.AverageCurvature(0, 0.25*lbR0, 0, Nz);
            const double theta    = LB.CalculateContactAngle(sigma_sl, sigma_sv, sigma_lv);

            // TODO Calculate pressures
            // const double pL       = LB.Pressure(   0,    0,    0);
            // const double pV       = LB.Pressure(Nx-1, Ny-1, Nz-1);
            // const double dp       = pL - pV;

            // TODO Calculate interface tensions (via curvature and pressure)
            //const double sigma_laplace = (std::abs(kappa) > DBL_EPSILON) ? dp/kappa : 0;

            // Calculate volumes
            std::array<double,2> lbVolumes = LB.CalculateFluidVolumes();
            const double volumeL  = lbVolumes[0]*dx*dx*dx;
            const double volumeV  = lbVolumes[1]*dx*dx*dx;

            // Calculate masses
            const double massF = LB.CalculateFluidMass()[0];

            // Calculate densities
#ifdef MPI_PARALLEL
            double locDensityVapor  = (MPI_RANK==MPI_SIZE-1) ? LB.Density(Phase, OPSettings.Nx, Ny, (Nz == 1)?0:1, 0) : 0;
            double locDensityLiquid = (MPI_RANK==0) ? LB.Density(Phase,0,0,0,0) : 0;
            double DensityVapor  = 0;
            double DensityLiquid = 0;
            MPI_Allreduce(&locDensityVapor , &DensityVapor , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&locDensityLiquid, &DensityLiquid, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
            double DensityVapor  = LB.Density(Phase, Nx, Ny, (Nz == 1)?0:1, 0);
            double DensityLiquid = LB.Density(Phase,0,0,0,0);
#endif

            // Sort forces in ascending order to achieve better numerical accuracy
            std::sort(WallForces.begin(), WallForces.end(), [](double i, double j) {return std::abs(i) < std::abs(j);});
            const size_t N = WallForces.size() > 0 ? WallForces.size() : 1;

            // Calculate average Force on wall
            long double AverageWallForce = 0.0;
            for (auto WallForce: WallForces) AverageWallForce += WallForce;
            AverageWallForce /= N;

            // Calculate standard deviation of Force on wall
            long double stdDevWallForce = 0.0;
            for (auto WallForce: WallForces)
            {
                stdDevWallForce += std::pow(WallForce - AverageWallForce,2);
            }
            stdDevWallForce /= N;
            stdDevWallForce = std::sqrt(stdDevWallForce);
            WallForces.clear();

            const double RMid       = LB.CalculateRadius(OPSettings, 0, 0);
            const double WallForce  = Wall.Force[1];
            const double F_ex       = ThreeDim ? -0.25*Pi*RMid*sigma_lv : -sigma_lv*dx;
            if (TimeRescaled < 0)
            if (std::abs(F_ex) > 100*DBL_EPSILON)
            {
                 ErrorForce = std::abs(AverageWallForce-F_ex)/std::abs(F_ex);
            }

            // Calculate wall distance between the walls
            //const double WallSpeed         = Wall.Vcm[1];
            const double WallGap           = (2.0*Wall.Rcm[1]*dx - R0);
            const double WallGapRescaled   = (2.0*Wall.Rcm[1]*dx - R0)/R0;

            double WallGapRescaledEx = 2.0;
            if (TimeRescaled>0)
            {
                if (ThreeDim) WallGapRescaledEx = 2.0 -Pi*TimeRescaled*TimeRescaled;
                else          WallGapRescaledEx = 2.0 - 2*TimeRescaled*TimeRescaled;
            }

            if (TimeRescaled>0)
            {
                ErrorGap  = std::abs(WallGapRescaled-WallGapRescaledEx);
                ErrorGap /= WallGapRescaledEx;
            }

            std::array<std::stringstream,2> line;
            line[1] << std::scientific << std::setprecision(16);
            Info::WriteTimeStep(RTC.tStep, RTC.nSteps);
            Info::WriteWithLog(line, RTC.tStep, "Time [s]", Time);
            Info::WriteWithLog(line, RTC.tStep, "Time Rescaled", TimeRescaled);
            Info::Write("");
            Info::WriteWithLog(line, RTC.tStep, "Volume of Vapor[m^3]", volumeV);
            Info::WriteWithLog(line, RTC.tStep, "Volume of Liquid [m^3]", volumeL);
            Info::WriteWithLog(line, RTC.tStep, "Volume of Wall [m^3]", volumeW);
            Info::Write("");
            Info::WriteWithLog(line, RTC.tStep, "Mass of Fluid [kg]", massF);
            Info::WriteWithLog(line, RTC.tStep, "Mass of Wall [kg]", SolidMass);
            Info::Write("");
            Info::WriteWithLog(line, RTC.tStep, "Desity of Vapor [kg/m^3]", DensityVapor);
            Info::WriteWithLog(line, RTC.tStep, "Desity of Liquid [kg/m^3]", DensityLiquid);
            Info::WriteWithLog(line, RTC.tStep, "Desity of Wall [kg/m^3]", SolidDensity);
            Info::Write("");
            Info::WriteWithLog(line, RTC.tStep, "Surface Tension (integral) [kg/s^2]", sigma_lv);
            //Info::WriteWigthLog("Surface tension (dp/kappa) [kg/s^2]", sigma_laplace);
            Info::WriteWithLog(line, RTC.tStep, "Expected Surface Tension [kg/s^2]", LB.SurfaceTension[0]);
            Info::Write("");
            Info::WriteWithLog(line, RTC.tStep, "Contact Angle [deg]", theta);
            Info::WriteWithLog(line, RTC.tStep, "Wetting Parameter [%]", LB.Wetting[FluidPhaseIdx][WallPhaseIdx]);
            Info::Write("");
            Info::WriteWithLog(line, RTC.tStep, "Capillary Diameter [m]", 2*RMid);
            Info::WriteWithLog(line, RTC.tStep, "Capillary Diameter Rescaled", 2*RMid/R0);
            Info::WriteWithLog(line, RTC.tStep, "Curvature [1/m]", kappa);
            Info::WriteWithLog(line, RTC.tStep, "Curvature Rescaled", kappa*R0);
            Info::Write("");
            Info::WriteWithLog(line, RTC.tStep, "Force on Wall [kg*m/s^2]", WallForce);
            Info::WriteWithLog(line, RTC.tStep, "Average Force on Wall [kg*m/s^2]", AverageWallForce);
            Info::WriteWithLog(line, RTC.tStep, "Expected Force on Wall [kg*m/s^2]", F_ex);
            Info::Write("-");
            Info::WriteWithLog(line, RTC.tStep, "Relative Force Error", ErrorForce);
            Info::Write("");
            Info::WriteWithLog(line, RTC.tStep, "Distance Between Walls [m]", WallGap);
            Info::WriteWithLog(line, RTC.tStep, "Expected Distance Between Walls [m]", (Time > 0) ? LB.SurfaceTension[0]*dx/(SolidMass/lbR0)*Time*Time : WallGap);
            Info::Write("");
            Info::WriteWithLog(line, RTC.tStep, "Rescaled Distance Between Walls", WallGapRescaled);
            Info::WriteWithLog(line, RTC.tStep, "Expected Rescaled Distance between Walls", WallGapRescaledEx);
            Info::Write("-");
            Info::WriteWithLog(line, RTC.tStep, "Relative Distance Error", ErrorGap);
            Info::Write("");
            Info::WriteLineToLogfile(log, line, RTC.tStep);

            if (RTC.tStep < tAdvect) tc = std::sqrt(SolidMass/sigma_lv);
            if (WallGapRescaled < 1.25) break;
        }

        // Solve lattice Boltzmann
        LB.Solve(Phase, Vel, BC);

        // Let walls move
        if (AdvectSolids and RTC.tStep > tAdvect)
        {
            // Allow only uni-axial motion
            Phase.FieldsStatistics[WallIdx].Force.setX(0);
            Phase.FieldsStatistics[WallIdx].Force.setZ(0);
            Phase.FieldsStatistics[WallIdx].Torque = {0,0,0};

            // Calculate wall velocities
            InteractionSolidFluid::CalculateSolidVelocities(Phase, Vel, BC, OPSettings, RTC.dt);

            // Move solid walls
            Phase.Advect(ADHR, Vel, BC, LB, RTC.dt, RTC.tStep);
        }
        else if (RTC.tStep > R0*R0/dx/dx)
        {

            if (EnforceDiameter)
            {
                const double DMid = 2*LB.CalculateRadius(OPSettings, 0, 0);
                LB.EnforceLiquidDiameter(BC, DMid, 2*R0);
            }

            if (EnforceVolume)
            {
                std::array<double,2> volume = LB.CalculateFluidVolumes();
                const double lbVCurrent = volume[0];
                const double lbVGoal    = ThreeDim ? lbR0*lbR0*lbR0 : lbR0*lbR0;
                LB.EnforceVolume(BC, lbVCurrent, lbVGoal);
            }

            if (EnforceCurvature)
            {
                const double kappa  = LB.AverageCurvature(0, 0.25*lbR0, 0, Nz);
                LB.EnforceCurvature(LB.Wetting[FluidPhaseIdx][WallPhaseIdx], kappa, ThreeDim ? 1.0/R0 : 0.0);
            }
        }
    }// End of time loop

    std::ofstream resultsSim ("Results.sim");
    if (resultsSim.is_open())
    {
        resultsSim << 2 << "\n";
        resultsSim << "ErrorWallForce "    << ErrorForce << " " << ErrorForce << "\n";
        resultsSim << "ErrorWallDistance " << ErrorGap   << " " << ErrorGap;
        resultsSim.flush();
        resultsSim.close();
    }

#ifdef MPI_PARALLEL
    }
    MPI_Finalize();
#endif

    return 0;
}

