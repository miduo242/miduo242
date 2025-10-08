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

#include "BoundaryConditions.h"
#include "AdvectionHR/AdvectionHR.h"
#include "DoubleObstacle.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "FluidDynamics/InteractionSolidFluid.h"
#include "FluidDynamics/InteractionSolidSolid.h"
#include "GrainInfo.h"
#include "Info.h"
#include "Initializations.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "VTK.h"
#include "Velocities.h"
#include <memory>
#include <iomanip>

bool SortAbs (double i, double j) {return std::abs(i) < std::abs(j);}

using namespace openphase;

class LocalLBM : public FlowSolverLBM                                           ///< Modification of FLowSolverLBM which are only used in this example
{
 public:
    LocalLBM(Settings& locSettings, double in_dt): FlowSolverLBM(locSettings, in_dt){};

    double rho_l = 1.8141;                                                      ///<  Initial liquid density
    double rho_v = 0.2;                                                         ///<  Vapor gas density
    double FluidMass;

    double CalculateLiquidVolume(void);
    void EnforceLiquidConcentration(const double& CCurrent,
            const double& CGoal, const int& InitialRadius);

    void InitializeSingle();                                                    ///<  Initializes a single density
    void InitializeSphere(const double Radius,
            const int i0, const int j0, const int k0);                          ///<  Initializes a single sphere
    void InitializeFinalize(PhaseField& Phase, Velocities& Vel,
            BoundaryConditions& BC);
};

void WriteDistanceStatisticToFile(std::ofstream& file, const int tStep,
        PhaseField& Phase, const unsigned int N,
        double& OldDistanceAverage,
        const char sep = ',');

void DetermineSolidParticlePositions(
        std::vector<double>& i0, std::vector<double>& j0, std::vector<double>& k0,
        const Settings& OPSettings, const int NSolids, const double RadiusSolid,
        const double RadiusLiquid);

int main(int argc, char *argv[])
{

//#ifdef DEBUG
    // Detects floating point errors on Linux machines
    //feenableexcept(FE_INEXACT);  // NOTE: THIS DETECTS ROUNDING ERRORS
    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID);    // domain error occurred
    feenableexcept(FE_OVERFLOW);   // the result was too large
    feenableexcept(FE_UNDERFLOW);  // the result was too small
//#endif

    std::string InputFileName = "ProjectInput.opi";

    Info::WriteBlankLine();
    Info::WriteLineInsert("Lattice Boltzmann Capillary Bridge Benchmark");
    Info::WriteStandard("Source", InputFileName);

    std::fstream inp(InputFileName, std::ios::in);
    if (!inp)
    {
        std::stringstream message;
        message << "File \"" << InputFileName << "\" could not be opened";
        throw std::runtime_error(message.str());
    };
    std::stringstream inp_data;
    inp_data << inp.rdbuf();
    inp.close();


    int moduleLocation = UserInterface::FindModuleLocation(inp_data, "SinterLiquidState");

    //const bool UseInterfaceDiffusion = false; //Note does note work with advection yet!
    const bool   UseDoubleObstacle   = UserInterface::ReadParameterB(inp_data, moduleLocation, std::string("bUseDoubleObstacle"));//true;
    const bool   UseVolumeFix        = UserInterface::ReadParameterB(inp_data, moduleLocation, std::string("bUseVolumeFix"));//true;
    const double LiquidConcentration = UserInterface::ReadParameterD(inp_data, moduleLocation, std::string("dLiquidConcentration")); //0.20;  // Concentration of liquid
    const double SolidDensity        = UserInterface::ReadParameterD(inp_data, moduleLocation, std::string("dSolidDensity"));//4.00;   //Lattice units
    const int    InitialRadius       = UserInterface::ReadParameterD(inp_data, moduleLocation, std::string("iInitialRadius"));//10    // Initial radius of solids
    const int    NSolids             = UserInterface::ReadParameterD(inp_data, moduleLocation, std::string("iNSolids"));//3;     // Number of solid particles
    const int    iDistance           = UserInterface::ReadParameterD(inp_data, moduleLocation, std::string("iDistance"));//10;    // Initial distance between particles
    const int    tAdvect             = UserInterface::ReadParameterD(inp_data, moduleLocation, std::string("iTAdvect"));//2000;  // time to start advection

    Settings       OPSettings (DefaultInputFileName);
    RunTimeControl RTC        (OPSettings);

    AdvectionHR         ADHR  (OPSettings);
    BoundaryConditions  BC    (OPSettings);
    DoubleObstacle      DO    (OPSettings);
    InterfaceProperties IP    (OPSettings);
    LocalLBM            LB    (OPSettings, RTC.dt);
    PhaseField          Phase (OPSettings);
    Velocities          Vel   (OPSettings);

    // Set initial conditions
    if(RTC.Restart)
    {
        Phase.Read(BC, RTC.tStart);
        LB   .Read(BC, RTC.tStart);
    }
    else
    {
        RTC.tStart = 0;

        // Radius of liquid and solid spheres
        const double RadiusSolid  = InitialRadius;
        const double RadiusLiquid = (1.0+3.0*LiquidConcentration)*RadiusSolid;

        // Radius of liquid and solid spheres
        Info::WriteLine("-");
        Info::WriteStandard("Liquid radius", RadiusLiquid);
        Info::WriteStandard("Solid radius", RadiusSolid);

        // Set position of particles
        std::vector<double> i0;
        std::vector<double> j0;
        std::vector<double> k0;

        i0.resize(NSolids,0);
        j0.resize(NSolids,0);
        k0.resize(NSolids,0);

        // Initialize vapor phase
        LB.InitializeSingle();
        const int fluid = Initializations::Single(Phase, 0, BC, OPSettings);
        Phase.FieldsStatistics[fluid].State = AggregateStates::Liquid;

        // If the distance is to small, a single liquid bubble is initialed
        DetermineSolidParticlePositions(i0, j0, k0, OPSettings, NSolids,
                RadiusSolid, RadiusSolid + iDistance/2.0);

        // Initialize liquid container
        LB.InitializeSphere(RadiusLiquid,
                OPSettings.Nx/2.0, OPSettings.Ny/2.0, OPSettings.Nz/2.0);

        for(int n = 0; n < NSolids; n++)
        {
            // Initialize solid core
            const unsigned int idx = Initializations::Sphere(Phase, 1,
                    RadiusSolid, i0[n], j0[n], k0[n], BC, OPSettings);

            // Set solid properties
            Phase.FieldsStatistics[idx].Density  = SolidDensity;
            Phase.FieldsStatistics[idx].Mobile = true;
            Phase.FieldsStatistics[idx].State    = AggregateStates::Solid;
        }

        // Finalise initialization
        Phase.SetBoundaryConditions(BC);
        LB.InitializeFinalize(Phase, Vel, BC);
        LB.SetBoundaryConditions(BC);
    }

    // Calculate initial center of mass
    InteractionSolidFluid::CollectGrainsStatistics(Phase, BC, OPSettings);

    // Declare and initialize statistic files
    std::ofstream LogDistance(DefaultTextDir + "stat_distance.csv", std::ofstream::out|std::ofstream::app);
    std::ofstream Log(DefaultTextDir + "stat.csv", std::ofstream::out|std::ofstream::app);
    std::ofstream* LogSolid;
    LogSolid = new std::ofstream [NSolids];
    for(int n = 0; n < NSolids; n++)
    {
        LogSolid[n].open(DefaultTextDir + "stat_solid" + std::to_string(n) + ".csv", std::ofstream::out|std::ofstream::app);
        LogSolid[n] << std::scientific << std::endl;
        LogSolid[n] << std::setprecision(16) << std::endl;
    }

    // Write header
    Log <<"tStep,LiquidConcentration,LiquidVolume,SolidVolume,FluidMass,"
        <<"SolidInterfceEnergy,ObstacleNodes,AvForce,StdDevForce" << std::endl;
    Log << std::scientific << std::endl;
    Log << std::setprecision(16) << std::endl;

    //Calculate reference volume
    std::vector<double> RefVolume;
    InteractionSolidSolid::SetRefVolume(Phase, RefVolume);

    // Start time loop
    const double dt = RTC.dt;
    double OldDistance = 0;
    std::vector<double> Forces;
    for(int tStep = RTC.tStart; tStep <= RTC.nSteps; ++tStep)
    {
        //  Output to VTK file
        if (!(tStep%RTC.tFileWrite))
        {
            Phase.WriteVTK(tStep, OPSettings);
            LB   .WriteVTK(tStep, Phase, OPSettings);
            LB .lbWriteVTK(tStep, Phase, OPSettings);
        }

        //  Output raw data to file
        if (!(tStep%RTC.tRestartWrite))
        {
            Phase.Write(tStep);
            LB.Write(tStep);
        }

        Forces.push_back(Phase.FieldsStatistics[1].Force.abs());

        //  Output to screen
        if(!(tStep%RTC.tScreenWrite))
        {
            Info::WriteTimeStep(tStep, RTC.nSteps);

            double VolumeSolids = 0;
            int    grainNr      = 0;

            // Print out statistics
            dVector3 MomentumSolids = {0.0,0.0,0.0};
            for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
            if (Phase.FieldsStatistics[idx].State == AggregateStates::Solid)
            {
                VolumeSolids += Phase.FieldsStatistics[idx].Volume;

                const double   dx       = OPSettings.dx;
                const double   Mass     = Phase.FieldsStatistics[idx].Volume*
                                          Phase.FieldsStatistics[idx].Density*dx*dx*dx;
                const dVector3 Momentum = Phase.FieldsStatistics[idx].Vcm*Mass;
                MomentumSolids += Momentum;

                Info::WriteSimple("Grain " + std::to_string(grainNr));
                Info::WriteLine("-");
                Info::Write("Density",  Phase.FieldsStatistics[idx].Density);
                Info::Write("Position", Phase.FieldsStatistics[idx].Rcm);
                Info::Write("Velocity", Phase.FieldsStatistics[idx].Vcm);
                Info::Write("Force",    Phase.FieldsStatistics[idx].Force);
                Info::Write("Torque",   Phase.FieldsStatistics[idx].Torque);
                Info::WriteBlankLine();

                Phase.FieldsStatistics[idx].WriteTable(LogSolid[grainNr], tStep);
                grainNr++;
            }

            const double SolidIntEnergy = DO.Energy(Phase,IP);
            const double VolumeLiquid   = LB.CalculateLiquidVolume();
            const double TempCLiquid    = VolumeLiquid/(VolumeLiquid + VolumeSolids);
            const double FluidMass      = LB.CalculateFluidMass()[0];

            const dVector3 FluidMomentum = LB.CalculateFluidMomentum(Phase)[0];

            const size_t ObstacleNodes = LB.CountObstacleNodes();


            // Sort forces in ascending order to achieve better numerical accuracy
            std::sort(Forces.begin(),Forces.end(), SortAbs);
            const size_t N = Forces.size() > 0 ? Forces.size() : 1;

            // Calculate average Force on wall
            long double AverageForce = 0.0;
            for (auto Force: Forces) AverageForce += Force;
            AverageForce /= N;

            // Calculate standard deviation of Force on wall
            long double stdDevForce = 0.0;
            for (auto Force: Forces)
            {
                stdDevForce += std::pow(Force - AverageForce,2);
            }
            stdDevForce /= N;
            stdDevForce = std::sqrt(stdDevForce);
            Forces.clear();

            Info::Write("Solid volume ",          VolumeSolids);
            Info::Write("Solid obstacle nodes ",  ObstacleNodes);
            Info::Write("Solid momentum",         MomentumSolids);
            Info::Write("Solid interface energy", SolidIntEnergy);
            Info::WriteBlankLine();
            Info::Write("Liquid volume",          VolumeLiquid);
            Info::Write("Liquid density",         LB.rho_l * LB.dRho);
            Info::Write("Liquid concentration",   TempCLiquid);
            Info::WriteBlankLine();
            Info::Write("Fluid mass",             FluidMass);
            Info::Write("Fluid momentum",         FluidMomentum);
            Info::WriteBlankLine();
            Info::Write("Average force",          AverageForce);
            Info::Write("stdDev of force",        stdDevForce);
            Info::WriteBlankLine();

            Log << tStep          << ","
                << TempCLiquid    << ","
                << VolumeLiquid   << ","
                << VolumeSolids   << ","
                << FluidMass      << ","
                << SolidIntEnergy << ","
                << ObstacleNodes  << ","
                << AverageForce   << ","
                << stdDevForce    << ","
                << std::endl;

            WriteDistanceStatisticToFile(LogDistance, tStep, Phase, NSolids, OldDistance);
        }

        // Enforce Liquid concentration
        if (tStep > 1*tAdvect/10 and tStep < 8*tAdvect/10)
        {
            double VolumeSolids = 0;
            for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
            if (Phase.FieldsStatistics[idx].State == AggregateStates::Solid)
            {
                VolumeSolids += Phase.FieldsStatistics[idx].Volume;
            }

            const double VolumeLiquid = LB.CalculateLiquidVolume();
            const double TempCLiquid  = VolumeLiquid/(VolumeLiquid + VolumeSolids);

            LB.EnforceLiquidConcentration(TempCLiquid, LiquidConcentration, InitialRadius);
        }

        IP.Set(Phase);
        LB.Solve(Phase, Vel, BC);

        // Let solid move
        if (tStep > tAdvect)
        {
            // Calculate solid/fluid and solid/solid interaction
            //InteractionSolidSolid::Calculate(Phase, LB, dt);
            for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
            {
                Phase.FieldsStatistics[idx].Torque = {0,0,0}; // Torque may destroy simulation
                Phase.FieldsStatistics[idx].Force[1] = 0.0;
                Phase.FieldsStatistics[idx].Force[2] = 0.0;
            }
            InteractionSolidFluid::CalculateSolidVelocities(Phase, Vel, BC, OPSettings, dt);
            LB.EnforceSolidMomentum(Phase,{0.0,0.0,0.0});
            Phase.Advect(ADHR, Vel, BC, LB, RTC.dt, RTC.tStep);

            // Use phase field to restore interface
            if(UseDoubleObstacle and !(tStep%10))
            {
                DO.CalculatePhaseFieldIncrements(Phase, IP);
                Phase.MergeIncrements(BC, dt, true);
            }

            // Applies volume fix
            size_t iterations = 0;
            while (UseVolumeFix and InteractionSolidSolid::VolumeError(Phase, RefVolume, 1) > 1.0e-4)
            {
                IP.Set(Phase);
                InteractionSolidSolid::PreserveVolume(Phase, RefVolume, 1, dt);
                Phase.MergeIncrements(BC, dt);
                ++iterations;
                if (iterations > 100)
                {
                    std::cout << "Max Iterations in VolumeControl!" << std::endl;
                    break;
                }
            }
        }
    }// End of time loop
    return 0;
}

void DetermineSolidParticlePositions(
        std::vector<double>& i0, std::vector<double>& j0, std::vector<double>& k0,
        const Settings& OPSettings, const int NSolids, const double RadiusSolid,
        const double RadiusLiquid)
{
    const double offset = std::max(RadiusLiquid, RadiusSolid+OPSettings.iWidth/2.0);
    if (NSolids == 1)
    {
        i0[0] = OPSettings.Nx/2.0;
        j0[0] = OPSettings.Ny/2.0;
        k0[0] = OPSettings.Nz/2.0;
    }
    else if (NSolids == 2)
    {
        i0[0] = OPSettings.Nx/2.0 - offset;
        j0[0] = OPSettings.Ny/2.0;
        k0[0] = OPSettings.Nz/2.0;

        i0[1] = OPSettings.Nx/2.0 + offset;
        j0[1] = OPSettings.Ny/2.0;
        k0[1] = OPSettings.Nz/2.0;
    }
    else if (NSolids == 3)
    {
        const double h  = std::sqrt(3.0) * offset;
        const double O  = tan(Pi/180.0*30) * offset;       // distance from base line to center of equilateral triangle

        i0[0] = OPSettings.Nx/2.0 - offset;
        j0[0] = OPSettings.Ny/2.0 - O;
        k0[0] = OPSettings.Nz/2.0;

        i0[1] = OPSettings.Nx/2.0 + offset;
        j0[1] = OPSettings.Ny/2.0 - O;
        k0[1] = OPSettings.Nz/2.0;

        i0[2] = OPSettings.Nx/2.0;
        j0[2] = OPSettings.Ny/2.0 - O + h;
        k0[2] = OPSettings.Nz/2.0;
    }
    else if (NSolids == 4)
    {
        const double O1 = tan(Pi/180.0*30) * offset;       // distance from base line to center of equilateral triangle
        const double h1 = std::sqrt(3.0)/2.0 * 2.0*offset; // height of equilateral triangle
        const double h2 = std::sqrt(6.0)/3.0 * 2.0*offset; // median of tetrahedron

        i0[0] = OPSettings.Nx/2.0 - offset;
        j0[0] = OPSettings.Ny/2.0 - O1;
        k0[0] = OPSettings.Nz/2.0 - h2/3.0;

        i0[1] = OPSettings.Nx/2.0 + offset;
        j0[1] = OPSettings.Ny/2.0 - O1;
        k0[1] = OPSettings.Nz/2.0 - h2/3.0;

        i0[2] = OPSettings.Nx/2.0;
        j0[2] = OPSettings.Ny/2.0 - O1 + h1;
        k0[2] = OPSettings.Nz/2.0 - h2/3.0;

        i0[3] = OPSettings.Nx/2.0;
        j0[3] = OPSettings.Ny/2.0;
        k0[3] = OPSettings.Nz/2.0 + h2*2.0/3.0;
    }
    else
    {
        //TODO Random distribution
    }
}

void WriteDistanceStatisticToFile(std::ofstream& file, const int tStep,
        PhaseField& Phase, const unsigned int N, double& OldDistanceAverage,
        const char sep)
{
    file << std::scientific << std::endl;
    file << std::setprecision(16) << std::endl;
    if (tStep == 0)
    {
        file << "tStep" << sep;
        for(size_t n =     0; n < Phase.FieldsStatistics.size(); ++n)
        for(size_t m = n + 1; m < Phase.FieldsStatistics.size(); ++m)
        {
            Grain grainA = Phase.FieldsStatistics[n];
            Grain grainB = Phase.FieldsStatistics[m];

            if (grainA.State == AggregateStates::Solid and grainB.State == AggregateStates::Solid)
            {
                file << "R" << n << "_" << m << sep;
            }
        }
        file << "RAverage" << sep << "dRAverage" << std::endl;
    }
    else
    {
        // Calculate distances
        double DistanceAverage = 0.0;
        file  << tStep << sep;
        for(size_t n =     0; n < Phase.FieldsStatistics.size(); ++n)
        for(size_t m = n + 1; m < Phase.FieldsStatistics.size(); ++m)
        {
            Grain grainA = Phase.FieldsStatistics[n];
            Grain grainB = Phase.FieldsStatistics[m];

            if (grainA.State == AggregateStates::Solid and grainB.State == AggregateStates::Solid)
            {
                const double buffer = (grainA.Rcm - grainB.Rcm).abs();
                DistanceAverage += buffer;
                file << buffer << sep;
            }
        }
        DistanceAverage /= (N*(N+1)/2-N);
        file << DistanceAverage << sep;

        if (tStep > 1) file << 0 << std::endl;
        else file << DistanceAverage - OldDistanceAverage << std::endl;

        OldDistanceAverage = DistanceAverage;
    }
}

void LocalLBM::EnforceLiquidConcentration(const double& CCurrent, const double& CGoal, const int& InitialRadius)
{
    if (std::abs(CCurrent-CGoal) > 0.0001)
    {
        const double DeltaDensity = -0.1/InitialRadius*(CGoal-CCurrent)*rho_l;

        // Fix density
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
        if (not Obstacle(i,j,k))
        {
            const double DistLiquid = abs(DensityWetting(i,j,k)({0}) - rho_l);
            const double DistVapor  = abs(DensityWetting(i,j,k)({0}) - rho_v);

            if (DistLiquid < DistVapor)
            {
                lbPopulations(i,j,k)({0}) = lbPopulations(i,j,k)({0}) -
                    FlowSolverLBM::EquilibriumDistribution(DensityWetting(i,j,k)({0})/dRho + DeltaDensity, lbWeights,{0.,0.,0.}) +
                    FlowSolverLBM::EquilibriumDistribution(DensityWetting(i,j,k)({0})/dRho, lbWeights, {0.,0.,0.});
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        // Update entire fluid mass in order to ensure mass conservation
        FluidMass = CalculateFluidMass()[0];
    }
}

double LocalLBM::CalculateLiquidVolume(void)
{
    int VolumeLiquid   = 0;
    int VolumeVapor    = 0;
    //int VolumeObstacle = 0;

    //double MassLiquid = 0;
    //double MassVapor  = 0;

    double New_rho_l = 0;
    double New_rho_v = 1;

    STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0)
    if (not Obstacle(i,j,k))
    {
        const double DistLiquid = abs(DensityWetting(i,j,k)({0})/dRho - rho_l);
        const double DistVapor  = abs(DensityWetting(i,j,k)({0})/dRho - rho_v);

        if (DistLiquid < DistVapor)
        {
            ++VolumeLiquid;
            if (DensityWetting(i,j,k)({0}) > New_rho_l)
                New_rho_l = DensityWetting(i,j,k)({0})/dRho;
            //MassLiquid += DensityWetting(i,j,k)({0})/dRho;
        }
        else
        {
            ++VolumeVapor;
            if (DensityWetting(i,j,k)({0}) < New_rho_v)
                New_rho_v = DensityWetting(i,j,k)({0})/dRho;
            //MassVapor += DensityWetting(i,j,k)({0})/dRho;
        }
    }
    //else VolumeObstacle++;
    STORAGE_LOOP_END

    //Update liquid and vapor density
    //rho_l = MassLiquid/VolumeLiquid;
    //rho_v = MassVapor/VolumeVapor;

    rho_l = New_rho_l;
    rho_v = New_rho_v;

    return VolumeLiquid;
}

void LocalLBM::InitializeSingle()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    {
        Obstacle(i,j,k) = 0;
        DensityWetting(i,j,k)({0}) = rho_v*dRho;
        MomentumDensity(i,j,k)({0}) = {0.,0.,0.};
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LocalLBM::InitializeSphere(const double Radius,
        const int i0, const int j0, const int k0)
{
    const double D = 3;
    const int Bcells = DensityWetting.Bcells();
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,Bcells,)
    {
        const int ii = i - i0;
        const int jj = j - j0;
        const int kk = k - k0;
        const int rr = sqrt(ii*ii + jj*jj + kk*kk);

        DensityWetting(i,j,k)({0}) = std::max((0.5*(rho_l+rho_v) -0.5*(rho_l-rho_v)*tanh(2.*(rr-Radius)/(D)))*dRho, DensityWetting(i,j,k)({0}));
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LocalLBM::InitializeFinalize(PhaseField& Phase, Velocities& Vel,
        BoundaryConditions& BC)
{
    // Initialize lattice Boltzmann populations for liquid phase
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,DensityWetting.Bcells(),)
    {
        lbPopulations(i,j,k)({0}) = FlowSolverLBM::EquilibriumDistribution(DensityWetting(i,j,k)({0})/dRho, lbWeights, MomentumDensity(i,j,k)({0})/dm);
        ForceDensity(i,j,k)({0}).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    DetectObstacles(Phase);
    SetObstacleNodes(Phase,Vel);

}
