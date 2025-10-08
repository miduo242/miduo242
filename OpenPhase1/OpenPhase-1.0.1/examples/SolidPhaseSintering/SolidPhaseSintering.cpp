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
#include "SolidPhaseSintering.h"
#include "GrandPotential/Density.h"
#include "GrandPotential/Solver.h"
#include "Initializations.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Tools/TimeInfo.h"
#include <iostream>

namespace op = openphase;

#ifdef MPI_PARALLEL
   int MPI_RANK;
   int MPI_SIZE;
#endif

int main(int argc, char *argv[])
{
#ifdef MPI_PARALLEL
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);
    std::cout << "MPI_RANK: " << MPI_RANK << "/" << MPI_SIZE << "\n";
#endif

    std::string InputFileName = op::DefaultInputFileName;
    if (argc > 1) InputFileName = argv[1];

    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID  );  // domain error occurred

    op::Settings Settings(InputFileName);

    op::BoundaryConditions       BC    (Settings, InputFileName);
    op::DoubleObstacle           DO    (Settings, InputFileName);
    op::GrandPotential::Solver   GPS   (Settings, InputFileName);
    op::GrandPotential::Density  omega (Settings, InputFileName);
    op::InterfaceProperties      IP    (Settings, InputFileName);
    op::PhaseField               Phase (Settings, InputFileName);
    op::RunTimeControl           RTC   (Settings, InputFileName);
    op::TimeInfo                 Timer (Settings, "Execution Time Statistics");

    double Temperature = 600; // [K] TODO use temperature class
    int GravityDirection = 1; // TODO move this to GrandPotential::Density::ReadInput

    omega.Set(Temperature,GPS.ChemicalPotential,GravityDirection);

    if (RTC.RestartPossible()) // Restart Simulation
    {
        Phase.Read(BC, false);
        if (RTC.tStart != 0) GPS.Read(BC, Phase, omega);
        else GPS.SetInitial(Phase, omega, BC);
    }
    else // Initialize Simulation
    {
        RTC.tStart = 0;
        double MeanRadius, StdRadius;
        int X0, XN, Y0, YN, Z0, ZN;
        ReadInitializationInuptParameters(Settings, InputFileName, MeanRadius,
                StdRadius, X0, XN, Y0, YN, Z0, ZN);
        auto PhaseIdx = [X0, XN, Settings](double i, double j, double k)
        {
            int idx = 1;
            double step = (XN-X0)/(Settings.Nphases-1);
            for (size_t n = 2; n < Settings.Nphases; n++)
                if (i >= X0 + (n-1)*step) idx++; 
            return idx;
        };

        op::Initializations::Single(Phase,0,BC,Settings);
        op::Initializations::FillRectangularWithSpheres(
                Phase, BC, Settings, PhaseIdx, MeanRadius, StdRadius,
                X0, XN, Y0, YN, Z0, ZN);
        GPS.SetInitial(Phase, omega, BC);
    }

    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        if(RTC.WriteToScreen())
        {
            op::Info::WriteTimeStep(RTC);
            DoSinteringDiagnostics(Settings, Phase, DO, IP, omega, GPS, RTC);
            Phase.FieldsStatistics.WriteTableVolumes(RTC.tStep);
            Timer.PrintWallClockSummary();
        }

        if (RTC.WriteVTK())
        {
            Phase.WriteVTK (RTC.tStep, Settings, false, 8);
            GPS  .WriteVTK (RTC.tStep, Settings, Phase, omega, 8);
        }

        if(RTC.WriteRawData())
        {
            RTC  .Write();
            Phase.Write();
            GPS  .Write();
        }
                                                                                Timer.SetStart();
        IP.Set(Phase);                                                          Timer.SetTimeStamp("IP.Set(Phase)");
        DO.CalculatePhaseFieldIncrements(Phase, IP);                            Timer.SetTimeStamp("DO.CalculatePhaseFieldIncrements");
        GPS.Solve(Phase,omega,BC,IP,RTC.dt);                                    Timer.SetTimeStamp("GPS.Solve");
        Phase.MergeIncrements(BC,RTC.dt);                                       Timer.SetTimeStamp("Phase.MergeIncrements");
    }
    return EXIT_SUCCESS;
}
