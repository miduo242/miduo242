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
#include "DoubleObstacle.h"
#include "Initializations.h"
#include "InterfaceDiffusion.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "Tools/TimeInfo.h"

namespace op = openphase;

int main()
{
    std::string InputFileName = op::DefaultInputFileName;

    // Detects floating point errors on Linux machines
    //feenableexcept(FE_INEXACT);  // NOTE: THIS DETECTS ROUNDING ERRORS
    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID);    // domain error occurred
    feenableexcept(FE_OVERFLOW);   // the result was too large
    feenableexcept(FE_UNDERFLOW);  // the result was too small

    op::BoundaryConditions  BC;
    op::DoubleObstacle      DO;
    op::InterfaceDiffusion  ID;
    op::InterfaceProperties IP;
    op::PhaseField          Phase;
    op::RunTimeControl      RTC;
    op::Settings            OPSettings;
    op::TimeInfo            Timer;

    OPSettings.Initialize();
    OPSettings.ReadInput(InputFileName);

    BC.Initialize(OPSettings);
    BC.ReadInput(InputFileName);
    //ES.Setup_MPI(OPSettings, BC);

    DO   .Initialize(OPSettings);
    ID   .Initialize(OPSettings);
    IP   .Initialize(OPSettings);
    Phase.Initialize(OPSettings);
    RTC  .Initialize(OPSettings);
    Timer.Initialize(OPSettings, "Timer");

    RTC  .ReadInput(InputFileName);
    DO   .ReadInput(InputFileName);
    ID   .ReadInput(InputFileName);
    IP   .ReadInput(InputFileName);
    Phase.ReadInput(InputFileName);

    const unsigned int Nx = OPSettings.Nx;
    const unsigned int Ny = OPSettings.Ny;
    const unsigned int Nz = OPSettings.Nz;

    int index[3] = {};
    // Set initial geometry of the phases
    index[0] = op::Initializations::Single(Phase, 0, BC, OPSettings);
    index[1] = op::Initializations::Sphere(Phase, 1, Ny/4, Nx/2, Ny/2, Nz/2, BC,
            OPSettings);
    index[2] = op::Initializations::SectionalPlane(Phase,  2,
            {Nx/2.,Ny/2.,Nz/2.}, {0,1,0}, BC, OPSettings);


    const double Volume20 = Phase.FieldsStatistics[index[1]].Volume;

    double ElapsedTime = 0.0;
    //----------------------------- The Time Loop ----------------------------//
    for(int tStep = RTC.tStart; tStep <= RTC.nSteps; tStep++)
    {
        // Output to files
        if (!(tStep%RTC.tFileWrite))
        {
            Phase.WriteVTK(tStep, OPSettings);
        }

        // Output to screens
        if(!(tStep%RTC.tScreenWrite))
        {
            const double AvInterfaceEnergie = DO.Energy(Phase,IP);
            const double Volume1 = Phase.FieldsStatistics[index[0]].Volume;
            const double Volume2 = Phase.FieldsStatistics[index[1]].Volume;
            const double Volume3 = Phase.FieldsStatistics[index[2]].Volume;

            double Error = std::abs(Volume2 - Volume20)/Volume20;

            ElapsedTime += RTC.dt;
            op::Info::WriteTimeStep(tStep, RTC.nSteps);
            op::Info::Write("Elapsed time     [s]", ElapsedTime);
            op::Info::Write("Interface energy [J]", AvInterfaceEnergie);
            op::Info::Write("Volume phase 1", Volume1);
            op::Info::Write("Volume phase 2", Volume2);
            op::Info::Write("Volume phase 3", Volume3);
            op::Info::WriteBlankLine();
            op::Info::Write("Relative error", Error);
            op::Info::WriteBlankLine();

            const double MaxError = 50*DBL_EPSILON;
            if ((Error > MaxError) or (tStep == RTC.nSteps))
            {
                if (Error < MaxError)
                {
                    op::Info::WriteBlankLine();
                    op::Info::WriteLine("_");
                    op::Info::WriteBlankLine();
                    op::Info::WriteSimple("Benchmark successfully completed");
                    op::Info::WriteLine("_");
                    op::Info::WriteBlankLine();
                }

                // Write simulation results
                std::ofstream resultsSim ("Results.sim");
                if (resultsSim.is_open())
                {
                    resultsSim << 2 << "\n";
                    resultsSim << "RelativeVolumeError " << Error << " "
                        << MaxError << "\n";
                    resultsSim.flush();
                    resultsSim.close();
                }
                break;
            }
        }


        // Actual calculation of the interface diffusion
        Timer.SetStart();
        IP.Set(Phase);
        Timer.SetTimeStamp("Sigma.Set()");
        ID.CalculatePhaseFieldIncrements(Phase, IP);
        Timer.SetTimeStamp("ID.Calculate()");
        Phase.MergeIncrements(BC, RTC.dt, false);
        Timer.SetTimeStamp("Phase.MergeIncrements()");

        if(!(tStep%RTC.tScreenWrite)) Timer.PrintWallClockSummary();
    }
    return 0;
}
