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

#include "Settings.h"
#include "RunTimeControl.h"
#include "InterfaceProperties.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "Composition.h"
#include "Temperature.h"
#include "EquilibriumPartitionDiffusionBinary.h"
#include "BoundaryConditions.h"
#include "Initializations.h"
#include "Tools/TimeInfo.h"
#include "Tools/MicrostructureAnalysis.h"

using namespace std;
using namespace openphase;

/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
{
    Settings                            OPSettings;
    OPSettings.ReadInput();

    RunTimeControl                      RTC(OPSettings);
    PhaseField                          Phi(OPSettings);
    DoubleObstacle                      DO(OPSettings);
    InterfaceProperties                 IP(OPSettings);
    EquilibriumPartitionDiffusionBinary DF(OPSettings);
    Composition                         Cx(OPSettings);
    Temperature                         Tx(OPSettings);
    DrivingForce                        dG(OPSettings);
    BoundaryConditions                  BC(OPSettings);
    TimeInfo                            Timer(OPSettings, "Execution Time Statistics");

    if(RTC.Restart)
    {
        cout << "Restart data being read!" << endl;
        cout << "Restart time step: " << RTC.tStart << endl;

        Phi.Read(BC, RTC.tStart);
        Cx.Read(BC, RTC.tStart);
        Tx.Read(BC, RTC.tStart);

        cout << "Done reading restart parameters!" << endl;
    }
    else
    {
        int idx0 = Initializations::Single(Phi, 0, BC, OPSettings);
        Phi.FieldsStatistics[idx0].State = AggregateStates::Liquid;

        int idx1 = Phi.PlantGrainNucleus(1, OPSettings.Nx/2, OPSettings.Ny/2, OPSettings.Nz/2);
        //int idx1 = Initializations::Sphere(Phi, 1, 5, 0, 0, 0, BC, OPSettings);
        Phi.FieldsStatistics[idx1].State = AggregateStates::Solid;
        //EulerAngles angle({0.0, 0.0, 0.0},ZXZ);
        //Phi.FieldsStatistics[idx1].Orientation = angle.getQuaternion();

        Cx.SetInitialMoleFractions(Phi);
        Tx.SetInitial(BC);
        DF.SetDiffusionCoefficients(Phi, Tx);
    }

    cout << "Initialization stage done!" << endl;
    cout << "Entering the Time Loop!!!" << endl;


    RTC.dt = DF.ReportMaximumTimeStep();
    cout << "Time step: " << RTC.dt << endl;

    //-------------- The Time Loop -------------//
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        Timer.SetStart();

        IP.Set(Phi);
        DF.CalculateInterfaceMobility(Phi, Cx, Tx, BC, IP);
        Timer.SetTimeStamp("Calculate/Set interface properties");

        DO.CalculatePhaseFieldIncrements(Phi, IP);
        Timer.SetTimeStamp("Get PsiDot");

        DF.GetDrivingForce(Phi, Cx, Tx, dG);
        Timer.SetTimeStamp("Chemical driving force");

        dG.Average(Phi, BC);
        Timer.SetTimeStamp("Driving force average");

        if (RTC.WriteVTK())
        {
            dG.WriteVTKforPhases(RTC.tStep, OPSettings, Phi);
        }
        Timer.SetTimeStamp("Write driving force");

        dG.MergePhaseFieldIncrements(Phi, IP);
        Timer.SetTimeStamp("Driving force merge to Psi");

        Phi.NormalizeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("Normalize Psi");

        DF.Solve(Phi, Cx, Tx, BC, RTC.dt);
        Timer.SetTimeStamp("Solve diffusion");

        Tx.Set(BC, Phi, RTC.dt);

        Timer.SetTimeStamp("Set temperature");
        Phi.MergeIncrements(BC, RTC.dt);

        Timer.SetTimeStamp("Merge phase fields");

        //  Output to file
        if (RTC.WriteVTK())
        {
            // Write data in VTK format
            Phi.WriteVTK(RTC.tStep,OPSettings);
            Cx.WriteVTK(RTC.tStep,OPSettings);
            Tx.WriteVTK(RTC.tStep,OPSettings);
            Cx.WriteStatistics(RTC.tStep, RTC.dt);
            IP.WriteVTK(Phi,RTC.tStep);
        }
        if (RTC.WriteRawData())
        {
            // Write raw data
            Phi.Write(RTC.tStep);
            Cx.Write(RTC.tStep);
            Tx.Write(RTC.tStep);
        }
        Timer.SetTimeStamp("File Output");

        //  Output to screen
        if(RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message = Info::GetStandard("Interface energy density", to_string(I_En));
            Info::WriteTimeStep(RTC, message);

            dG.PrintDiagnostics();
            Phi.PrintVolumeFractions();
            Timer.PrintWallClockSummary();
        }
    } //end time loop

    return 0;
}
