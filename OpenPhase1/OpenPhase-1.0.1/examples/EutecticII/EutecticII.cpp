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

using namespace std;
using namespace openphase;
/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
{
    Settings                        OPSettings;
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

    cout << "Initialization stage! Done!" << endl;

    if(RTC.Restart)
    {
        cout << "Restart data being read!";

        Phi.Read(BC, RTC.tStart);
        Cx.Read(BC, RTC.tStart);
        Tx.Read(BC, RTC.tStart);
        cout << " Done!" << endl;
    }
    else
    {
        RTC.tStart = -1;
        double iRadius = 0.25*OPSettings.Nx;
        // single lamella
        Initializations::Fractional(Phi, 0, 1, 0.5*iRadius, BC, OPSettings);
        Initializations::Sphere(Phi, 2, iRadius, 0.5*(OPSettings.Nx-1), 0.5*(OPSettings.Ny-1), 0, BC, OPSettings);

        // two lamellae
        /*Initializations::Fractional(Phi, 0, 1, iRadius/2.0, BC, OPSettings);
        Initializations::Sphere(Phi, 2, iRadius/2, (OPSettings.Nx)/4, (OPSettings.Ny)/2, 0, BC, OPSettings);
        Initializations::Sphere(Phi, 2, iRadius/2, (OPSettings.Nx)*3/4, (OPSettings.Ny)/2, 0, BC, OPSettings);
        */
        Cx.SetInitialMoleFractions(Phi);
        Tx.SetInitial(BC);
    }

    cout << "Entering the Time Loop!!!" << endl;
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
Timer.SetStart();
        IP.Set(Phi);
Timer.SetTimeStamp("IP.Set");
        DO.CalculatePhaseFieldIncrements(Phi, IP);
Timer.SetTimeStamp("DO.CalculatePhaseFieldIncrements");
        DF.GetDrivingForce(Phi, Cx, Tx, dG);
Timer.SetTimeStamp("DF.GetDrivingForce");
        dG.Average(Phi, BC);
Timer.SetTimeStamp("dG.Average");
        //if (RTC.WriteVTK()) dG.WriteVTKforPhases(Phi,tStep);
        dG.MergePhaseFieldIncrements(Phi, IP);
Timer.SetTimeStamp("dG.MergePhaseFieldIncrements");
        Phi.NormalizeIncrements(BC, RTC.dt);
Timer.SetTimeStamp("Phi.NormalizeIncrements");
        DF.Solve(Phi, Cx, Tx, BC, RTC.dt);
Timer.SetTimeStamp("DF.Solve");
        Tx.Set(BC, Phi, RTC.dt);
Timer.SetTimeStamp("Tx.Set");
        Phi.MergeIncrements(BC, RTC.dt);
Timer.SetTimeStamp("Phi.MergeIncrements");

        /*if (Phi.Fields((OPSettings.Nx)/4, (OPSettings.Ny)/2, (OPSettings.Nz)/2)[0] <= 0.4)
        {// Moving frame
            Phi.MoveFrame(0,0,1, BC);
            Tx.MoveFrame(0,0,1, BC);
            Cx.MoveFrame(0,0,1, BC);
        }*/
        if (RTC.WriteVTK())
        {// Output to file in VTK format
            Phi.WriteVTK(RTC.tStep, OPSettings);
            Cx.WriteVTK(RTC.tStep, OPSettings);
            Cx.WriteStatistics(RTC.tStep, RTC.dt);
            //Tx.WriteVTK(tStep);
        }

        if (RTC.WriteRawData())
        {// Output to file of raw data
            Phi.Write(RTC.tStep);
            Cx.Write(RTC.tStep);
            Tx.Write(RTC.tStep);
        }
Timer.SetTimeStamp("File output");
        if (RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message  = Info::GetStandard("Interface energy density", to_string(I_En));
            Info::WriteTimeStep(RTC, message);

            dG.PrintDiagnostics();
            Phi.PrintPFVolumes();
            Timer.PrintWallClockSummary();
        }
    } //end time loop
    return 0;
}
