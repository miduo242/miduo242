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

    if (RTC.Restart)
    {
        cout << "Restart data being read!";

        Phi.Read(BC, RTC.tStart);
        Cx.Read(BC, RTC.tStart);
        Tx.Read(BC, RTC.tStart);
        cout << " Done!" << endl;
    }
    else
    {
        int iRadius = 6;
        Initializations::Fractional(Phi, 0, 1, iRadius / 2, BC, OPSettings);
        Initializations::Sphere(Phi, 2, iRadius, 0, 0, 0, BC, OPSettings);
        Cx.SetInitialMoleFractions(Phi);
        Tx.SetInitial(BC);
        cout << "Initialization stage done!" << endl;
    }

    cout << "Entering the Time Loop!!!" << endl;

    //-------------- The Time Loop -------------//

    for (RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        IP.Set(Phi);
        DO.CalculatePhaseFieldIncrements(Phi, IP);

        DF.GetDrivingForce(Phi, Cx, Tx, dG);
        dG.Average(Phi, BC);
        dG.MergePhaseFieldIncrements(Phi, IP);
        Phi.NormalizeIncrements(BC, RTC.dt);

        DF.Solve(Phi, Cx, Tx, BC, RTC.dt);
        Tx.Set(BC, Phi, RTC.dt);
        Phi.MergeIncrements(BC, RTC.dt);

        if (RTC.WriteVTK())
        {//  Output to file in VTK format
            Phi.WriteVTK(RTC.tStep, OPSettings);
            Cx.WriteVTK(RTC.tStep, OPSettings);
            Tx.WriteVTK(RTC.tStep, OPSettings);
            Cx.WriteStatistics(RTC.tStep, RTC.dt);
        }

        if (RTC.WriteRawData())
        {//  Output to file of raw data
            Phi.Write(RTC.tStep);
            Cx.Write(RTC.tStep);
            Tx.Write(RTC.tStep);
        }

        if (RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);

            std::string message  = Info::GetStandard("Interface energy density", to_string(I_En));
                        message += Info::GetStandard("Temperature [K]"         , to_string(Tx(0,0,0)));
            Info::WriteTimeStep(RTC, message);

            dG.PrintDiagnostics();
            Phi.PrintVolumeFractions();
        }
    } //end time loop

    return 0;
}
