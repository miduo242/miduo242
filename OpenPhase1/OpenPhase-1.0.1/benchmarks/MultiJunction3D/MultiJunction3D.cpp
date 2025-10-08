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
#include "Initializations.h"
#include "BoundaryConditions.h"

using namespace std;
using namespace openphase;

/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
{
    Settings                    OPSettings;
    OPSettings.ReadInput();

    RunTimeControl              RTC(OPSettings);
    PhaseField                  Phi(OPSettings);
    DoubleObstacle              DO(OPSettings);
    InterfaceProperties         IP(OPSettings);
    BoundaryConditions          BC(OPSettings);

    Initializations::Young4(Phi, 0, BC, OPSettings);

    int x1 = (OPSettings.Nx)/2;
    int y1 = (OPSettings.Ny)/2;
    int z1 = (OPSettings.Nz)/2;

    cout << "Entering the Time Loop!!!" << endl;

    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        IP.Set(Phi);
        DO.CalculatePhaseFieldIncrements(Phi, IP);
        Phi.NormalizeIncrements(BC, RTC.dt);
        Phi.MergeIncrements(BC, RTC.dt);

        // Keeping the quadruple junction in the center of the simulation box
        double x2 = 0;
        double y2 = 0;
        double z2 = 0;
        double count = 1.0;
        for(int i = 0;i < OPSettings.Nx; ++i)
        for(int j = 0;j < OPSettings.Ny; ++j)
        for(int k = 0;k < OPSettings.Nz; ++k)
        if(Phi.Fields(i,j,k).size() == 4)
        {
            count ++;
            x2 += i;
            y2 += j;
            z2 += k;
        }
        if(count)
        {
            x2 /= count;
            y2 /= count;
            z2 /= count;
        }

        double tempdx = x2 - x1;
        double tempdy = y2 - y1;
        double tempdz = z2 - z1;

        int dx = 0;
        int dy = 0;
        int dz = 0;
        if (tempdx >= 1) dx = 1;
        if (tempdy >= 1) dy = 1;
        if (tempdz >= 1) dz = 1;
        if (tempdx <= -1) dx = -1;
        if (tempdy <= -1) dy = -1;
        if (tempdz <= -1) dz = -1;
        if (abs(dx) + abs(dy) + abs(dz) != 0)
        {
            Phi.MoveFrame(dx, dy, dz, BC);
        }

        if (RTC.WriteVTK())
        {
           //  Output to file
           Phi.WriteVTK(RTC.tStep,OPSettings);
        }

        if (RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message  = Info::GetStandard("Interface energy density", to_string(I_En));
            Info::WriteTimeStep(RTC, message);

            //  For control (values at quadruple junction)
            Phi.PrintPointStatistics(x1,y1,z1);
            Phi.PrintPointStatistics(int(x2),int(y2),int(z2));
        }
    } //end of time loop
    return 0;
}
