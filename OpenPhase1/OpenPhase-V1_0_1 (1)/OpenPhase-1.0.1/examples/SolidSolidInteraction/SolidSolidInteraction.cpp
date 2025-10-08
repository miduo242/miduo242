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
#include "InterfaceProperties.h"
#include "RunTimeControl.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "Composition.h"
#include "Temperature.h"
#include "Initializations.h"
#include "FluidDynamics/InteractionSolidFluid.h"
#include "FluidDynamics/InteractionSolidSolid.h"
#include "Velocities.h"
#include "BoundaryConditions.h"
#include "Tools/TimeInfo.h"
#include <iostream>

using namespace std;
using namespace openphase;

/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
  {
    feenableexcept(FE_INVALID | FE_OVERFLOW);
    Settings                        OPSettings;
    OPSettings.ReadInput();

    RunTimeControl                  RTC(OPSettings);
    PhaseField                      Phi(OPSettings);
    DoubleObstacle                  DO(OPSettings);
    InterfaceProperties             IP(OPSettings);
    Velocities                      Vel(OPSettings);
    BoundaryConditions              BC(OPSettings);
    TimeInfo                        Timer;

    Timer.Initialize(OPSettings, "Execution Time Statistics");


    double dt = RTC.dt;
    int sol1 = 0;
    int sol2 = 0;
    //int sol3 = 0;
    vector<int> systemPhFields;
    ignore_result(system("mkdir VTK/Density"));
    if(RTC.Restart)
    {
        cout << "Restart data being read!";

        Phi.Read(BC, RTC.tStart);
        cout << " Done!" << endl;
    }
    else
    {
        RTC.tStart = 0;
        int ad = Initializations::Single(Phi, 0, BC, OPSettings);
        sol1 = Initializations::Sphere(Phi,1,15,30,50,40,BC, OPSettings);
        sol2 = Initializations::Sphere(Phi,1,15,70,50,60,BC, OPSettings);
        Phi.FieldsStatistics[sol1].Mobile = true;
        Phi.FieldsStatistics[sol1].Density = 1.0;
        Phi.FieldsStatistics[sol2].Density = 1.0;
        Phi.FieldsStatistics[sol2].Mobile = true;

        dVector3 initVel;
        initVel.set_to_zero();
        initVel[0] =  0.1;
        initVel[1] =  0.0;
        initVel[2] =  0.0;
        Phi.FieldsStatistics[sol1].Vcm = initVel;
        initVel[0] = -0.1;
        initVel[1] = -0.0;
        initVel[2] = -0.0;
        Phi.FieldsStatistics[sol2].Vcm = initVel;
    }

    cout << "Initialization stage done!" << endl;

    cout << "Entering the Time Loop!!!" << endl;
	std::vector<double> RefVolume;
    InteractionSolidSolid::SetRefVolume(Phi, RefVolume);
    //-------------- The Time Loop -------------//

    for(int tStep = RTC.tStart; tStep < RTC.nSteps+1; tStep++)
    {
        DO.CalculatePhaseFieldIncrements(Phi, IP);
        Phi.NormalizeIncrements(BC, RTC.dt);
        Phi.MergeIncrements(BC, dt);
        int iterations = 0;
        while (InteractionSolidSolid::VolumeError(Phi, RefVolume, 1) > 1.0e-4)
        {
            InteractionSolidSolid::PreserveVolume(Phi, RefVolume, 1, dt);
            Phi.MergeIncrements(BC, dt);
            ++iterations;
            if (iterations > 10)
            {
                std::cout << "Max Iterations in VolumeControl!" << std::endl;
                break;
            }
        }
        for(unsigned int idx = 0; idx < Phi.FieldsStatistics.size(); idx++)
        {
            Phi.FieldsStatistics[idx].Force.set_to_zero();
            Phi.FieldsStatistics[idx].Torque.set_to_zero();
        }
        InteractionSolidSolid::Calculate(Phi, RTC.dt, 4, 3, 0.01, 0.0);
        InteractionSolidFluid::CalculateSolidVelocities(Phi, Vel, BC, OPSettings, RTC.dt);
        InteractionSolidSolid::AdvectSolid(Phi, BC, RTC.dt);
        Phi.NormalizeIncrements(BC, RTC.dt);
        Phi.MergeIncrements(BC, RTC.dt);
        if (!(tStep%RTC.tFileWrite))
        {
            Phi.WriteVTK(tStep, OPSettings);
        }
        cout << "Time step: " << tStep << endl;
        time_t rawtime;
        time(&rawtime);

    } //end time loop
    return 0;
}
