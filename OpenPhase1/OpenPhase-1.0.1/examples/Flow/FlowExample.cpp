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
#include "Initializations.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"
#include "TextOutput.h"
#include "Tools/TimeInfo.h"
#include "Temperature.h"
#include "Composition.h"
#include "InterfaceProperties.h"
#include "DoubleObstacle.h"
#include "DrivingForce.h"
#include "EquilibriumPartitionDiffusionBinary.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "Velocities.h"
#include "AdvectionHR/AdvectionHR.h"
#include "BoundaryConditions.h"

using namespace std;
using namespace openphase;

/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
{
    string InputFile = "ProjectInput.opi";

    Settings                            OPSettings(InputFile);
    RunTimeControl                      RTC(OPSettings,InputFile);
    PhaseField                          Phi(OPSettings);
    BoundaryConditions                  BC(OPSettings,InputFile);
    Temperature                         Tx(OPSettings,InputFile);
    Composition                         Cx(OPSettings,InputFile);
    InterfaceProperties                 IP(OPSettings,InputFile);
    DoubleObstacle                      DO(OPSettings);
    DrivingForce                        dG(OPSettings,InputFile);
    EquilibriumPartitionDiffusionBinary DF(OPSettings,InputFile);
    FlowSolverLBM                       FL(OPSettings, RTC.dt, InputFile);
    Velocities                          Vel(OPSettings);
    AdvectionHR                         ADHR (OPSettings);
    TimeInfo                            Timer;

    int x1 = OPSettings.Nx/2;
    int y1 = OPSettings.Ny/2;
    int z1 = OPSettings.Nz/2;

    if(RTC.Restart)
    {
        cout << "Restart data being read!";
        Phi.Read(BC, RTC.tStart);
        Cx.Read(BC, RTC.tStart);
        cout << " Done!" << endl;
    }
    else
    {
        // Single liquid phase starting configuration
        size_t idx1 = Initializations::Single(Phi, 0, BC, OPSettings);
        // Adding the nucleus in the middle of the simulation domain
        size_t idx2 = Initializations::Sphere(Phi, 1, 5.0, x1, y1, z1, BC, OPSettings);
        // Orienting the nucleus as desired
        //EulerAngles locAngle({Pi/2.0, 0.0*Pi/8.0, Pi/8.0}, XYZ);
        EulerAngles locAngle({5.0*Pi/8.0, 3.0*Pi/8.0, Pi/8.0}, XYZ);
        Phi.FieldsStatistics[idx2].Orientation = locAngle.getQuaternion();

        // Liquid channel surrounded by alpha phase
        /*std::vector<size_t> ids = Initializations::TwoWalls(Phi, 0, 1, 10,BC,OPSettings);*/
    }

    FL.SetUniformVelocity(BC, dVector3::ZeroVector());
    Cx.SetInitialMoleFractions(Phi);
    Tx.SetInitial(BC);
    DF.SetDiffusionCoefficients(Phi, Tx);
    cout << "Initialization stage done!" << endl;
    cout << "Entering the Time Loop!!!" << endl;

    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        Timer.SetStart();

        IP.Set(Phi, Tx);
        Timer.SetTimeStamp("Interface properties");

        DO.CalculatePhaseFieldIncrements(Phi, IP);
        Timer.SetTimeStamp("DO.CalculatePhaseFieldIncrements");

        DF.GetDrivingForce(Phi, Cx, Tx, dG);
        Timer.SetTimeStamp("DF.dG");

        dG.Average(Phi, BC);
        Timer.SetTimeStamp("dG.Average");

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

        FL.Solve(Phi, Cx, Vel, BC);
        Vel.SetBoundaryConditions(BC);
        Timer.SetTimeStamp("FL.Solve");

        DF.CalculatePhaseConcentrations(Phi, Tx, Cx);
        Cx.Advect(ADHR, Vel, Phi, BC, RTC.dt, RTC.tStep);
        DF.CalculatePhaseConcentrations(Phi, Tx, Cx);
        Timer.SetTimeStamp("Adv.Advect Cx");

        //  Output to file
        if (RTC.WriteVTK())
        {
            // Write data in VTK format
            Phi.WriteVTK(RTC.tStep,OPSettings);
            Cx.WriteVTK(RTC.tStep, OPSettings);
            Tx.WriteVTK(RTC.tStep, OPSettings);
            Vel.WriteVTK(RTC.tStep, OPSettings);
            //Cx.WriteStatistics(tStep, OPSettings.dt);
            FL.WriteVTK(RTC.tStep, Phi, OPSettings);
        }

        // Write restart output
        if(RTC.WriteRawData())
        {
            // Write raw data
            Phi.Write(RTC.tStep);
            Cx.Write(RTC.tStep);
            Tx.Write(RTC.tStep);
            FL.Write(RTC.tStep);
        }

        //  Output to screen
        if(RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message = Info::GetStandard("Interface energy density", to_string(I_En));
            Info::WriteTimeStep(RTC, message);

            //  Statistics
            Phi.PrintPointStatistics(x1,y1,z1);
            Cx.PrintPointStatistics(x1,y1,z1);
            Tx.PrintPointStatistics(x1,y1,z1);
            Vel.PrintPointStatistics(x1,y1,z1);
            dG.PrintDiagnostics();
            Phi.PrintPFVolumes();
        }
    } //end time loop
    return 0;
}
