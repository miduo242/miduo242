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
#include "Nucleation.h"
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
    Nucleation                          Nuc(OPSettings);

    cout << "Initialization stage done!" << endl;

    if(RTC.Restart)
    {
        cout << "Restart data being read!";
        Phi.Read(BC, RTC.tStart);
        Cx.Read(BC, RTC.tStart);
        Tx.Read(BC, RTC.tStart);
        Nuc.Read(RTC.tStart);
        cout << " Done!" << endl;
    }
    else
    {
        // Single liquid phase starting configuration
        size_t idx1 = Initializations::Single(Phi, 0, BC, OPSettings);
        // Adding the nucleus in the middle of the simulation domain
        size_t idx2 = Phi.PlantGrainNucleus(1, (OPSettings.Nx)/2,(OPSettings.Ny)/2,(OPSettings.Nz)/2);
        // Orienting the nucleus as desired
        //EulerAngles locAngle({Pi/2.0, 0.0*Pi/8.0, Pi/8.0}, XYZ);
        EulerAngles locAngle({5.0*Pi/8.0, 3.0*Pi/8.0, Pi/8.0}, XYZ);
        Phi.FieldsStatistics[idx2].Orientation = locAngle.getQuaternion();

        // Liquid channel surrounded by alpha phase
        /*std::vector<size_t> ids = Initializations::TwoWalls(Phi, 0, 1, 10,BC,OPSettings);
        Phi.FieldsStatistics[ids[0]].State = Liquid;
        Phi.FieldsStatistics[ids[1]].State = Solid;*/

        // Liquid pocket surrounded by alpha phase
        /*size_t idx1 = Initializations::Single(Phi, 1, BC, OPSettings);
        Phi.FieldsStatistics[idx1].State = AggregateStates::Solid;
        //size_t idx2 = Initializations::Sphere(Phi, 0, OPSettings.Nz/2 - 4, (OPSettings.Nx)/2,(OPSettings.Ny)/2,(OPSettings.Nz)/2, BC, OPSettings);
        size_t idx2 = Initializations::Ellipsoid(Phi, 0, OPSettings.Nx/2 - 5, OPSettings.Ny/2 - 5, OPSettings.Nz/2 - 5, (OPSettings.Nx)/2, (OPSettings.Ny)/2, (OPSettings.Nz)/2, BC, OPSettings);
        Phi.FieldsStatistics[idx2].State = AggregateStates::Liquid;*/

        // Setting initial values for composition and temperature
        Cx.SetInitialMoleFractions(Phi);
        Tx.SetInitial(BC);
        DF.SetDiffusionCoefficients(Phi, Tx);
    }

    fstream tt_out("TimeTemperature.dat", ios::out);
    tt_out << "Time    dt    Temperature" << endl;
    tt_out.close();

    cout << "Entering the Time Loop!!!" << endl;

    //-------------- The Time Loop -------------//
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        Nuc.GenerateNucleationSites(Phi, Tx);
        Nuc.PlantNuclei(Phi, RTC.tStep);
        IP.Set(Phi,Tx);
        DF.CalculateInterfaceMobility(Phi,Cx,Tx,BC,IP);
        DO.CalculatePhaseFieldIncrements(Phi, IP);
        DF.GetDrivingForce(Phi, Cx, Tx, dG);
        dG.Average(Phi, BC);
        Nuc.CheckNuclei(Phi, IP, dG, RTC.tStep);

        if (RTC.WriteVTK())
        {//  Write data in VTK format
            dG.WriteVTK(RTC.tStep, OPSettings, 1, 0);
            dG.WriteVTKforPhases(RTC.tStep, OPSettings, Phi);
        }

        dG.MergePhaseFieldIncrements(Phi, IP);
        Phi.NormalizeIncrements(BC, RTC.dt);
        DF.Solve(Phi, Cx, Tx, BC, RTC.dt);
        Phi.MergeIncrements(BC, RTC.dt);
        //Tx.Set(BC, Phi, 6.2e8, 1.773e6, 1, dt);
        Tx.Set(BC, Phi, RTC.dt);

        if (RTC.WriteRawData())
        {// Write raw data
            Phi.Write(RTC.tStep);
            Cx.Write(RTC.tStep);
            Tx.Write(RTC.tStep);
            Nuc.Write(RTC.tStep);
        }
        if (RTC.WriteVTK())
        {// Write data in VTK format
            Phi.WriteVTK(RTC.tStep,OPSettings);
            Phi.WriteLaplacianVTK(RTC.tStep, OPSettings, 1);

            Cx.WriteVTK(RTC.tStep,OPSettings);
            Cx.WriteStatistics(RTC.tStep, RTC.dt);
        }
        if(RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message = Info::GetStandard("Interface energy density", to_string(I_En));
            Info::WriteTimeStep(RTC, message);

            dG.PrintDiagnostics();
            Phi.PrintVolumeFractions();

            fstream tt_out("TimeTemperature.dat", ios::out | ios::app);
            tt_out << RTC.SimulationTime <<  "   " << RTC.dt <<  "   " << Tx(1,1,1) << endl;
            tt_out.close();
        }
    } //end time loop
    return 0;
}
