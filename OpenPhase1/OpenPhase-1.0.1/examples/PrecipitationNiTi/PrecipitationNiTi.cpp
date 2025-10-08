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

#include "EquilibriumPartitionDiffusionBinary.h"
#include "Settings.h"
#include "RunTimeControl.h"
#include "InterfaceProperties.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "Composition.h"
#include "Temperature.h"
#include "Initializations.h"
#include "Mechanics/ElasticProperties.h"
#include "Mechanics/ElasticitySolverSpectral.h"
#include "BoundaryConditions.h"
#include "Nucleation.h"
#include "Info.h"
#include "Tools/TimeInfo.h"

using namespace std;
using namespace openphase;
/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
{
#ifdef MPI_PARALLEL
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);
    {
#endif

    Settings                         OPSettings;
    OPSettings.ReadInput();

    RunTimeControl                   RTC(OPSettings);
    PhaseField                       Phi(OPSettings);
    DoubleObstacle                   DO(OPSettings);
    InterfaceProperties              IP(OPSettings);
    EquilibriumPartitionDiffusionBinary DF(OPSettings);
    Composition                      Cx(OPSettings);
    Temperature                      Tx(OPSettings);
    DrivingForce                     dG(OPSettings);
    ElasticProperties                EP(OPSettings);
    ElasticitySolverSpectral         ES(OPSettings);
    BoundaryConditions               BC(OPSettings);
    Nucleation                       Nuc(OPSettings);

    TimeInfo                         Timer;
    Timer.Initialize(OPSettings, "Execution Time Statistics");

    cout << "Initialization stage! Done!" << endl;

    //int x1 = OPSettings.TotalNx/2;
    //int y1 = OPSettings.Ny/2;
    //int z1 = OPSettings.Nz/2;

    int matrixIdx = Initializations::Single(Phi, 0, BC, OPSettings);
    EulerAngles ph1({0,0,0},XYZ);
    Phi.FieldsStatistics[matrixIdx].Orientation = ph1.getQuaternion();
    // First particle nucleation:
    //int particleIdx = Phi.PlantGrainNucleus(1, x1, y1, z1);
    //Phi.FieldsStatistics[particleIdx].Variant = 0;
    //Phi.FieldsStatistics[particleIdx].Orientation = ph1.getQuaternion();

    // Second particle nucleation:
    //int particleIdx2 = Phi.PlantGrainNucleus(1, x1, y1, z1+5);
    //Phi.FieldsStatistics[particleIdx2].Orientation = ph1.getQuaternion();
    //Phi.FieldsStatistics[particleIdx2].Variant = 1;

    Cx.SetInitialMoleFractions(Phi);
    EP.SetGrainsProperties(Phi);

    //EP.SetEffectiveElasticConstants(Phi,Cx);
    EP.SetEffectiveElasticConstants(Phi);
    EP.SetEffectiveTransformationStretches(Phi);

    Tx.SetInitial(BC);

    /*fstream output_file;
    output_file.open(DefaultTextDir + "ParticleSize.txt", fstream::out);

    output_file << "tStep" << "\t" << "sim_time" << "\t"
                << "part_diameter" << "\t" << "part_height" << endl;
    output_file.close();*/

    cout << "Entering the Time Loop!!!" << endl;

    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
    Timer.SetStart();

        Nuc.GenerateNucleationSites(Phi, Tx);
        Nuc.PlantNuclei(Phi, RTC.tStep);
    Timer.SetTimeStamp("Plant nuclei");

        IP.Set(Phi);
    Timer.SetTimeStamp("Set interface properties");

        DO.CalculatePhaseFieldIncrements(Phi, IP);
        //DO.CalculatePhaseFieldIncrements(Phi, IP, dG);
    Timer.SetTimeStamp("Get PhiDot");

        EP.SetGrainsProperties(Phi);
        EP.SetEffectiveElasticConstants(Phi/*,Cx*/);
        EP.SetEffectiveTransformationStretches(Phi);
        EP.CalculateDrivingForce(Phi, dG);
    Timer.SetTimeStamp("Elastic Driving Force");

        DF.GetDrivingForce(Phi, Cx, Tx, dG);
    Timer.SetTimeStamp("Chemical Driving Force");

        dG.Average(Phi, BC);
    Timer.SetTimeStamp("Driving Force Average");

        Nuc.CheckNuclei(Phi, IP, dG, RTC.tStep);
    Timer.SetTimeStamp("Check nuclei");

        if (RTC.WriteVTK())//  Output to file
        {
            dG.WriteVTKforPhases(RTC.tStep, OPSettings, Phi);
        }
    Timer.SetTimeStamp("Driving Force VTK Output");
        dG.MergePhaseFieldIncrements(Phi, IP);
        Timer.SetTimeStamp("Merge driving force");

        Phi.NormalizeIncrements(BC, RTC.dt);
    Timer.SetTimeStamp("Merge Velocities");

        ES.Solve(EP, BC, RTC.dt);
    Timer.SetTimeStamp("Elastic Solver");

        //EP.CalculateChemicalPotentialContribution(Phi, DF);
        DF.Solve(Phi, Cx, Tx, BC, RTC.dt);
    Timer.SetTimeStamp("Diffusion Solver");

        Phi.MergeIncrements(BC, RTC.dt);
    Timer.SetTimeStamp("Merge PhaseFields");

        Tx.Set(BC, Phi, RTC.dt);
    Timer.SetTimeStamp("Set Temperature");

        if (RTC.WriteRawData())
        {
            Phi.Write(RTC.tStep);
            Cx.Write(RTC.tStep);
            Tx.Write(RTC.tStep);
        }
        if (RTC.WriteVTK())
        {
            Phi.WriteVTK(RTC.tStep,OPSettings);
            Cx.WriteVTK(RTC.tStep,OPSettings);
            Cx.WriteStatistics(RTC.tStep, RTC.dt);
            EP.WriteStressesVTK(RTC.tStep,OPSettings);
            EP.WriteTotalStrainsVTK(RTC.tStep,OPSettings);
            EP.WriteEigenStrainsVTK(RTC.tStep,OPSettings);
            EP.WriteEffectiveElasticConstantsVTK(RTC.tStep,OPSettings);
        }
    Timer.SetTimeStamp("File Output");
        if (RTC.WriteToScreen())
        {
            double E_En = EP.AverageEnergyDensity();
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message  = Info::GetStandard("Interface energy density", to_string(I_En));
                        message += Info::GetStandard("Elastic energy density", to_string(E_En));
            Info::WriteTimeStep(RTC, message);

            //Phi.PrintPointStatistics(x1,y1,z1);
            //Cx.PrintPointStatistics(x1,y1,z1);
            //Tx.PrintPointStatistics(x1,y1,z1);
            //EP.PrintPointStatistics(x1,y1,z1);
            dG.PrintDiagnostics();

            /*double part_diameter = 0.0;
            double part_height = 0.0;
            // collect data over line and calculate the particle's diameter and height:
            for(int i = 1; i < OPSettings.Nx+1; ++i)
            {
                part_diameter += Phi.Fields(i, OPSettings.Ny/2+1, OPSettings.Nz/2+1)[3];
                part_height   += Phi.Fields(OPSettings.Nx/2+1, OPSettings.Ny/2+1, i)[3];
            }

            fstream output_file;
            output_file.open(DefaultTextDir + "ParticleSize.txt", fstream::out | fstream::app);
            output_file << RTC.tStep << "\t" << sim_time << "\t"
                        << part_diameter << "\t" << part_height << endl;
            output_file.close();*/

            //Timer.PrintWallClockSummary();
            Phi.PrintVolumeFractions();
        }
    } //end time loop
    simulation_end();

#ifdef MPI_PARALLEL
    }
    MPI_Finalize ();
#endif
    return 0;
}
