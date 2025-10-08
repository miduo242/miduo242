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
#include "BoundaryConditions.h"
#include "InterfaceProperties.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "Initializations.h"
#include "Mechanics/ElasticProperties.h"
#include "Mechanics/ElasticitySolverSpectral.h"
#include "UserDrivingForce.h"
#include "Temperature.h"
#include "Nucleation.h"
#include "TextOutput.h"

using namespace std;
using namespace openphase;

void WriteTemp_VMStress(ElasticProperties& EP, Temperature& Tx, int tStep, double RealTime, string OutFile);

int main(int argc, char *argv[])
{

#ifdef MPI_PARALLEL
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);
    {
#endif

    string MartensiteVolumeFractionFile = "MartensiteVolumeFraction.dat";
    remove(MartensiteVolumeFractionFile.c_str());

    string InputFile = "ProjectInput.opi";

    Settings                    OPSettings(InputFile);
    RunTimeControl              RTC(OPSettings,InputFile);
    BoundaryConditions          BC(OPSettings,InputFile);
    PhaseField                  Phi(OPSettings);
    DoubleObstacle              DO(OPSettings);
    InterfaceProperties         IP(OPSettings,InputFile);
    ElasticProperties           EP(OPSettings,InputFile);
    ElasticitySolverSpectral    ES(OPSettings);
    Temperature                 Tx(OPSettings,InputFile);
    Nucleation                  Nuc(OPSettings,InputFile);
    DrivingForce                dG(OPSettings,InputFile);
    UserDrivingForce            UDF(OPSettings,InputFile);
    Crystallography             CR(OPSettings, InputFile);

    string TempFile      = OPSettings.TextDir + "TemperatureFile.txt";
    string RSSFile       = OPSettings.TextDir + "RSSFile.txt";

    if(RTC.Restart)
    {
        cout << "Restart data being read! " << endl;
        Phi.Read(BC, RTC.tStart);
        Tx.Read(BC, RTC.tStart);
        Nuc.Read(RTC.tStart);
        EP.Read(BC, RTC.tStart);

        EP.SetGrainsProperties(Phi);

        cout << "Restart: solving mechanical problem! " << endl;
        EP.SetEffectiveElasticConstants(Phi);
        EP.SetEffectiveTransformationStretches(Phi);
        ES.Solve(EP,  BC, RTC.dt);
        RTC.SimulationTime = RTC.dt * RTC.tStart;
        cout << "Restart Done!" << endl;
    }
    else
    {
        Initializations::Single(Phi, 0, BC, OPSettings);
        //==========================================================================
        dMatrix3x3 locRot ({0.983877, -0.0648783, 0.166667,
                            0.052973,  0.995782,  0.0749147,
                           -0.170824, -0.064878,  0.983163});

        EP.PhaseElasticConstants[1].rotate(locRot);
        EP.PhaseTransformationStretches[1] = EP.PhaseTransformationStretches[1].rotatedU(locRot);

        Tx.SetInitial(BC);
        EP.SetGrainsProperties(Phi);
    }

    //---------------------------------------------------------------------------------------------------------------------------------------------
        cout << "Entering the Time Loop!!!" << endl;
    //---------------------------------------------------------------------------------------------------------------------------------------------

    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
    	Nuc.GenerateNucleationSites(Phi, Tx);
    	Nuc.PlantNuclei(Phi, RTC.tStep);

        EP.SetGrainsProperties(Phi);
        IP.Set(Phi);

        DO.CalculatePhaseFieldIncrements(Phi, IP);
        UDF.SetDrivingForce(Phi, dG, Tx);

        EP.SetEffectiveElasticConstants(Phi);
        EP.SetEffectiveTransformationStretches(Phi);
        EP.CalculateDrivingForce(Phi, dG);

        dG.Average(Phi, BC);
        Nuc.CheckNuclei(Phi, IP, dG, RTC.tStep);

        if (RTC.WriteVTK())
        {
            dG.WriteVTKforPhases(RTC.tStep, OPSettings, Phi);
        }

        dG.MergePhaseFieldIncrements(Phi, IP);
        Phi.NormalizeIncrements(BC, RTC.dt);

        if(Tx(0,0,0) > 273)
        {
            Tx.Set(BC, Phi, RTC.dt);
        }
        Phi.MergeIncrements(BC, RTC.dt);
        EP.SetGrainsProperties(Phi);
        EP.SetEffectiveElasticConstants(Phi);
        EP.SetEffectiveTransformationStretches(Phi);
        ES.Solve(EP, BC, RTC.dt);

        if (RTC.WriteToScreen())
        {/***********************************************/
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            double E_En = EP.AverageEnergyDensity();
            double T_En = E_En + I_En;

            //  Output to screen
            Info::WriteStandard("Simulation Time",RTC.SimulationTime);

            std::string message  = Info::GetStandard("Interface energy density", to_string(I_En));
                        message += Info::GetStandard("Elastic energy density", to_string(E_En));
                        message += Info::GetStandard("Total energy density", to_string(T_En));
                        message += Info::GetStandard("Temperature [K]", to_string(Tx(0,0,0)));
            Info::WriteTimeStep(RTC, message);

            TextOutput::WriteMultipleValues({"Temperature","Volume Fraction"}, { Tx(0,0,0),
                    1.0 - Phi.FieldsStatistics[0].Volume/(OPSettings.TotalNx*OPSettings.Ny*OPSettings.Nz)}, MartensiteVolumeFractionFile, RTC.tStep);
            dG.PrintDiagnostics();
            Phi.PrintVolumeFractions();
        }

        if (RTC.WriteVTK())
        {
            Phi.WriteVTK(RTC.tStep, OPSettings);
            EP.WriteCauchyStressesVTK(RTC.tStep, OPSettings);
            Phi.WriteDistortedVTK(RTC.tStep, OPSettings, EP);
        }
        WriteTemp_VMStress(EP, Tx, RTC.tStep, RTC.SimulationTime, TempFile);

        if (RTC.WriteRawData())
        {
            Phi.Write(RTC.tStep);
            EP.Write(RTC.tStep);
            Nuc.Write(RTC.tStep);
            Tx.Write(RTC.tStep);
        }
    }
#ifdef MPI_PARALLEL
    }
    MPI_Finalize ();
#endif
    return 0;
}

void WriteTemp_VMStress(ElasticProperties& EP, Temperature& Tx, int tStep, double RealTime, string FileName)
{
    double avgStress = 0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.Stresses,0,reduction (+:avgStress))
    {
        avgStress += EP.CauchyStress(i,j,k).Mises();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    double tmpStress = avgStress;
    MPI_Allreduce(&tmpStress, &avgStress, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    avgStress /= (EP.TotalNx*EP.TotalNy*EP.TotalNz);
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
    {
#endif
    ofstream output_file;
    if (tStep == 0)
    {
        output_file.open(FileName.c_str(), ios::out);
        output_file << "time" << "\t\t"
                    << "Temperature" << "\t\t"
                    << "VonMises";
        output_file << endl;
    }
    output_file.open(FileName.c_str(), ios::app);
    output_file << RealTime << "\t\t"
                << Tx(0,0,0) << "\t\t"
                << avgStress
                << endl;
    output_file.close();
#ifdef MPI_PARALLEL
    }
#endif
}

