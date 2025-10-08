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
#include "Settings.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "Composition.h"
#include "Temperature.h"
#include "Initializations.h"
#include "InterfaceProperties.h"
#include "Mechanics/ElasticProperties.h"
#include "Mechanics/ElasticitySolverSpectral.h"
#include "EquilibriumPartitionDiffusionBinary.h"
#include "TextOutput.h"
#include "BoundaryConditions.h"
#include "RunTimeControl.h"
#include "Tools/TimeInfo.h"
#include "Tools.h"
#include "VTK.h"

using namespace std;
using namespace openphase;

enum class Structure{Matrix, Precipitate, Multiple};

class LocSuperalloys
{
public:
   struct FileName
   {
       string StressesFile;
       string StrainsFile;
   } Files;
   Structure Microstructure;
    LocSuperalloys(){};
    LocSuperalloys(Settings& OPSettings)
    {
        ReadInput(OPSettings);
    };
    void ReadInput(Settings& OPSettings);
    void CreateMicrostructure(PhaseField& Phi, BoundaryConditions& BC, Settings& OPSettings, Composition& Cx);
    std::string TextDir;
};

void WriteEnergy(double EE, double IE, RunTimeControl& RTC, double time, string OutFile);

/****************************************************** <<< The Main >>> ****************************************************************/
int main(int argc, char *argv[])
{
#ifdef MPI_PARALLEL
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);
    {
#endif
    string InputFile = "ProjectInput.opi";

    Settings                        OPSettings(InputFile);
    RunTimeControl                  RTC(OPSettings,InputFile);
    Composition                     Cx(OPSettings,InputFile);
    BoundaryConditions              BC(OPSettings,InputFile);
    PhaseField                      Phi(OPSettings);
    DoubleObstacle                  DO(OPSettings);
    InterfaceProperties             IP(OPSettings,InputFile);
    EquilibriumPartitionDiffusionBinary DF(OPSettings);
    Temperature                     Tx(OPSettings);
    DrivingForce                    dG(OPSettings);
    ElasticProperties               EP(OPSettings);
    ElasticitySolverSpectral        ES(OPSettings);
    LocSuperalloys                  LS(OPSettings);
    TimeInfo                        Timer(OPSettings, "Execution Time Statistics");

    Info::WriteLineInsert("Initialization stage done!");

    if(RTC.Restart)
    {
        Info::WriteLineInsert("Restart data being read! ");
        Phi.Read(BC, RTC.tStart);
        Cx.Read(BC, RTC.tStart);
        Tx.Read(BC, RTC.tStart);
        EP.Read(BC,RTC.tStart);
        EP.SetGrainsProperties(Phi);
        EP.SetEffectiveElasticConstants(Phi);
        EP.SetEffectiveTransformationStretches(Phi);
        ES.Solve(EP, BC, RTC.dt);
        RTC.tStart -= 1;
        Info::WriteLineInsert("Done","+");
    }
    else
    {
        Tx.SetInitial(BC);
        RTC.tStart = -1;
        LS.CreateMicrostructure(Phi,BC,OPSettings,Cx);
        EP.SetGrainsProperties(Phi);
        Cx.SetInitialMoleFractions(Phi);
    }
    for(size_t n = 0; n < Phi.FieldsStatistics.size(); n++)
    {
        double Q1 = 0.0 * Pi/180.0;
        double Q2 = 0.0 * Pi/180.0;
        double Q3 = 0.0 * Pi/180.0;
        EulerAngles locAngles({Q1, Q2, Q3}, XYZ);
        Phi.FieldsStatistics[n].Orientation = locAngles.getQuaternion();
    }
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Info::WriteLineInsert("Entering the Time Loop!!!");
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        double dt = RTC.dt;
        dG.Clear();
Timer.SetStart();
        IP.Set(Phi);
Timer.SetTimeStamp("Set Interface Properties");
        DO.CalculatePhaseFieldIncrements(Phi, IP);
Timer.SetTimeStamp("Calculate PhaseField Increments");
        DF.GetDrivingForce(Phi, Cx, Tx, dG);
Timer.SetTimeStamp("Chemical Driving Force");
        EP.SetEffectiveElasticConstants(Phi);
        EP.SetEffectiveTransformationStretches(Phi);
Timer.SetTimeStamp("Elastic Properties I");
        EP.CalculateDrivingForce(Phi, dG);
Timer.SetTimeStamp("Elastic Driving Force");
        dG.Average(Phi, BC);
Timer.SetTimeStamp("Driving Force Average");
        dG.MergePhaseFieldIncrements(Phi, IP);
Timer.SetTimeStamp("Merge PhaseField Increments");
        Phi.NormalizeIncrements(BC, dt);
Timer.SetTimeStamp("Normalize PhaseField Increments");
        EP.SetEffectiveElasticConstants(Phi);
        EP.SetEffectiveTransformationStretches(Phi);
Timer.SetTimeStamp("Elastic Properties II");
        ES.Solve(EP, BC, dt);
Timer.SetTimeStamp("Elastic Solver");
        DF.SetDiffusionCoefficients(Phi, Tx);
        DF.Solve(Phi, Cx, Tx, BC, RTC.dt, false);
Timer.SetTimeStamp("Diffusion Solver");
        Phi.MergeIncrements(BC, dt);
Timer.SetTimeStamp("MergeIncrements");
        Tx.Set(BC,Phi,RTC.dt);
Timer.SetTimeStamp("Set Temperature");
        if (RTC.WriteRawData())
        {
            // Write rawdata
            Phi.Write(RTC.tStep);
            Cx.Write(RTC.tStep);
            Tx.Write(RTC.tStep);
            EP.Write(RTC.tStep);
        }
        if (RTC.WriteVTK())//  Output to file
        {
            Phi.WriteVTK(RTC.tStep,OPSettings,false,3);
            Cx.WriteVTK(RTC.tStep,OPSettings);
            EP.WriteStressesVTK(RTC.tStep,OPSettings,3);
            Phi.WriteAverageVolume(RTC.tStep, 1);
            Cx.WriteStatistics(RTC.tStep, dt);
        }
Timer.SetTimeStamp("File Output");
        if (RTC.WriteToScreen())
        {
            double E_En = EP.AverageEnergyDensity();
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            double T_En = E_En + I_En;


            std::string message  = Info::GetStandard("Elastic energy density", to_string(E_En));
                        message += Info::GetStandard("Interface energy density", to_string(I_En));
                        message += Info::GetStandard("Total energy density", to_string(T_En));
            Info::WriteTimeStep(RTC, message);

            Phi.PrintVolumeFractions();
            Info::WriteLine();
            dG.PrintDiagnostics();

            Info::WriteLine("-");
            Info::WriteLine("|");
            Info::WriteLine("-");
        } // print if statement

//--------------------------------------------|| Write to File ||-------------------------------------------
            TextOutput::AverageStress(EP, LS.Files.StressesFile, RTC.tStep);
            TextOutput::AverageStrain(EP, LS.Files.StrainsFile, RTC.tStep);
//----------------------------------------------------------------------------------------------------------
    } //end time loop

    simulation_end();
#ifdef MPI_PARALLEL
    }
    MPI_Finalize ();
#endif
    return 0;
}
void LocSuperalloys::ReadInput(Settings& OPSettings)
{
    string thisclassname = "LocSuperalloys";
    Info::WriteBlankLine();
    Info::WriteLineInsert(thisclassname);
    fstream inpF(DefaultInputFileName.c_str(), ios::in);
    if (!inpF)
    {
        Info::WriteExit("File \"" + DefaultInputFileName + "\" could not be opened",
                thisclassname, "ReadInput");
        exit(1);
    };
    stringstream inp;
    inp << inpF.rdbuf();
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
    string tmp         = UserInterface::ReadParameterK(inp,moduleLocation, "Microstructure");
    TextDir = OPSettings.TextDir;

    if(tmp == "MATRIX")      Microstructure = Structure::Matrix;
    if(tmp == "PRECIPITATE") Microstructure = Structure::Precipitate;
    if(tmp == "MULTIPLE")    Microstructure = Structure::Multiple;

    Files.StressesFile      = OPSettings.TextDir + "StressesFile.txt";
    Files.StrainsFile       = OPSettings.TextDir + "StrainsFile.txt";

    remove(Files.StressesFile.c_str());
    remove(Files.StrainsFile.c_str());

}
void LocSuperalloys::CreateMicrostructure(PhaseField& Phi, BoundaryConditions& BC, Settings& OPSettings, Composition& Cx)
{
    switch (Microstructure)
    {
        case Structure::Matrix:
        {
            Initializations::Single(Phi, 0, BC, OPSettings);
            break;
        }
        case Structure::Precipitate:
        {
            Initializations::Single(Phi, 0, BC, OPSettings);
            Initializations::Sphere(Phi, 1, 5, (OPSettings.Nx+1)/2, (OPSettings.Ny+1)/2, (OPSettings.Nz+1)/2, BC, OPSettings);
            break;
        }
        case Structure::Multiple:
        {
            Initializations::Single(Phi, 0, BC, OPSettings);
            Info::WriteLineInsert("Planting particles on a regular grid with random offset! ");
            switch (Phi.Resolution)
            {
                case Resolutions::Single:
                {
                    Initializations::QuasiRandomNuclei(Phi, OPSettings,1,8,-1,5);
                    break;
                }
                case Resolutions::Double:
                {
                    Initializations::QuasiRandomNuclei(Phi, OPSettings,1,6,-1,5);
                    break;
                }
            }
            Info::WriteLineInsert("Done");
            break;
        }
    }
}
void WriteEnergy(double EE, double IE, RunTimeControl& RTC,  double time, string OutFile)
{
    ofstream output_file;
    stringstream FileName;
    FileName<<OutFile.c_str()<<".txt";
    if (time == 0 or (RTC.Restart and time == RTC.tStart+1))
    {
        remove(FileName.str().c_str());
        output_file.open(FileName.str().c_str(), ios::out);
        output_file << setw(12) << left << "time"
                    << setw(20) << left << "Elastic Energy"
                    << setw(20) << left << "Interface Energy"
                    << setw(20) << left << "Total Energy";
        output_file<< endl;
        output_file.close();
    }
    output_file.open(FileName.str().c_str(), ios::app);
    output_file << setw(12) << left << time
                << setw(20) << left << EE
                << setw(20) << left << IE
                << setw(20) << left << IE + EE;
    output_file<< endl;
    output_file.close();
}
