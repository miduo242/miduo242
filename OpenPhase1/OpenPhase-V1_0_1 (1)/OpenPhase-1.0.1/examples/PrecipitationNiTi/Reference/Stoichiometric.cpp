#include "Types.h"
#include "Settings.h"
#include "InterfaceEnergy.h"
#include "InterfaceMobility.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "InterfaceField.h"
#include "DrivingForce.h"
#include "Composition.h"
#include "Temperature.h"
#include "EquilibriumPartitionDiffusion.h"
#include "Initializations.h"
#include "ElasticProperties.h"
#include "ElasticitySteinbach.h"
#include "SpectralElasticSolver.h"
#include "BoundaryConditions.h"
#include "Info.h"
#include "TimeInfo.h"

using namespace std;
using namespace openphase;
/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
{
    Settings                        OPSettings;
    PhaseField                      Phi;
    InterfaceField                  Psi;
    DoubleObstacle                  DO;
    IdenticalInterfaceMobility      Mu;
    IdenticalInterfaceEnergy        Sigma;
    EquilibriumPartitionDiffusion   DF;
    Composition                     Cx;
    Temperature                     Tx;
    DrivingForce                    dG;
    ElasticitySteinbach             EM;
    ElasticProperties               EP;
    SpectralElasticSolver           ES;  
    BoundaryConditions              BC;
    TimeInfo                        Timer;

    OPSettings.ReadInput(ProjectFileName);

    Mu.Initialize(OPSettings);
    Sigma.Initialize(OPSettings);
    Phi.Initialize(OPSettings);
    Psi.Initialize(OPSettings);
    dG.Initialize(OPSettings);
    DO.Initialize(OPSettings);
    DF.Initialize(OPSettings);
    Tx.Initialize(OPSettings);
    Cx.Initialize(OPSettings);
    EP.Initialize(OPSettings);
    ES.Initialize(OPSettings);
    EM.Initialize(OPSettings); 
    BC.Initialize(OPSettings);
    Timer.Initialize(OPSettings, "Execution Time Statistics");

    DF.ReadInput(ProjectDir + DiffusivityInputFileName);
    Cx.ReadInput(ProjectDir + DiffusivityInputFileName);
    Tx.ReadInput(ProjectDir + DiffusivityInputFileName);
    EP.ReadInput(ProjectDir + ElasticityInputFileName);
    BC.ReadInput(ProjectDir + BoundaryConditionsInputFileName);

    cout << "Initialization stage! Done!" << endl;

    OPSettings.tStart = -1;

    int x1 = (OPSettings.Nx+1)/2;
    int y1 = (OPSettings.Ny+1)/2;
    int z1 = (OPSettings.Nz+1)/2;
    double sim_time = 0.0;

    Initializations::Single(Phi, 0, BC, OPSettings);
    Phi.PlantGrainNucleus(OPSettings.iPhiIndex, x1, y1, z1, BC);
    
    DF.SetInitialComposition(Phi, Cx);

    EP.SetGrainsProperties(Phi);

    Tx.SetInitial(BC);
    
    fstream output_file;
    output_file.open("ParticleSize.txt", fstream::out);

    output_file << "tStep" << "\t" << "sim_time" << "\t"
                << "part_diameter" << "\t" << "part_height" << endl;
    output_file.close();

    cout << "Entering the Time Loop!!!" << endl;
    for(int tStep = OPSettings.tStart+1; tStep < OPSettings.nSteps+1; tStep++)
    {
    Timer.SetStart();

        sim_time += OPSettings.dt;

        Sigma.Calculate(Phi);
        Mu.Calculate(Phi);
        double I_En = 0.0;

        if(!(tStep % 1))
        {
            EM.SetEffectiveElasticConstants(Phi, EP);
            EM.SetEffectiveEigenStrains(Phi, EP);
            /*int iterations = */ES.Solve(EP, BC, 1.0e-6, 10000.0, 100);
            //cout << "ES.Solve() Iterations: " << iterations << endl;
        }
    Timer.SetTimeStamp("Elastic Solver");
        if (!(tStep%OPSettings.tScreenWrite)) I_En = DO.Energy(Phi, Sigma);

        DO.GetPsiDot(Phi, Sigma, Mu, Psi);
    Timer.SetTimeStamp("Get PhiDot");
        EM.SetEffectiveElasticConstants(Phi, EP);
        EM.SetEffectiveEigenStrains(Phi, EP);
    Timer.SetTimeStamp("Elastic Properties");
        EM.GetDrivingForce(Phi, EP, dG);
        if (!(tStep%OPSettings.tFileWrite))//  Output to file
            dG.WriteVTK(tStep,0,3);
    Timer.SetTimeStamp("Elastic Driving Force");
        DF.GetDrivingForce(Phi, Cx, Tx, dG);
        if (!(tStep%OPSettings.tFileWrite))//  Output to file
            dG.WriteVTK(tStep+1,0,3);
    Timer.SetTimeStamp("Chemical Driving Force");
        dG.Average(Phi, BC);
    Timer.SetTimeStamp("Driving Force Average");
        dG.MergePsi(Phi, Sigma, Mu, Psi);
        Psi.Normalize(Phi, BC, OPSettings.dt);
    Timer.SetTimeStamp("Merge Velocities");
        DF.Solve(Phi, Psi, Cx, Tx, BC, OPSettings.dt);
    Timer.SetTimeStamp("Diffusion Solver");
        Psi.Merge(Phi, BC, OPSettings.dt);
    Timer.SetTimeStamp("Merge PhaseFields");
        Tx.Set(BC, OPSettings.dt);
    Timer.SetTimeStamp("Set Temperature");
        if (!(tStep%OPSettings.tRestartWrite))
        {
            // Write raw data
            Phi.Write(tStep);
            Cx.Write(tStep);
            Tx.Write(tStep);
        }
        if (!(tStep%OPSettings.tFileWrite))
        {
            // Write data in VTK format
            Phi.WriteVTK(tStep);
            Cx.WriteVTK(tStep);
            Cx.WriteStatistics(tStep, OPSettings.dt);
            EP.WriteStressesVTK(tStep);
            //dG.WriteVTK(tStep,0,3);
            EP.WriteStrainsVTK(tStep);
            //EP.WriteDisplacementsVTK(tStep);
        }
    Timer.SetTimeStamp("File Output");
        if (!(tStep%OPSettings.tScreenWrite))
        {
            time_t rawtime;
            time(&rawtime);
            //  Output to screen
            cout << "+++++++++++++++++++++++++++++++++++++++++\n"
                 << "Time Step:          " << tStep << "\n"
                 << "Wall Clock Time:    " << ctime(&rawtime) << "\n";
            double E_En = EM.Energy(EP);
            double T_En = E_En + I_En;
            cout << "=========================================\n"
                 << "Elastic Energy:     " << E_En << "\n"
                 << "Interface Energy:   " << I_En << "\n"
                 << "-----------------------------------------\n"
                 << "Total Energy:       " << T_En << "\n"
                 << "=========================================\n" << endl;

            Phi.PrintPointStatistics(x1,y1,z1);
            Cx.PrintPointStatistics(x1,y1,z1);
            Tx.PrintPointStatistics(x1,y1,z1);
            EP.PrintPointStatistics(x1,y1,z1);
            dG.PrintDiagnostics();

            double part_diameter = 0.;
            double part_height = 0.;
            // collect data over line and calculation of the particle's diameter and height:
            for(int i = 1; i < OPSettings.Nx+1; ++i)
            {
                part_diameter += Phi.Fields(i, OPSettings.Ny/2+1, OPSettings.Nz/2+1)[3];
                part_height   += Phi.Fields(OPSettings.Nx/2+1, OPSettings.Ny/2+1, i)[3];
            }

            fstream output_file;
            output_file.open("ParticleSize.txt", fstream::out | fstream::app);
            output_file << tStep << "\t" << sim_time << "\t"
                        << part_diameter << "\t" << part_height << endl;
            output_file.close();

            Timer.PrintWallClockSummary();
            Phi.PrintPFVolumes();
        }
    } //end time loop
   return 0;
}
