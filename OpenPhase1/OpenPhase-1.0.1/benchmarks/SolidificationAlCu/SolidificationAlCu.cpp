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
#include "BoundaryConditions.h"
#include "Initializations.h"
#include "Tools/TimeInfo.h"
#include "Tools/MicrostructureAnalysis.h"
#include "EquilibriumPartitionDiffusionBinary.h"

using namespace std;
using namespace openphase;

void WriteVector(std::vector<double>& Output, Settings& locStettings, int tStep);
void ReadVector(std::vector<double>& Input, Settings& locStettings, int tStep);
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
    BoundaryConditions                  BC_C(OPSettings);
    TimeInfo                            Timer(OPSettings, "Execution Time Statistics");

    BC_C.BCNZ = BoundaryConditionTypes::Fixed; // Composition boundary condition at far end Z boundary

    // 1D diffusion domain extension
    int NzEXT = 4*OPSettings.Nz;
    vector<double> Cz_1D(NzEXT+2);
    vector<double> delta_Cz_1D(NzEXT+2, 0.0);

    // End of 1D diffusion domain extension
    ofstream fit_out;
    ofstream tip_velocity;

    if(RTC.Restart)
    {
        cout << "Restart data being read!" << endl;
        cout << "Restart time step: " << RTC.tStart << endl;
        Phi.Read(BC, RTC.tStart);
        Cx.Read(BC_C, RTC.tStart);
        Tx.Read(BC, RTC.tStart);
        ReadVector(Cz_1D, OPSettings, RTC.tStart);
        DF.SetDiffusionCoefficients(Phi, Tx);
        for(int i = -OPSettings.dNx; i < OPSettings.Nx+OPSettings.dNx; i++)
        for(int j = -OPSettings.dNy; j < OPSettings.Ny+OPSettings.dNy; j++)
        {
            Cx.MoleFractionsTotal(i,j,OPSettings.Nz)({0}) = Cz_1D[1];
            Cx.MoleFractions(i,j,OPSettings.Nz)({0,0}) = Cz_1D[1];
        }
        cout << "Done reading restart parameters!" << endl;
    }
    else
    {
        ignore_result(system("rm -rf DendriteProperties"));
        ignore_result(system("mkdir DendriteProperties"));

        tip_velocity.open("tStep-tip_velocity-tip_temperature.dat", ios::out);
        tip_velocity << "tStep \t tip_velocity \t tip_temperature" << endl;
        tip_velocity.close();

        Initializations::Single(Phi, 0, BC, OPSettings);
        Phi.PlantGrainNucleus(1, 0, 0, 0);
        //int idx1 = Initializations::Sphere(Phi, 1, 5, 0, 0, 0, BC, OPSettings);
        //EulerAngles angle({0.0, 0.0, 0.0},ZXZ);
        //Phi.FieldsStatistics[idx1].Orientation = angle.getQuaternion();
        Cx.SetInitialMoleFractions(Phi);
        Tx.SetInitial(BC);
        DF.SetDiffusionCoefficients(Phi, Tx);

        // 1D diffusion domain extension
        double Cz_ave = 0.0;
        for(int i = 0; i < OPSettings.Nx; i++)
        for(int j = 0; j < OPSettings.Ny; j++)
        {
            Cz_ave += Cx.MoleFractionsTotal(i,j,OPSettings.Nz-1)({0});
        }
        Cz_ave /= double(OPSettings.Nx*OPSettings.Ny);
        for(int k = 0; k < NzEXT+2; k++)
        {
            Cz_1D[k] = Cz_ave;
        }
        // End of 1D diffusion domain extension
    }
    cout << "Initialization stage done!" << endl;

    int x1 = 0;
    int y1 = 0;
    int z1 = (OPSettings.Nz)/3 + OPSettings.iWidth;

    // Evaluating the dendrite tip and trunk radius
    double old_tip_pos = 0.0;
    double old_time = 0.0;
    // End of evaluating the dendrite tip and trunk radius

    //RTC.dt = DF.ReportMaximumTimeStep();
    cout << "Time step: " << RTC.dt << endl;

    cout << "Entering the Time Loop!!!" << endl;
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        Timer.SetStart();
        IP.Set(Phi);
        DF.CalculateInterfaceMobility(Phi, Cx, Tx, BC, IP);
        Timer.SetTimeStamp("Set IP");
        DO.CalculatePhaseFieldIncrements(Phi, IP);
        Timer.SetTimeStamp("Get PsiDot");
        DF.GetDrivingForce(Phi, Cx, Tx, dG);
        Timer.SetTimeStamp("Chemical Driving Force");
        dG.Average(Phi, BC);
        Timer.SetTimeStamp("Driving Force Average");
        if (RTC.WriteVTK()) dG.WriteVTK(RTC.tStep,OPSettings, 1, 0);
        Timer.SetTimeStamp("Driving Force WriteVTK");
        dG.MergePhaseFieldIncrements(Phi, IP);
        Timer.SetTimeStamp("Driving Force merge to Psi");
        Phi.NormalizeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("Psi Normalize");
        DF.Solve(Phi, Cx, Tx, BC_C, RTC.dt);
        // 1D diffusion domain extension
        {
            int Nx = OPSettings.Nx;
            int Ny = OPSettings.Ny;
            int Nz = OPSettings.Nz;
            double dt = RTC.dt;
            double dx_2 = 1.0/(OPSettings.dx*OPSettings.dx);
            Cz_1D[0] = 0.0;
            for(int i = 0; i < Nx; i++)
            for(int j = 0; j < Ny; j++)
            {
                Cz_1D[0] += Cx.MoleFractionsTotal(i,j,Nz-1)({0});
            }
            Cz_1D[0] /= double(Nx*Ny);
            for(int k = 1; k < NzEXT; k++)
            {
                delta_Cz_1D[k] = DF.DC(Nx-1,Ny-1,Nz-1)({0})*(Cz_1D[k-1] - 2.0*Cz_1D[k] + Cz_1D[k+1])*dx_2;
            }
            for(int k = 1; k < NzEXT; k++)
            {
                Cz_1D[k] += delta_Cz_1D[k]*dt;
            }
            for(int i = -OPSettings.dNx; i < Nx+OPSettings.dNx; i++)
            for(int j = -OPSettings.dNy; j < Ny+OPSettings.dNy; j++)
            {
                Cx.MoleFractionsTotal(i,j,OPSettings.Nz)({0}) = Cz_1D[1];
                Cx.MoleFractions(i,j,OPSettings.Nz)({0,0}) = Cz_1D[1];
            }
        }
        // End of 1D diffusion domain extension
        Timer.SetTimeStamp("Solve diffusion");
        Tx.Set(BC, Phi, RTC.dt);
        Timer.SetTimeStamp("Set temperature");
        Phi.MergeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("Merge Phase Fields");
        if (Phi.Fields(0, 0, OPSettings.Nz/3)[0] <= 0.4)
        {
            Phi.MoveFrame(0,0,1, BC);
            Tx.MoveFrame(0,0,1, BC);
            Cx.MoveFrame(0,0,1, BC_C);
            for(int k = 0; k < NzEXT+1; k++)
            {
                Cz_1D[k] = Cz_1D[k+1];
            }
            old_tip_pos -= 1;
        }
        // Evaluating the dendrite tip and trunk radius
        {
            double dx        = OPSettings.dx;
            double dt        = RTC.dt;
            double delta     = 1.0e-4;
            int    tip_pos   = old_tip_pos;
            double tip_vel   = 0.0;
            bool   tip_found = false;
            for(int z = (OPSettings.Nz)/4; z < OPSettings.Nz - 1; z++)
            if(Phi.Fields(0, 0, z)[0] > 0.5 - delta and
               Phi.Fields(0, 0, z)[0] < 0.5 + delta and
               z != old_tip_pos)
            {
                tip_pos = z;
                tip_vel = (tip_pos - old_tip_pos)*dx/(RTC.tStep*dt - old_time);
                old_tip_pos = z;
                old_time = RTC.tStep*dt;
                tip_found = true;
                tip_velocity.open("tStep-tip_velocity-tip_temperature.dat", ios::app);
                tip_velocity << RTC.tStep << ", " << tip_vel << ", " << Tx.at(0,0,z) << endl;
                tip_velocity.close();
                break;
            }
            if (tip_found)
            {
                z1 = tip_pos;
                vector<double> parabola (tip_pos, 0.0);
                vector<double>   radius (tip_pos, 0.0);
                vector<double>   radius2(tip_pos, 0.0);
                vector<double> tip_position(tip_pos, 0.0);
                for(int z = tip_pos-1; z >= 0; z--)
                {
                    parabola[z] = MicrostructureAnalysis::FindValuePosition(Phi,1,0.5,(dVector3){0,0,(double)z},(dVector3){1,0,0})[0];
                }
                for(int z = tip_pos-1; z > 0; z--)
                {
                    double xx1 = -parabola[z] - 1;
                    double xx2 = 0.0;
                    double xx3 = parabola[z];
                    double zz1 = z;
                    double zz2 = tip_pos;
                    double zz3 = z;
                    radius[z] = dx/(2.0*((-(xx2*zz1) + xx3*zz1 + xx1*zz2 - xx3*zz2 - xx1*zz3 + xx2*zz3)/
                                        ((xx1 - xx2)*(xx2 - xx3)*(xx1 - xx3))));
                    xx1 = parabola[z-1];
                    xx3 = parabola[z];
                    zz1 = z-1;
                    zz3 = z;
                    radius2[z] = dx/(2.0*(zz3 - zz1)/(pow(xx1,2) - pow(xx3,2)));
                    tip_position[z] = -((pow(xx3,2)*zz1 - pow(xx1,2)*zz3)/(pow(xx1,2) - pow(xx3,2)));
                }
                fstream parabola_out(string("DendriteProperties/DendriteContour_") + to_string(RTC.tStep)+ string(".txt"), ios::out);
                parabola_out << tip_pos << " " << 0 << endl;
                for(int z = tip_pos-1; z >= 0; z--)
                {
                    parabola_out << z << " " << parabola[z] << endl;
                }
                parabola_out.close();
                fstream curvature_radius_out(string("DendriteProperties/DendriteTipRadius_") + to_string(RTC.tStep)+ string(".txt"), ios::out);
                curvature_radius_out << "TipPosition" << " " << "TipRadius"<< " " << "TipPosition2" << " "<< "TipRadius2" << " " << "dTipRadius_dz" << " " << "d2TipRadius_dz2" << " " << "d2TipRadius_dz2-5point" << endl;
                curvature_radius_out << tip_pos << " " << 0 << " " << tip_pos << " " << 0 << " " << 0 << " " << 0 << endl;
                for(int z = tip_pos-1; z > 1; z--)
                {
                    curvature_radius_out << z << " " << radius[z]<< " " << tip_position[z] <<  " " << radius2[z] << " " << radius[z-1] - radius[z+1] << " " << radius[z-1]+radius[z+1]-2.0*radius[z] << endl;
                }
                curvature_radius_out.close();
            }
        }
        // End of evaluating the dendrite tip and trunk curvature
        Timer.SetTimeStamp("Moving Frame");
        //  Output to file
        if (RTC.WriteVTK())
        {
            // Write data in VTK format
            Phi.WriteVTK(RTC.tStep,OPSettings);
            Cx.WriteVTK(RTC.tStep,OPSettings);
            Tx.WriteVTK(RTC.tStep,OPSettings);
            Cx.WriteStatistics(RTC.tStep, RTC.dt);
            IP.WriteVTK(Phi,RTC.tStep);
        }
        if (RTC.WriteRawData())
        {
            // Write raw data
            Phi.Write(RTC.tStep);
            Cx.Write(RTC.tStep);
            Tx.Write(RTC.tStep);
            WriteVector(Cz_1D, OPSettings, RTC.tStep);
        }
        //Timer.SetTimeStamp("File Output");
        //  Output to screen
        if(RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message  = Info::GetStandard("Interface energy density", to_string(I_En));
            Info::WriteTimeStep(RTC, message);
            // Statistics
            Phi.PrintPointStatistics(x1,y1,z1);
            Cx.PrintPointStatistics(x1,y1,z1);
            Tx.PrintPointStatistics(x1,y1,z1);
            dG.PrintDiagnostics();
            Phi.PrintPFVolumes();
            Timer.PrintWallClockSummary();
        }
    } //end time loop
    return 0;
}
void WriteVector(std::vector<double>& Output, Settings& locStettings, int tStep)
{
    string FileName = UserInterface::MakeFileName(locStettings.RawDataDir, "Composition1D_" , tStep, ".dat");
    ofstream out(FileName.c_str(), ios::out | ios::binary);
    for(size_t n = 0; n < Output.size(); n++)
    {
        double out_value = Output[n];
        out.write(reinterpret_cast<char*>(&out_value), sizeof(double));
    }
};
void ReadVector(std::vector<double>& Input, Settings& locStettings, int tStep)
{
    string FileName = UserInterface::MakeFileName(locStettings.RawDataDir, "Composition1D_" , tStep, ".dat");
    ifstream inp(FileName.c_str(), ios::in | ios::binary);
    for(size_t n = 0; n < Input.size(); n++)
    {
        double in_value;
        inp.read(reinterpret_cast<char*>(&in_value), sizeof(double));
        Input[n] = in_value;
    }
};
