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
#include "Initializations.h"
#include "Mechanics/ElasticProperties.h"
#include "Mechanics/ElasticitySolverSpectral.h"
#include "BoundaryConditions.h"
#include "Tools/TimeInfo.h"
#include "PhysicalConstants.h"

using namespace std;
using namespace openphase;

void SetDiffusionCoefficientsCx(PhaseField& Phase, Temperature& Tx,
                                  EquilibriumPartitionDiffusionBinary& DF,
                                  Composition& Elements)                        ///<  Sets concentration dependent diffusion coefficients in each point
{
    size_t Comp = 0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DF.DC,DF.DC.Bcells(),)
    {
        double invRT  = 1.0/(PhysicalConstants::R*Tx(i,j,k));
        double T = Tx(i,j,k);
        for (size_t n = 0; n < Phase.Nphases; ++n)
        {
            if (n == 0)
            {
                double locCx = Elements.MoleFractions(i,j,k)({0,Comp});
                DF.DC(i,j,k)({n}) = 4.45e-7*(1.0 + (locCx/(100.0 - locCx)) *
                                            (1.0 - (locCx/(100.0 - locCx))) * (8339.9/T)) *
                exp(-((1.0/T) - 2.221e-4)*(17767.0 - (locCx/(100.0 - locCx)) * 26436.0));
            }
            else
            {
                DF.DC(i, j, k)({n}) = DF.DC0[n] * exp(-DF.AE[n] * invRT);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ScaleDiffusionCoefficientsStress(PhaseField& Phase, Temperature& Tx,
                                  EquilibriumPartitionDiffusionBinary& DF,
                                  ElasticProperties& EP, double A)              ///<  Scales diffusion coefficient according to its stress dependence
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DF.DC,DF.DC.Bcells(),)
    {
        double invRT  = 1.0/(PhysicalConstants::R*Tx(i,j,k));

        for (size_t n = 0; n < Phase.Nphases; ++n)
        {
            if (n == 0)
            {
                DF.DC(i,j,k)({n}) *=

                (1.0 - (A*EP.Stresses(i,j,k).Pressure()/150e6))*
                exp(-30000.0*100.0*EP.Stresses(i,j,k).Pressure()*invRT);

                // sets minimum diffusion coefficient
                DF.DC(i,j,k)({n}) = max(3.3e-14, DF.DC(i,j,k)({n}));
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ScaleDiffusionCoefficientsInterface(PhaseField& Phase,
                                         EquilibriumPartitionDiffusionBinary& DF,
                                         double DInterface)                     ///<  Scales diffusion coefficients in the interface by the factor "DInterface"
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DF.DC,DF.DC.Bcells(),)
    {
        if (Phase.Interface(i,j,k))
        for (size_t n = 0; n < Phase.Nphases; n++)
        {
            DF.DC(i,j,k)({n}) *= DInterface;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void SetEffectiveShearFreeElasticConstants(PhaseField& Phase, ElasticProperties& EP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.EffectiveElasticConstants,EP.EffectiveElasticConstants.Bcells(),)
    {
        if(Phase.Interface(i,j,k))
        {
            EP.EffectiveElasticConstants(i,j,k).set_to_zero();
            dMatrix6x6 TempCompliances;
            TempCompliances.set_to_zero();
            bool gamma_present = false;
            bool alpha_present = false;
            bool cementite_present = false;
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                TempCompliances += (EP.Compliances[alpha->index] * alpha->value);
                if(Phase.FieldsStatistics[alpha->index].Phase == 0) gamma_present = true;
                if(Phase.FieldsStatistics[alpha->index].Phase == 1) alpha_present = true;
                if(Phase.FieldsStatistics[alpha->index].Phase == 2) cementite_present = true;
            }
            if(gamma_present and alpha_present)

            {
                TempCompliances(3, 3) = 8.6e-11;
                TempCompliances(4, 4) = 8.6e-11;
                TempCompliances(5, 5) = 8.6e-11;
                TempCompliances(0, 1) = -TempCompliances(3,3)/2.0 + TempCompliances(0,0);
                TempCompliances(0, 2) = -TempCompliances(3,3)/2.0 + TempCompliances(0,0);
                TempCompliances(1, 2) = -TempCompliances(3,3)/2.0 + TempCompliances(0,0);
                TempCompliances(1, 0) = -TempCompliances(3,3)/2.0 + TempCompliances(0,0);
                TempCompliances(2, 0) = -TempCompliances(3,3)/2.0 + TempCompliances(0,0);
                TempCompliances(2, 1) = -TempCompliances(3,3)/2.0 + TempCompliances(0,0);
            }
            if(gamma_present and cementite_present)
            {
                TempCompliances(3, 3) = 10.42e-11;
                TempCompliances(4, 4) = 10.42e-11;
                TempCompliances(5, 5) = 10.42e-11;
                TempCompliances(0, 1) = -TempCompliances(3, 3)/2.0 + TempCompliances(0, 0);
                TempCompliances(0, 2) = -TempCompliances(3, 3)/2.0 + TempCompliances(0, 0);
                TempCompliances(1, 2) = -TempCompliances(3, 3)/2.0 + TempCompliances(0, 0);
                TempCompliances(1, 0) = -TempCompliances(3, 3)/2.0 + TempCompliances(0, 0);
                TempCompliances(2, 0) = -TempCompliances(3, 3)/2.0 + TempCompliances(0, 0);
                TempCompliances(2, 1) = -TempCompliances(3, 3)/2.0 + TempCompliances(0, 0);
            }
            EP.EffectiveElasticConstants(i,j,k) = TempCompliances.inverted();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
{
    Settings                        OPSettings;
    OPSettings.ReadInput();

    RunTimeControl                      RTC(OPSettings);
    PhaseField                          Phi(OPSettings);
    DoubleObstacle                      DO(OPSettings);
    InterfaceProperties                 IP(OPSettings);
    EquilibriumPartitionDiffusionBinary DF(OPSettings);
    Composition                         Cx(OPSettings);
    Temperature                         Tx(OPSettings);
    DrivingForce                        dG(OPSettings);
    ElasticProperties                   EP(OPSettings);
    BoundaryConditions                  BC(OPSettings);
    ElasticitySolverSpectral            ES(OPSettings);
    TimeInfo                            Timer(OPSettings, "Execution Time Statistics");

    cout << "Initialization stage! Done!" << endl;

    if(RTC.Restart)
    {
        cout << "Restart data being read!";
        Phi.Read(BC, RTC.tStart);
        Cx.Read(BC, RTC.tStart);
        Tx.Read(BC, RTC.tStart);
        cout << " Done!" << endl;
    }
    else
    {
        Initializations::Single(Phi, 0, BC, OPSettings);
        Initializations::Zlayer(Phi, 1, OPSettings.Nz/2, OPSettings.Nx/10, BC, OPSettings);
        Initializations::Sphere(Phi, 2, OPSettings.Nx/12,
                                        OPSettings.Nx/2,
                                        OPSettings.Ny/2,
                                        OPSettings.Nz/2, BC, OPSettings);
        Cx.SetInitialMoleFractions(Phi);
        Tx.SetInitial(BC);
        EP.SetGrainsProperties(Phi);
    }

    int x1 = (OPSettings.Nx+1)/2;
    int y1 = (OPSettings.Ny+1)/2;
    int z1 = 13;

    double tStepOld = RTC.tStart;
    double AusteniteVolumeOld = Phi.FieldsStatistics[0].Volume;

    ofstream GrowthVelocity("GrowthVelocity.dat", ios::out);
    GrowthVelocity << "time(s) \t Velocity(m/s)" << endl;
    GrowthVelocity.close();

    cout << "Entering the Time Loop!!!" << endl;

    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        Timer.SetStart();
        IP.Set(Phi, Tx);
        Timer.SetTimeStamp("Set Interface Properties");
        DO.CalculatePhaseFieldIncrements(Phi, IP);
        Timer.SetTimeStamp("Get PhiDot");
        DF.GetDrivingForce(Phi, Cx, Tx, dG);
        Timer.SetTimeStamp("Chemical Driving Force");
        EP.SetEffectiveElasticConstants(Phi);
        EP.SetEffectiveTransformationStretches(Phi, Cx);
        Timer.SetTimeStamp("Elastic Properties");
        EP.CalculateChemicalPotentialContribution(Phi, DF);
        Timer.SetTimeStamp("Elastic Chemical Potential");
        EP.CalculateDrivingForce(Phi, dG);
        Timer.SetTimeStamp("Elastic Driving Force");
        if(!(RTC.tStep % 10))
        {
            ES.Solve(EP, BC, RTC.dt);
        }
        Timer.SetTimeStamp("Elastic Solver");
        dG.Average(Phi, BC);
        Timer.SetTimeStamp("Driving Force Average");
        dG.MergePhaseFieldIncrements(Phi, IP);
        Phi.NormalizeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("Merge Velocities");

        //Setting up modified diffusion coefficient and elastic constants

/*1a*/  DF.SetDiffusionCoefficients(Phi, Tx);
/*1b*/  //SetDiffusionCoefficientsCx(Phi, Tx, DF, Cx);
/*2*/   //ScaleDiffusionCoefficientsInterface(Phi, DF, 5);
/*3*/   //ScaleDiffusionCoefficientsStress(Phi, Tx, DF, EP, 10);
/*4*/   //SetEffectiveShearFreeElasticConstants(Phi, EP);

        // End setting up modified diffusion coefficient and elastic constants

        DF.Solve(Phi, Cx, Tx, BC, RTC.dt, false);
        Timer.SetTimeStamp("Diffusion Solver");
        Phi.MergeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("Merge PhaseFields");

        if (Phi.Fields((OPSettings.Nx+1)/2, (OPSettings.Ny+1)/2, (OPSettings.Nz+1)/3)[0] <= 0.4)
        {
            double deltaV = Phi.FieldsStatistics[0].Volume - AusteniteVolumeOld;
            double deltaL = deltaV*OPSettings.dx/(OPSettings.Nx*OPSettings.Ny);

            double vel = -0.5*deltaL/((RTC.tStep - tStepOld)*RTC.dt);

            GrowthVelocity.open("GrowthVelocity.dat", ios::app);
            GrowthVelocity << RTC.tStep*RTC.dt << "\t" << vel << endl;
            GrowthVelocity.close();

            Phi.ConsumePlane(0,0,1, (OPSettings.Nx+1)/2, (OPSettings.Ny+1)/2, (OPSettings.Nz+1)/2, BC);
            Tx.ConsumePlane(0,0,1, (OPSettings.Nx+1)/2, (OPSettings.Ny+1)/2, (OPSettings.Nz+1)/2, BC);
            Cx.ConsumePlane(0,0,1, (OPSettings.Nx+1)/2, (OPSettings.Ny+1)/2, (OPSettings.Nz+1)/2, BC);

            AusteniteVolumeOld = Phi.FieldsStatistics[0].Volume;
            tStepOld = RTC.tStep;
        }

        Tx.Set(BC, Phi, RTC.dt);
        Timer.SetTimeStamp("Set Temperature");
        if (RTC.WriteVTK())
        {
            Phi.WriteVTK(RTC.tStep, OPSettings);
            Cx.WriteVTK(RTC.tStep, OPSettings);
            Cx.WriteStatistics(RTC.tStep, RTC.dt);
            EP.WriteStressesVTK(RTC.tStep, OPSettings);
            EP.WriteTotalStrainsVTK(RTC.tStep, OPSettings);
        }
        if (RTC.WriteRawData())
        {
            Phi.Write(RTC.tStep);
            Cx.Write(RTC.tStep);
            Tx.Write(RTC.tStep);
        }
        Timer.SetTimeStamp("File Output");
        if (RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            double E_En = EP.AverageEnergyDensity();
            double T_En = E_En + I_En;
            std::string message  = Info::GetStandard("Interface energy density", to_string(I_En));
                        message += Info::GetStandard("Elastic energy density", to_string(E_En));
                        message += Info::GetStandard("Total energy density", to_string(T_En));

            Info::WriteTimeStep(RTC, message);

            Phi.PrintPointStatistics(x1,y1,z1);

            //Cx.PrintPointStatistics(x1,y1,z1);
            //Tx.PrintPointStatistics(x1,y1,z1);
            //EP.PrintPointStatistics(x1,y1,z1);
            dG.PrintDiagnostics();
            Phi.PrintPFVolumes();
            //Timer.PrintWallClockSummary();
        }
    } //end time loop

   return 0;
}
