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
#include "FluidDynamics/InteractionSolidFluid.h"
#include "Velocities.h"
#include "AdvectionHR/AdvectionHR.h"
#include "BoundaryConditions.h"

using namespace std;
using namespace openphase;

void SetInletVelocityX0(FlowSolverLBM& FL, const BoundaryConditions& BC, const PhaseField& Phi);

/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
{
#ifdef MPI_PARALLEL

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);

#endif
    //feenableexcept(FE_DIVBYZERO);  // Devision by zero
    //feenableexcept(FE_INVALID);    // domain error occurred
    //feenableexcept(FE_OVERFLOW);   // the result was too large
    //feenableexcept(FE_UNDERFLOW);  // the result was too small

    string InputFile = "ProjectInput.opi";

    Settings                            OPSettings(InputFile);
    RunTimeControl                      RTC(OPSettings,InputFile);
    PhaseField                          Phi(OPSettings);
    BoundaryConditions                  BC(OPSettings,InputFile);
    Temperature                         Tx(OPSettings,InputFile);
    Composition                         Cx(OPSettings,InputFile);
    FlowSolverLBM                       FL(OPSettings, RTC.dt, InputFile);
    Velocities                          Vel(OPSettings);

    if(RTC.Restart)
    {
        //cout << "Restart data being read!";
        Phi.Read(BC, RTC.tStart);
        //Cx.Read(BC, RTC.tStart);
        //Tx.SetInitial(OPSettings, BC, Phi, 0);
        //Tx.Read(RTC.tStart);
        //cout << " Done!" << endl;
    }
    else
    {
        size_t idx0 = Initializations::Single(Phi, 0, BC, OPSettings);
        //vector<size_t> idx01 =  Initializations::TwoWalls( Phi,  0,  1,  2, BC,  OPSettings);

        double Nz=OPSettings.TotalNz-1.0;
        double Ny=OPSettings.TotalNy-1.0;

        double Xshift=0.0;
        double ND = Nz/6.0;
        double Xline1=OPSettings.TotalNx*(8.0/28.0)+Xshift;
        double Xline2=Xline1+sqrt(2.0)*ND;
        double Xline3=Xline2+sqrt(2.0)*ND;
        double Xline4=Xline3+sqrt(2.0)*ND;

        double z11=ND;
        double z21=Nz/2.0;
        double z31=Nz-ND;

        double y11=ND;
        double y21=Ny/2.0;
        double y31=Ny-ND;

        double z12=Nz/2.0-3.0*ND;
        double z22=Nz/2.0-ND;
        double z32=Nz/2.0+ND;
        double z42=Nz/2.0+3.0*ND;

        double y12=0.0;
        double y22=2.0*ND;
        double y32=4.0*ND;
        double y42=6.0*ND;

        size_t idx1  = Initializations::Sphere(Phi, 1, ND,  Xline1, y21, z21, BC,  OPSettings);
        size_t idx2  = Initializations::Sphere(Phi, 1, ND,  Xline1, y21, z11, BC,  OPSettings);
        size_t idx3  = Initializations::Sphere(Phi, 1, ND,  Xline1, y21, z31, BC,  OPSettings);

        size_t idx4  = Initializations::Sphere(Phi, 1, ND,  Xline1, y11, z21, BC,  OPSettings);
        size_t idx5  = Initializations::Sphere(Phi, 1, ND,  Xline1, y11, z11, BC,  OPSettings);
        size_t idx6  = Initializations::Sphere(Phi, 1, ND,  Xline1, y11, z31, BC,  OPSettings);

        size_t idx7  = Initializations::Sphere(Phi, 1, ND,  Xline1, y31, z21, BC,  OPSettings);
        size_t idx8  = Initializations::Sphere(Phi, 1, ND,  Xline1, y31, z11, BC,  OPSettings);
        size_t idx9  = Initializations::Sphere(Phi, 1, ND,  Xline1, y31, z31, BC,  OPSettings);



        size_t idx10  = Initializations::Sphere(Phi, 1, ND,  Xline2, y12, z12, BC,  OPSettings);
        size_t idx11  = Initializations::Sphere(Phi, 1, ND,  Xline2, y12, z22, BC,  OPSettings);
        size_t idx12  = Initializations::Sphere(Phi, 1, ND,  Xline2, y12, z32, BC,  OPSettings);
        size_t idx13  = Initializations::Sphere(Phi, 1, ND,  Xline2, y12, z42, BC,  OPSettings);

        size_t idx14  = Initializations::Sphere(Phi, 1, ND,  Xline2, y22, z12, BC,  OPSettings);
        size_t idx15  = Initializations::Sphere(Phi, 1, ND,  Xline2, y22, z22, BC,  OPSettings);
        size_t idx16  = Initializations::Sphere(Phi, 1, ND,  Xline2, y22, z32, BC,  OPSettings);
        size_t idx17  = Initializations::Sphere(Phi, 1, ND,  Xline2, y22, z42, BC,  OPSettings);

        size_t idx18  = Initializations::Sphere(Phi, 1, ND,  Xline2, y32, z12, BC,  OPSettings);
        size_t idx19  = Initializations::Sphere(Phi, 1, ND,  Xline2, y32, z22, BC,  OPSettings);
        size_t idx20  = Initializations::Sphere(Phi, 1, ND,  Xline2, y32, z32, BC,  OPSettings);
        size_t idx21  = Initializations::Sphere(Phi, 1, ND,  Xline2, y32, z42, BC,  OPSettings);

        size_t idx22  = Initializations::Sphere(Phi, 1, ND,  Xline2, y42, z12, BC,  OPSettings);
        size_t idx23  = Initializations::Sphere(Phi, 1, ND,  Xline2, y42, z22, BC,  OPSettings);
        size_t idx24  = Initializations::Sphere(Phi, 1, ND,  Xline2, y42, z32, BC,  OPSettings);
        size_t idx25  = Initializations::Sphere(Phi, 1, ND,  Xline2, y42, z42, BC,  OPSettings);



        size_t idx26  = Initializations::Sphere(Phi, 1, ND,  Xline3, y21, z21, BC,  OPSettings);
        size_t idx27  = Initializations::Sphere(Phi, 1, ND,  Xline3, y21, z11, BC,  OPSettings);
        size_t idx28  = Initializations::Sphere(Phi, 1, ND,  Xline3, y21, z31, BC,  OPSettings);

        size_t idx29  = Initializations::Sphere(Phi, 1, ND,  Xline3, y11, z21, BC,  OPSettings);
        size_t idx30  = Initializations::Sphere(Phi, 1, ND,  Xline3, y11, z11, BC,  OPSettings);
        size_t idx31  = Initializations::Sphere(Phi, 1, ND,  Xline3, y11, z31, BC,  OPSettings);

        size_t idx32  = Initializations::Sphere(Phi, 1, ND,  Xline3, y31, z21, BC,  OPSettings);
        size_t idx33  = Initializations::Sphere(Phi, 1, ND,  Xline3, y31, z11, BC,  OPSettings);
        size_t idx34  = Initializations::Sphere(Phi, 1, ND,  Xline3, y31, z31, BC,  OPSettings);


        size_t idx35  = Initializations::Sphere(Phi, 1, ND,  Xline4, y12, z12, BC,  OPSettings);
        size_t idx36  = Initializations::Sphere(Phi, 1, ND,  Xline4, y12, z22, BC,  OPSettings);
        size_t idx37  = Initializations::Sphere(Phi, 1, ND,  Xline4, y12, z32, BC,  OPSettings);
        size_t idx38  = Initializations::Sphere(Phi, 1, ND,  Xline4, y12, z42, BC,  OPSettings);

        size_t idx39  = Initializations::Sphere(Phi, 1, ND,  Xline4, y22, z12, BC,  OPSettings);
        size_t idx40  = Initializations::Sphere(Phi, 1, ND,  Xline4, y22, z22, BC,  OPSettings);
        size_t idx41  = Initializations::Sphere(Phi, 1, ND,  Xline4, y22, z32, BC,  OPSettings);
        size_t idx42  = Initializations::Sphere(Phi, 1, ND,  Xline4, y22, z42, BC,  OPSettings);

        size_t idx43  = Initializations::Sphere(Phi, 1, ND,  Xline4, y32, z12, BC,  OPSettings);
        size_t idx44  = Initializations::Sphere(Phi, 1, ND,  Xline4, y32, z22, BC,  OPSettings);
        size_t idx45  = Initializations::Sphere(Phi, 1, ND,  Xline4, y32, z32, BC,  OPSettings);
        size_t idx46  = Initializations::Sphere(Phi, 1, ND,  Xline4, y32, z42, BC,  OPSettings);

        size_t idx47  = Initializations::Sphere(Phi, 1, ND,  Xline4, y42, z12, BC,  OPSettings);
        size_t idx48  = Initializations::Sphere(Phi, 1, ND,  Xline4, y42, z22, BC,  OPSettings);
        size_t idx49  = Initializations::Sphere(Phi, 1, ND,  Xline4, y42, z32, BC,  OPSettings);
        size_t idx50  = Initializations::Sphere(Phi, 1, ND,  Xline4, y42, z42, BC,  OPSettings);

        Cx.SetInitialMoleFractions(Phi);
        Cx.SetBoundaryConditions(BC);

        Tx.SetInitial(BC);
        Tx.SetBoundaryConditions(BC);
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx.Tx, Tx.Tx.Bcells(),)
    {
        Tx.TxOld(i,j,k)=Tx.Tx(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Tx.SetBoundaryConditions(BC);

    double cp = 1000.0;                       ///< specific heat capacity
    Cx.SpecificHeatCapacityMix = cp;

    double Ts=273;
    double mus=1.68e-05; //it's for air
    double S=110.5;
    double R=PhysicalConstants::R;     ///<  Universal gas constant  j/mol k

    double AW_air;
    AW_air=(0.21*2.0*Cx.PT.GetData("O").AtomicWeight+0.79*2.0*Cx.PT.GetData("N").AtomicWeight)/1000.0;
    Cx.AtomicWeightMixture = AW_air;
    double Rm=R/AW_air;                 ///< gas constant for air j/kg k

    double Tref=(Tx.TBC0X+Tx.TSphere)/2.0;
    double Rhoref=FL.Pth0/(Rm*Tref);
    double muref = mus*(Tref/Ts)*sqrt(Tref/Ts)*(Ts+S)/(Tref+S);
    double nuref=muref/Rhoref;
    double Kref=cp*muref/Tx.Pr;

    double Tin1=Tx.TBC0X;
    double Rhoin1=FL.Pth0/(Rm*Tin1);
    double mucold = mus*(Tin1/Ts)*sqrt(Tin1/Ts)*(Ts+S)/(Tin1+S);
    double rhocold= FL.Pth0/(Rm*Tin1);
    double nucold = mucold/rhocold;
    double kcold=nucold*rhocold*Kref/(Rhoref*nuref); //if cp is constant

    double ND = OPSettings.TotalNz/6.0;

    FL.CalculateDensityTC(Tx, Cx, BC);

    for (size_t idx = 0; idx < Phi.FieldsStatistics.size(); idx++)
    {
        Phi.FieldsStatistics[idx].Force  = {0,0,0};
        Phi.FieldsStatistics[idx].Torque = {0,0,0};
        Phi.FieldsStatistics[idx].Acm    = {0,0,0};
        Phi.FieldsStatistics[idx].aAcc   = {0,0,0};
    }


    FL.DetectObstacles(Phi);
    FL.SetObstacleNodes(Phi, Vel);
    FL.SetInitialPopulationsTC(BC);
    FL.SetBoundaryConditions(BC);
    FL.CalculateHydrodynamicPressureAndMomentum(Tx, Cx, Vel);
    FL.CalculateFluidVelocities(Vel, Phi, BC);


    fstream Data;
    Data.open(OPSettings.TextDir+"Data.txt", ios_base::out);
    // writing data
    Data<<"Time step                           = ";
    Data<<RTC.dt<<endl;
    Data<<"Grid Spacing                        = ";
    Data<<OPSettings.dx<<endl;
    Data<<"dx/dt                               = ";
    Data<<OPSettings.dx/RTC.dt<<endl;
    Data<<"System Size in X Direction          = ";
    Data<<OPSettings.TotalNx<<endl;
    Data<<"System Size in Y Direction          = ";
    Data<<OPSettings.TotalNy<<endl;
    Data<<"System Size in Z Direction          = ";
    Data<<OPSettings.TotalNz<<endl;
    Data<<"Length in Z Direction              = ";
    Data<<OPSettings.TotalNz*OPSettings.dx<<endl;
    Data<<"Total Size in the Domain            = ";
    Data<<OPSettings.TotalNx*OPSettings.TotalNy*OPSettings.TotalNz/1000000.0<<" million"<<endl;
    Data<<"Minimum Kinematic Viscosity (m2/s)  = ";
    Data<<nucold<<endl;
    Data<<"Mean Kinematic Viscosity           = ";
    Data<<nuref<<endl;
    Data<<"Minimum LB Kinematic Viscosity      = ";
    Data<<nucold*RTC.dt/(OPSettings.dx*OPSettings.dx)<<endl;
    Data<<"Relaxation Time in LB               = ";
    Data<<0.5+3.0*nucold*RTC.dt/(OPSettings.dx*OPSettings.dx)<<endl;
    Data<<"Thermal Conductivity               = ";
    Data<<kcold<<endl;
    Data<<"Initial density for 273 k          = ";
    Data<<rhocold<<endl;
    Data<<"Reference density                  = ";
    Data<<Rhoref<<endl;
    Data<<"Reference dynamic viscosity        = ";
    Data<<muref<<endl;
    Data<<"specific Heat capacity (j/kgK)     = ";
    Data<<cp<<endl;
    Data<<"1/2 * dx^2/dt                      = ";
    Data<<(0.5*OPSettings.dx*OPSettings.dx/RTC.dt)<<endl;
    Data<<"Prandtl Number                     = ";
    Data<< Tx.Pr <<endl;
    Data<<"Rayleigh number                    = ";
    Data<< Tx.Pr*(-FL.GA[0])*(Rhoref*Rhoref)*(Tx.TSphere-Tx.TBC0X)*pow(2.0*ND*FL.dx,3)/(Tref*muref*muref) <<endl;
    Data<<"Initial Pressure(Pa)               = ";
    Data<< FL.Pth0 <<endl;
    Data<<"Initial Temperature(K)             = ";
    Data<< Tx.T0 <<endl;
    Data<<"Reynolds Number                    = ";
    Data<< FL.U0X*(2.0*ND*FL.dx)/nucold <<endl;
    Data<<"Inlet Velocity(m/s)                = ";
    Data<< FL.U0X <<endl;
    Data<<"Diameter of cylinder(m)            = ";
    Data<<2.0*ND*FL.dx<<endl;

    //end of writing data

#ifdef MPI_PARALLEL
if(MPI_RANK==0)
{

    cout << "Initialization stage done!" << endl;
    cout << "Entering the Time Loop!!!" << endl;
}
#endif

    double sumall=0.0;
    double NosavingT=1.0;
    std::clock_t c_start = std::clock();
    double sumtime=0.0;
    fstream timedata;
    timedata.open(OPSettings.TextDir+"timedata.txt",ios_base::out);

    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        FL.SetDivVelZeroNearObst();

        FL.CollisionTC(Tx, Cx, Vel, Phi);
        FL.SetBoundaryConditions(BC);

        SetInletVelocityX0(FL, BC, Phi);

        FL.Propagation(Phi,BC);
        FL.CalculateDensityTC(Tx, Cx, BC);
        FL.ApplyForces(Phi, Vel, Cx, Tx);

        FL.CalculateHydrodynamicPressureAndMomentum(Tx, Cx, Vel);
        FL.CalculateFluidVelocities(Vel, Phi, BC);

        //  Output to file
        if (RTC.WriteVTK())
        {
            // Write data in VTK format
            //Phi.WriteVTK(RTC.tStep,OPSettings);
            //Cx.WriteVTK(RTC.tStep, OPSettings);
            //Tx.WriteVTK(RTC.tStep, OPSettings);
            //Vel.WriteVTK(RTC.tStep, OPSettings);
            FL.WriteVTK(RTC.tStep, Phi, OPSettings);
        }
        //  Output to screen


        if(RTC.WriteToScreen())
        {

        #ifdef MPI_PARALLEL
            if(MPI_RANK==0)
            {
            cout << "Marching in time : "<<RTC.tStep << endl;
            cout << "====================================================================" << endl;

            std::clock_t c_end = std::clock();
            double time_elapsed_ms1 = 1000.0 *(c_end-c_start) / CLOCKS_PER_SEC;
            sumtime = time_elapsed_ms1;

            timedata<<"time step= ";
            timedata<<RTC.tStep<<", ";
            timedata<<"real time= ";
            timedata<<sumtime/1000.0<<" s"<<endl;
            }
        #endif

        }
    } //end time loop

    #ifdef MPI_PARALLEL
        MPI_Finalize ();
    #endif
    return 0;
}

void SetInletVelocityX0(FlowSolverLBM& FL, const BoundaryConditions& BC, const PhaseField& Phi)
{

#ifdef MPI_PARALLEL
{
    if(MPI_CART_RANK[0]==0)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.DensityWetting,FL.DensityWetting.Bcells(),)
        {
            if(i==0)
            {
                if(!FL.Obstacle(i,j,k-1) and !FL.Obstacle(i,j,k+1))
                {
                    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
                    {
                        if(FL.Do_ThermalComp)
                        {
                        	for(int ii = -FL.dNx; ii <= FL.dNx; ++ii)
                            {
                                for(int jj = -FL.dNy; jj <= FL.dNy; ++jj)
                                {
                                    for(int kk = -FL.dNz; kk <= FL.dNz; ++kk)
                                    {
                                        double lbUx1=0.0;
                                        double lbUy1=0.0;
                                        double lbUz1=0.0;
                                        for (auto it = Phi.Fields(i,j,k).cbegin(); it != Phi.Fields(i,j,k).cend(); it++)
                                        {
                                        	if (Phi.FieldsStatistics[it->index].State == AggregateStates::Liquid or
                                        				Phi.FieldsStatistics[it->index].State == AggregateStates::Gas)
                                        	{
                                        		lbUx1=FL.U0X*FL.dt/FL.dx * it->value;
                                        		lbUy1=FL.U0Y*FL.dt/FL.dx * it->value;
                                        		lbUz1=FL.U0Z*FL.dt/FL.dx * it->value;
                                        	}
                                        }
                                        double cu = ii*lbUx1 + jj*lbUy1 + kk*lbUz1;
                                        if(ii==1 and kk==0)
                                        {
                                        	
                                            FL.lbPopulations(i-ii,j-jj,k-kk)({n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k+kk)({n})(-ii,-jj,-kk)
                                                                        +2.0*FL.lbWeights[ii+1][jj+1][kk+1] * FL.DensityWetting(i-ii,j,k)({n})/FL.dRho * cu;
                                        }
                                        if(ii==1 and kk==-1)
                                        {
                                        	FL.lbPopulations(i-ii,j-jj,k-kk)({n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k+kk)({n})(-ii,-jj,-kk)
							            			+2.0*FL.lbWeights[ii+1][jj+1][kk+1] * FL.DensityWetting(i-ii,j,k)({n})/FL.dRho * cu;
                                        }

                                        if(ii==1 and kk==1)
                                        {
                                        	FL.lbPopulations(i-ii,j-jj,k-kk)({n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k+kk)({n})(-ii,-jj,-kk)
							            			+2.0*FL.lbWeights[ii+1][jj+1][kk+1] * FL.DensityWetting(i-ii,j,k)({n})/FL.dRho * cu;

                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}
#endif

#ifndef MPI_PARALLEL
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.DensityWetting,FL.DensityWetting.Bcells(),)
        {
            if(i==0)
            {
                if(!FL.Obstacle(i,j,k-1) and !FL.Obstacle(i,j,k+1))
                {
                    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
                    {
                        if(FL.Do_ThermalComp)
                        {
                        	for(int ii = -FL.dNx; ii <= FL.dNx; ++ii)
                            {
                                for(int jj = -FL.dNy; jj <= FL.dNy; ++jj)
                                {
                                    for(int kk = -FL.dNz; kk <= FL.dNz; ++kk)
                                    {
                                        double lbUx1=0.0;
                                        double lbUy1=0.0;
                                        double lbUz1=0.0;
                                        for (auto it = Phi.Fields(i,j,k).cbegin(); it != Phi.Fields(i,j,k).cend(); it++)
                                        {
                                        	if (Phi.FieldsStatistics[it->index].State == AggregateStates::Liquid or
                                        				Phi.FieldsStatistics[it->index].State == AggregateStates::Gas)
                                        	{
                                        		lbUx1=FL.U0X*FL.dt/FL.dx * it->value;
                                        		lbUy1=FL.U0Y*FL.dt/FL.dx * it->value;
                                        		lbUz1=FL.U0Z*FL.dt/FL.dx * it->value;
                                        	}
                                        }
                                        double cu = ii*lbUx1 + jj*lbUy1 + kk*lbUz1;
                                        if(ii==1 and kk==0)
                                        {
                                        	
                                            FL.lbPopulations(i-ii,j-jj,k-kk)({n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k+kk)({n})(-ii,-jj,-kk)
                                                                        +2.0*FL.lbWeights[ii+1][jj+1][kk+1] * FL.DensityWetting(i-ii,j,k)({n})/FL.dRho * cu;
                                        }
                                        if(ii==1 and kk==-1)
                                        {
                                        	FL.lbPopulations(i-ii,j-jj,k-kk)({n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k+kk)({n})(-ii,-jj,-kk)
							            			+2.0*FL.lbWeights[ii+1][jj+1][kk+1] * FL.DensityWetting(i-ii,j,k)({n})/FL.dRho * cu;
                                        }

                                        if(ii==1 and kk==1)
                                        {
                                        	FL.lbPopulations(i-ii,j-jj,k-kk)({n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k+kk)({n})(-ii,-jj,-kk)
							            			+2.0*FL.lbWeights[ii+1][jj+1][kk+1] * FL.DensityWetting(i-ii,j,k)({n})/FL.dRho * cu;

                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
}
#endif

}
