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
#include "HeatDiffusion.h"
#include "HeatSources.h"

using namespace std;
using namespace openphase;

void SolveAdvectionWENO5_DiffusedWall( Temperature& Tx, BoundaryConditions& BC, Velocities& Vel, double dt);
void SolveDiffusionCD2_DiffusedWall_Cons(PhaseField& Phase, Temperature& Tx, HeatDiffusion& HD, BoundaryConditions& BC,  FlowSolverLBM& FL, double dt);
void SetTemperatureForSpherebyPhaseIndex(PhaseField& Phase, Temperature& Tx, int List[], int len);
void SetInletVelocityX0(FlowSolverLBM& FL, const BoundaryConditions& BC, const PhaseField& Phi);
void SetCornerCorrectionX0(FlowSolverLBM& FL, const BoundaryConditions& BC, const PhaseField& Phi);
std::vector<double> CalculateHeatTransferOfGrainLocallyIR2L(PhaseField& Phase, Temperature& Tx, 
                        FlowSolverLBM& FL, HeatDiffusion& HD, Composition& Cx, Settings& locSettings,
											           const BoundaryConditions& BC, int GrainIndex);

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
    HeatDiffusion						HD(OPSettings,InputFile);
    //AdvectionHR                         ADHR (OPSettings);
    //TimeInfo                            Timer;

    double Diameter = 0.01;
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
	    //size_t idx1 = Initializations::Single(Phi, 0, BC, OPSettings);
	    vector<size_t> idx01 =  Initializations::TwoWalls( Phi,  0,  1,  15, BC,  OPSettings);

        //double Xshift=10.0/FL.dx/1000.;
        double Xshift=20.0;

        double firstGap = 0.002;
        double firstColumn = Diameter/2.0 + firstGap;
        double Gap_Zdir=1.23*Diameter;
        double Gap_Xdir=Gap_Zdir * cos(3.14/6.0);
        double DistinK=Gap_Zdir/FL.dx;

        double Xline1=firstColumn/FL.dx+Xshift;
        double Xline2=(Gap_Xdir+firstColumn)/FL.dx+Xshift;
        double Xline3=(2.0*Gap_Xdir+firstColumn)/FL.dx+Xshift;
        double ND = Diameter/FL.dx/2.0;

        double Sref=OPSettings.TotalNz/2.0;
     
        //second column
	    size_t idx2  = Initializations::Sphere(Phi, 1, ND,  Xline2, 0.0, Sref, BC,  OPSettings);
	    size_t idx3  = Initializations::Sphere(Phi, 1, ND,  Xline2, 0.0, Sref+DistinK, BC,  OPSettings);
	    size_t idx4  = Initializations::Sphere(Phi, 1, ND,  Xline2, 0.0, Sref+2.0*DistinK, BC,  OPSettings);
	    size_t idx5  = Initializations::Sphere(Phi, 1, ND,  Xline2, 0.0, Sref+3.0*DistinK, BC,  OPSettings);
	    size_t idx6  = Initializations::Sphere(Phi, 1, ND,  Xline2, 0.0, Sref-DistinK, BC,  OPSettings);
	    size_t idx7  = Initializations::Sphere(Phi, 1, ND,  Xline2, 0.0, Sref-2.0*DistinK, BC,  OPSettings);
	    size_t idx8  = Initializations::Sphere(Phi, 1, ND,  Xline2, 0.0, Sref-3.0*DistinK, BC,  OPSettings);


	    //first column
	    size_t idx9  = Initializations::Sphere(Phi, 1, ND,  Xline1, 0.0, Sref+DistinK/2.0, BC,  OPSettings);
	    size_t idx10 = Initializations::Sphere(Phi, 1, ND,  Xline1, 0.0, Sref+DistinK/2.0+DistinK, BC,  OPSettings);
	    size_t idx11 = Initializations::Sphere(Phi, 1, ND,  Xline1, 0.0, Sref+DistinK/2.0+2.0*DistinK, BC,  OPSettings);
	    size_t idx12 = Initializations::Sphere(Phi, 1, ND,  Xline1, 0.0, Sref-DistinK/2.0, BC,  OPSettings);
	    size_t idx13 = Initializations::Sphere(Phi, 1, ND,  Xline1, 0.0, Sref-DistinK/2.0-DistinK, BC,  OPSettings);
	    size_t idx14 = Initializations::Sphere(Phi, 1, ND,  Xline1, 0.0, Sref-DistinK/2.0-2.0*DistinK, BC,  OPSettings);
        
	    //third column
        size_t idx15 = Initializations::Sphere(Phi, 1, ND,  Xline3, 0.0, Sref+DistinK/2.0, BC,  OPSettings);
	    size_t idx16 = Initializations::Sphere(Phi, 1, ND,  Xline3, 0.0, Sref+DistinK/2.0+DistinK, BC,  OPSettings);
	    size_t idx17 = Initializations::Sphere(Phi, 1, ND,  Xline3, 0.0, Sref+DistinK/2.0+2.0*DistinK, BC,  OPSettings);
	    size_t idx18 = Initializations::Sphere(Phi, 1, ND,  Xline3, 0.0, Sref-DistinK/2.0, BC,  OPSettings);
	    size_t idx19 = Initializations::Sphere(Phi, 1, ND,  Xline3, 0.0, Sref-DistinK/2.0-DistinK, BC,  OPSettings);
	    size_t idx20 = Initializations::Sphere(Phi, 1, ND,  Xline3, 0.0, Sref-DistinK/2.0-2.0*DistinK, BC,  OPSettings);

         
         
	 	double k0w1 = Sref-7.5/FL.dx/1000.0;
	 	double k1w1 = Sref-5.5/FL.dx/1000.0;
	 	double kw0  = (k0w1 + k1w1)/2.0;
	 	double NT1 = k1w1 - k0w1;
	 	double k0w2 = Sref+5.5/FL.dx/1000.0;
	 	double k1w2 = Sref+7.5/FL.dx/1000.0;
	 	double kw1  = (k0w2 + k1w2)/2.0;
	 	double NT2 = k1w2 - k0w2;
	    size_t idnozzlewall1=Initializations::Rectangular(Phi,  1, Xshift,  0.0, NT1, Xshift/2.0-OPSettings.iWidth/2.0-1.0,  0.0,  kw0, BC, OPSettings);
	    size_t idnozzlewall2=Initializations::Rectangular(Phi,  1, Xshift,  0.0, NT2, Xshift/2.0-OPSettings.iWidth/2.0-1.0,  0.0,  kw1, BC, OPSettings);
         
         
        Cx.SetInitialMoleFractions(Phi);
        Cx.SetBoundaryConditions(BC);

        Tx.SetInitial(BC);
        Tx.SetBoundaryConditions(BC);
     }

    constexpr int len=3;
    int List[len]={2,9,12};
    SetTemperatureForSpherebyPhaseIndex(Phi, Tx, List, len);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN (i,j,k,Tx.Tx, Tx.Tx.Bcells(),)
    {
    	Tx.TxOld(i,j,k) =Tx.Tx(i,j,k) ;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Tx.SetBoundaryConditions(BC);

    double cp  = 1000.0;                       ///< specific heat capacity
    Cx.SpecificHeatCapacityMix = cp;

	double Ts=273;
	double mus=1.68e-05; //it's for air
	double S=110.5;
    double R= PhysicalConstants::R;     ///<  Universal gas constant  j/mol k

    double AW_air;
    AW_air=(0.21*2.0*Cx.PT.GetData("O").AtomicWeight+0.79*2.0*Cx.PT.GetData("N").AtomicWeight)/1000.0;
    Cx.AtomicWeightMixture = AW_air;
    double Rm=R/AW_air;				    ///< gas constant for air j/kg k

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

    vector<bool>checkgrain(len, false);
    vector<int> igc;
    int ncg=0;

	OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx.Tx,0,)
    {
        for (auto it = Phi.Fields(i,j,k).cbegin(); it != Phi.Fields(i,j,k).cend(); ++it)
        {
            int GrainlocIndex = it -> index;
            for (int ig = 0; ig < len; ig++)
            {
            	if(!checkgrain[ig])
            	{
            		if(List[ig]==GrainlocIndex)
            		{
                    	int PhaseIndex= Phi.FieldsStatistics[GrainlocIndex].Phase;
                        if(Phi.Fractions(i,j,k)[PhaseIndex]>0.5)
                        {
                			checkgrain[ig] = true;
                			ncg++;
                			igc.push_back(List[ig]);
                        }
            		}
            	}
            }
        }
    }
	OMP_PARALLEL_STORAGE_LOOP_END

	vector <fstream> HeatTransferdata(ncg);
	vector <fstream> NusseltNumber(ncg);
	for(int ist=0; ist<ncg; ist++)
	{
		string sg = to_string(igc[ist]);
		int id=0;
		#ifdef MPI_PARALLEL
	    	id = MPI_RANK;
		#endif
	    string sr = to_string(id);
	    string FileName = OPSettings.TextDir+"HeatData_Grain#"+sg+"_Proc#"+sr+".txt";
	    string FileName2 = OPSettings.TextDir+"NusseltNumber_Grain#"+sg+"_Proc#"+sr+".txt";
	    HeatTransferdata[ist].open(FileName,ios_base::out);
	    NusseltNumber[ist].open(FileName2,ios_base::out);
	}

	double ND = 0.01/FL.dx/2.0;
	double NuCoeff;
	if(Tx.TSphere != Tx.TBC0X)
	{
		NuCoeff	= Diameter / (Kref * (Tx.TSphere  - Tx.TBC0X ) );
	}
	else
	{
		NuCoeff=0.0;
	}
	vector<vector<double>> HeatTransfer(len);

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
	Data<<"Time step                  	       = ";
	Data<<RTC.dt<<endl;
	Data<<"Grid Spacing               	       = ";
	Data<<OPSettings.dx<<endl;
	Data<<"dx/dt                     	       = ";
	Data<<OPSettings.dx/RTC.dt<<endl;
	Data<<"System Size in X Direction 	       = ";
	Data<<OPSettings.TotalNx<<endl;
	Data<<"System Size in Y Direction 	       = ";
	Data<<OPSettings.TotalNy<<endl;
	Data<<"System Size in Z Direction 	       = ";
    Data<<OPSettings.TotalNz<<endl;
    Data<<"Length in Z Direction              = ";
    Data<<OPSettings.TotalNz*OPSettings.dx<<endl;
    Data<<"Total Size in the Domain   	       = ";
    Data<<OPSettings.TotalNx*OPSettings.TotalNy*OPSettings.TotalNz<<endl;
    Data<<"Minimum Kinematic Viscosity        = ";
    Data<<nucold<<endl;
    Data<<"Mean Kinematic Viscosity           = ";
    Data<<nuref<<endl;
    Data<<"Minimum LB Kinematic Viscosity 	   = ";
    Data<<nucold*RTC.dt/(OPSettings.dx*OPSettings.dx)<<endl;
    Data<<"Relaxation Time in LB        	   = ";
    Data<<0.5+3.0*nucold*RTC.dt/(OPSettings.dx*OPSettings.dx)<<endl;
    Data<<"Thermal Conductivity               = ";
    Data<<kcold<<endl;
    Data<<"Initial density for 273 k          = ";
    Data<<rhocold<<endl;
    Data<<"Reference density                  = ";
    Data<<Rhoref<<endl;
    Data<<"Reference dynamic viscosity        = ";
    Data<<muref<<endl;
    Data<<"specific Heat capacity             = ";
    Data<<cp<<endl;
    Data<<"1/2 * dx^2/dt                      = ";
    Data<<(0.5*OPSettings.dx*OPSettings.dx/RTC.dt)<<endl;
    Data<<"Prandtl Number                     = ";
    Data<< Tx.Pr <<endl;
    Data<<"Rayleigh number                    = ";
    Data<< Tx.Pr*(-FL.GA[0])*(Rhoref*Rhoref)*(Tx.TSphere-Tx.TBC0X)*pow(2.0*ND*FL.dx,3)/(Tref*muref*muref) <<endl;
    Data<<"Initial Pressure                   = ";
    Data<< FL.Pth0 <<endl;
    Data<<"Initial Temperature                = ";
    Data<< Tx.T0 <<endl;
    Data<<"Reynolds Number                    = ";
    Data<< FL.U0X*(2.0*ND*FL.dx)/nucold <<endl;
    Data<<"Inlet Velocity                     = ";
    Data<< FL.U0X <<endl;
    Data<<"Diameter of cylinder               = ";
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
    	HD.SetEffectiveProperties(Phi, Tx, FL);

        SolveAdvectionWENO5_DiffusedWall(Tx, BC, Vel, RTC.dt);
        SolveDiffusionCD2_DiffusedWall_Cons(Phi, Tx, HD, BC, FL, RTC.dt);
        Tx.SetBoundaryConditions(BC);
        SetTemperatureForSpherebyPhaseIndex(Phi, Tx, List, len);

        FL.SetDivVelZeroNearObst();

        FL.CollisionTC(Tx, Cx, Vel, Phi);
        FL.SetBoundaryConditions(BC);

        SetInletVelocityX0(FL, BC, Phi);
        SetCornerCorrectionX0(FL, BC, Phi);

        FL.Propagation(Phi,BC);
        FL.CalculateDensityTC(Tx, Cx, BC);
        FL.ApplyForces(Phi, Vel, Cx, Tx);

        FL.CalculateHydrodynamicPressureAndMomentum(Tx, Cx, Vel);
        FL.CalculateFluidVelocities(Vel, Phi, BC);

	    OMP_PARALLEL_STORAGE_LOOP_BEGIN (i,j,k,Tx.Tx,Tx.Tx.Bcells(),)
	    {
	    	Tx.TxOld(i,j,k) = Tx.Tx(i,j,k) ;
	    }
	    OMP_PARALLEL_STORAGE_LOOP_END

        //Timer.SetTimeStamp("FL.Solve");
        //  Output to file
        if (RTC.WriteVTK())
        {
            // Write data in VTK format
            //Phi.WriteVTK(RTC.tStep,OPSettings);
            //Cx.WriteVTK(RTC.tStep, OPSettings);
            Tx.WriteVTK(RTC.tStep, OPSettings);
            //Vel.WriteVTK(RTC.tStep, OPSettings);
            FL.WriteVTK(RTC.tStep, Phi, OPSettings);
        }
        //  Output to screen


        if(RTC.WriteToScreen())
        {
        	for (int ig=0; ig<len; ig++)
        	{
        		vector<double> HeatTransferLoc = CalculateHeatTransferOfGrainLocallyIR2L(Phi, Tx, FL, HD, Cx,OPSettings, BC, List[ig]);
        	    if(HeatTransferLoc[2]!=0.0)
        	    {
        	        for(int ist=0; ist<ncg; ist++)
        	        {
        	        	if(igc[ist]==List[ig])
        	        	{
        	    	    	HeatTransferdata[ist]<<RTC.tStep<<" ";
        	    	    	HeatTransferdata[ist]<<HeatTransferLoc[2]<<" ";
        	    	    	HeatTransferdata[ist]<<HeatTransferLoc[3]<<endl;

        	    	    	NusseltNumber[ist]<<RTC.tStep<<" ";
        	    	    	NusseltNumber[ist]<<HeatTransferLoc[2] * NuCoeff<<" ";
        	    	    	NusseltNumber[ist]<<HeatTransferLoc[3]<<endl;

        	        	}
        	        }
        	    }
        	}

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

void SolveAdvectionWENO5_DiffusedWall( Temperature& Tx, BoundaryConditions& BC, Velocities& Vel, double dt)
{

    double gama1m = 1.0/10.0;
    double gama2m=  3.0/5.0;
    double gama3m = 3.0/10.0;

    double gama1p = 3.0/10.0;
    double gama2p = 3.0/5.0;
    double gama3p = 1.0/10.0;
    double ep = Tx.dx * Tx.dx;
    //double ep = 1.0e-6;

    std::vector<int> dir(3);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx.Tx,0,)
    {
    	for (int direction=0; direction<3; ++direction)
    	{
    		dir[0] = (direction == 0)?1:0;
    		dir[1] = (direction == 1)?1:0;
    	    dir[2] = (direction == 2)?1:0;

            double C=0.0;

            if(direction == 0 && Tx.dNx == 1)
            {
                C=1.0;
            }
            else if(direction == 1 && Tx.dNy == 1)
            {
                C=1.0;
            }
            else if(direction == 2 && Tx.dNz == 1)
            {
                C=1.0;
            }
            if(C==1)
            {
    	        double q   = Tx.TxOld(i,j,k);
    	        double qp  = Tx.TxOld(i+dir[0]   ,j+dir[1]   ,k+dir[2]  );
    	        double qpp = Tx.TxOld(i+2*dir[0] ,j+2*dir[1] ,k+2*dir[2]);
    	        double qm  = Tx.TxOld(i-dir[0]   ,j-dir[1]   ,k-dir[2]  );
    	        double qmm = Tx.TxOld(i-2*dir[0] ,j-2*dir[1] ,k-2*dir[2]);
    	        double v   = Vel.Average(i,j,k)[direction];

    	        double beta1 = 13/12*(qmm-2.0*qm+q)*(qmm-2.0*qm+q)+1.0/4.0*(qmm-4.0*qm+3.0*q)*(qmm-4.0*qm+3.0*q);
    	        double beta2 = 13/12*(qm-2.0*q+qp)*(qm-2.0*q+qp)+1.0/4.0*(qm-qp)*(qm-qp);
    	        double beta3 = 13/12*(q-2.0*qp+qpp)*(q-2.0*qp+qpp)+1.0/4.0*(3.0*q-4.0*qp+qpp)*(3.0*q-4.0*qp+qpp);

    	        double Wk1m=gama1m*(1.0  +  (abs(beta1-beta3)/(ep*ep+beta1)) *(abs(beta1-beta3)/(ep*ep+beta1))  );  //weno-Z
    	        double Wk2m=gama2m*(1.0  +  (abs(beta1-beta3)/(ep*ep+beta2)) *(abs(beta1-beta3)/(ep*ep+beta2))  );
    	        double Wk3m=gama3m*(1.0  +  (abs(beta1-beta3)/(ep*ep+beta3)) *(abs(beta1-beta3)/(ep*ep+beta3))  );

    	        double Wk1p=gama1p*(1.0  +  (abs(beta1-beta3)/(ep*ep+beta1)) *(abs(beta1-beta3)/(ep*ep+beta1))  );
    	        double Wk2p=gama2p*(1.0  +  (abs(beta1-beta3)/(ep*ep+beta2)) *(abs(beta1-beta3)/(ep*ep+beta2))  );
    	        double Wk3p=gama3p*(1.0  +  (abs(beta1-beta3)/(ep*ep+beta3)) *(abs(beta1-beta3)/(ep*ep+beta3))  );

    	        double Wi1m=Wk1m/(Wk1m+Wk2m+Wk3m);
    	        double Wi2m=Wk2m/(Wk1m+Wk2m+Wk3m);
    	        double Wi3m=Wk3m/(Wk1m+Wk2m+Wk3m);

    	        double Wi1p=Wk1p/(Wk1p+Wk2p+Wk3p);
    	        double Wi2p=Wk2p/(Wk1p+Wk2p+Wk3p);
    	        double Wi3p=Wk3p/(Wk1p+Wk2p+Wk3p);

    	        double fhp1 = 1.0/3.0*qmm-7.0/6.0*qm+11.0/6.0*q;
    	        double fhp2 =-1.0/6.0*qm+5.0/6.0*q+1.0/3.0*qp;
    	        double fhp3 = 1.0/3.0*q+5.0/6.0*qp-1.0/6.0*qpp;

    	        double fhm1 =-1.0/6.0*qmm+5.0/6.0*qm+1.0/3.0*q;
    	        double fhm2 = 1.0/3.0*qm+5.0/6.0*q-1.0/6.0*qp;
    	        double fhm3 = 11.0/6.0*q-7.0/6.0*qp+1.0/3.0*qpp;

    	        double fhp = Wi1m*fhp1 + Wi2m*fhp2 +Wi3m*fhp3;
    	        double fhm = Wi1p*fhm1 + Wi2p*fhm2 +Wi3p*fhm3;

    	        Tx.Tx(i, j, k) += - dt/Tx.dx * v * (fhp-fhm);
            }
    	}
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void SolveDiffusionCD2_DiffusedWall_Cons(PhaseField& Phase, Temperature& Tx, HeatDiffusion& HD, BoundaryConditions& BC,  FlowSolverLBM& FL, double dt)
{
    std::vector<int> dir(3);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx.Tx,0,)
    {
    	double cp     = HD.EffectiveHeatCapacity(i,j,k);
    	double kc     = HD.EffectiveThermalConductivity(i,j,k);
        double Tc     = Tx.TxOld(i,j,k);
        double Diff = 0.0;

    	for (int direction=0; direction<3; ++direction)
    	{
    		dir[0] = (direction == 0)?1:0;
    		dir[1] = (direction == 1)?1:0;
    	    dir[2] = (direction == 2)?1:0;

    	    double ke    = 0.5*(HD.EffectiveThermalConductivity(i+dir[0],j+dir[1],k+dir[2])+HD.EffectiveThermalConductivity(i,j,k));
    	    double kw    = 0.5*(HD.EffectiveThermalConductivity(i,j,k)+HD.EffectiveThermalConductivity(i-dir[0],j-dir[1],k-dir[2]));

            double Tp     = Tx.TxOld(i+dir[0],j+dir[1],k+dir[2]);
    	    double Tm     = Tx.TxOld(i-dir[0],j-dir[1],k-dir[2]);

            double C=0.0;

            if(direction == 0 && Tx.dNx == 1)
            {
                C=1.0;
            }
            else if(direction == 1 && Tx.dNy == 1)
            {
                C=1.0;
            }
            else if(direction == 2 && Tx.dNz == 1)
            {
                C=1.0;
            }

            if(C==1)
            {
            	//if(Phase.Fields(i,j,k).flag!=1)
            	//{
            	//	Diff += 1.0/cp * dt/(Tx.dx*Tx.dx)* (  ke*(Tp-Tc) - kw*(Tc-Tm) );
            	//}
            	//else
            	//{
            		Diff +=  kc/cp * (Tp +Tm - 2.0* Tc)/(Tx.dx*Tx.dx) ;
            	//}
            }
        }

    	Tx(i, j, k) +=  dt * Diff;
    	FL.DivergenceVel(i,j,k)({0})= Diff / Tx.TxOld(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void SetTemperatureForSpherebyPhaseIndex(PhaseField& Phase, Temperature& Tx, int List[], int len)
{
	OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx.Tx,0,)
    {
        for ( int ilen = 0; ilen < len; ilen++)
        {
            for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
            {
                int GrinIndex = it -> index;
                int PhaseIndex= Phase.FieldsStatistics[it->index].Phase;
                if(GrinIndex==List[ilen])
                {
                    if(Phase.Fractions(i,j,k)[PhaseIndex]>=0.95)
                    {
				        Tx.Tx(i,j,k) = Tx.TSphere;
                    }
                }
           }
        }
    }
	OMP_PARALLEL_STORAGE_LOOP_END
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

void SetCornerCorrectionX0(FlowSolverLBM& FL, const BoundaryConditions& BC, const PhaseField& Phi)
{

#ifdef MPI_PARALLEL
    if(MPI_CART_RANK[0]==0)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.DensityWetting,FL.DensityWetting.Bcells(),)
        {
            if(i==0)
            {
                if(!FL.Obstacle(i,j,k) and FL.Obstacle(i,j,k-1))   //left corners
                {
                    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
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

                                    if(ii==1 and kk==1)
                                    {
                                    	FL.lbPopulations(i,j,k)({n})(-ii,-jj,-kk)=FL.lbPopulations(i+ii,j+jj,k+kk)({n})(-ii,-jj,-kk);
                                    }
                                    if(ii==1 and kk==-1)
                                    {
                                    	FL.lbPopulations(i-ii,j-jj,k-kk)({n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k-kk)({n})(-ii,-jj,kk);
                                    	FL.lbPopulations(i,j,k)({n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k-kk)({n})(-ii,-jj,kk);
                                    }
                                }
                            }
                        }

                    }
                }
                if(!FL.Obstacle(i,j,k) and FL.Obstacle(i,j,k+1))   //right corners
                {
                    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
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
                                    	FL.lbPopulations(i,j,k)({n})(-ii,-jj,-kk)=FL.lbPopulations(i+ii,j+jj,k+kk)({n})(-ii,-jj,-kk);
                                    }
                                    if(ii==1 and kk==1)
                                    {
                                    	FL.lbPopulations(i-ii,j-jj,k-kk)({n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k-kk)({n})(-ii,-jj,kk);
                                    	FL.lbPopulations(i,j,k)({n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k-kk)({n})(-ii,-jj,kk);
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

#ifndef MPI_PARALLEL
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.DensityWetting,FL.DensityWetting.Bcells(),)
        {
            if(i==0)
            {
                if(!FL.Obstacle(i,j,k) and FL.Obstacle(i,j,k-1))   //left corners
                {
                    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
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

                                    if(ii==1 and kk==1)
                                    {
                                    	FL.lbPopulations(i,j,k)({n})(-ii,-jj,-kk)=FL.lbPopulations(i+ii,j+jj,k+kk)({n})(-ii,-jj,-kk);
                                    }
                                    if(ii==1 and kk==-1)
                                    {
                                    	FL.lbPopulations(i-ii,j-jj,k-kk)({n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k-kk)({n})(-ii,-jj,kk);
                                    	FL.lbPopulations(i,j,k)({n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k-kk)({n})(-ii,-jj,kk);
                                    }
                                }
                            }
                        }

                    }
                }
                if(!FL.Obstacle(i,j,k) and FL.Obstacle(i,j,k+1))   //right corners
                {
                    for (size_t n = 0; n < FL.N_Fluid_Comp; ++n)
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
                                    	FL.lbPopulations(i,j,k)({n})(-ii,-jj,-kk)=FL.lbPopulations(i+ii,j+jj,k+kk)({n})(-ii,-jj,-kk);
                                    }
                                    if(ii==1 and kk==1)
                                    {
                                    	FL.lbPopulations(i-ii,j-jj,k-kk)({n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k-kk)({n})(-ii,-jj,kk);
                                    	FL.lbPopulations(i,j,k)({n})(ii,jj,kk)=FL.lbPopulations(i+ii,j+jj,k-kk)({n})(-ii,-jj,kk);
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

std::vector<double> CalculateHeatTransferOfGrainLocallyIR2L(PhaseField& Phase, Temperature& Tx, FlowSolverLBM& FL, HeatDiffusion& HD, Composition& Cx, Settings& locSettings,
											const BoundaryConditions& BC, int GrainIndex)
{
	double Q=0.0;
	double NoP=0.0;
	double Rank=0.0;
	#ifdef MPI_PARALLEL
	Rank=MPI_RANK;
	#endif
	double GI=GrainIndex;

	OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx.Tx,0,)
    {
        bool doneOnce=false;
     	for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
        {
            int GrainlocIndex = it -> index;

            if(GI==GrainlocIndex)
            {
            	int PhaseIndex = Phase.FieldsStatistics[GrainlocIndex].Phase;
            	if(Phase.Fractions(i,j,k)[PhaseIndex]>0.0)
            	{
            		int PhaseIndex = Phase.FieldsStatistics[GrainlocIndex].Phase;
            		for(int ii = -Tx.dNx; ii <= Tx.dNx; ++ii)
            		for(int jj = -Tx.dNy; jj <= Tx.dNy; ++jj)
            		for(int kk = -Tx.dNz; kk <= Tx.dNz; ++kk)
            		{
            			if(Phase.Fractions(i+ii,j+jj,k+kk)[PhaseIndex]==0.0)
            			{
                            if(!doneOnce)
            				{
                                size_t PhaseGasIndex;
                                for (auto ito = Phase.Fields(i+ii,j+jj,k+kk).cbegin(); ito != Phase.Fields(i+ii,j+jj,k+kk).cend(); ito++)
                                {
                                    if (Phase.FieldsStatistics[ito->index].State == AggregateStates::Liquid or
                                    Phase.FieldsStatistics[ito->index].State == AggregateStates::Gas)
                                    {
                                        PhaseGasIndex = ito->index ;
                                    }
                                }
                                dVector3 Norm = Phase.Normal(i,j,k, GrainlocIndex, PhaseGasIndex);
                                double m = Norm[0]/Norm[2];
                                double Alpha = atan2(Norm[0],Norm[2]);

                                //double dist = 1.5*dx;
                                double dist = 1.5;
                                double xp = double(i+Phase.OffsetX) + dist*sin(Alpha);
                                double yp = double(j+Phase.OffsetY);
                                double zp = double(k+Phase.OffsetZ) + dist*cos(Alpha);

                                double x1 = int(xp);
                                double x2 = x1+1.0;
                                double z1 = int(zp);
                                double z2 = z1 +1.0;

                                double Tpnew = (x2-xp)/(x2-x1)*( (z2-zp)/(z2-z1)* Tx.TxOld(x1-Phase.OffsetX,yp-Phase.OffsetY,z1-Phase.OffsetZ)
                                                + (zp-z1)/(z2-z1)* Tx.TxOld(x1-Phase.OffsetX,yp-Phase.OffsetY,z2-Phase.OffsetZ) )
                                                +(xp-x1)/(x2-x1)*( (z2-zp)/(z2-z1)* Tx.TxOld(x2-Phase.OffsetX,yp-Phase.OffsetY,z1-Phase.OffsetZ)
                                                + (zp-z1)/(z2-z1)* Tx.TxOld(x2-Phase.OffsetX,yp-Phase.OffsetY,z2-Phase.OffsetZ) );
                                NoP++;
                                double kpnew = (x2-xp)/(x2-x1)*( (z2-zp)/(z2-z1)* HD.EffectiveThermalConductivity(x1-Phase.OffsetX,yp-Phase.OffsetY,z1-Phase.OffsetZ)
                                                + (zp-z1)/(z2-z1)* HD.EffectiveThermalConductivity(x1-Phase.OffsetX,yp-Phase.OffsetY,z2-Phase.OffsetZ) )
                                                +(xp-x1)/(x2-x1)*( (z2-zp)/(z2-z1)* HD.EffectiveThermalConductivity(x2-Phase.OffsetX,yp-Phase.OffsetY,z1-Phase.OffsetZ)
                                                + (zp-z1)/(z2-z1)* HD.EffectiveThermalConductivity(x2-Phase.OffsetX,yp-Phase.OffsetY,z2-Phase.OffsetZ) );

                                double dist_2 = 3.0;
                                double xp_2 = double(i+Phase.OffsetX) + dist_2*sin(Alpha);
                                double yp_2 = double(j+Phase.OffsetY);
                                double zp_2 = double(k+Phase.OffsetZ) + dist_2*cos(Alpha);

                                double x1_2 = int(xp_2);
                                double x2_2 = x1_2+1.0;
                                double z1_2 = int(zp_2);
                                double z2_2 = z1_2 +1.0;

                                double Tpnew_2 = (x2_2-xp_2)/(x2_2-x1_2)*( (z2_2-zp_2)/(z2_2-z1_2)* Tx.TxOld(x1_2-Phase.OffsetX,yp_2-Phase.OffsetY,z1_2-Phase.OffsetZ)
                                                + (zp_2-z1_2)/(z2_2-z1_2)* Tx.TxOld(x1_2-Phase.OffsetX,yp_2-Phase.OffsetY,z2_2-Phase.OffsetZ) )
                                                +(xp_2-x1_2)/(x2_2-x1_2)*( (z2_2-zp_2)/(z2_2-z1_2)* Tx.TxOld(x2_2-Phase.OffsetX,yp_2-Phase.OffsetY,z1_2-Phase.OffsetZ)
                                                + (zp_2-z1_2)/(z2_2-z1_2)* Tx.TxOld(x2_2-Phase.OffsetX,yp_2-Phase.OffsetY,z2_2-Phase.OffsetZ) );

                                double kpnew_2 = (x2_2-xp_2)/(x2_2-x1_2)*( (z2_2-zp_2)/(z2_2-z1_2)* HD.EffectiveThermalConductivity(x1_2-Phase.OffsetX,yp_2-Phase.OffsetY,z1_2-Phase.OffsetZ)
                                        + (zp_2-z1_2)/(z2_2-z1_2)* HD.EffectiveThermalConductivity(x1_2-Phase.OffsetX,yp_2-Phase.OffsetY,z2_2-Phase.OffsetZ) )
                                        +(xp_2-x1_2)/(x2_2-x1_2)*( (z2_2-zp_2)/(z2_2-z1_2)* HD.EffectiveThermalConductivity(x2_2-Phase.OffsetX,yp_2-Phase.OffsetY,z1_2-Phase.OffsetZ)
                                        + (zp_2-z1_2)/(z2_2-z1_2)* HD.EffectiveThermalConductivity(x2_2-Phase.OffsetX,yp_2-Phase.OffsetY,z2_2-Phase.OffsetZ) );

                                double kpn = 0.5*(kpnew+kpnew_2);

                                Q += kpn *( Tpnew - Tpnew_2 )/((dist_2-dist)*Tx.dx) ;
                                doneOnce=true;


                            }
            			}
            		}
            	}
            }
        }
    }
	OMP_PARALLEL_STORAGE_LOOP_END
	if(NoP!=0) Q /= NoP;
	return vector<double>({Rank,GI,Q,NoP});
}
