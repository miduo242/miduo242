/*
 *   This file is part of the OpenPhase (R) software library.
 *
 *   Copyright (c) 2009-2020 Ruhr-Universitaet Bochum,
 *                 Universitaetsstrasse 150, D-44801 Bochum, Germany
 *             AND 2018-2020 OpenPhase Solutions GmbH,
 *                 Wasserstrasse 494, D-44795 Bochum, Germany.
 *
 *    All rights reserved.
 *
 *
 *    DEVELOPMENT VERSION, DO NOT PUBLISH OR DISTRIBUTE.
 *
 *
 *   OpenPhase (R) is a joint development of Interdisciplinary Centre for
 *   Advanced Materials Simulation (ICAMS), Ruhr University Bochum
 *   and OpenPhase Solutions GmbH.
 *
 *   File created :   2020
 *   Main contributors :   Samad Vakili; Oleg Shchyglo;
 *
 */


#include "Settings.h"
#include "RunTimeControl.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"
#include "Initializations.h"
#include "Composition.h"
#include "Temperature.h"
#include "AdvectionHR/AdvectionHR.h"
#include "FluidDynamics/D3Q27.h"
#include "Velocities.h"
#include "Tools/CSVParser.h"
#include "VTK.h"
#include "Tools/TimeInfo.h"
#include "Info.h"
#include "Base/Includes.h"

using namespace std;
using namespace openphase;

class MPFSolver: public OPObject
{
 public:
    MPFSolver(){};
    MPFSolver(const Settings& locSettings, unsigned int boundary = 1)
    {
        Initialize(locSettings, boundary);
    };
    void Initialize(const Settings& locSettings, int boundary = 2);
    void ReadInput(string InputFileName);
    int Nx;
    int Ny;
    int Nz;
    int dNx;
    int dNy;
    int dNz;
    int Nphases;
    double dx;
    double dt;
    double Eta;

    int    Nb;                                                                  ///< Initial number of bubbles
    double Rmin;                                                                ///< Minimum initial bubble radius
    double Rmax;                                                                ///< Maximum initial bubble radius

    double p_crit;                                                              ///< Critical pressure for bubbles merging
    double p_factor;                                                            ///< Critical pressure magnifying factor
    double mu;
    double specific_v;
    double Cs2_Bub;
    dVector3 Uinit;
    double IntEng;
    double IntEngC;
    double ReactiveMassDensity;                                                 ///< Density of the gas forming agent in liquid [kg/m^3]
    double ActiveLayer;                                                         ///< Thickness of the active (gas producing) layer around the bubble [m]

    dMatrixNxN Wsqr;
    dMatrixNxN gamma;
    dVectorN   Mass;                                                            ///< Mass of the individual bubbles
    dVectorN   f;                                                               ///< Bulk free energy density
    dVectorN   Rho0;                                                            ///< Initial density of phases

    Storage3D<NodeA,0> RawIntF;
    dMatrix3x3 VelocityGradients(Velocities& Vel, const int i, const int j, const int k) const;  ///< Calculates local velocity gradient tensor
    dMatrix3x3 AverageVelocityGradients;                                        ///< Average velocity gradient in the simulation domain
    void CalculateAverageVelocityGradients(Velocities& Vel);                               ///< Calculates average velocity gradient in the simulation domain
    void MergeBubblesInContact(TimeInfo& Timer, BoundaryConditions& BC, PhaseField& Phi, int tStep);
    void DensityUpdate(TimeInfo& Timer, BoundaryConditions& BC, PhaseField& Phi, const int tStep);
    void PrintLocalData(TimeInfo& Timer, PhaseField& Phi) const;
    dVector3 AdvectVelocity(Velocities& Vel, const int i, const int j, const int k);             ///< This is changed and it is different than the one in the library
    dVector3 ViscousForce(Velocities& Vel, const int i, const int j, const int k);
    dVector3 InterfacialForce(TimeInfo& Timer, BoundaryConditions& BC,PhaseField& Phi, const int i, const int j, const int k);
    void SolveVelocity_PhaseField(TimeInfo& Timer, BoundaryConditions& BC, PhaseField& Phi, Velocities& Vel);         // bulk force (hydrostatic pressure force)
    void InterfaceEnergyCurv(TimeInfo& Timer, PhaseField& Phi);
    void WriteVTK(const int tStep, Settings& locSettings, PhaseField& Phi, Velocities& Vel,
                  const int precision = 16) const;
    void WriteSupplementary(const int tStep) const;
    void ReadSupplementary(const int tStep);
    string RawDataDir;
};

/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
{
    Settings           OPSettings;
    OPSettings.ReadInput();

    RunTimeControl     RTC(OPSettings);
    PhaseField         Phi(OPSettings);
    Velocities         Vel(OPSettings);
    BoundaryConditions BC(OPSettings);
    TimeInfo           Timer(OPSettings, "Execution Time Statistics");
    MPFSolver          MPFS;

    MPFS.Initialize(OPSettings);
    MPFS.ReadInput("ProjectInput.opi");
    MPFS.dt = RTC.dt;

    cout << "Initialization stage! Done!" << endl;

    if (RTC.Restart)
    {
        Info::WriteBlankLine();
        cout << "Restart data being read!" << endl;
        Phi.Read(BC, RTC.tStart);
        Vel.Read(BC, RTC.tStart);
        MPFS.ReadSupplementary(RTC.tStart);
        MPFS.DensityUpdate(Timer, BC, Phi, RTC.tStart);
        cout << "Done reading restart parameters!" << endl;
    }
    else
    {
        // Set initial geometry of the phases
        int LiquidPFIndex = Initializations::Single(Phi, 0, BC, OPSettings);
        //random bubble nucleation
        Initializations::FillGrainWithSpheres(Phi,LiquidPFIndex,1,
                MPFS.Rmin, MPFS.Rmax, BC, OPSettings, MPFS.Nb);
        Phi.WriteVTK(0, OPSettings);

        //Single bubble initialization
        //Initializations::Sphere(Phi, 1, iRadius, (OPSettings.Nx+1)/2, (OPSettings.Ny+1)/2, (OPSettings.Nz+1)/2, BC, OPSettings);

        //Two bubble initialization
        //Initializations::Sphere(Phi, 1, iRadius, (OPSettings.Nx+1)/4, (OPSettings.Ny+1)/2, (OPSettings.Nz+1)/2, BC, OPSettings);
        //Initializations::Sphere(Phi, 1, iRadius, 3*(OPSettings.Nx+1)/4, (OPSettings.Ny+1)/2, (OPSettings.Nz+1)/2, BC, OPSettings);

        MPFS.Mass.Allocate(Phi.FieldsStatistics.size());
        MPFS.f.Allocate(Phi.FieldsStatistics.size());

        for (size_t i = 0; i < Phi.FieldsStatistics.size(); i++)
        {
            double Vol   = MPFS.specific_v * Phi.FieldsStatistics[i].Volume;
            MPFS.Mass[i] = MPFS.Rho0[Phi.FieldsStatistics[i].Phase]*Vol;
        }
        double rho_ref = 0.9*MPFS.Rho0[1];

        for (size_t n = 0; n < Phi.FieldsStatistics.size(); n++)
        {
            size_t pIndex = Phi.FieldsStatistics[n].Phase;
            if(pIndex == 1)
            {
                Phi.FieldsStatistics[n].Density = 0.1*MPFS.Rho0[1]*(rand()%1000)/1000.0 + rho_ref;
            }
            else
            {
                Phi.FieldsStatistics[n].Density = MPFS.Rho0[0];
            }
        }

        MPFS.f[0] = -(6.4*MPFS.Rho0[0]/(6000.0-MPFS.Rho0[0]) - 0.00000075*MPFS.Rho0[0]);

        Vel.SetAverage(MPFS.Uinit);
        MPFS.DensityUpdate(Timer, BC, Phi, 0);
        Vel.SetBoundaryConditions(BC);
    }

    /**************************************************************start from here*/

    Info::WriteSimple("Starting simulation...");

    ignore_result(system("if [ -d obs ] ; then rm -rf obs; fi"));
    ignore_result(system("mkdir obs"));

    ofstream out_time;
    out_time.open("obs/time_terms.dat");

    //-------------- The Time Loop -------------//
    for (RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        Timer.SetStart();

        MPFS.InterfaceEnergyCurv(Timer, Phi);
        MPFS.SolveVelocity_PhaseField(Timer, BC, Phi, Vel);
        MPFS.MergeBubblesInContact(Timer, BC, Phi, RTC.tStep);
        MPFS.DensityUpdate(Timer, BC, Phi, RTC.tStep);
        Timer.SetTimeStamp("WritingTimeFile");

        // Output to files
        if (RTC.WriteVTK())
        {
            Phi.WriteVTK(RTC.tStep, OPSettings);
            Timer.SetTimeStamp("VTKPhi");
            Vel.WriteVTK(RTC.tStep, OPSettings);
            Timer.SetTimeStamp("VTKVel");
            MPFS.WriteVTK(RTC.tStep, OPSettings, Phi, Vel);
            Timer.SetTimeStamp("VTKLocal");
        }

        // Write restart output
        if (RTC.WriteRawData())
        {
            if (RTC.tStep == 20000 )  RTC.tFileWrite = 5000;
            if (RTC.tStep == 60000 )  RTC.tFileWrite = 10000;
            if (RTC.tStep == 100000)  RTC.tFileWrite = 20000;
            if (RTC.tStep == 200000)  RTC.tFileWrite = 50000;
            // Output PhaseField in binary format
            Phi.Write(RTC.tStep);
            Vel.Write(RTC.tStep);
            MPFS.WriteSupplementary(RTC.tStep);
        }

        // Output to screens
        if (RTC.WriteToScreen())
        {
            Info::WriteTimeStep(RTC.tStep, RTC.nSteps);
            //  Statistics
            //MPFS.PrintLocalData(Timer, Phi);
        }

        out_time << RTC.tStep << " ";
        for (size_t n = 0; n < Phi.FieldsStatistics.size(); n++)
        {
            out_time << Phi.FieldsStatistics[n].Density << " "
                     << Phi.FieldsStatistics[n].Volume*MPFS.specific_v << " "
                     << MPFS.f[n] << " ";
        }
        out_time << MPFS.IntEng << " " << MPFS.IntEngC << endl;
    }
    out_time.close();

    return 0;
}

void MPFSolver::Initialize(const Settings& locSettings, int boundary)
{
    thisclassname = "MPFSolver";

    Nx      = locSettings.Nx;
    Ny      = locSettings.Ny;
    Nz      = locSettings.Nz;
    dNx     = locSettings.dNx;
    dNy     = locSettings.dNy;
    dNz     = locSettings.dNz;
    dx      = locSettings.dx;
    Eta     = locSettings.Eta;

    Nphases = locSettings.Nphases;
    RawDataDir = locSettings.RawDataDir;

    specific_v = pow(dx, dNx + dNy + dNz);

    RawIntF.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, 3);

    Rho0.Allocate(Nphases);
    Mass.Allocate(Nphases);
    f.Allocate(Nphases);
    Wsqr.Allocate(Nphases);                                                     //interface energy coefficients
    gamma.Allocate(Nphases);

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void MPFSolver::ReadInput(string InputFileName)
{
    Info::WriteLineInsert("Multiphase Flow settings");
    Info::WriteStandard("Source", InputFileName.c_str());

    fstream inpF(InputFileName.c_str(), ios::in);

    if (!inpF)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };
    stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();

    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    Nb       = UserInterface::ReadParameterI(inp, moduleLocation, "Nb");
    Rmin     = UserInterface::ReadParameterD(inp, moduleLocation, "Rmin");
    Rmax     = UserInterface::ReadParameterD(inp, moduleLocation, "Rmax");

    mu       = UserInterface::ReadParameterD(inp, moduleLocation, "Mu");
    Cs2_Bub  = UserInterface::ReadParameterD(inp, moduleLocation, "Cs2_Bub");
    ReactiveMassDensity = UserInterface::ReadParameterD(inp, moduleLocation, "ReMass");
    ActiveLayer = UserInterface::ReadParameterD(inp, moduleLocation, "ALayer");
    p_factor    = UserInterface::ReadParameterD(inp, moduleLocation, "Mfactor");

    for (int alpha = 0; alpha < Nphases; alpha++)
    {
        stringstream counter;
        counter << alpha;
        Rho0[alpha] = UserInterface::ReadParameterD(inp, moduleLocation, string("Rho0_") + counter.str(), true, 0.0);
    }

    for (int alpha = 0; alpha < Nphases; alpha++)
    for (int beta = alpha; beta < Nphases; beta++)
    {
        stringstream counter;
        counter << alpha << "_" << beta;
        Wsqr(alpha,beta) = UserInterface::ReadParameterD(inp, moduleLocation, string("Wsqr_") + counter.str(), true, 0.0);
        Wsqr(beta,alpha) = Wsqr(alpha,beta);
        gamma(alpha,beta) = Pi*Pi*Wsqr(alpha,beta)/(Eta*Eta);
        gamma(beta,alpha) = gamma(alpha,beta);
    }

    p_crit = (dNx+dNy+dNz - 1.0)*Wsqr(0,1)*Pi*Pi/(8.0*Eta*Eta); //-2.0*f[0] is from average pressure
    Uinit[0] = UserInterface::ReadParameterD(inp, moduleLocation, string("u"));
    Uinit[1] = UserInterface::ReadParameterD(inp, moduleLocation, string("v"));
    Uinit[2] = UserInterface::ReadParameterD(inp, moduleLocation, string("w"));

    Info::WriteLine();
    Info::WriteBlankLine();
}

dMatrix3x3 MPFSolver::VelocityGradients(Velocities& Vel, const int i, const int j, const int k) const
{
    dMatrix3x3 dVel;
    double factor = 0.5/dx;

    dVel(0,0) = (Vel.Average(i+1,j,k)[0] - Vel.Average(i-1,j,k)[0])*factor;
    dVel(0,1) = (Vel.Average(i,j+1,k)[0] - Vel.Average(i,j-1,k)[0])*factor;
    dVel(0,2) = (Vel.Average(i,j,k+1)[0] - Vel.Average(i,j,k-1)[0])*factor;

    dVel(1,0) = (Vel.Average(i+1,j,k)[1] - Vel.Average(i-1,j,k)[1])*factor;
    dVel(1,1) = (Vel.Average(i,j+1,k)[1] - Vel.Average(i,j-1,k)[1])*factor;
    dVel(1,2) = (Vel.Average(i,j,k+1)[1] - Vel.Average(i,j,k-1)[1])*factor;

    dVel(2,0) = (Vel.Average(i+1,j,k)[2] - Vel.Average(i-1,j,k)[2])*factor;
    dVel(2,1) = (Vel.Average(i,j+1,k)[2] - Vel.Average(i,j-1,k)[2])*factor;
    dVel(2,2) = (Vel.Average(i,j,k+1)[2] - Vel.Average(i,j,k-1)[2])*factor;

    return dVel;
}

void MPFSolver::CalculateAverageVelocityGradients(Velocities& Vel)
{
    dMatrix3x3 dVelAVG;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Vel.Average,0,reduction(dMatrix3x3SUM:dVelAVG))
    {
        dVelAVG += VelocityGradients(Vel,i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    AverageVelocityGradients = dVelAVG/(Nx*Ny*Nz);
}

void MPFSolver::MergeBubblesInContact(TimeInfo& Timer, BoundaryConditions& BC ,PhaseField& Phi, int tStep)
{
    int Nthreads = omp_get_max_threads();
    vector<NodeAB> Overlap(Nthreads);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phi.Fields,0,)
    {
        int thread_num = omp_get_thread_num();
        if(Phi.Fields(i,j,k).flag)
        for (auto alpha = Phi.Fields(i,j,k).cbegin();
                  alpha != Phi.Fields(i,j,k).cend() - 1; ++alpha)
        for (auto  beta = alpha + 1;
                   beta != Phi.Fields(i,j,k).cend(); beta++)
        if(Phi.FieldsStatistics[alpha->index].Phase != 0 and
           Phi.FieldsStatistics[ beta->index].Phase != 0)                       // if phases of both phase-fields are right
        {
            Overlap[thread_num].set_sym1(alpha->index, beta->index, 1.0);       // phase-fields overlap is detected

            if (fabs(f[alpha->index] + f[beta->index]) > p_crit*p_factor)          // if pressure from both bubbles is above critical
            {
                Overlap[thread_num].set_sym2(alpha->index, beta->index, 1.0);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    for (int it = 1; it < Nthreads; it++)
    {
        Overlap[0] += Overlap[it];
    }

    for (auto it = Overlap[0].begin(); it < Overlap[0].end(); it++)
    {
        if (it->value1 == 1.0 and it->value2 == 1.0)                            // if both conditions (overlap and pressure) hold -> merge bubbles
        {
            cout << "Selective combine phase-fields: phase " << it->indexB
                 << " is merged into phase " << it->indexA
                 << ". Time step: " << tStep << endl;
            Phi.SelectiveCombinePhaseFields(BC,it->indexA,it->indexB);
            Mass[it->indexA] += Mass[it->indexB];
            Mass[it->indexB] = 0;
        }
    }

    Timer.SetTimeStamp("Merge");
}

void MPFSolver::DensityUpdate(TimeInfo& Timer, BoundaryConditions& BC, PhaseField& Phi, const int tStep)
{
    double ReactiveMassReduction = 0.0;
    double LiquidVolume = 0.0;

    for (size_t Np = 0; Np < Phi.FieldsStatistics.size(); Np++)
    {
        if(Phi.FieldsStatistics[Np].Phase == 1)
        {
            switch(dNx+dNy+dNz)
            {
                case 2: //2D
                {
                    double locRadius = sqrt(Phi.FieldsStatistics[Np].Volume/Pi)*dx;
                    double locSurfaceVolume = 2.0*Pi*locRadius*ActiveLayer;
                    double MassIncrease = locSurfaceVolume*ReactiveMassDensity;
                    ReactiveMassReduction += MassIncrease;
                    Mass[Np] += MassIncrease;
                    //cout << locSurfaceVolume*ReactiveMassDensity << endl;
                    break;
                }
                case 3: //3D
                {
                    double locRadius = pow(0.75*Phi.FieldsStatistics[Np].Volume/Pi,1.0/3.0)*dx;
                    double locSurfaceVolume = 4.0*Pi*locRadius*locRadius*ActiveLayer;
                    double MassIncrease = locSurfaceVolume*ReactiveMassDensity;
                    ReactiveMassReduction += MassIncrease;
                    Mass[Np] += MassIncrease;
                    break;
                }
            }
        }
        else
        {
            LiquidVolume += Phi.FieldsStatistics[Np].Volume;
        }
    }

    ReactiveMassDensity -= ReactiveMassReduction/(LiquidVolume*dx*dx*dx);

    for (size_t n = 0; n < Phi.FieldsStatistics.size(); n++)
    if(Phi.FieldsStatistics[n].Exist)
    {
        if(Phi.FieldsStatistics[n].Phase == 1)
        {
            double vol = specific_v * Phi.FieldsStatistics[n].Volume;
            Phi.FieldsStatistics[n].Density = Mass[n] / vol;
            f[n] = -Cs2_Bub*Phi.FieldsStatistics[n].Density;
        }
        else
        {
            f[n] = -(6.4*Phi.FieldsStatistics[n].Density/(6000.0 - Phi.FieldsStatistics[n].Density) - 0.00000075*Phi.FieldsStatistics[n].Density);
        }
    }
    Timer.SetTimeStamp("DensityUpdate");
}

void MPFSolver::PrintLocalData(TimeInfo& Timer, PhaseField& Phi) const
{
    for(size_t idx = 0; idx < Phi.FieldsStatistics.size(); idx++)
    {
        if(Phi.FieldsStatistics[idx].Volume > 0.0)
        {
            Info::WriteStandardNarrow( "Phase-Field", std::to_string(idx));
            Info::WriteStandardNarrow( "Density", std::to_string(Phi.FieldsStatistics[idx].Density));
            Info::WriteStandardNarrow( "Volume", std::to_string(Phi.FieldsStatistics[idx].Volume*specific_v));
            //Info::WriteStandardNarrow( "dVavg", AverageVelocityGradients.print());
        }
    }
    Timer.SetTimeStamp("PrintLocalData");
}

dVector3 MPFSolver::AdvectVelocity(Velocities& Vel, const int i, const int j, const int k)
{
    const double dxHalf_Inv = 0.5/dx;
    dVector3 locVec;

    for (auto dir = 0; dir < 3; dir++)
    {
        locVec[dir] = dxHalf_Inv*(
                ( Vel.Average(i,j,k)[0] - fabs(Vel.Average(i,j,k)[0]) )*(  Vel.Average(i+1,j,k)[dir] - Vel.Average(i,j,k)[dir] ) +
                ( Vel.Average(i,j,k)[0] + fabs(Vel.Average(i,j,k)[0]) )*( -Vel.Average(i-1,j,k)[dir] + Vel.Average(i,j,k)[dir] ) +
                ( Vel.Average(i,j,k)[1] - fabs(Vel.Average(i,j,k)[1]) )*(  Vel.Average(i,j+1,k)[dir] - Vel.Average(i,j,k)[dir] ) +
                ( Vel.Average(i,j,k)[1] + fabs(Vel.Average(i,j,k)[1]) )*( -Vel.Average(i,j-1,k)[dir] + Vel.Average(i,j,k)[dir] ) +
                ( Vel.Average(i,j,k)[2] - fabs(Vel.Average(i,j,k)[2]) )*(  Vel.Average(i,j,k+1)[dir] - Vel.Average(i,j,k)[dir] ) +
                ( Vel.Average(i,j,k)[2] + fabs(Vel.Average(i,j,k)[2]) )*( -Vel.Average(i,j,k-1)[dir] + Vel.Average(i,j,k)[dir] ) );
    }
    return locVec;
}

dVector3 MPFSolver::ViscousForce(Velocities& Vel, const int i, const int j, const int k)
{
    dVector3 locVec;
    double dx_sqr = dx*dx;

    for (int ii = -1; ii <= +1; ++ii)
    for (int jj = -1; jj <= +1; ++jj)
    for (int kk = -1; kk <= +1; ++kk)
    {
        locVec += Vel.Average(i+ii,j+jj,k+kk)*LBStencil3D[ii+1][jj+1][kk+1];
    }

    locVec = (locVec - Vel.Average(i,j,k)) * 6.0/dx_sqr;

    locVec[0] += (
                   Vel.Average(i+1,j,k)[0] + Vel.Average(i-1,j,k)[0] - 2.0*Vel.Average(i,j,k)[0] +
            0.25*( Vel.Average(i+1,j+1,k)[1]-Vel.Average(i+1,j-1,k)[1]-Vel.Average(i-1,j+1,k)[1]+Vel.Average(i-1,j-1,k)[1] ) +
            0.25*( Vel.Average(i+1,j,k+1)[2]-Vel.Average(i+1,j,k-1)[2]-Vel.Average(i-1,j,k+1)[2]+Vel.Average(i-1,j,k-1)[2] ) )/dx_sqr;

    locVec[1] += (
            0.25*( Vel.Average(i+1,j+1,k)[0]-Vel.Average(i+1,j-1,k)[0]-Vel.Average(i-1,j+1,k)[0]+Vel.Average(i-1,j-1,k)[0] ) +
                   Vel.Average(i,j+1,k)[1] + Vel.Average(i,j-1,k)[1] - 2.0*Vel.Average(i,j,k)[1] +
            0.25*( Vel.Average(i,j+1,k+1)[2]-Vel.Average(i,j+1,k-1)[2]-Vel.Average(i,j-1,k+1)[2]+Vel.Average(i,j-1,k-1)[2] ) )/dx_sqr;

    locVec[2] += (
            0.25*( Vel.Average(i+1,j,k+1)[0]-Vel.Average(i+1,j,k-1)[0]-Vel.Average(i-1,j,k+1)[0]+Vel.Average(i-1,j,k-1)[0] ) +
            0.25*( Vel.Average(i,j+1,k+1)[1]-Vel.Average(i,j+1,k-1)[1]-Vel.Average(i,j-1,k+1)[1]+Vel.Average(i,j-1,k-1)[1] ) +
                   Vel.Average(i,j,k+1)[2] + Vel.Average(i,j,k-1)[2] - 2.0*Vel.Average(i,j,k)[2] )/dx_sqr;

    return locVec;
}

dVector3 MPFSolver::InterfacialForce(TimeInfo& Timer, BoundaryConditions& BC,PhaseField& Phi, const int i, const int j, const int k)
{
    dVector3 locVec;

    double pre_f = 8.0/Pi;

    if (Phi.Fields(i,j,k).flag)
    for (auto alpha = Phi.Fields(i,j,k).cbegin();
              alpha != Phi.Fields(i,j,k).cend()-1; alpha++)
    for (auto beta = alpha+1;
              beta != Phi.Fields(i,j,k).cend(); beta++)
    {
        int  pIndexA = Phi.FieldsStatistics[alpha->index].Phase;
        int  pIndexB = Phi.FieldsStatistics[ beta->index].Phase;

        double phi_a = f[beta->index]*pre_f*sqrt(alpha->value*beta->value) -
                       f[alpha->index]*pre_f*sqrt(alpha->value*beta->value) -
                       0.5*Wsqr(pIndexA,pIndexB)*(beta->laplacian - alpha->laplacian) -
                       0.5*gamma(pIndexA,pIndexB)*(beta->value - alpha->value);

        if(Phi.Fields(i,j,k).size() > 2)
        for(auto ksi = Phi.Fields(i,j,k).cbegin();
                 ksi != Phi.Fields(i,j,k).cend(); ++ksi)
        if((ksi != alpha) && (ksi != beta))
        {
            int  pIndexG = Phi.FieldsStatistics[ ksi->index].Phase;

            phi_a += - 0.5*(Wsqr(pIndexA,pIndexG) - Wsqr(pIndexB,pIndexG))*ksi->laplacian -
                       0.5*(gamma(pIndexA,pIndexG) - gamma(pIndexB,pIndexG))*ksi->value;
        }

        phi_a *= 1.0/Phi.Fields(i,j,k).size();

        RawIntF(i,j,k).add_value(alpha->index, phi_a);
        RawIntF(i,j,k).add_value(beta->index, -phi_a);

        locVec += (alpha->gradient - beta->gradient)*phi_a;
    }
    return locVec;
}

void MPFSolver::SolveVelocity_PhaseField(TimeInfo& Timer, BoundaryConditions& BC, PhaseField& Phi,
                                         Velocities& Vel)
{
    CalculateAverageVelocityGradients(Vel);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phi.Fields,0,)
    {
        RawIntF(i,j,k).clear();

        dVector3 locVec_Adv = AdvectVelocity(Vel,i,j,k),
                 locVec_Vis = ViscousForce(Vel,i,j,k),
                 locVec_Int = InterfacialForce(Timer,BC,Phi,i,j,k);

        double RhoTotal = 0.0;
        for (auto alpha = Phi.Fields(i,j,k).cbegin();
                  alpha != Phi.Fields(i,j,k).cend(); alpha++)
        {
            RhoTotal += Phi.FieldsStatistics[alpha->index].Density*alpha->value;
        }

        if (RhoTotal < DBL_EPSILON)
        {
            cout << "*********     Error: Density is zero!     **********"<< endl;
            cout << i << " " << j << " " << k << " " << RhoTotal << endl;
            exit(1);
        }
        double dt_rho = dt/RhoTotal;

        Vel.Average(i,j,k)[0] += -dt*locVec_Adv[0] + dt_rho*( mu*locVec_Vis[0] - locVec_Int[0] );
        Vel.Average(i,j,k)[1] += -dt*locVec_Adv[1] + dt_rho*( mu*locVec_Vis[1] - locVec_Int[1] );
        Vel.Average(i,j,k)[2] += -dt*locVec_Adv[2] + dt_rho*( mu*locVec_Vis[2] - locVec_Int[2] );
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phi.Fields,0,)
    if (Phi.Fields(i,j,k).flag)
    for (auto alpha = Phi.Fields(i,j,k).begin();
              alpha != Phi.Fields(i,j,k).end(); alpha++)
    {
        alpha->value += dt*( -(Vel.Average(i,j,k)*alpha->gradient)
                       + RawIntF(i,j,k).get_value(alpha->index));
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Phi.Finalize(BC);
    Vel.SetBoundaryConditions(BC);
    Timer.SetTimeStamp("SolvePhaseField");
}

void MPFSolver::InterfaceEnergyCurv(TimeInfo& Timer, PhaseField& Phi)
{
    double sum[2] = {0.0, 0.0};
    int i = Nx/2;
    int j = Ny/2;
    for (int k = Nz/2+1; k < Nz; k++)
    {
        sum[0] += Wsqr(0,1)*Phi.Fields(i,j,k).get_gradient(1)[2] * Phi.Fields(i,j,k).get_gradient(1)[2]*dx;
        sum[1] += Wsqr(0,1)*Phi.Fields(i,j,k).get_gradient(1)[2] * Phi.Fields(i,j,k).get_gradient(1)[2]/(k-Nz/2);
    }
    IntEng  = sum[0];
    IntEngC = sum[1];
    Timer.SetTimeStamp("InterfaceEnergy");
}

void MPFSolver::WriteVTK(const int tStep, Settings& locSettings, PhaseField& Phi, Velocities& Vel, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = UserInterface::MakeFileName(locSettings.VTKDir, thisclassname+"_", tStep, ".vts");

    ListOfFields.push_back((VTK::Field_t) {"RawInterfacial", [this](int i,int j,int k){return RawIntF(i,j,k)[1];}});
    ListOfFields.push_back((VTK::Field_t) {"VelocityGradients_", [&Vel,this](int i,int j,int k){return VelocityGradients(Vel,i, j, k);}});
    ListOfFields.push_back((VTK::Field_t) {"Density",        [this, &Phi](int i,int j,int k)
    {
        double RhoTotal = 0.0;
        for (auto alpha = Phi.Fields(i,j,k).cbegin();
                  alpha != Phi.Fields(i,j,k).cend(); alpha++)
        {
            RhoTotal += Phi.FieldsStatistics[alpha->index].Density*alpha->value;
        }
        return RhoTotal;
    }});

    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void MPFSolver::WriteSupplementary(const int tStep) const
{
    std::string FileName = UserInterface::MakeFileName(RawDataDir,"SUP_",tStep,".dat");

    std::ofstream out(FileName.c_str(), std::ios::out | std::ios::binary);

    if (!out)
    {
        std::cout << "File \"" << FileName << "\" could not be created! Terminating!!!" << std::endl;
        exit(1);
    }

    double value = f[0];
    out.write(reinterpret_cast<char*>(&value), sizeof(double));
    for (size_t alpha = 0; alpha < Mass.size(); alpha++)
    {
        double valueR = Mass[alpha];
        out.write(reinterpret_cast<char*>(&valueR), sizeof(double));
    }
    out.close();
}

void MPFSolver::ReadSupplementary(const int tStep)
{
    std::string FileName = UserInterface::MakeFileName(RawDataDir,"SUP_", tStep, ".dat");

    std::fstream inp(FileName.c_str(), std::ios::in | std::ios::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be opened", thisclassname, "Read");
        exit(1);
    }

    inp.read(reinterpret_cast<char*>(&f[0]), sizeof(double));
    for (size_t alpha = 0; alpha < Mass.size(); alpha++)
    {
        inp.read(reinterpret_cast<char*>(&Mass[alpha]), sizeof(double));
    }
    Info::WriteStandardNarrow( "Supplementary data", "Binary Input Read");
}
