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
 *
 *   File created :   2011
 *   Main contributors :   Efim Borukhovich; Oleg Shchyglo
 *
 */

#include "ElectricalPotential.h"
#include "Composition.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Info.h"
#include "Base/UserInterface.h"
#include "BoundaryConditions.h"
#include "VTK.h"

namespace openphase
{

using namespace std;
ElectricalPotential::ElectricalPotential(Settings& locSettings,
                                         const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}
void ElectricalPotential::Initialize(Settings& locSettings)
{
    thisclassname = "ElectricalPotential";

    EpsInv  = 1.0/8.85418781762e-12;
    Faraday = 1.602e-19*6.022e23;

    CathodeCharge = 0.0;
    AnodeCharge   = 0.0;

    Nx      = locSettings.Nx;
    Ny      = locSettings.Ny;
    Nz      = locSettings.Nz;

    dNx     = locSettings.dNx;
    dNy     = locSettings.dNy;
    dNz     = locSettings.dNz;

    Size    = Nx*Ny*Nz;
    Nz2     = locSettings.Nz/2+1;
    ftSize  = Nx*Ny*Nz2;
    dx      = locSettings.dx;
    Ncomp   = locSettings.Ncomp;
    Nphases = locSettings.Nphases;

    DPi_nX = 2.0*Pi/double(Nx);
    DPi_nY = 2.0*Pi/double(Ny);
    DPi_nZ = 2.0*Pi/double(Nz);

    Q[0] = new double[ftSize];
    Q[1] = new double[ftSize];
    Q[2] = new double[ftSize];
    QXYZ();

    ftRHS       = new complex<double>[ftSize];
    ftPotential = new complex<double>[ftSize];
    RHS         = new double[Size];
    rlPotential = new double[Size];

    ForwardPlan  = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, RHS,reinterpret_cast<fftw_complex*> (ftRHS),FFTW_ESTIMATE);
    BackwardPlan = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, reinterpret_cast<fftw_complex*> (ftPotential), rlPotential, FFTW_ESTIMATE);

    VTKDir = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;

//    PotentialDataDir     = "PotentialData/";
//    int ignore = system(string("mkdir " + PotentialDataDir).c_str());
    size_t Bcells = locSettings.Bcells;
    Potential.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
//    Rho.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);

//    MolarVolume.Allocate(nPhases, nComp);
//    NuRef.Allocate(nPhases, nComp);
//    RhoRef.Allocate(nPhases, nComp);
//    Lattice.Allocate(nPhases, nComp);

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void ElectricalPotential::ReadInput(string InputFileName)
{
    // Raphael: NOTE use PhysicalConstants.h if necessary
    //const double N_Avogadro = 6.02214129e23;                        ///< Avogadro constant [mol^-1]
    //const double e_charge   = 1.602176565e-19;                      ///< elementar charge [C]

    Info::WriteLineInsert("ElectricalPotential input");
    Info::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };

    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);

    inp.close();
}

void ElectricalPotential::ReadInput(stringstream& inp)
{
    const double N_Avogadro = 6.02214129e23;                        ///< Avogadro constant [mol^-1]
    const double e_charge   = 1.602176565e-19;                      ///< elementar charge [C]

    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    Nelectrodes = UserInterface::ReadParameterI(inp, moduleLocation, "Nelectrodes");

    electrodes.resize(Nelectrodes);
    for(int i = 0; i < Nelectrodes; i++)
    {
        electrodes[i].index = i+1;
        electrodes[i].charge = 0;
        electrodes[i].interfaceVolume = 0.0;
        for(int dir = 0; dir < 3; dir++)
        {
            electrodes[i].position[dir] = 0;
        }
        stringstream converter;
        converter << "CntElectr" << i+1;
        electrodes[i].counterElectrode = UserInterface::ReadParameterI(inp, moduleLocation, converter.str())-1;
    }

    for (int q = 0 ; q < Ncomp; q++)
    {
        stringstream converter;
        converter << "Z" << q;
        ElementarCharges.push_back(UserInterface::ReadParameterI(inp, moduleLocation, converter.str()));
        MolarCharge.push_back(ElementarCharges[q]*N_Avogadro*e_charge);
    }

    Info::WriteLine();
    Info::WriteBlankLine();
}

//void ElectricalPotential::CalculateChargeDensity(MolarDensity& Nu, Diffusion& DF)
//{
//    #pragma omp parallel //OMP BEGIN
//    {
//        int nThreads = omp_get_num_threads();
//        int myThread = omp_get_thread_num();
//
//        int BegR=(myThread*(Nx+2))/nThreads;
//        int EndR=((myThread+1)*(Nx+2))/nThreads;
//
//        for (int i = BegR; i < EndR; i++)
//        for (int j =    0; j < Ny+2; j++)
//        for (int k =    0; k < Nz+2; k++)
//        {
//            Rho(i,j,k) = 0.;
//
//            double C[2];
//            C[0]=DF.Cx(i,j,k)/100.;
//            C[1]=1.-C[0];
//
//            for(int q = 0; q < nComp; q++)
//                Rho(i,j,k) += C[q] * Nu.Nu(i,j,k) * MolarCharge[q];
//        }//loop over the domain
//    }//OMP END
//}

//void ElectricalPotential::SetElectrodeCharges(PhaseField& Phase, Composition& Cx)
//{
//    CalculateElectrodeInterfaceVolumes(Phase);
//}

inline double ElectricalPotential::ChargeDensity(int i, int j, int k,
                                                 const PhaseField& Phase,
                                                 const Composition& Cx)
{
//    for(int eIdx = 0; eIdx < Nelectrodes; eIdx++)
//    {
//        if(i == electrodes[eIdx].position[0] and
//           j == electrodes[eIdx].position[1] and
//           k == electrodes[eIdx].position[2])
//        {
//            return electrodes[eIdx].charge;
//        }
//    }
    if(Phase.Interface(i,j,k))
    {
        int Index = -1;
        for(int eIdx = 0; eIdx < Nelectrodes; eIdx++)
        {
            if(Phase.Fields(i,j,k)[eIdx+1] > 0) Index=eIdx;                     ///assumes that there is no interface between electrodes
        }

//        cout<<"charge at "<<i<<" "<<j<<" "<<k<<" is: "
//            <<electrodes[Index].charge/**Phase.Fields(i,j,k)[Index]/
//                       electrodes[Index].interfaceVolume +
//                       Cx.phase(i,j,k,0)*Phase.Fields(i,j,k)[0]*Faraday/(13.*8.472e-5)*/
//           <<endl;
        if(Index<0)
        {
            stringstream message;
            message << "Electrode Index not found!" << endl;
            Info::WriteExit(message.str(),thisclassname, "ChargeDensity()");
            exit(13);
        }
        return electrodes[Index].charge*Phase.Fields(i,j,k)[Index+1]/
               electrodes[Index].interfaceVolume +
               Cx.MoleFractions(i,j,k)({0,0})*Phase.Fields(i,j,k)[0]*Faraday/(13.0*8.472e-5);
    }
    else
    {
        return Cx.MoleFractions(i,j,k)({0,0})*Faraday/(13.0*8.472e-5);                        ///13 atoms in a molecule of PC, 8.472e-5 mol/m^3
    }
}

void ElectricalPotential::CalculateElectrodeInterfaceVolumes(const PhaseField& Phase) //TODO: parallelize
{
    for(int eIdx = 0; eIdx < Nelectrodes; eIdx++)
    {
        electrodes[eIdx].interfaceVolume = 0;
    }//electrodes loop
    for(int eIdx = 0; eIdx < Nelectrodes; eIdx++)
    {
        for(int i = 0; i < Nx; i++)
        for(int j = 0; j < Ny; j++)
        for(int k = 0; k < Nz; k++)
        if(Phase.Interface(i,j,k))
        {
            electrodes[eIdx].interfaceVolume += Phase.Fields(i,j,k)[eIdx+1];
        }
    }//electrodes loop
}

//inline double ElectricalPotential::ChargeDensity(int i, int j, int k, MolarDensity& Nu, Composition& Cx)
//{
//    double result = 0.;
//
//    double C[2];
//    C[0]=Cx.total(i,j,k)/100.;
//    C[1]=1.-C[0];
//
//    for(int q = 0; q < Ncomp; q++)
//        result += C[q] * Nu.Nu(i,j,k) * MolarCharge[q];
//    return result;
//}

void ElectricalPotential::QXYZ(void)
{
    ///copied from the spectral solver
    for(int i = 0; i < Nx ; i++)
    for(int j = 0; j < Ny ; j++)
    for(int k = 0; k < Nz2; k++)
    {
         int XYZ = k + Nz2*(j + Ny*i);

         Q[0][XYZ] = DPi_nX*(i*(i <= Nx/2) - (Nx-i)*(i > Nx/2))/dx;
         Q[1][XYZ] = DPi_nY*(j*(j <= Ny/2) - (Ny-j)*(j > Ny/2))/dx;
         Q[2][XYZ] = DPi_nZ*(k*(k <= Nz/2) - (Nz-k)*(k > Nz/2))/dx;
    }
}

void ElectricalPotential::Solve(PhaseField& Phase, Composition& Cx, BoundaryConditions& BC)
{
    CalculateElectrodeInterfaceVolumes(Phase);
    for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
    for (int k = 0; k < Nz; k++)
    {
        RHS[k + Nz*(j + Ny*i)] = ChargeDensity(i, j, k, Phase, Cx);//Rho(i+1,j+1,k+1); //sin(i*6.28/Nx);
    }

    fftw_execute(ForwardPlan);

    for(int XYZ = 0; XYZ < ftSize; XYZ++)
    {
        /// ^phi = - ^rho / (i*q*i*q)
        double Q_sqr = Q[0][XYZ]*Q[0][XYZ] + Q[1][XYZ]*Q[1][XYZ] + Q[2][XYZ]*Q[2][XYZ];
        if (Q_sqr == 0)
        {
            ftPotential[XYZ] = 0.0;
        }
        else
        {
            ftPotential[XYZ] = ftRHS[XYZ] / Q_sqr * EpsInv;
        }

    }

    fftw_execute(BackwardPlan);

    for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
    for (int k = 0; k < Nz; k++)
    {
        Potential(i,j,k) = rlPotential[k + Nz*(j + Ny*i)]/Size;
    }

    ///Boundary Conditions (periodic only):
    BC.SetX(Potential);
    BC.SetY(Potential);
    BC.SetZ(Potential);
    /*
    for (int i=1; i < Nx+1; i++)
    for (int k=1; k < Nz+1; k++)
    {
        Potential(i,    0, k) = Potential(i, Ny, k);
        Potential(i, Ny+1, k) = Potential(i,  1, k);
    }

    for (int j=1; j < Ny+1; j++)
    for (int k=1; k < Nz+1; k++)
    {
        Potential(    0, j, k) = Potential( Nx, j, k);
        Potential( Nx+1, j, k) = Potential(  1, j, k);
    }

    for (int i=0; i<Nx+2; i++)
    for (int j=0; j<Ny+2; j++)
    {
        Potential(i, j,    0) = Potential(i, j, Nz);
        Potential(i, j, Nz+1) = Potential(i, j,  1);
    }

    //  Edges
    for (int i=1; i < Nx+1; i++)
    {
        Potential(i,    0,    0) = Potential(i, Ny, Nz);
        Potential(i, Ny+1,    0) = Potential(i,  1, Nz);
        Potential(i,    0, Nz+1) = Potential(i, Ny,  1);
        Potential(i, Ny+1, Nz+1) = Potential(i,  1,  1);
    }

    for (int j=1; j < Ny+1; j++)
    {
        Potential(0,    j,    0) = Potential(Nx, j, Nz);
        Potential(Nx+1, j,    0) = Potential(1,  j, Nz);
        Potential(0,    j, Nz+1) = Potential(Nx, j,  1);
        Potential(Nx+1, j, Nz+1) = Potential(1,  j,  1);
    }

    for (int k=1; k < Nz+1; k++)
    {
        Potential(0,    0,    k) = Potential(Nx, Ny, k);
        Potential(Nx+1, Ny+1, k) = Potential(1,   1, k);
        Potential(0,    Ny+1, k) = Potential(Nx,  1, k);
        Potential(Nx+1,    0, k) = Potential(1,  Ny, k);
    }

    //  Corners
    Potential(0,    0,    0) = Potential(Nx, Ny, Nz);
    Potential(Nx+1, 0,    0) = Potential(1,  Ny, Nz);
    Potential(0, Ny+1,    0) = Potential(Nx,  1, Nz);
    Potential(0,    0, Nz+1) = Potential(Nx, Ny,  1);
    Potential(Nx+1, Ny+1, 0) = Potential(1,   1, Nz);
    Potential(Nx+1, 0, Nz+1) = Potential(1,  Ny,  1);
    Potential(0, Ny+1, Nz+1) = Potential(Nx,  1,  1);
    Potential(Nx+1, Ny+1, Nz+1) = Potential(1, 1, 1);
    */
    // End Boundary Conditions
}

void ElectricalPotential::WriteChargeDensityVTK(const int tStep, const Settings& locSettings, const PhaseField& Phase, const Composition& Cx)
{
    CalculateElectrodeInterfaceVolumes(Phase);

    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = UserInterface::MakeFileName(locSettings.VTKDir, "ChargeDensity_", tStep, ".vts");
    ListOfFields.push_back((VTK::Field_t) {"ChargeDensity", [this,&Phase,&Cx](int i,int j,int k){return ChargeDensity(i, j, k, Phase, Cx)/*Rho(i,j,k)*/;}});
    VTK::Write(Filename, locSettings, ListOfFields);

}  //  WriteVTK

void ElectricalPotential::WritePotentialVTK(const int tStep, const Settings& locSettings) const
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = UserInterface::MakeFileName(locSettings.VTKDir, "ElectricalPotential_", tStep, ".vts");
    ListOfFields.push_back((VTK::Field_t) {"ElPot", [this](int i,int j,int k){return Potential(i,j,k);}});
    VTK::Write(Filename, locSettings, ListOfFields);
}  //  WriteVTK

}
