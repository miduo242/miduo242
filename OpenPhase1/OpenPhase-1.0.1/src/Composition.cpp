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
 *   File created :   2012
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Matthias Stratmann
 *
 */

#include "Composition.h"
#include "Mechanics/ElasticProperties.h"
#include "Settings.h"
#include "Info.h"
#include "BoundaryConditions.h"
#include "Base/UserInterface.h"
#include "VTK.h"
#include "Velocities.h"
#include "AdvectionHR/AdvectionHR.h"
#include "Chemistry/PeriodicTable.h"
#include "H5Interface.h"

namespace openphase
{
using namespace std;
/*************************************************************************/

void Composition::Initialize(Settings& locSettings)
{
    thisclassname = "Composition";

    TotalNx = locSettings.TotalNx;
    OffsetX = locSettings.OffsetX;
    TotalNy = locSettings.TotalNy;
    OffsetY = locSettings.OffsetY;
    TotalNz = locSettings.TotalNz;
    OffsetZ = locSettings.OffsetZ;

    Nx  = locSettings.Nx;
    Ny  = locSettings.Ny;
    Nz  = locSettings.Nz;

    dNx = locSettings.dNx;
    dNy = locSettings.dNy;
    dNz = locSettings.dNz;

    dx      = locSettings.dx;
    Ncomp   = locSettings.Ncomp;
    Nphases = locSettings.Nphases;

    ElementNames = locSettings.ElementNames;
    PhaseNames = locSettings.PhaseNames;

    int Bcells = locSettings.Bcells;

    MoleFractions.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Nphases, Ncomp}, Bcells);

    Norm.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Nphases, Ncomp}, Bcells);
    MoleFractionsDot.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Nphases, Ncomp}, Bcells);

    MoleFractionsTotal.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Ncomp}, Bcells);

    Limiting.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);

    Initial.Allocate({Nphases, Ncomp});

    MoleFractionsAverage.Allocate({{Nphases, Ncomp}});
    MoleFractionsTotalAverage.Allocate({{Ncomp}});
    MoleFractionsInterfaceAverage.Allocate({{Nphases, Nphases, Ncomp}});
    MoleFractionsMIN.Allocate({{Nphases, Ncomp}});
    MoleFractionsMAX.Allocate({{Nphases, Ncomp}});

    TotInitial.resize(Ncomp);

    VTKDir = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;
    TextDir = locSettings.TextDir;

    PT.Initialize();

    AtStart = true;
    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void Composition::ReadInput(const string InputFileName)
{
    fstream inp(InputFileName.c_str(), ios::in);
    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput");
        exit(1);
    };

    Info::WriteLine();
    Info::WriteLineInsert("Composition properties");
    Info::WriteStandard("Source", InputFileName);

    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);

    inp.close();
    Info::WriteLine();
}

void Composition::ReadInput(stringstream& inp)
{
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    for(size_t n = 0; n < Ncomp; n++)
    for(size_t alpha = 0; alpha < Nphases; alpha++)
    {
        stringstream converter;
        converter << alpha << string("_") << ElementNames[n];
        string counter = converter.str();

        Initial({alpha, n}) = UserInterface::ReadParameterD(inp, moduleLocation, string("C0_") + converter.str());
        MoleFractionsMIN({alpha, n}) = UserInterface::ReadParameterD(inp, moduleLocation, string("CMIN_") + converter.str());
        MoleFractionsMAX({alpha, n}) = UserInterface::ReadParameterD(inp, moduleLocation, string("CMAX_") + converter.str());
    }

    Info::WriteLine();
    Info::WriteBlankLine();
}

void Composition::SetInitialMoleFractions(PhaseField& Phi, int mode)
{
    /* Simple composition set up, uniform component fractions per phase. */

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,MoleFractionsTotal, MoleFractionsTotal.Bcells(),)
    {
        MoleFractions(i,j,k) = Initial;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    CalculateTotalMoleFractions(Phi);
    CollectStatistics(Phi);
    CalculateTotalMolarVolume();

    Info::WriteLine("Total molar volume: " + to_string(TotalMolarVolume));
}


Tensor<double, 1> Composition::WeightFractionsTotal(int x, int y, int z) const
{
    /**This function calculates weight fractions of all elements in a given
     * grid cell.*/

    Tensor<double, 1> result({Ncomp});

    double div = 0.0;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        double locWeight = MoleFractionsTotal(x,y,z)({comp})*PT.GetData(ElementNames[comp]).AtomicWeight;
        div += locWeight;
        result({comp}) = locWeight;
    }
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        result({comp}) /= div;
    }
    return result;
}

void Composition::CalculateTotalMoleFractions(PhaseField& Phi)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,MoleFractionsTotal,0,)
    {
        MoleFractionsTotal(i,j,k).set_to_zero();

        for(size_t comp = 0; comp < Ncomp; comp++)
        for(size_t alpha = 0; alpha < Nphases; alpha++)
        {
            MoleFractionsTotal(i,j,k)({comp}) += Phi.Fractions(i,j,k)({alpha})*
                                             MoleFractions(i,j,k)({alpha,comp});
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Composition::CalculateTotalMolarVolume(void)
{
    TotalMolarVolume = 0.0;

    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        TotalMolarVolume += PT.GetData(ElementNames[comp]).MolarVolume
                            *MoleFractionsTotalAverage[comp];
    }
}

void Composition::CalculateMoleFractionsTotalAverage(void)
{
    Tensor<double, 1> locMoleFractionsTotalAverage({Ncomp});
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(x,y,z,MoleFractionsTotal,0,reduction(TensorD1Sum: locMoleFractionsTotalAverage))
    {
        locMoleFractionsTotalAverage += MoleFractionsTotal(x,y,z);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    MoleFractionsTotalAverage = locMoleFractionsTotalAverage;

#ifdef MPI_PARALLEL
    MPI_Allreduce(locMoleFractionsTotalAverage.data(), MoleFractionsTotalAverage.data(), Ncomp, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    MoleFractionsTotalAverage /= double(TotalNx*TotalNy*TotalNz);
}

void Composition::CalculateMoleFractionsAverage(PhaseField& Phase)
{
    Tensor<double, 3> locMoleFractionsInterfaceAverage({Nphases,Nphases,Ncomp});
    Tensor<double, 2> locMoleFractionsAverage({Nphases,Ncomp});

    Tensor<double, 2> NpointsAB({Nphases,Nphases});
    Tensor<double, 1> NpointsA({Nphases});
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(x,y,z,MoleFractions,0,reduction(TensorD2Sum: locMoleFractionsAverage) reduction(TensorD3Sum: locMoleFractionsInterfaceAverage) reduction(TensorD1Sum: NpointsA) reduction(TensorD2Sum: NpointsAB))
    {
        for(size_t alpha = 0; alpha < Nphases; alpha++)
        {
            if(Phase.Fractions(x,y,z)({alpha}) != 0.0)
            {
                NpointsA({alpha}) += 1;
                for(size_t comp = 0; comp < Ncomp; comp++)
                {
                    locMoleFractionsAverage({alpha, comp}) += MoleFractions(x,y,z)({alpha, comp});
                }
                for(size_t beta = 0; beta < Nphases; beta++)
                {
                    if(Phase.Fractions(x,y,z)({beta}) != 0.0)
                    {
                        NpointsAB({alpha, beta}) += 1;
                        for(size_t comp = 0; comp < Ncomp; comp++)
                        {
                            locMoleFractionsInterfaceAverage({alpha, beta, comp}) += Phase.Fractions(x,y,z)({alpha})*MoleFractions(x,y,z)({alpha, comp}) +
                                                                                     Phase.Fractions(x,y,z)({ beta})*MoleFractions(x,y,z)({ beta, comp});
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    MoleFractionsAverage = locMoleFractionsAverage;
    MoleFractionsInterfaceAverage = locMoleFractionsInterfaceAverage;

#ifdef MPI_PARALLEL
    Tensor<double, 1> tmpNpointsA = NpointsA;
    Tensor<double, 2> tmpNpointsAB = NpointsAB;
    MPI_Allreduce(tmpNpointsA.data(), NpointsA.data(), Nphases, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(tmpNpointsAB.data(), NpointsAB.data(), Nphases*Nphases, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(locMoleFractionsAverage.data(), MoleFractionsAverage.data(), Nphases*Ncomp, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(locMoleFractionsInterfaceAverage.data(), MoleFractionsInterfaceAverage.data(), Nphases*Nphases*Ncomp, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    for(size_t alpha = 0; alpha < Nphases; alpha++)
    {
        if(NpointsA({alpha}) != 0)
        {
            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                MoleFractionsAverage({alpha, comp}) /= NpointsA({alpha});
            }
        }
        else
        {
            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                MoleFractionsAverage({alpha, comp}) = MoleFractionsTotalAverage({comp});
            }
        }
    }

    for(size_t alpha = 0; alpha < Nphases; alpha++)
    for(size_t  beta = 0;  beta < Nphases; beta++)
    {
        if(NpointsAB({alpha, beta}) != 0)
        {
            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                MoleFractionsInterfaceAverage({alpha, beta, comp}) /= NpointsAB({alpha, beta});
            }
        }
        else
        {
            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                MoleFractionsInterfaceAverage({alpha, beta, comp}) = 0.5*(MoleFractionsAverage({alpha, comp}) + MoleFractionsAverage({beta, comp}));
            }
        }
    }
}

void Composition::Remesh(const int newNx, const int newNy, const int newNz, const BoundaryConditions& BC)
{
    MoleFractions.Remesh(newNx, newNy, newNz);
    MoleFractionsTotal.Remesh(newNx, newNy, newNz);
    Limiting.Reallocate(newNx, newNy, newNz);
    Norm.Reallocate(newNx, newNy, newNz);

    Nx = newNx;
    Ny = newNy;
    Nz = newNz;
    SetBoundaryConditions(BC);
    Info::WriteStandard(thisclassname, "Remeshed");
}

void Composition::MoveFrame(const int dx, const int dy, const int dz,
                            const BoundaryConditions& BC)
{
    int xBeg = (dx >= 0) + (dx < 0)*(Nx) - 1;
    int xEnd = (dx >= 0)*(Nx) + (dx < 0) - 1;
    int xInc = 1 - 2*(dx < 0);

    int yBeg = (dy >= 0) + (dy < 0)*(Ny) - 1;
    int yEnd = (dy >= 0)*(Ny) + (dy < 0) - 1;
    int yInc = 1 - 2*(dy < 0);

    int zBeg = (dz >= 0) + (dz < 0)*(Nz) - 1;
    int zEnd = (dz >= 0)*(Nz) + (dz < 0) - 1;
    int zInc = 1 - 2*(dz < 0);

    for(int i = xBeg; ((dx >= 0) and (i <= xEnd)) or ((dx < 0) and (i >= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j <= yEnd)) or ((dy < 0) and (j >= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k <= zEnd)) or ((dz < 0) and (k >= zEnd)); k += zInc)
    {
        MoleFractions(i,j,k) = MoleFractions(i + dx*dNx, j + dy*dNy, k + dz*dNz);
        MoleFractionsTotal(i,j,k) = MoleFractionsTotal(i + dx*dNx, j + dy*dNy, k + dz*dNz);
    }
    SetBoundaryConditions(BC);
    Info::WriteStandard(thisclassname, "Frame moved");
}

void Composition::ConsumePlane(const int dx, const int dy, const int dz,
                               const int x, const int y, const int z,
                               const BoundaryConditions& BC)
{
    int xBeg = (dx >= 0) + (dx < 0)*(Nx) - 1;
    int xEnd = (dx >= 0)*(Nx) + (dx < 0) - 1;
    int xInc = 1 - 2*(dx < 0);

    int yBeg = (dy >= 0) + (dy < 0)*(Ny) - 1;
    int yEnd = (dy >= 0)*(Ny) + (dy < 0) - 1;
    int yInc = 1 - 2*(dy < 0);

    int zBeg = (dz >= 0) + (dz < 0)*(Nz) - 1;
    int zEnd = (dz >= 0)*(Nz) + (dz < 0) - 1;
    int zInc = 1 - 2*(dz < 0);

    for(int i = xBeg; ((dx >= 0) and (i <= xEnd)) or ((dx < 0) and (i >= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j <= yEnd)) or ((dy < 0) and (j >= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k <= zEnd)) or ((dz < 0) and (k >= zEnd)); k += zInc)
    if((i-x)*dx+(j-y)*dy+(k-z)*dz >= 0)
    {
        MoleFractions(i,j,k) = MoleFractions(i + dx*dNx, j + dy*dNy, k + dz*dNz);
        MoleFractionsTotal(i,j,k) = MoleFractionsTotal(i + dx*dNx, j + dy*dNy, k + dz*dNz);
    }
    xBeg = (dx >= 0)*(Nx) + (dx < 0) - 1;
    xEnd = (dx >= 0) + (dx < 0)*(Nx) - 1;
    xInc = 2*(dx < 0) - 1;

    yBeg = (dy >= 0)*(Ny) + (dy < 0) - 1;
    yEnd = (dy >= 0) + (dy < 0)*(Ny) - 1;
    yInc = 2*(dy < 0) - 1;

    zBeg = (dz >= 0)*(Nz) + (dz < 0) - 1;
    zEnd = (dz >= 0) + (dz < 0)*(Nz) - 1;
    zInc = 2*(dz < 0) - 1;

    for(int i = xBeg; ((dx >= 0) and (i >= xEnd)) or ((dx < 0) and (i <= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j >= yEnd)) or ((dy < 0) and (j <= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k >= zEnd)) or ((dz < 0) and (k <= zEnd)); k += zInc)
    if((i-x)*dx+(j-y)*dy+(k-z)*dz < 0)
    {
        MoleFractions(i,j,k) = MoleFractions(i - dx*dNx, j - dy*dNy, k - dz*dNz);
        MoleFractionsTotal(i,j,k) = MoleFractionsTotal(i - dx*dNx, j - dy*dNy, k - dz*dNz);
    }
    SetBoundaryConditions(BC);
    Info::WriteStandard(thisclassname, "Plane consumed");
}

void Composition::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(MoleFractions);
    BC.SetY(MoleFractions);
    BC.SetZ(MoleFractions);

    BC.SetX(MoleFractionsTotal);
    BC.SetY(MoleFractionsTotal);
    BC.SetZ(MoleFractionsTotal);
}

void Composition::SetLimitsBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(Norm);
    BC.SetY(Norm);
    BC.SetZ(Norm);

    BC.SetXFlags(Limiting);
    BC.SetYFlags(Limiting);
    BC.SetZFlags(Limiting);
}

void Composition::Write(int tStep)
{

#ifdef MPI_PARALLEL
    string FileName = UserInterface::MakeFileName(RawDataDir, thisclassname + "_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = UserInterface::MakeFileName(RawDataDir, thisclassname + "_", tStep, ".dat");
#endif

    ofstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        std::stringstream message;
        message << "File \"" << FileName << "\" could not be created";
        Info::WriteExit(message.str(), thisclassname, "Write()");
#ifdef DEBUG
        throw std::runtime_error(message.str());
#else
        std::exit(EXIT_FAILURE);
#endif
    };

    int tmp = Nx;
    out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
    tmp = Ny;
    out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
    tmp = Nz;
    out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
    size_t tmp2 = Nphases;
    out.write(reinterpret_cast<char*>(&tmp2), sizeof(size_t));
    tmp2 = Ncomp;
    out.write(reinterpret_cast<char*>(&tmp2), sizeof(size_t));

    STORAGE_LOOP_BEGIN(i,j,k,MoleFractions,0)
        out.write(reinterpret_cast<char*>(MoleFractions(i,j,k).data()), MoleFractions(i,j,k).size()*sizeof(double));
    STORAGE_LOOP_END
    STORAGE_LOOP_BEGIN(i,j,k,MoleFractionsTotal,0)
        out.write(reinterpret_cast<char*>(MoleFractionsTotal(i,j,k).data()), MoleFractionsTotal(i,j,k).size()*sizeof(double));
    STORAGE_LOOP_END
}

bool Composition::Read(BoundaryConditions& BC, int tStep)
{

#ifdef MPI_PARALLEL
    string FileName = UserInterface::MakeFileName(RawDataDir, thisclassname + "_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = UserInterface::MakeFileName(RawDataDir, thisclassname + "_", tStep, ".dat");
#endif


    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        std::stringstream message;
        message << "File \"" << FileName << "\" could not be opened";
        Info::WriteWarning(message.str(), thisclassname, "Read()");
        return false;
    }

    int locNx = Nx;
    int locNy = Ny;
    int locNz = Nz;
    size_t locNphases = Nphases;
    size_t locNcomp = Ncomp;
    inp.read(reinterpret_cast<char*>(&locNx), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNy), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNz), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNphases), sizeof(size_t));
    inp.read(reinterpret_cast<char*>(&locNcomp), sizeof(size_t));
    if(locNx != Nx or locNy != Ny or locNz != Nz or locNphases != Nphases or locNcomp != Ncomp)
    {
        stringstream message;
        message << "Inconsistent system dimensions!\n"
                << "Input data dimensions: (" << locNx << ", " << locNy << ", " << locNz << ") grid points.\n"
                << "Input data Nphases: " << locNphases << "\n"
                << "Input data Ncomp: " << locNcomp << "\n"
                << "Required data dimensions: (" << Nx << ", " << Ny << ", " << Nz << ") grid points.\n"
                << "Required data Nphases: " << Nphases << "\n"
                << "Required data Ncomp: " << Ncomp << "\n";
        Info::WriteWarning(message.str(), thisclassname, "Read()");
        return false;
    }

    STORAGE_LOOP_BEGIN(i,j,k,MoleFractions,0)
        inp.read(reinterpret_cast<char*>(MoleFractions(i,j,k).data()), MoleFractions(i,j,k).size()*sizeof(double));
    STORAGE_LOOP_END
    STORAGE_LOOP_BEGIN(i,j,k,MoleFractionsTotal,0)
        inp.read(reinterpret_cast<char*>(MoleFractionsTotal(i,j,k).data()), MoleFractionsTotal(i,j,k).size()*sizeof(double));
    STORAGE_LOOP_END
    STORAGE_LOOP_BEGIN(i,j,k,MoleFractionsTotal,0)
        for(size_t n = 0; n < Ncomp; n++)
        {
            TotInitial[n] += MoleFractionsTotal(i, j, k)({n});
        }
    STORAGE_LOOP_END
    for(size_t n = 0; n < Ncomp; n++)
    {
        TotInitial[n] /= double(Nx*Ny*Nz);
    }
    SetBoundaryConditions(BC);

    // Calculation MoleFractionsTotalAverage is need for CalculateTotalMolarVolume()
    CalculateMoleFractionsTotalAverage();
    CalculateTotalMolarVolume();
    Info::WriteStandard(thisclassname, "Binary input loaded");
    return true;
}


void Composition::WriteH5(int tStep, H5Interface& H5)
{
    #ifdef H5OP
    std::vector<double> dbuffer;
    dbuffer.push_back(Nx);
    dbuffer.push_back(Ny);
    dbuffer.push_back(Nz);
    dbuffer.push_back(Nphases);
    dbuffer.push_back(Ncomp);
    H5.WriteCheckPoint(tStep, "CxDomain", dbuffer);
    dbuffer.clear();
    dbuffer = MoleFractions.pack();
    H5.WriteCheckPoint(tStep, "MoleFractions", dbuffer);
    dbuffer.clear();
    dbuffer = MoleFractionsTotal.pack();
    H5.WriteCheckPoint(tStep, "MoleFractionsTotal", dbuffer);
    #else
    std::cerr << "OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"" << std::endl;
    exit(5);
    #endif
}

bool Composition::ReadH5(int tStep, H5Interface& H5)
{
    #ifdef H5OP
    std::vector<double> dbuffer;
    H5.ReadCheckPoint(tStep, "CxDomain", dbuffer);
    int locNx = dbuffer[0];
    int locNy = dbuffer[1];
    int locNz = dbuffer[2];
    size_t locNphases = dbuffer[3];
    size_t locNcomp = dbuffer[4];
    dbuffer.clear();
    if(locNx != Nx or locNy != Ny or locNz != Nz or locNphases != Nphases or locNcomp != Ncomp)
    {
        stringstream message;
        message << "Inconsistent system dimensions!\n"
                << "Input data dimensions: (" << locNx << ", " << locNy << ", " << locNz << ") grid points.\n"
                << "Input data Nphases: " << locNphases << "\n"
                << "Input data Ncomp: " << locNcomp << "\n"
                << "Required data dimensions: (" << Nx << ", " << Ny << ", " << Nz << ") grid points.\n"
                << "Required data Nphases: " << Nphases << "\n"
                << "Required data Ncomp: " << Ncomp << "\n";
        Info::WriteWarning(message.str(), thisclassname, "Read()");
        return false;
    }
    H5.ReadCheckPoint(tStep, "MoleFractions", dbuffer);
    MoleFractions.unpack(dbuffer);
    dbuffer.clear();
    H5.ReadCheckPoint(tStep, "MoleFractionsTotal", dbuffer);
    MoleFractionsTotal.unpack(dbuffer);
    return true;
    #else
    std::cerr << "OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"" << std::endl;
    exit(5);
    #endif
}


void Composition::ReadSubset(BoundaryConditions& BC,
                             string FileName,
                             int offsetX,
                             int offsetY,
                             int offsetZ,
                             bool legacy_format,
                             int inpNx,
                             int inpNy,
                             int inpNz)
{
    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be opened", thisclassname, "ReadSubset()");
        exit(1);
    };
    int locNx = Nx;
    int locNy = Ny;
    int locNz = Nz;
    size_t locNphases = Nphases;
    size_t locNcomp = Ncomp;

    if(legacy_format)
    {
        locNx = inpNx;
        locNy = inpNy;
        locNz = inpNz;
    }
    else
    {
        inp.read(reinterpret_cast<char*>(&locNx), sizeof(int));
        inp.read(reinterpret_cast<char*>(&locNy), sizeof(int));
        inp.read(reinterpret_cast<char*>(&locNz), sizeof(int));
        inp.read(reinterpret_cast<char*>(&locNphases), sizeof(int));
        inp.read(reinterpret_cast<char*>(&locNcomp), sizeof(int));

    }
    if(Nx + offsetX > locNx or
       Ny + offsetY > locNy or
       Nz + offsetZ > locNz or
       locNphases != Nphases or
       locNcomp != Ncomp)
    {
        stringstream message;
        message << "Inconsistent system dimensions!\n"
                << "Input data dimensions: (" << locNx << ", " << locNy << ", " << locNz << ") grid points.\n"
                << "Input data Nphases: " << locNphases << "\n"
                << "Input data Ncomp: " << locNcomp << "\n"
                << "Required data dimensions: (" << Nx << ", " << Ny << ", " << Nz << ") grid points.\n"
                << "Required data Nphases: " << Nphases << "\n"
                << "Required data Ncomp: " << Ncomp << "\n";
        Info::WriteExit(message.str(), thisclassname, "ReadSubset()");
        exit(1);
    }

    for(int i = 0; i < locNx; i++)
    for(int j = 0; j < locNy; j++)
    for(int k = 0; k < locNz; k++)
    {
        Tensor<double, 2> locPhaseComp;
        locPhaseComp.Allocate({Nphases, Ncomp});
        inp.read(reinterpret_cast<char*>(locPhaseComp.data()), locPhaseComp.size()*sizeof(double));

        if((i >= offsetX and j >= offsetY and k >= offsetZ) and
           (i < Nx + offsetX and j < Ny + offsetY and k < Nz + offsetZ))
        {
            int x = i - offsetX;
            int y = j - offsetY;
            int z = k - offsetZ;
            MoleFractions(x,y,z) = locPhaseComp;
        }
    }

    for(int i = 0; i < locNx; i++)
    for(int j = 0; j < locNy; j++)
    for(int k = 0; k < locNz; k++)
    {
        Tensor<double, 1> locTotalComp;
        locTotalComp.Allocate({Ncomp});
        inp.read(reinterpret_cast<char*>(locTotalComp.data()), locTotalComp.size()*sizeof(double));

        if((i >= offsetX and j >= offsetY and k >= offsetZ) and
           (i < Nx + offsetX and j < Ny + offsetY and k < Nz + offsetZ))
        {
            int x = i - offsetX;
            int y = j - offsetY;
            int z = k - offsetZ;
            MoleFractionsTotal(x,y,z) = locTotalComp;
        }
    }

    STORAGE_LOOP_BEGIN(i,j,k,MoleFractionsTotal,0)
        for(size_t n = 0; n < Ncomp; n++)
        {
            TotInitial[n] += MoleFractionsTotal(i, j, k)({n});
        }
    STORAGE_LOOP_END
    for(size_t n = 0; n < Ncomp; n++)
    {
        TotInitial[n] /= double(Nx*Ny*Nz);
    }

    SetBoundaryConditions(BC);
    Info::WriteStandard(thisclassname, "Binary input loaded");
}

void Composition::WriteDistortedVTK(int tStep, const Settings& locSettings, const ElasticProperties& EP)
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        const std::string nameComp1 = "WeightFractionsTotal_" + ElementNames[comp];
        ListOfFields.push_back((VTK::Field_t) {nameComp1, [this, comp](int i,int j,int k){return WeightFractionsTotal(i,j,k)({comp});}});
        const std::string nameComp2 = "MoleFractionsTotal_" + ElementNames[comp];
        ListOfFields.push_back((VTK::Field_t) {nameComp2, [this, comp](int i,int j,int k){return MoleFractionsTotal(i,j,k)({comp});}});
        for (size_t alpha = 0; alpha < Nphases; ++alpha)
        {
            const std::string namePhase = "MoleFractionsPhase_" + ElementNames[comp] + "(" + to_string(alpha) + ")";
            ListOfFields.push_back((VTK::Field_t) {namePhase, [this, comp, alpha](int i,int j,int k){return MoleFractions(i,j,k)({alpha,comp});}});
        }
    }
    std::string Filename = UserInterface::MakeFileName(VTKDir, thisclassname + "Distorted_", tStep, ".vts");

    VTK::WriteDistorted(Filename, locSettings, EP, ListOfFields);
}

void Composition::WriteVTK(const int tStep, const Settings& locSettings, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        const std::string nameComp1 = "WeightFractionsTotal_" + ElementNames[comp];
        ListOfFields.push_back((VTK::Field_t) {nameComp1, [this, comp](int i,int j,int k){return WeightFractionsTotal(i,j,k)({comp});}});
        const std::string nameComp2 = "MoleFractionsTotal_" + ElementNames[comp];
        ListOfFields.push_back((VTK::Field_t) {nameComp2, [this, comp](int i,int j,int k){return MoleFractionsTotal(i,j,k)({comp});}});
        for (size_t alpha = 0; alpha < Nphases; ++alpha)
        {
            const std::string namePhase = "MoleFractionsPhase_" + ElementNames[comp] + "(" + to_string(alpha) + ")";
            ListOfFields.push_back((VTK::Field_t) {namePhase, [this, comp, alpha](int i,int j,int k){return MoleFractions(i,j,k)({alpha,comp});}});
        }
    }
    std::string Filename = UserInterface::MakeFileName(VTKDir, thisclassname + '_', tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}

void Composition::PrintPointStatistics(int x, int y, int z)
{
    cout << "Components:\t";
    for(size_t comp = 0; comp < Ncomp; comp++) cout << ElementNames[comp] << "\t";
    cout << endl;
    cout << "Total:\t\t";
    for(size_t comp = 0; comp < Ncomp; comp++) cout << MoleFractionsTotal(x, y, z)({comp}) << "\t";
    cout << endl;
    cout << "Phase:  " << endl;

    for (size_t n = 0; n < Nphases; n++)
    {
        cout << "Phase " << n << ":\t";
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            cout << MoleFractions(x,y,z)({n, comp}) << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

void Composition::WriteStatistics(int tStep, double dt)
{
    vector<double> total(Ncomp, 0);
    vector<double> deviation(Ncomp, 0);
    double sim_time = tStep*dt;

    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        for(int i = 0; i < Nx; i++)
        for(int j = 0; j < Ny; j++)
        for(int k = 0; k < Nz; k++)
        {
            total[comp] += MoleFractionsTotal(i,j,k)({comp});
        }
        total[comp] /= double(Nx*Ny*Nz);
        if (AtStart)
        {
            TotInitial[comp] = total[comp];
        }
        deviation[comp] = TotInitial[comp] - total[comp];
    }

    AtStart = false;

    ofstream output_file;
    if (tStep == 0)
    {
        output_file.open(TextDir + "CompositionStatistics.opd", ios::out);
        output_file << "#   sim_time" << "   ";
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            output_file << "total_" << ElementNames[comp] << "   "
                        << "deviation_" << ElementNames[comp] << "   ";
        }
        output_file<< endl;
        output_file.close();
    }

    output_file.open(TextDir + "CompositionStatistics.opd", ios::app);
    output_file << sim_time << "   ";
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        output_file << MoleFractionsTotalAverage({comp})  << "   " << deviation[comp] << "   ";
    }
    output_file << endl;
    output_file.close();
}

void Composition::Advect(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi, const BoundaryConditions& BC, const double dt, const double tStep)
{
    Adv.AdvectField(MoleFractions, Vel, BC, dx, dt, tStep);
    Adv.AdvectField(MoleFractionsTotal, Vel, BC, dx, dt, tStep);
    //CalculateTotalMoleFractions(Phi);
}

}// namespace openphase
