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
 *   Main contributors :   Oleg Shchyglo; Raphael Schiedung; Marvin Tegeler;
 *                         Matthias Stratmann
 *
 */

#include "AdvectionHR/AdvectionHR.h"
#include "Base/UserInterface.h"
#include "BoundaryConditions.h"
#include "Composition.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Temperature.h"
#include "VTK.h"
#include "Velocities.h"

namespace openphase
{
using namespace std;

Temperature::Temperature(Settings& locSettings, const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void Temperature::Initialize(Settings& locSettings)
{
    thisclassname = "Temperature";

    TotalNx = locSettings.TotalNx;
    OffsetX = locSettings.OffsetX;
    TotalNy = locSettings.TotalNy;
    OffsetY = locSettings.OffsetY;
    TotalNz = locSettings.TotalNz;
    OffsetZ = locSettings.OffsetZ;

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;

    dNx = locSettings.dNx;
    dNy = locSettings.dNy;
    dNz = locSettings.dNz;

    Nphases = locSettings.Nphases;
    dx = locSettings.dx;

    Tmin = 0.0;
    Tmax = 0.0;
    Tavg = 0.0;

    Tneg = false;

    ExtensionsActive = false;

    size_t Bcells = locSettings.Bcells;
    Tx   .Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
    TxOld.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);

    dT_dr.set_to_zero();
    r0.set_to_zero();

    VTKDir = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void Temperature::ReadInput(const string InputFileName)
{
    Info::WriteLineInsert("Temperature input");
    Info::WriteStandard("Source", InputFileName.c_str());

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };

    Info::WriteBlankLine();
    Info::WriteLineInsert("Temperature properties");
    Info::WriteStandard("Source", InputFileName.c_str());

    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);

    inp.close();
}

void Temperature::ReadInput(stringstream& inp)
{
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    T0 = UserInterface::ReadParameterD(inp, moduleLocation, string("T0"));

    r0[0] = UserInterface::ReadParameterD(inp, moduleLocation, string("R0X"), false, 0.0);
    r0[1] = UserInterface::ReadParameterD(inp, moduleLocation, string("R0Y"), false, 0.0);
    r0[2] = UserInterface::ReadParameterD(inp, moduleLocation, string("R0Z"), false, 0.0);

    dT_dr[0] = UserInterface::ReadParameterD(inp, moduleLocation, string("DT_DRX"), false, 0.0);
    dT_dr[1] = UserInterface::ReadParameterD(inp, moduleLocation, string("DT_DRY"), false, 0.0);
    dT_dr[2] = UserInterface::ReadParameterD(inp, moduleLocation, string("DT_DRZ"), false, 0.0);

    dT_dt = UserInterface::ReadParameterD(inp, moduleLocation, string("DT_Dt"), false, 0.0);

    Tneg = UserInterface::ReadParameterB(inp, moduleLocation,"Tneg", false, false);

    TBC0X = UserInterface::ReadParameterD(inp, moduleLocation, string("TBC0X"), false, 0.0);
    TBCNX = UserInterface::ReadParameterD(inp, moduleLocation, string("TBCNX"), false, 0.0);
    TBC0Y = UserInterface::ReadParameterD(inp, moduleLocation, string("TBC0Y"), false, 0.0);
    TBCNY = UserInterface::ReadParameterD(inp, moduleLocation, string("TBCNY"), false, 0.0);
    TBC0Z = UserInterface::ReadParameterD(inp, moduleLocation, string("TBC0Z"), false, 0.0);
    TBCNZ = UserInterface::ReadParameterD(inp, moduleLocation, string("TBCNZ"), false, 0.0);

    TSphere = UserInterface::ReadParameterD(inp, moduleLocation, string("TSphere"), false, 0.0);
    Pr      = UserInterface::ReadParameterD(inp, moduleLocation, string("Pr"), false, 0.0);

    int X0_size = UserInterface::ReadParameterI(inp, moduleLocation, "Extension_X0", false, 0);
    int XN_size = UserInterface::ReadParameterI(inp, moduleLocation, "Extension_XN", false, 0);
    int Y0_size = UserInterface::ReadParameterI(inp, moduleLocation, "Extension_Y0", false, 0);
    int YN_size = UserInterface::ReadParameterI(inp, moduleLocation, "Extension_YN", false, 0);
    int Z0_size = UserInterface::ReadParameterI(inp, moduleLocation, "Extension_Z0", false, 0);
    int ZN_size = UserInterface::ReadParameterI(inp, moduleLocation, "Extension_ZN", false, 0);

    if(X0_size > 0) {ExtensionX0.Initialize(X0_size,(iVector3){-1, 0, 0}); ExtensionsActive = true;};
    if(XN_size > 0) {ExtensionXN.Initialize(XN_size,(iVector3){ 1, 0, 0}); ExtensionsActive = true;};
    if(Y0_size > 0) {ExtensionY0.Initialize(Y0_size,(iVector3){ 0,-1, 0}); ExtensionsActive = true;};
    if(YN_size > 0) {ExtensionYN.Initialize(YN_size,(iVector3){ 0, 1, 0}); ExtensionsActive = true;};
    if(Z0_size > 0) {ExtensionZ0.Initialize(Z0_size,(iVector3){ 0, 0,-1}); ExtensionsActive = true;};
    if(ZN_size > 0) {ExtensionZN.Initialize(ZN_size,(iVector3){ 0, 0, 1}); ExtensionsActive = true;};

    Info::WriteLine();
    Info::WriteBlankLine();
}

void Temperature::Set(const BoundaryConditions& BC, const PhaseField& Phase, const double dt)
{
    if(dT_dt != 0.0)
    {
        /* Applies equal temperature increment in all points. */

        IncrementWithValue(dt*dT_dt);
        SetMinMaxAvg();
        SetBoundaryConditions(BC);
    }
}

void Temperature::SetBoundaryConditions(const BoundaryConditions& BC)
{
    if (dNx > 0) BC.SetX(Tx);
    if (dNy > 0) BC.SetY(Tx);
    if (dNz > 0) BC.SetZ(Tx);
}

void Temperature::MoveFrame(const int dx, const int dy, const int dz,
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
         Tx(i,j,k) = Tx(i + dx*dNx, j + dy*dNy, k + dz*dNz);
    }
    SetBoundaryConditions(BC);

    Info::WriteStandard(thisclassname, "Frame moved");
}

void Temperature::ConsumePlane(const int dx, const int dy, const int dz,
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
        Tx(i,j,k) = Tx(i + dx*dNx, j + dy*dNy, k + dz*dNz);
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
        Tx(i,j,k) = Tx(i - dx*dNx, j - dy*dNy, k - dz*dNz);
    }
    SetBoundaryConditions(BC);

    Info::WriteStandard(thisclassname, "Plane consumed");
}

void Temperature::Write(const int tStep) //should be const
{
#ifdef MPI_PARALLEL
    string FileName = UserInterface::MakeFileName(RawDataDir,thisclassname+"_"+ std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = UserInterface::MakeFileName(RawDataDir,thisclassname+"_", tStep, ".dat");
#endif
    ofstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be created", thisclassname, "Write()");
        exit(1);
    };

    out.write(reinterpret_cast<char*>(&Tx[0]), Tx.tot_size()*sizeof(double));
    out.close();

    if(ExtensionsActive)
    {
        /* Write Extensions */
#ifdef MPI_PARALLEL
        string FileName2 = UserInterface::MakeFileName(RawDataDir,"TemperatureExtensions_"+ std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
        string FileName2 = UserInterface::MakeFileName(RawDataDir,"TemperatureExtensions_", tStep, ".dat");
#endif

        ofstream out2(FileName2.c_str(), ios::out | ios::binary);

        if (!out2)
        {
            Info::WriteExit("File \"" + FileName2 + "\" could not be created", thisclassname, "Write()");
            exit(1);
        };

        if(ExtensionX0.isActive()) ExtensionX0.write(out2);
        if(ExtensionXN.isActive()) ExtensionXN.write(out2);
        if(ExtensionY0.isActive()) ExtensionY0.write(out2);
        if(ExtensionYN.isActive()) ExtensionYN.write(out2);
        if(ExtensionZ0.isActive()) ExtensionZ0.write(out2);
        if(ExtensionZN.isActive()) ExtensionZN.write(out2);
    }
}

bool Temperature::Read(BoundaryConditions& BC, const int tStep)
{
#ifdef MPI_PARALLEL
    string FileName = UserInterface::MakeFileName(RawDataDir,thisclassname+"_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = UserInterface::MakeFileName(RawDataDir,thisclassname+"_", tStep, ".dat");
#endif

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteWarning("File \"" + FileName + "\" could not be opened", thisclassname, "Read()");
        return false;
    };

    inp.read(reinterpret_cast<char*>(&Tx[0]), Tx.tot_size()*sizeof(double));
    inp.close();

    SetMinMaxAvg();
    SetBoundaryConditions(BC);

    if(ExtensionsActive)
    {
        /* Read Extensions */
#ifdef MPI_PARALLEL
        string FileName2 = UserInterface::MakeFileName(RawDataDir,"TemperatureExtensions_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
        string FileName2 = UserInterface::MakeFileName(RawDataDir,"TemperatureExtensions_", tStep, ".dat");
#endif

        ifstream inp2(FileName2.c_str(), ios::in | ios::binary);

        if (!inp2)
        {
            Info::WriteWarning("File \"" + FileName2 + "\" could not be opened", thisclassname, "Read()");
            return false;
        };

        if(ExtensionX0.isActive()) ExtensionX0.read(inp2);
        if(ExtensionXN.isActive()) ExtensionXN.read(inp2);
        if(ExtensionY0.isActive()) ExtensionY0.read(inp2);
        if(ExtensionYN.isActive()) ExtensionYN.read(inp2);
        if(ExtensionZ0.isActive()) ExtensionZ0.read(inp2);
        if(ExtensionZN.isActive()) ExtensionZN.read(inp2);
    }
    Info::WriteStandard(thisclassname, "Binary input loaded");

    return true;
}

void Temperature::WriteH5(const int tStep, H5Interface& H5)
{
    #ifdef H5OP
    std::vector<double> dbuffer;
    dbuffer.push_back(Nx);
    dbuffer.push_back(Ny);
    dbuffer.push_back(Nz);
    H5.WriteCheckPoint(tStep, "TxDomain", dbuffer);
    dbuffer.clear();

    dbuffer = Tx.pack();
    H5.WriteCheckPoint(tStep, "Tx", dbuffer);


    if(ExtensionX0.isActive()) {
        H5.WriteCheckPoint(tStep, "TxExX0", ExtensionX0.Data);
    }
    if(ExtensionXN.isActive()){
        H5.WriteCheckPoint(tStep, "TxExXN", ExtensionXN.Data);
    }
    if(ExtensionY0.isActive()) {
        H5.WriteCheckPoint(tStep, "TxExY0", ExtensionY0.Data);
    }
    if(ExtensionYN.isActive()) {
        H5.WriteCheckPoint(tStep, "TxExYN", ExtensionYN.Data);
    }
    if(ExtensionZ0.isActive()) {
        H5.WriteCheckPoint(tStep, "TxExZ0", ExtensionZ0.Data);
    }
    if(ExtensionZN.isActive()) {
        H5.WriteCheckPoint(tStep, "TxExZN", ExtensionZN.Data);
    }

    #else
    std::cerr << "OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"" << std::endl;
    exit(5);
    #endif
}

bool Temperature::ReadH5(int tStep, H5Interface& H5)
{
    #ifdef H5OP
    std::vector<double> dbuffer;
    H5.ReadCheckPoint(tStep, "PFDomain", dbuffer);
    int locNx = dbuffer[0];
    int locNy = dbuffer[1];
    int locNz = dbuffer[2];
    dbuffer.clear();
    if(locNx != Nx or locNy != Ny or locNz != Nz)
    {
        stringstream message;
        message << "Inconsistent system dimensions!\n"
                << "Input data dimensions: (" << locNx << ", " << locNy << ", " << locNz << ") grid points.\n"
                << "Required data dimensions: (" << Nx << ", " << Ny << ", " << Nz << ") grid points.\n";
        Info::WriteWarning(message.str(), thisclassname, "Read()");
        return false;
    }

    H5.ReadCheckPoint(tStep, "Tx", dbuffer);
    Tx.unpack(dbuffer);

    if(ExtensionX0.isActive()) {
        H5.ReadCheckPoint(tStep, "TxExX0", ExtensionX0.Data);
    }
    if(ExtensionXN.isActive()){
        H5.ReadCheckPoint(tStep, "TxExXN", ExtensionXN.Data);
    }
    if(ExtensionY0.isActive()) {
        H5.ReadCheckPoint(tStep, "TxExY0", ExtensionY0.Data);
    }
    if(ExtensionYN.isActive()) {
        H5.ReadCheckPoint(tStep, "TxExYN", ExtensionYN.Data);
    }
    if(ExtensionZ0.isActive()) {
        H5.ReadCheckPoint(tStep, "TxExZ0", ExtensionZ0.Data);
    }
    if(ExtensionZN.isActive()) {
        H5.ReadCheckPoint(tStep, "TxExZN", ExtensionZN.Data);
    }

    return true;
    #else
    std::cerr << "OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"" << std::endl;
    exit(5);
    #endif
}

void Temperature::Remesh(const int newNx, const int newNy, const int newNz,
                                                    const BoundaryConditions& BC)
{
    Tx.Remesh(newNx, newNy, newNz);

    if(TxDot.IsAllocated())
    {
        TxDot.Reallocate(newNx, newNy, newNz);
    }

    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    SetBoundaryConditions(BC);
    Info::WriteStandard(thisclassname, "Remeshed");
}

void Temperature::SetToValue(double value)
{
    if(value < 0.0 and !Tneg)
    {
        stringstream message;
        message << "Negative temperature detected!\n"
                << "If it is an intended behavior use <$Tneg : Yes>\n"
                << "in the temperature input to allow negative temperatures.\n";
        Info::WriteExit(message.str(), thisclassname, "SetToValue()");
        exit(1);
    }
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,)
    {
        Tx(i,j,k) = value;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Temperature::SetToValueWithGrad(double value, dVector3 dT_dr, dVector3 r0)
{
    bool limit = false;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,Tx.Bcells(),reduction(or:limit))
    {
        Tx(i, j, k) = value + (dT_dr[0] * (i - r0[0]) + dT_dr[1] * (j - r0[1]) + dT_dr[2] * (k - r0[2]))*dx;

        if(Tx(i,j,k) < 0.0)
        {
            limit = true;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    if(limit and !Tneg)
    {
        stringstream message;
        message << "Negative temperature detected!\n"
                << "If it is an intended behavior use <$Tneg : Yes>\n"
                << "in the temperature input to allow negative temperatures.\n";
        Info::WriteExit(message.str(), thisclassname, "IncrementWithValue()");
        exit(1);
    }
}

void Temperature::IncrementWithValue(double value)
{
    bool limit = false;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,reduction(or:limit))
    {
        Tx(i,j,k) += value;

        if(Tx(i,j,k) < 0.0)
        {
            limit = true;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    if(limit and !Tneg)
    {
        stringstream message;
        message << "Negative temperature detected!\n"
                << "If it is an intended behavior use <$Tneg : Yes>\n"
                << "in the temperature input to allow negative temperatures.\n";
        Info::WriteExit(message.str(), thisclassname, "IncrementWithValue()");
        exit(1);
    }
}

void Temperature::SetMinMaxAvg()
{
    double locTmin = 6000.0;
    double locTmax =    0.0;
    double locTavg =    0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,reduction(min:locTmin) reduction(max:locTmax) reduction(+:locTavg))
    {
        locTmin = min(locTmin,Tx(i,j,k));
        locTmax = max(locTmax,Tx(i,j,k));
        locTavg += Tx(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Tmin = locTmin;
    Tmax = locTmax;
    Tavg = locTavg;

#ifdef MPI_PARALLEL
    MPI_Allreduce(&locTmin, &Tmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&locTmax, &Tmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&locTavg, &Tavg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    Tavg /= TotalNx*TotalNy*TotalNz;
}

void Temperature::SetInitial(const BoundaryConditions& BC)
{
    bool limit = false;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,Tx.Bcells(),reduction(or:limit))
    {
        Tx(i,j,k) = T0 + (dT_dr[0]*(i - r0[0]) +
                          dT_dr[1]*(j - r0[1]) +
                          dT_dr[2]*(k - r0[2]))*dx;

        if(Tx(i,j,k) < 0.0)
        {
            limit = true;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    //NX0
    if (TBC0X > 0.0)
    for(long int i = -Tx.BcellsX(); i < 0; i++)
    for(long int j = -Tx.BcellsY(); j < Tx.sizeY() + Tx.BcellsY(); j++)
    for(long int k = -Tx.BcellsZ(); k < Tx.sizeZ() + Tx.BcellsZ(); k++)
    {
        Tx(i,j,k) = TBC0X;
    }

    //NY0
    if (TBC0Y > 0.0)
    for(long int i = -Tx.BcellsX(); i < Tx.sizeX() + Tx.BcellsX(); i++)
    for(long int j = -Tx.BcellsY(); j < 0; j++)
    for(long int k = -Tx.BcellsZ(); k < Tx.sizeZ() + Tx.BcellsZ(); k++)
    {
        Tx(i,j,k) = TBC0Y;
    }

    //NZ0
    if (TBC0Z > 0.0)
    for(long int i = -Tx.BcellsX(); i < Tx.sizeX() + Tx.BcellsX(); i++)
    for(long int j = -Tx.BcellsY(); j < Tx.sizeY() + Tx.BcellsY(); j++)
    for(long int k = -Tx.BcellsZ(); k < 0; k++)
    {
        Tx(i,j,k) = TBC0Z;
    }

    //NXN
    if (TBCNX > 0.0)
    for(long int i =    Tx.sizeX(); i < Tx.sizeX() + Tx.BcellsX(); i++)
    for(long int j = -Tx.BcellsY(); j < Tx.sizeY() + Tx.BcellsY(); j++)
    for(long int k = -Tx.BcellsZ(); k < Tx.sizeZ() + Tx.BcellsZ(); k++)
    {
        Tx(i,j,k) = TBCNX;
    }

    //NYN
    if (TBCNY > 0.0)
    for(long int i = -Tx.BcellsX(); i < Tx.sizeX() + Tx.BcellsX(); i++)
    for(long int j =    Tx.sizeY(); j < Tx.sizeY() + Tx.BcellsY(); j++)
    for(long int k = -Tx.BcellsZ(); k < Tx.sizeZ() + Tx.BcellsZ(); k++)
    {
        Tx(i,j,k) = TBCNY;
    }

    //NZN
    if (TBCNZ > 0.0)
    for(long int i = -Tx.BcellsX(); i < Tx.sizeX() + Tx.BcellsX(); i++)
    for(long int j = -Tx.BcellsY(); j < Tx.sizeY() + Tx.BcellsY(); j++)
    for(long int k =    Tx.sizeZ(); k < Tx.sizeZ() + Tx.BcellsZ(); k++)
    {
        Tx(i,j,k) = TBCNZ;
    }

#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            BC.CommunicateX(Tx);
            BC.CommunicateY(Tx);
            BC.CommunicateZ(Tx);
        }
        else BC.Communicate(Tx);
#endif

    if(limit and !Tneg)
    {
        stringstream message;
        message << "Negative temperature detected!\n"
                << "If it is an intended behavior use <$Tneg : Yes>\n"
                << "in the temperature input to allow negative temperatures\n"
                << "or adjust initial temperature input parameters.\n";
        Info::WriteExit(message.str(), thisclassname, "SetInitial()");
        exit(1);
    }

    SetMinMaxAvg();
    SetBoundaryConditions(BC);

    if(ExtensionX0.isActive()) {ExtensionX0.SetInitial(*this);}
    if(ExtensionXN.isActive()) {ExtensionXN.SetInitial(*this);}
    if(ExtensionY0.isActive()) {ExtensionY0.SetInitial(*this);}
    if(ExtensionYN.isActive()) {ExtensionYN.SetInitial(*this);}
    if(ExtensionZ0.isActive()) {ExtensionZ0.SetInitial(*this);}
    if(ExtensionZN.isActive()) {ExtensionZN.SetInitial(*this);}
}

void Temperature::PrintPointStatistics(const int x, const int y, const int z) const
{
    Info::WriteStandard(thisclassname, "T(" + std::to_string(x) + ","
                                            + std::to_string(y) + ","
                                            + std::to_string(z) + ") = "
                                            + std::to_string(Tx(x, y, z)));
}
void Temperature::PrintStatistics() const
{
    std::string message  = "\n";
                message += Info::GetStandard("Tavg = ", std::to_string(Tavg));
                message += Info::GetStandard("Tmin = ", std::to_string(Tmin));
                message += Info::GetStandard("Tmax = ", std::to_string(Tmax));
    Info::WriteStandard(thisclassname, message);

}
void Temperature::WriteVTK(const int tStep, Settings& locSettings) const
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = UserInterface::MakeFileName(locSettings.VTKDir, thisclassname+"_", tStep, ".vts");
    ListOfFields.push_back((VTK::Field_t) {"T", [this](int i,int j,int k){return Tx(i,j,k);}});
    VTK::Write(Filename, locSettings, ListOfFields);
}
void Temperature::WriteGradientVTK(const int tStep, Settings& locSettings) const
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = UserInterface::MakeFileName(locSettings.VTKDir, thisclassname+"Gradient_", tStep, ".vts");
    ListOfFields.push_back((VTK::Field_t) {"dT_dx", [this](int i,int j,int k){return (dVector3){0.5*(Tx(i+1,j,k)-Tx(i-1,j,k))/dx,0.5*(Tx(i,j+1,k)-Tx(i,j-1,k))/dx,0.5*(Tx(i,j,k+1)-Tx(i,j,k-1))/dx};}});
    VTK::Write(Filename, locSettings, ListOfFields);
}

Temperature& Temperature::operator=(const Temperature& rhs)
{
    // protect against invalid self-assignment and copy of uninitialized object
    if (this != &rhs and rhs.thisclassname == "Temperature")
    {
        thisclassname = rhs.thisclassname;

        TotalNx = rhs.TotalNx;
        OffsetX = rhs.OffsetX;
        TotalNy = rhs.TotalNy;
        OffsetY = rhs.OffsetY;
        TotalNz = rhs.TotalNz;
        OffsetZ = rhs.OffsetZ;

        Nx = rhs.Nx;
        Ny = rhs.Ny;
        Nz = rhs.Nz;

        dNx = rhs.dNx;
        dNy = rhs.dNy;
        dNz = rhs.dNz;

        Nphases  = rhs.Nphases;

        dx = rhs.dx;

        dT_dr = rhs.dT_dr;
        r0 = rhs.r0;
        T0 = rhs.T0;
        dT_dt = rhs.dT_dt;

        Tmin = rhs.Tmin;
        Tmax = rhs.Tmax;
        Tavg = rhs.Tavg;
        Tneg = rhs.Tneg;

        Tx = rhs.Tx;
        TxDot = rhs.TxDot;

        VTKDir = rhs.VTKDir;
        RawDataDir = rhs.RawDataDir;

        ExtensionX0 = rhs.ExtensionX0;
        ExtensionXN = rhs.ExtensionXN;
        ExtensionY0 = rhs.ExtensionY0;
        ExtensionYN = rhs.ExtensionYN;
        ExtensionZ0 = rhs.ExtensionZ0;
        ExtensionZN = rhs.ExtensionZN;

        ExtensionsActive = rhs.ExtensionsActive;
    }
    return *this;
}

void Temperature::Advect(AdvectionHR& Adv, const Velocities& Vel,
                         PhaseField& Phi, const BoundaryConditions& BC,
                         const double dt, const double tStep)
{
    Adv.AdvectField(Tx, Vel, BC, dx, dt, tStep);
}

double Temperature::Average(void)
{
    double result = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(x,y,z,Tx,0,reduction(+:result))
    {
        result += Tx(x,y,z);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    result /= (Nx*Ny*Nz);

    return result;
}

double Temperature::Average(PhaseField& Phi, size_t n, size_t m)
{
    double result = 0.0;
    int interfacepoints = 0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(x,y,z,Tx,0,reduction(+:result,interfacepoints))
    {
        if(Phi.Interface(x,y,z))
        {
            bool valid = false;

            for (auto alpha = Phi.Fields(x,y,z).cbegin();
                      alpha != Phi.Fields(x,y,z).cend(); ++alpha)
            for (auto beta  = Phi.Fields(x,y,z).cbegin();
                      beta  != Phi.Fields(x,y,z).cend(); ++beta)
            {
                size_t PIdxA = Phi.FieldsStatistics[alpha->index].Phase;
                size_t PIdxB = Phi.FieldsStatistics[beta->index].Phase;

                if((PIdxA == n and PIdxB == m) or (PIdxA == m and PIdxB == n))
                {
                    valid = true;
                }
            }

            if(valid)
            {
                interfacepoints++;
                result += Tx(x,y,z);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    int rinterfacepoints = interfacepoints;
    MPI_Allreduce(&rinterfacepoints, &interfacepoints, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    double rresult = result;
    MPI_Allreduce(&rresult, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    if(interfacepoints > 0)
    {
        result /= double(interfacepoints);
    }
    else
    {
        result = Tavg;
    }

    return result;
}

}// namespace openphase
