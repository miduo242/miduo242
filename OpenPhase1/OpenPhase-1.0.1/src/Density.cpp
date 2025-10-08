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
 *   File created :   2014
 *   Main contributors :   Oleg Shchyglo;
 *
 */

#include "Density.h"
#include "Settings.h"
#include "Info.h"
#include "PhaseField.h"
#include "Composition.h"
#include "Temperature.h"
#include "BoundaryConditions.h"
#include "Base/UserInterface.h"
#include "VTK.h"

namespace openphase
{
using namespace std;

Density::Density(Settings& locSettings, const string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}
void Density::Initialize(Settings& locSettings)
{
    thisclassname = "Density";

    TotalNx = locSettings.TotalNx;
    OffsetX = locSettings.OffsetX;
    TotalNy = locSettings.TotalNy;
    OffsetY = locSettings.OffsetY;
    TotalNz = locSettings.TotalNz;
    OffsetZ = locSettings.OffsetZ;

    Nx       = locSettings.Nx;
    Ny       = locSettings.Ny;
    Nz       = locSettings.Nz;

    dNx      = locSettings.dNx;
    dNy      = locSettings.dNy;
    dNz      = locSettings.dNz;

    dx      = locSettings.dx;
    Nphases = locSettings.Nphases;

    size_t Bcells = locSettings.Bcells;
    Phase.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Nphases}, Bcells);
    Total.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);

    Initial.Allocate({Nphases});

    VTKDir = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void Density::ReadInput(const string InputFileName)
{
    Info::WriteLineInsert("Density input");
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

void Density::ReadInput(stringstream& inp)
{
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    for(size_t alpha = 0; alpha != Nphases; alpha++)
    {
        stringstream converter;
        converter << string("Rho0_") << alpha;

        Initial({alpha}) = UserInterface::ReadParameterD(inp, moduleLocation, converter.str());
    }

    Info::WriteLine();
    Info::WriteBlankLine();
}

void Density::Set(PhaseField& PF, Composition& Cx, Temperature& Tx)
{
	size_t Ncomp = Cx.Ncomp;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Total,0,)
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            Total(i,j,k) = 0.0;
        }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Total,0,)
        if(PF.Fields(i,j,k).flag)
        {
            for(size_t comp = 0; comp < Ncomp; comp++)
            for(size_t alpha = 0; alpha != Nphases; alpha++)
            {
                Total(i,j,k) += PF.Fractions(i,j,k)({alpha})*
                                Cx.MoleFractions(i,j,k)({alpha, comp});
            }
        }
        else
        {
        	size_t alpha = PF.FieldsStatistics[PF.Fields(i,j,k).front().index].Phase;
            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                Total(i,j,k) = Cx.MoleFractions(i,j,k)({alpha, comp});
            }
        }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Density::SetInitial(PhaseField& PF)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Total,0,)
    {
        Total(i,j,k) = 0.0;
        Tensor<double,1> locFractions = PF.Fractions(i,j,k);
        for(size_t alpha = 0; alpha != Nphases; alpha++)
        {
            Phase(i,j,k)({alpha}) = Initial({alpha});
            Total(i,j,k) += locFractions({alpha})*Initial({alpha});
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Density::Remesh(int newNx, int newNy, int newNz, BoundaryConditions& BC)
{
    Total.Remesh(newNx, newNy, newNz);
    Phase.Remesh(newNx, newNy, newNz);

    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    SetBoundaryConditions(BC);
    Info::WriteStandard(thisclassname, "Remeshed");
}

void Density::MoveFrame(int dx, int dy, int dz, BoundaryConditions& BC)
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

    SetBoundaryConditions(BC);

    for(int i = xBeg; ((dx >= 0) and (i <= xEnd)) or ((dx < 0) and (i >= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j <= yEnd)) or ((dy < 0) and (j >= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k <= zEnd)) or ((dz < 0) and (k >= zEnd)); k += zInc)
    {
        Phase(i,j,k) = Phase(i + dx, j + dy, k + dz);
        Total(i, j, k) = Total(i + dx, j + dy, k + dz);
    }

    SetBoundaryConditions(BC);

    Info::WriteStandard(thisclassname, "Frame moved.");
}

void Density::SetBoundaryConditions(BoundaryConditions& BC)
{
    BC.SetX(Phase);
    BC.SetY(Phase);
    BC.SetZ(Phase);

    BC.SetX(Total);
    BC.SetY(Total);
    BC.SetZ(Total);
}

void Density::Write(const int tStep)
{
#ifdef MPI_PARALLEL
    string FileName = UserInterface::MakeFileName(RawDataDir, thisclassname+"_"+std::to_string(MPI_RANK)+"_", tStep, ".dat");
#else
    string FileName = UserInterface::MakeFileName(RawDataDir, thisclassname+"_", tStep, ".dat");
#endif

    ofstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be created", thisclassname, "Write()");
        exit(1);
    };
    STORAGE_LOOP_BEGIN(i,j,k,Phase,Phase.Bcells())
        out.write(reinterpret_cast<char*>(Phase(i,j,k).data()), Phase(i,j,k).size()*sizeof(double));
    STORAGE_LOOP_END
    STORAGE_LOOP_BEGIN(i,j,k,Total,Total.Bcells())
        out.write(reinterpret_cast<char*>(&Total(i,j,k)), sizeof(double));
    STORAGE_LOOP_END
}

void Density::Read(const int tStep)
{
#ifdef MPI_PARALLEL
    string FileName = UserInterface::MakeFileName(RawDataDir, thisclassname+"_"+std::to_string(MPI_RANK)+"_", tStep, ".dat");
#else
    string FileName = UserInterface::MakeFileName(RawDataDir, thisclassname+"_", tStep, ".dat");
#endif

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be opened", thisclassname, "Read");
        exit(1);
    };

    STORAGE_LOOP_BEGIN(i,j,k,Phase,Phase.Bcells())
        inp.read(reinterpret_cast<char*>(Phase(i,j,k).data()), Phase(i,j,k).size()*sizeof(double));
    STORAGE_LOOP_END
    STORAGE_LOOP_BEGIN(i,j,k,Total,Total.Bcells())
        inp.read(reinterpret_cast<char*>(&Total(i,j,k)), sizeof(double));
    STORAGE_LOOP_END
    Info::WriteStandard(thisclassname, "Binary Input Read Successfully.");
}

void Density::WriteVTK(int tStep, const Settings& locSettings, const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = UserInterface::MakeFileName(locSettings.VTKDir, thisclassname+"_", tStep, ".vts");

    ListOfFields.push_back((VTK::Field_t){"Total",  [this](int i,int j,int k){return Total(i,j,k);}});
    for(size_t n = 0; n < Nphases; n++)
    {
        ListOfFields.push_back((VTK::Field_t){"Phase_" + std::to_string(n), [n,this](int i,int j,int k){return Phase(i,j,k)({n});}});
    }
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void Density::PrintPointStatistics(int x, int y, int z)
{
    cout << "Density:\t";
    cout << "Total:\t" << Total(x, y, z) << endl;

    for (size_t n = 0; n < Nphases; n++)
    {
        cout << "Phase " << n << ":\t" << Phase(x,y,z)({n}) << endl;
    }
}

Density& Density::operator=(const Density& rhs)
{
    // protect against invalid self-assignment and copy of unitialized object
    if (this != &rhs and rhs.thisclassname != "")
    {
        thisclassname = rhs.thisclassname;

        TotalNx = rhs.TotalNx;
        OffsetX = rhs.OffsetX;
        TotalNy = rhs.TotalNy;
        OffsetY = rhs.OffsetY;
        TotalNz = rhs.TotalNz;
        OffsetZ = rhs.OffsetZ;

        Nx      = rhs.Nx;
        Ny      = rhs.Ny;
        Nz      = rhs.Nz;

        dNx     = rhs.dNx;
        dNy     = rhs.dNy;
        dNz     = rhs.dNz;

        dx      = rhs.dx;

        Nphases = rhs.Nphases;

        if (Phase.IsNotAllocated())
        {
            Phase.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Nphases}, rhs.Phase.Bcells());
            Total.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, rhs.Total.Bcells());
            Initial.Allocate({Nphases});
        }
        else if (not Phase.IsSize(Nx, Ny, Nz))
        {
            Phase.Reallocate(Nx, Ny, Nz);
            Total.Reallocate(Nx, Ny, Nz);
            Initial.Reallocate({Nphases});
        }

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase,Phase.Bcells(),)
            Phase(i,j,k) = rhs.Phase(i,j,k);
            Total(i,j,k) = rhs.Total(i,j,k);
        OMP_PARALLEL_STORAGE_LOOP_END

        for (size_t n = 0; n < Nphases; n++)
        {
            Initial({n}) = rhs.Initial({n});
        }
    }
    return *this;
}

}// namespace openphase
