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
 *   File created :   2019
 *   Main contributors :   Marvin Tegeler; Raphael Schiedung
 *
 */

#include "Base/SolveLinearSystem.h"
#include "Base/SparseMatrix.h"
#include "Base/Tensor.h"
#include "Base/UserInterface.h"
#include "Info.h"
#include "Magnetism/LinearMagneticSolver.h"
#include "Mechanics/ElasticProperties.h"
#include "PhaseField.h"
#include "PhysicalConstants.h"
#include "VTK.h"

#include <algorithm>

namespace openphase
{

LinearMagneticSolver::LinearMagneticSolver(const Settings& locSettings, const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void LinearMagneticSolver::Initialize(const Settings& locSettings)
{
    thisclassname = "LinearMagneticSolver";

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;

    OffsetX = locSettings.OffsetX;

    dNx = locSettings.dNx;
    dNy = locSettings.dNy;
    dNz = locSettings.dNz;

    Nphases = locSettings.Nphases;

    size_t N = Nx*Ny*Nz;

    size_t Bcells = std::max(locSettings.Bcells, 2ul);
    phi.Allocate  (Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
    chi.Allocate  (Nx, Ny, Nz, dNx, dNy, dNz, Bcells);

    A.Allocate(N,7);
    b.resize(N);
    x.resize(N);

    dx = locSettings.dx;

    VTKDir     = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;
    TextDir    = locSettings.TextDir;

    initialized = true;
    Info::Write(thisclassname, "Initialized");
}

void LinearMagneticSolver::ReadInput(const std::string InputFileName)
{
    std::fstream inp(InputFileName.c_str(), std::ios::in);
    if (!inp)
    {
        Info::WriteExit("File " + InputFileName + " could not be opened",
                        thisclassname, "ReadInput()");
        exit(16);
    }

    Info::WriteBlankLine();
    Info::WriteLineInsert("LinearMagneticSolver");
    Info::WriteStandard("Source", InputFileName);

    std::stringstream data;
    data << inp.rdbuf();

    ReadInput(data);

    inp.close();
}

void LinearMagneticSolver::ReadInput(std::stringstream& inp)
{
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    H0x     = UserInterface::ReadParameterD(inp, moduleLocation,"H0X",    false, 0.0);
    H0y     = UserInterface::ReadParameterD(inp, moduleLocation,"H0Y",    false, 0.0);
    H0z     = UserInterface::ReadParameterD(inp, moduleLocation,"H0Z",    false, 0.0);
    dH0x_dx = UserInterface::ReadParameterD(inp, moduleLocation,"H0XX",   false, 0.0);
    dH0y_dx = UserInterface::ReadParameterD(inp, moduleLocation,"H0YX",   false, 0.0);
    dH0z_dx = UserInterface::ReadParameterD(inp, moduleLocation,"H0ZX",   false, 0.0);
    dH0x_dy = UserInterface::ReadParameterD(inp, moduleLocation,"H0XY",   false, 0.0);
    dH0y_dy = UserInterface::ReadParameterD(inp, moduleLocation,"H0YY",   false, 0.0);
    dH0z_dy = UserInterface::ReadParameterD(inp, moduleLocation,"H0ZY",   false, 0.0);
    dH0x_dz = UserInterface::ReadParameterD(inp, moduleLocation,"H0XZ",   false, 0.0);
    dH0y_dz = UserInterface::ReadParameterD(inp, moduleLocation,"H0YZ",   false, 0.0);
    dH0z_dz = UserInterface::ReadParameterD(inp, moduleLocation,"H0ZZ",   false, 0.0);
    //MaxResidual = UserInterface::ReadParameterD(inp, moduleLocation,"MaxRes", false, 1.0e-2);

    PhaseChi.resize(Nphases);
    for (size_t i = 0; i < Nphases; ++i)
    {
        std::stringstream converter;
        converter << i;
        std::string counter = converter.str();
        PhaseChi[i] = UserInterface::ReadParameterD(inp, moduleLocation, std::string("chi_") + counter);
    }

    Info::WriteLine();
    Info::WriteBlankLine();
}

void LinearMagneticSolver::SetEffectiveSusceptibility(
        const PhaseField& Phase, const BoundaryConditions& BC)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,chi,0)
    {
        chi(i,j,k) = 0.0;
        for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
        {
            chi(i,j,k) += PhaseChi[Phase.FieldsStatistics[it->index].Phase]*it->value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    BC.SetX(chi);
    BC.SetY(chi);
    BC.SetZ(chi);
}

void LinearMagneticSolver::Solve(const BoundaryConditions& BC, double MaxResidual)
{
#ifdef MPI_PARALLEL
    Info::WriteExit("MPI Parallelism not implemented !", "LinearMagneticSolver", "Solve");
    exit(EXIT_FAILURE);
#endif

    A.clear();
    for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j)
    for (int k = 0; k < Nz; ++k)
    {
        const double dchi_dx = (dNx) ? (chi(i+1,j,k)-chi(i-1,j,k))/2.0/dx : 0;
        const double dchi_dy = (dNy) ? (chi(i,j+1,k)-chi(i,j-1,k))/2.0/dx : 0;
        const double dchi_dz = (dNz) ? (chi(i,j,k+1)-chi(i,j,k-1))/2.0/dx : 0;

        const double mu_rel = (chi(i,j,k)+1);

        const double dphi0_dx = (dNx) ?  -H0x - (dH0x_dx*i + dH0x_dy*j + dH0x_dz*k)*dx : 0;
        const double dphi0_dy = (dNy) ?  -H0y - (dH0y_dx*i + dH0y_dy*j + dH0y_dz*k)*dx : 0;
        const double dphi0_dz = (dNz) ?  -H0z - (dH0z_dx*i + dH0z_dy*j + dH0z_dz*k)*dx : 0;
        const double Laplacephi0 = -dH0x_dx - dH0y_dy - dH0z_dz;

        x[idx(i,j,k)] = 0.0;
        b[idx(i,j,k)] = -dchi_dx*dphi0_dx - dchi_dy*dphi0_dy - dchi_dz*dphi0_dz - mu_rel*Laplacephi0;
        A(idx(i,j,k),idx(i,j,k)) = -6.0*mu_rel/dx/dx;

        const double Axp =  dchi_dx/2.0/dx + mu_rel/dx/dx;
        const double Axm = -dchi_dx/2.0/dx + mu_rel/dx/dx;
        const double Ayp =  dchi_dy/2.0/dx + mu_rel/dx/dx;
        const double Aym = -dchi_dy/2.0/dx + mu_rel/dx/dx;
        const double Azp =  dchi_dz/2.0/dx + mu_rel/dx/dx;
        const double Azm = -dchi_dz/2.0/dx + mu_rel/dx/dx;

        if (i+1 < Nx) A(idx(i,j,k),idx(i+1,j,k)) = Axp;
        else switch (BC.BCNX)
        {
            case BoundaryConditionTypes::NoFlux:   { A(idx(i,j,k),idx(i-1,j,k)) += Axp; break; }
            case BoundaryConditionTypes::Periodic: { A(idx(i,j,k),idx(  0,j,k)) += Axp; break; }
            //case Free:     { A(idx(i,j,k),idx(Nx-1,j,k)) += 2*Axp;
            //                 A(idx(i,j,k),idx(Nx-2,j,k)) -=   Axp; break; }
            case BoundaryConditionTypes::Fixed: break; // Fixes potential to zero which results in an arbitrary but constant magnetic field along the boundary
            default: { throw std::invalid_argument("Boundary Condition not implemented!"); }
        }
        if (i >= 1 ) A(idx(i,j,k),idx(i-1,j,k)) = Axm;
        else switch (BC.BC0X)
        {
            case BoundaryConditionTypes::NoFlux:   { A(idx(i,j,k),idx( i+1,j,k)) += Axm; break; }
            case BoundaryConditionTypes::Periodic: { A(idx(i,j,k),idx(Nx-1,j,k)) += Axm; break; }
            //case Free:     { A(idx(i,j,k),idx(   0,j,k)) += 2*Axm;
            //                 A(idx(i,j,k),idx(   1,j,k)) -=   Axm; break; }
            case BoundaryConditionTypes::Fixed: break;   // Fixes potential to zero which results in an arbitrary but constant magnetic field along the boundary
            default: { throw std::invalid_argument("Boundary Condition not implemented!"); }
        }
        if (j+1 < Ny) A(idx(i,j,k),idx(i,j+1,k)) = Ayp;
        else switch (BC.BCNY)
        {
            case BoundaryConditionTypes::NoFlux:   { A(idx(i,j,k),idx(i, j-1,k)) += Ayp; break; }
            case BoundaryConditionTypes::Periodic: { A(idx(i,j,k),idx(i,   0,k)) += Ayp; break; }
            case BoundaryConditionTypes::Fixed: break; // Fixes potential to zero which results in an arbitrary but constant magnetic field along the boundary
            default: { throw std::invalid_argument("Boundary Condition not implemented!"); }
        }
        if (j >= 1 ) A(idx(i,j,k),idx(i,j-1,k)) = Aym;
        else switch (BC.BC0Y)
        {
            case BoundaryConditionTypes::NoFlux:   { A(idx(i,j,k),idx(i, j+1,k)) += Aym; break; }
            case BoundaryConditionTypes::Periodic: { A(idx(i,j,k),idx(i,Ny-1,k)) += Aym; break; }
            case BoundaryConditionTypes::Fixed: break; // Fixes potential to zero which results in an arbitrary but constant magnetic field along the boundary

            default: { throw std::invalid_argument("Boundary Condition not implemented!"); }
        }
        if (k+1 < Nz) A(idx(i,j,k),idx(i,j,k+1)) = Azp;
        else switch (BC.BCNZ)
        {
            case BoundaryConditionTypes::NoFlux:   { A(idx(i,j,k),idx(i,j, k-1)) += Azp; break; }
            case BoundaryConditionTypes::Periodic: { A(idx(i,j,k),idx(i,j,   0)) += Azp; break; }
            case BoundaryConditionTypes::Fixed: break; // Fixes potential to zero which results in an arbitrary but constant magnetic field along the boundary
            default: { throw std::invalid_argument("Boundary Condition not implemented!"); }
        }
        if (k >= 1 ) A(idx(i,j,k),idx(i,j,k-1)) = Azm;
        else switch (BC.BC0Z)
        {
            case BoundaryConditionTypes::NoFlux:   { A(idx(i,j,k),idx(i,j, k+1)) += Azm; break; }
            case BoundaryConditionTypes::Periodic: { A(idx(i,j,k),idx(i,j,Nz-1)) += Azm; break; }
            case BoundaryConditionTypes::Fixed: break; // Fixes potential to zero which results in an arbitrary but constant magnetic field along the boundary
            default: { throw std::invalid_argument("Boundary Condition not implemented!"); }
        }
    }

    //SolveLinearSystem::Pivote(A,b);
    SolveLinearSystem::FastPivote(A,b);
    SolveLinearSystem::PreconditionJacobi(A,b);
    SolveLinearSystem::BiconjugateGradient(A,x,b,MaxResidual);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,phi,0)
    {
        phi(i,j,k) = x[idx(i,j,k)];
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    BC.SetX(phi);
    BC.SetY(phi);
    BC.SetZ(phi);
}

void LinearMagneticSolver::CalcForceDensity(ElasticProperties& EP,
        const BoundaryConditions& BC) const
{
    if(EP.ConsiderExternalForces)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.ForceDensity,0)
        {
            EP.ForceDensity(i,j,k) = ForceDensity(i,j,k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        BC.SetX(EP.ForceDensity);
        BC.SetY(EP.ForceDensity);
        BC.SetZ(EP.ForceDensity);
    }
    else
    {
        Info::WriteExit("ForceDensity Storage not allocated!", "LinearMagneticSolver", "CalcForceDensity");
        exit(EXIT_FAILURE);
    }
}

void LinearMagneticSolver::WriteVTK(const int tStep, const Settings& locSettings)
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t){"H"           , [this](int i,int j,int k){return H           (i, j, k);}});
    ListOfFields.push_back((VTK::Field_t){"B"           , [this](int i,int j,int k){return B           (i, j, k);}});
    ListOfFields.push_back((VTK::Field_t){"M"           , [this](int i,int j,int k){return M           (i, j, k);}});
    ListOfFields.push_back((VTK::Field_t){"mu"          , [this](int i,int j,int k){return mu          (i, j, k);}});
    ListOfFields.push_back((VTK::Field_t){"phi"         , [this](int i,int j,int k){return phi         (i, j, k);}});
    ListOfFields.push_back((VTK::Field_t){"chi"         , [this](int i,int j,int k){return chi         (i, j, k);}});
    ListOfFields.push_back((VTK::Field_t){"ForceDensity", [this](int i,int j,int k){return ForceDensity(i, j, k);}});

    std::string Filename = UserInterface::MakeFileName(VTKDir, thisclassname, tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}
}
