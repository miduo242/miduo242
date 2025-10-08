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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitri Medvedev
 *
 */

#include "Velocities.h"
#include "Info.h"
#include "Settings.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"
#include "Base/UserInterface.h"
#include "VTK.h"
#include "Mechanics/ElasticProperties.h"
#include "RunTimeControl.h"

namespace openphase
{
using namespace std;
/*************************************************************************/

Velocities::Velocities(Settings& locSettings)
{
    Initialize(locSettings);
}

void Velocities::Initialize(Settings& locSettings)
{
    thisclassname = "Velocities";

    Nx      = locSettings.Nx;
    Ny      = locSettings.Ny;
    Nz      = locSettings.Nz;

    dNx     = locSettings.dNx;
    dNy     = locSettings.dNy;
    dNz     = locSettings.dNz;

    Nphases = locSettings.Nphases;
    dx      = locSettings.dx;

    size_t Bcells = locSettings.Bcells;
    Phase.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Nphases}, Bcells);
    Average.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);

    Clear();

    VTKDir = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void Velocities::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetXVector(Average);
    BC.SetYVector(Average);
    BC.SetZVector(Average);

    BC.SetXVector(Phase);
    BC.SetYVector(Phase);
    BC.SetZVector(Phase);
}

void Velocities::Clear()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Average,Average.Bcells(),)
    {
        Average(i,j,k).set_to_zero();
        for(size_t alpha = 0; alpha < Nphases; ++alpha)
        {
            Phase(i, j, k)({alpha}).set_to_zero();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void Velocities::SetAverage(dVector3& value)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Average,Average.Bcells(),)
    {
        Average(i,j,k) = value;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Velocities::SetAverage(ElasticProperties& EP, double dt)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.Displacements,EP.Displacements.Bcells(),)
    {
        Average(i,j,k) = EP.Displacements(i,j,k) / dt;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Velocities::SetAllPhases(dVector3& value)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase,Phase.Bcells(),)
    for(size_t alpha = 0; alpha < Nphases; ++alpha)
    {
        Phase(i, j, k)({alpha}) = value;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Velocities::PrescribePhaseVelocities(PhaseField& Phi)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase,Phase.Bcells(),)
    for(size_t alpha = 0; alpha < Nphases; ++alpha)
    {
        Phase(i, j, k)({alpha}) = Average(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Velocities::CalculateAverage(const PhaseField& Phi)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Average,Average.Bcells(),)
    {
        Average(i,j,k).set_to_zero();

        if(Phi.Interface(i,j,k))
        {
            for(auto alpha = Phi.Fields(i,j,k).cbegin();
                     alpha != Phi.Fields(i,j,k).cend(); ++alpha)
            {
                size_t locPindex = Phi.FieldsStatistics[alpha->index].Phase;
                Average(i,j,k) += Phase(i, j, k)({locPindex})*alpha->value;
            }
        }
        else
        {
            size_t pInd = Phi.FieldsStatistics[Phi.Fields(i,j,k).front().index].Phase;
            Average(i,j,k) = Phase(i,j,k)({pInd});
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Velocities::WriteVTK(int tStep, Settings& locSettings) const
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = UserInterface::MakeFileName(locSettings.VTKDir, thisclassname+"_", tStep, ".vts");
    ListOfFields.push_back((VTK::Field_t){"VelocityAverage_",  [this](int i,int j,int k){return Average(i,j,k);}});

    for (size_t n = 0; n < Nphases; n++)
    {
        ListOfFields.push_back((VTK::Field_t){"VelocityPhase_" + std::to_string(n), [n,this](int i,int j,int k){return Phase(i,j,k)({n});}});
    }
    VTK::Write(Filename, locSettings, ListOfFields);
}

void Velocities::Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)
{
    Phase.Remesh(newNx, newNy, newNz);
    Average.Remesh(newNx, newNy, newNz);

    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    SetBoundaryConditions(BC);
    Info::WriteStandard(thisclassname, "Remeshed");
}

void Velocities::MoveFrame(int dx, int dy, int dz, BoundaryConditions& BC)
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
        Phase(i, j, k) = Phase(i + dx, j + dy, k + dz);
        Average(i, j, k) = Average(i + dx, j + dy, k + dz);
    }

    SetBoundaryConditions(BC);

    Info::WriteStandard(thisclassname, "Frame moved");

}

void Velocities::PrintPointStatistics(int x, int y, int z)
{
    stringstream message;
    message << "Point: " << x << " " << y << " " << z << endl;
    message << "Phase      Velocities " << endl;
    for (size_t alpha = 0; alpha < Nphases;  ++alpha)
    {
        message << alpha << ":    (" << Phase(x,y,z)({alpha})[0] << " "
                                     << Phase(x,y,z)({alpha})[1] << " "
                                     << Phase(x,y,z)({alpha})[2] << ")" << endl;
    }
    message << "Average Velocity: (" << Average(x,y,z)[0] << " "
                                     << Average(x,y,z)[1] << " "
                                     << Average(x,y,z)[2] << ")" << endl;
    Info::WriteSimple(message.str());
}

double Velocities::GetMaxVelocity()
{
    double MAXVelocity = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Average,0, reduction(max: MAXVelocity))
    {
        MAXVelocity = max(MAXVelocity, Average(i,j,k).abs());
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    double locMAXVelocity = MAXVelocity;
    MPI_Allreduce(&locMAXVelocity, &MAXVelocity, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    return MAXVelocity;
}

void Velocities::WriteStatistics(RunTimeControl& RTC)
{
    double MAXVelocity = GetMaxVelocity();

    ofstream output_file;
    if (RTC.tStep == 0)
    {
        output_file.open("VelocityStatistics.txt", ios::out);
        output_file << "tStep" << "\t\t\t" << "sim_time" << "\t\t\t"
                    << "max velocity" << endl;
        output_file.close();
    }

    output_file.open("VelocityStatistics.txt", ios::app);
    output_file << RTC.tStep << "\t\t\t" << RTC.SimulationTime << "\t\t\t" << MAXVelocity << endl;
    output_file.close();
}

Velocities& Velocities::operator= (const Velocities& rhs)
{
    // protect against invalid self-assignment and copy of unitialized object
    if (this != &rhs and rhs.thisclassname == "Velocities")
    {
        thisclassname = rhs.thisclassname;
        Nx = rhs.Nx;
        Ny = rhs.Ny;
        Nz = rhs.Nz;

        dNx = rhs.dNx;
        dNy = rhs.dNy;
        dNz = rhs.dNz;

        dx = rhs.dx;
        Nphases = rhs.Nphases;

        VTKDir = rhs.VTKDir;
        RawDataDir = rhs.RawDataDir;

        if (Phase.IsNotAllocated())
        {
            Phase.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Nphases}, rhs.Phase.Bcells());
            Average.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, rhs.Average.Bcells());
        }
        else if (not Phase.IsSize(rhs.Nx, rhs.Ny, rhs.Nz))
        {
            Phase.Reallocate(Nx, Ny, Nz);
            Average.Reallocate(Nx, Ny, Nz);
        }

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase,Phase.Bcells(),)
        {
            Average(i,j,k) = rhs.Average(i,j,k);
            for(size_t alpha = 0; alpha < Nphases; ++alpha)
            {
                Phase(i,j,k)({alpha}) = rhs.Phase(i,j,k)({alpha});
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    return *this;
}

void Velocities::Advect(BoundaryConditions& BC, double dt, AdvectionSchemes scheme)
{
    if(AverageDot.IsNotAllocated())
    {
        AverageDot.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, 1);
    }
    if(not AverageDot.IsSize(Nx, Ny, Nz))
    {
        AverageDot.Reallocate(Nx, Ny, Nz);
    }
    SetBoundaryConditions(BC);

    switch(scheme)
    {
        case AdvectionSchemes::Upwind:
        {
            const double dx2 = 0.5/dx;

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,AverageDot,0,)
            for (auto dir = 0; dir < 3; dir++)
            {
                AverageDot(i,j,k)[dir] = dx2 *
                     ((fabs(Average(i-1,j,k)[0]) + Average(i-1,j,k)[0])*Average(i-1,j,k)[dir] +
                      (fabs(Average(i,j-1,k)[1]) + Average(i,j-1,k)[1])*Average(i,j-1,k)[dir] +
                      (fabs(Average(i,j,k-1)[2]) + Average(i,j,k-1)[2])*Average(i,j,k-1)[dir] +
                      (fabs(Average(i+1,j,k)[0]) - Average(i+1,j,k)[0])*Average(i+1,j,k)[dir] +
                      (fabs(Average(i,j+1,k)[1]) - Average(i,j+1,k)[1])*Average(i,j+1,k)[dir] +
                      (fabs(Average(i,j,k+1)[2]) - Average(i,j,k+1)[2])*Average(i,j,k+1)[dir]) -
                      (fabs(Average(i,j,k)[0]) +
                       fabs(Average(i,j,k)[1]) +
                       fabs(Average(i,j,k)[2])) * Average(i, j, k)[dir]/dx;
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            break;
        }
        default:
        {
            Info::WriteExit("Wrong advection scheme selected",
                             thisclassname, "Velocities::Advect(BC, dt, scheme)");
            exit(13);
        }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,AverageDot,0,)
    for (int dir = 0; dir < 3; ++dir)
    {
        Average(i, j, k)[dir] += AverageDot(i, j, k)[dir]*dt;
        AverageDot(i, j, k)[dir] = 0.0;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    SetBoundaryConditions(BC);
}

void Velocities::Write(const int tStep) const
{
    std::string FileName = UserInterface::MakeFileName(RawDataDir,"Volocity_",tStep,".dat");

    std::ofstream out(FileName.c_str(), std::ios::out | std::ios::binary);

    if (!out)
    {
        std::stringstream message;
        message << "File \"" << FileName << "\" could not be created! Terminating!!!" << std::endl;
        throw std::runtime_error(message.str());
    };

    out.write(reinterpret_cast<const char*>(&Nx), sizeof(int));
    out.write(reinterpret_cast<const char*>(&Ny), sizeof(int));
    out.write(reinterpret_cast<const char*>(&Nz), sizeof(int));

    STORAGE_LOOP_BEGIN(i,j,k,Average,0)
    {
        out.write(reinterpret_cast<const char*>(&Average(i,j,k)[0]), sizeof(double));
        out.write(reinterpret_cast<const char*>(&Average(i,j,k)[1]), sizeof(double));
        out.write(reinterpret_cast<const char*>(&Average(i,j,k)[2]), sizeof(double));
    }
    STORAGE_LOOP_END

    STORAGE_LOOP_BEGIN(i,j,k,Phase,0)
    for(size_t n = 0; n < Nphases; n++)
    {
        out.write(reinterpret_cast<const char*>(&Phase(i,j,k)({n})[0]), sizeof(double));
        out.write(reinterpret_cast<const char*>(&Phase(i,j,k)({n})[1]), sizeof(double));
        out.write(reinterpret_cast<const char*>(&Phase(i,j,k)({n})[2]), sizeof(double));
    }
    STORAGE_LOOP_END

    out.close();
}

void Velocities::Read(BoundaryConditions& BC, const int tStep)
{
    std::string FileName = UserInterface::MakeFileName(RawDataDir,"Volocity_", tStep, ".dat");

    std::fstream inp(FileName.c_str(), std::ios::in | std::ios::binary);

    if (!inp)
    {
        std::stringstream message;
        message << "File \"" << FileName << "\" could not be opened";
        throw std::runtime_error(message.str());
    };

    int locNx = Nx;
    int locNy = Ny;
    int locNz = Nz;
    inp.read(reinterpret_cast<char*>(&Nx), sizeof(int));
    inp.read(reinterpret_cast<char*>(&Ny), sizeof(int));
    inp.read(reinterpret_cast<char*>(&Nz), sizeof(int));

    if(locNx != Nx or locNy != Ny or locNz != Nz)
    {
        std::stringstream message;
        message << "Inconsistent system dimensions!\n"
                << "Input data dimensions: (" << locNx
                << ", " << locNy << ", " << locNz << ") grid points.\n"
                << "Required data dimensions: (" << Nx
                << ", " << Ny << ", " << Nz << ") grid points.\n";
        throw std::runtime_error(message.str());
    }

    STORAGE_LOOP_BEGIN(i,j,k,Average,0)
    {
        inp.read(reinterpret_cast<char*>(&Average(i,j,k)[0]), sizeof(double));
        inp.read(reinterpret_cast<char*>(&Average(i,j,k)[1]), sizeof(double));
        inp.read(reinterpret_cast<char*>(&Average(i,j,k)[2]), sizeof(double));
    }
    STORAGE_LOOP_END

    STORAGE_LOOP_BEGIN(i,j,k,Phase,0)
    for(size_t n = 0; n < Nphases; n++)
    {
        inp.read(reinterpret_cast<char*>(&Phase(i,j,k)({n})[0]), sizeof(double));
        inp.read(reinterpret_cast<char*>(&Phase(i,j,k)({n})[1]), sizeof(double));
        inp.read(reinterpret_cast<char*>(&Phase(i,j,k)({n})[2]), sizeof(double));
    }
    STORAGE_LOOP_END

    inp.close();

    SetBoundaryConditions(BC);

    Info::WriteStandard(thisclassname, "Binary input loaded");
}


}// namespace openphase
