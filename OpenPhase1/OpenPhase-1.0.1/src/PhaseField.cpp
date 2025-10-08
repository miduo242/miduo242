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
 *   File created :   2009
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev;
 *                         Reza Darvishi Kamachali; Raphael Schiedung;
 *                         Johannes Goerler; Marvin Tegeler
 *
 */

#include "Base/UserInterface.h"
#include "BoundaryConditions.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "Info.h"
#include "Mechanics/ElasticProperties.h"
#include "PhaseField.h"
#include "Settings.h"
#include "VTK.h"
#include "Velocities.h"
#include "AdvectionHR/AdvectionHR.h"
#include "H5Interface.h"

namespace openphase
{
using namespace std;

void PhaseField::Initialize(Settings& locSettings)
{
    thisclassname = "PhaseField";
    Resolution = locSettings.Resolution;
    NucleationPresent = false;

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

    dx  = locSettings.dx;
    iWidth = locSettings.iWidth;
    Eta = locSettings.Eta;

    Nphases = locSettings.Nphases;
    PhaseNames = locSettings.PhaseNames;
    PhaseAggregateStates = locSettings.PhaseAggregateStates;

    ConsiderNucleusVolume = locSettings.ConsiderNucleusVolume;

    // extra cells needed for driving force averaging
    size_t Bcells = max(locSettings.Bcells, size_t(iWidth));

    if(Resolution == Resolutions::Double)
    {
        Bcells = max(locSettings.Bcells, size_t(iWidth/2 + 1));
    }
    Fields   .Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
    FieldsDot.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
    Fractions.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, { Nphases }, Bcells);

    FractionsTotal.resize(Nphases, 0.0);

    if(Resolution == Resolutions::Double)
    {
        FieldsDR   .Allocate((1+dNx)*Nx, (1+dNy)*Ny, (1+dNz)*Nz, dNx, dNy, dNz, Bcells*2);
        FieldsDotDR.Allocate((1+dNx)*Nx, (1+dNy)*Ny, (1+dNz)*Nz, dNx, dNy, dNz, Bcells*2);
    }
    switch(locSettings.ActiveDimensions())
    {
        case 1:
        {
            LStencil.Set(LaplacianStencil1D_3, dx,dNx,dNy,dNz);
            GStencil.Set(GradientStencil1D, dx,dNx,dNy,dNz);
            break;
        }
        case 2:
        {
            switch(locSettings.PhaseFieldLaplacianStencil)
            {
                case LaplacianStencils::Simple:
                {
                    LStencil.Set(LaplacianStencil2D_5, dx,dNx,dNy,dNz);
                    break;
                }
                case LaplacianStencils::Isotropic:
                {
                    LStencil.Set(LaplacianStencil2D_9, dx,dNx,dNy,dNz);
                    break;
                }
                case LaplacianStencils::LB:
                {
                    LStencil.Set(LaplacianStencil2D_LB, dx,dNx,dNy,dNz);
                    break;
                }
            }

            switch(locSettings.PhaseFieldGradientStencil)
            {
                case GradientStencils::Simple:
                {
                    GStencil.Set(GradientStencil1D, dx,dNx,dNy,dNz);
                    break;
                }
                case GradientStencils::Isotropic:
                {
                    GStencil.Set(GradientStencil2D, dx,dNx,dNy,dNz);
                    break;
                }
                case GradientStencils::LB:
                {
                    GStencil.Set(GradientStencil2D_LB, dx,dNx,dNy,dNz);
                    break;
                }
            }
            break;
        }
        case 3:
        {
            switch(locSettings.PhaseFieldLaplacianStencil)
            {
                case LaplacianStencils::Simple:
                {
                    LStencil.Set(LaplacianStencil3D_7, dx);
                    break;
                }
                case LaplacianStencils::Isotropic:
                {
                    LStencil.Set(LaplacianStencil3D_27a, dx);
                    break;
                }
                case LaplacianStencils::LB:
                {
                    LStencil.Set(LaplacianStencil3D_LB, dx);
                    break;
                }
            }

            switch(locSettings.PhaseFieldGradientStencil)
            {
                case GradientStencils::Simple:
                {
                    GStencil.Set(GradientStencil1D, dx);
                    break;
                }
                case GradientStencils::Isotropic:
                {
                    GStencil.Set(GradientStencil3D, dx);
                    break;
                }
                case GradientStencils::LB:
                {
                    GStencil.Set(GradientStencil3D_LB, dx);
                    break;
                }
            }
            break;
        }
    }
    RefVolume = Pi;

    if(Nx < iWidth) RefVolume *= (Nx+1)/2.0;
    else RefVolume *= 1.1*(iWidth);
    if(Ny < iWidth) RefVolume *= (Ny+1)/2.0;
    else RefVolume *= 1.1*(iWidth);
    if(Nz < iWidth) RefVolume *= (Nz+1)/2.0;
    else RefVolume *= 1.1*(iWidth);

    VTKDir = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;

    FieldsStatistics.VTKDir = VTKDir;
    FieldsStatistics.RawDataDir = RawDataDir;
    FieldsStatistics.TextDir = locSettings.TextDir;

    GrainPairLimits.content.resize(FieldsStatistics.size());

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void PhaseField::Clear()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,Fields.Bcells(),)
    {
        Fields(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    if(Resolution == Resolutions::Double)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,FieldsDR.Bcells(),)
        {
            FieldsDR(i,j,k).clear();
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

void PhaseField::CalculateFractions()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, Fields.Bcells(), )
    {
        for(size_t n = 0; n < Nphases; n++)
        {
            Fractions(i,j,k)({n}) = 0.0;
        }
        for (auto it = Fields(i,j,k).cbegin();
                  it != Fields(i,j,k).cend(); ++it)
        {
            size_t pIndex = FieldsStatistics[it->index].Phase;
            Fractions(i, j, k)({pIndex}) += it->value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    double totalVolume = double(TotalNx * TotalNy * TotalNz);

    for(size_t n = 0; n < Nphases; n++)
    {
        FractionsTotal[n] = 0.0;
    }
    for(size_t idx = 0; idx < FieldsStatistics.size(); idx++)
    {
        FractionsTotal[FieldsStatistics[idx].Phase] += FieldsStatistics[idx].Volume;
    }
    for(size_t n = 0; n < Nphases; n++)
    {
        FractionsTotal[n] /= totalVolume;
    }
}

void PhaseField::CalculateDerivatives(void)
{
    Storage3D< NodePF, 0 > FieldsCopy(Fields);

    const int offset = Fields.Bcells() - 1;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, offset, )
    {
        if (Fields(i, j, k).flag)
        {
            for (auto& Field : Fields(i, j, k))
            {
                Field.laplacian  = 0.0;
                Field.gradient.set_to_zero();
            }

            for (auto ls = LStencil.begin(); ls != LStencil.end(); ls++)
            {
                int ii = ls->di;
                int jj = ls->dj;
                int kk = ls->dk;

                for (auto it = FieldsCopy(i + ii, j + jj, k + kk).cbegin();
                          it != FieldsCopy(i + ii, j + jj, k + kk).cend(); ++it)
                if (it->value != 0.0)
                {
                    double laplacian = ls->weight * it->value;
                    Fields(i, j, k).add_laplacian(it->index, laplacian);
                }
            }
            for (auto gs = GStencil.cbegin(); gs != GStencil.cend(); ++gs)
            {
                int ii = gs->di;
                int jj = gs->dj;
                int kk = gs->dk;

                for (auto it = FieldsCopy(i + ii, j + jj, k + kk).cbegin();
                          it != FieldsCopy(i + ii, j + jj, k + kk).cend(); ++it)
                if (it->value != 0.0)
                {
                    double value_x = gs->weightX * it->value;
                    double value_y = gs->weightY * it->value;
                    double value_z = gs->weightZ * it->value;
                    Fields(i, j, k).add_gradient(it->index, (dVector3){value_x,value_y,value_z});
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::CalculateDerivativesDR(void)
{
    Storage3D< NodePF, 0 > FieldsDRCopy(FieldsDR);

    const int offset = FieldsDR.Bcells() - 1;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, FieldsDR, offset, )
    {
        if (FieldsDR(i, j, k).flag)
        {
            for (auto& Field : FieldsDR(i, j, k))
            {
                Field.laplacian  = 0.0;
                Field.gradient.set_to_zero();
            }

            for (auto ls = LStencil.begin(); ls != LStencil.end(); ls++)
            {
                int ii = ls->di;
                int jj = ls->dj;
                int kk = ls->dk;

                for (auto it = FieldsDRCopy(i + ii, j + jj, k + kk).cbegin();
                          it != FieldsDRCopy(i + ii, j + jj, k + kk).cend(); ++it)
                if (it->value != 0.0)
                {
                    double laplacian = 4.0*ls->weight * it->value;
                    FieldsDR(i, j, k).add_laplacian(it->index, laplacian);
                }
            }
            for (auto gs = GStencil.cbegin(); gs != GStencil.cend(); ++gs)
            {
                int ii = gs->di;
                int jj = gs->dj;
                int kk = gs->dk;

                for (auto it = FieldsDRCopy(i + ii, j + jj, k + kk).cbegin();
                          it != FieldsDRCopy(i + ii, j + jj, k + kk).cend(); ++it)
                if (it->value != 0.0)
                {
                    double value_x = 2.0*gs->weightX * it->value;
                    double value_y = 2.0*gs->weightY * it->value;
                    double value_z = 2.0*gs->weightZ * it->value;
                    FieldsDR(i, j, k).add_gradient(it->index, (dVector3){value_x,value_y,value_z});
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

dVector3 PhaseField::Normal(const int i, const int j, const int k, const size_t alpha, const size_t beta) const
{
    dVector3 locNormal = Fields(i,j,k).get_gradient(beta)*Fields(i,j,k).get_value(alpha) -
                         Fields(i,j,k).get_gradient(alpha)*Fields(i,j,k).get_value( beta);

    double Norm = sqrt(locNormal[0]*locNormal[0] +
                       locNormal[1]*locNormal[1] +
                       locNormal[2]*locNormal[2]);

    if (Norm > DBL_EPSILON)
    {
        Norm = 1.0/Norm;
        locNormal[0] *= Norm;
        locNormal[1] *= Norm;
        locNormal[2] *= Norm;
    }
    else
    {
        locNormal[0] = 0.0;
        locNormal[1] = 0.0;
        locNormal[2] = 0.0;
    }

    return locNormal;
}

dVector3 PhaseField::NormalDR(const int i, const int j, const int k, const size_t alpha, const size_t beta) const
{
    dVector3 locNormal = FieldsDR(i,j,k).get_gradient(beta)*FieldsDR(i,j,k).get_value(alpha) -
                         FieldsDR(i,j,k).get_gradient(alpha)*FieldsDR(i,j,k).get_value( beta);

    double Norm = sqrt(locNormal[0]*locNormal[0] +
                       locNormal[1]*locNormal[1] +
                       locNormal[2]*locNormal[2]);

    if (Norm > DBL_EPSILON)
    {
        Norm = 1.0/Norm;
        locNormal[0] *= Norm;
        locNormal[1] *= Norm;
        locNormal[2] *= Norm;
    }
    else
    {
        locNormal[0] = 0.0;
        locNormal[1] = 0.0;
        locNormal[2] = 0.0;
    }

    return locNormal;
}

NodeV3 PhaseField::Normals(const int i, const int j, const int k) const
{
    NodeV3 locNormals;

    for (auto alpha = Fields(i,j,k).cbegin();
              alpha != Fields(i,j,k).cend() - 1; ++alpha)
    for (auto beta = alpha + 1;
              beta != Fields(i,j,k).cend(); ++beta)
    {
        dVector3 value = beta->gradient*alpha->value -
                         alpha->gradient*beta->value;
        locNormals.set_asym(alpha->index, beta->index, value);
    }

    for (auto alpha = locNormals.begin();
              alpha != locNormals.end(); ++alpha)
    {
        double Norm = sqrt(alpha->X()*alpha->X() +
                           alpha->Y()*alpha->Y() +
                           alpha->Z()*alpha->Z());

        if (Norm > DBL_EPSILON)
        {
            Norm = 1.0/Norm;
            alpha->X() *= Norm;
            alpha->Y() *= Norm;
            alpha->Z() *= Norm;
        }
        else
        {
            alpha->X() = 0.0;
            alpha->Y() = 0.0;
            alpha->Z() = 0.0;
        }
    }
    return locNormals;
}

NodeV3 PhaseField::NormalsDR(const int i, const int j, const int k) const
{
    NodeV3 locNormals;

    for (auto alpha = FieldsDR(i,j,k).cbegin();
              alpha != FieldsDR(i,j,k).cend() - 1; ++alpha)
    for (auto beta = alpha + 1;
              beta != FieldsDR(i,j,k).cend(); ++beta)
    {
        dVector3 value = beta->gradient*alpha->value -
                         alpha->gradient*beta->value;
        locNormals.set_asym(alpha->index, beta->index, value);
    }

    for (auto alpha = locNormals.begin();
              alpha != locNormals.end(); ++alpha)
    {
        double Norm = sqrt(alpha->X()*alpha->X() +
                           alpha->Y()*alpha->Y() +
                           alpha->Z()*alpha->Z());

        if (Norm > DBL_EPSILON)
        {
            Norm = 1.0/Norm;
            alpha->X() *= Norm;
            alpha->Y() *= Norm;
            alpha->Z() *= Norm;
        }
        else
        {
            alpha->X() = 0.0;
            alpha->Y() = 0.0;
            alpha->Z() = 0.0;
        }
    }
    return locNormals;
}

NodeV3 PhaseField::NormalsPhase(const int i, const int j, const int k) const
{
    NodeV3 locNormals = Fields(i,j,k).get_gradients();

    for (auto alpha = locNormals.begin();
              alpha != locNormals.end(); ++alpha)
    {
        double Norm = sqrt(alpha->vector3[0]*alpha->vector3[0] +
                           alpha->vector3[1]*alpha->vector3[1] +
                           alpha->vector3[2]*alpha->vector3[2]);

        if (Norm > DBL_EPSILON)
        {
            Norm = 1.0/Norm;
            alpha->vector3[0] *= Norm;
            alpha->vector3[1] *= Norm;
            alpha->vector3[2] *= Norm;
        }
        else
        {
            alpha->vector3[0] = 0.0;
            alpha->vector3[1] = 0.0;
            alpha->vector3[2] = 0.0;
        }
    }
    return locNormals;
}

NodeV3 PhaseField::NormalsPhaseDR(const int i, const int j, const int k) const
{
    NodeV3 locNormals = FieldsDR(i,j,k).get_gradients();

    for (auto alpha = locNormals.begin();
              alpha != locNormals.end(); ++alpha)
    {
        double Norm = sqrt(alpha->X()*alpha->X() +
                           alpha->Y()*alpha->Y() +
                           alpha->Z()*alpha->Z());

        if (Norm > DBL_EPSILON)
        {
            Norm = 1.0/Norm;
            alpha->X() *= Norm;
            alpha->Y() *= Norm;
            alpha->Z() *= Norm;
        }
        else
        {
            alpha->X() = 0.0;
            alpha->Y() = 0.0;
            alpha->Z() = 0.0;
        }
    }
    return locNormals;
}

void PhaseField::CalculateVolumes(void)
{
    int Nthreads = 1;

    #ifdef _OPENMP
    Nthreads = omp_get_max_threads();
    #endif
    const size_t size = FieldsStatistics.size();
    vector<vector<double>> Volume(Nthreads);
    for(int t = 0; t < Nthreads; t++)
    {
        Volume[t].resize(size, 0.0);
    }
    // Calculate volumes in each OpenMP chunk
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        int thread = 0;

        #ifdef _OPENMP
        thread = omp_get_thread_num();
        #endif

        for(auto it = Fields(i,j,k).cbegin();
                 it != Fields(i,j,k).cend(); ++it)
        {
            Volume[thread][it->index] += it->value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    // Add volumes from different OpenMP chunks
    for(size_t idx = 0; idx < size; idx++)
    {
        FieldsStatistics[idx].Volume = 0.0;

        for(int t = 0; t < Nthreads; t++)
        {
            FieldsStatistics[idx].Volume += Volume[t][idx];
        }
    }

    // Update FieldsStatistics across MPI domains
#ifdef MPI_PARALLEL
    size_t loc_size = FieldsStatistics.size();
    size_t max_size = loc_size;

    MPI_Allreduce(&loc_size, &max_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    if (max_size > loc_size)
    {
        FieldsStatistics.Resize(max_size);
    }

    for(size_t idx = 0; idx < FieldsStatistics.size(); idx++)
    {
        double loc_volume = FieldsStatistics[idx].Volume;
        MPI_Allreduce(&loc_volume, &(FieldsStatistics[idx].Volume), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double loc_maxvolume = FieldsStatistics[idx].MAXVolume;
        MPI_Allreduce(&loc_maxvolume, &(FieldsStatistics[idx].MAXVolume), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        double loc_refvolume = FieldsStatistics[idx].RefVolume;
        MPI_Allreduce(&loc_refvolume, &(FieldsStatistics[idx].RefVolume), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        int loc_stage = FieldsStatistics[idx].Stage;
        MPI_Allreduce(&loc_stage, &(FieldsStatistics[idx].Stage), 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        size_t loc_variant = FieldsStatistics[idx].Variant;
        MPI_Allreduce(&loc_variant, &(FieldsStatistics[idx].Variant), 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        size_t loc_phase = FieldsStatistics[idx].Phase;
        MPI_Allreduce(&loc_phase, &(FieldsStatistics[idx].Phase), 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        //TODO: add other missing reductions
    }
#endif

    // Calculate MAXVolume and VolumeRatio for all phase fields
    for(size_t idx = 0; idx < FieldsStatistics.size(); idx++)
    {
        FieldsStatistics[idx].MAXVolume = max(FieldsStatistics[idx].Volume,
                                              FieldsStatistics[idx].MAXVolume);

        FieldsStatistics[idx].VolumeRatio = FieldsStatistics[idx].MAXVolume/
                                            FieldsStatistics[idx].RefVolume;
    }

    // Update statistics based on actual grains volume
    int NumberOfNuclei = 0;

    for(size_t idx = 0; idx < size; idx++)
    if(FieldsStatistics[idx].Exist && FieldsStatistics[idx].Volume <= 0.0)
    {
        FieldsStatistics[idx].Exist  = 0;
        FieldsStatistics[idx].Stage  = 0;
        FieldsStatistics[idx].Volume = 0.0;
        FieldsStatistics[idx].MAXVolume = 0.0;
        FieldsStatistics[idx].VolumeRatio = 1.0;
    }
    else if (FieldsStatistics[idx].Volume > 0.0)
    {
        FieldsStatistics[idx].Exist = 1;
        if(FieldsStatistics[idx].Stage == 2)
        {
            FieldsStatistics[idx].Stage = 1;
        }

        if(FieldsStatistics[idx].VolumeRatio > 1.0)
        {
            FieldsStatistics[idx].Stage = 0;
            FieldsStatistics[idx].VolumeRatio = 1.0;
        }
        NumberOfNuclei += FieldsStatistics[idx].Stage;
    }
    NucleationPresent = (NumberOfNuclei > 0);
}

void PhaseField::SetFlags(void)
{
    const int offset = Fields.Bcells()-1;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, offset,)
    {
        if(Fields(i,j,k).flag == 2)
        {
            for(int ii = -dNx; ii <= +dNx; ++ii)
            for(int jj = -dNy; jj <= +dNy; ++jj)
            for(int kk = -dNz; kk <= +dNz; ++kk)
            if(ii != 0 or jj != 0 or kk != 0)
            if(!(Fields(i+ii, j+jj, k+kk).flag))
            {
                Fields(i+ii, j+jj, k+kk).flag = 1;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::SetFlagsDR(void)
{
    const int offset = FieldsDR.Bcells()-1;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, FieldsDR, offset,)
    {
        if(FieldsDR(i,j,k).flag == 2)
        {
            for(int ii = -dNx; ii <= +dNx; ++ii)
            for(int jj = -dNy; jj <= +dNy; ++jj)
            for(int kk = -dNz; kk <= +dNz; ++kk)
            if(ii != 0 or jj != 0 or kk != 0)
            if(!(FieldsDR(i+ii, j+jj, k+kk).flag))
            {
                FieldsDR(i+ii, j+jj, k+kk).flag = 1;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::FinalizeSR(const BoundaryConditions& BC, bool finalize)
{
    if(finalize)
    {
        #ifdef MPI_PARALLEL
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,Fields.Bcells(),)
        #else
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
        #endif
        {
            if(Fields(i,j,k).flag)
            {
                Fields(i,j,k).finalize();
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    SetBoundaryConditions(BC);
    SetFlags();
    CalculateDerivatives();
    CalculateVolumes();
    CalculateFractions();
}

void PhaseField::FinalizeDR(const BoundaryConditions& BC, bool finalize)
{
    if(finalize)
    {
        #ifdef MPI_PARALLEL
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,FieldsDR.Bcells(),)
        #else
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,0,)
        #endif
        {
            if(FieldsDR(i,j,k).flag)
            {
                FieldsDR(i,j,k).finalize();
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    SetBoundaryConditionsDR(BC);
    SetFlagsDR();
    CalculateDerivativesDR();

    Coarsen();

    SetBoundaryConditions(BC);
    SetFlags();
    CalculateDerivatives();
    CalculateVolumes();
    CalculateFractions();
}

void PhaseField::FixSpreading(BoundaryConditions& BC, double cutoff)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        for (auto alpha = Fields(i,j,k).begin();
                  alpha != Fields(i,j,k).end(); alpha++)
        {
            if(alpha->value < cutoff and FieldsStatistics[alpha->index].Stage == 0)
            {
                alpha->value = 0.0;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Finalize(BC);
}

void PhaseField::Coarsen(void)
{
    double norm = 1.0/pow(2.0,dNx+dNy+dNz);
    long int fx = 1 + dNx;
    long int fy = 1 + dNy;
    long int fz = 1 + dNz;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        Fields(i,j,k).clear();

        for(int di = -dNx; di <= dNx; di+=2)
        for(int dj = -dNy; dj <= dNy; dj+=2)
        for(int dk = -dNz; dk <= dNz; dk+=2)
        {
            Fields(i,j,k) += FieldsDR(fx*i+(di+1)/2,fy*j+(dj+1)/2,fz*k+(dk+1)/2);
        }
        Fields(i,j,k) *= norm;
        Fields(i,j,k).finalize();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::CoarsenDot(void)
{
    double norm = 1.0/pow(2.0,dNx+dNy+dNz);
    long int fx = 1 + dNx;
    long int fy = 1 + dNy;
    long int fz = 1 + dNz;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDot,0,)
    {
        if(Fields(i,j,k).flag)
        {
            FieldsDot(i,j,k).clear();

            for(int di = -dNx; di <= dNx; di+=2)
            for(int dj = -dNy; dj <= dNy; dj+=2)
            for(int dk = -dNz; dk <= dNz; dk+=2)
            {
                FieldsDot(i,j,k).add_asym_pairs(FieldsDotDR(fx*i+(di+1)/2,fy*j+(dj+1)/2,fz*k+(dk+1)/2));
            }
            FieldsDot(i,j,k) *= norm;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::Refine(void)
{
    long int fx = 1 + dNx;
    long int fy = 1 + dNy;
    long int fz = 1 + dNz;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        //if(Fields(i,j,k).flag)
        for(int di = -dNx; di <= dNx; di+=2)
        for(int dj = -dNy; dj <= dNy; dj+=2)
        for(int dk = -dNz; dk <= dNz; dk+=2)
        {
            FieldsDR(fx*i+(di+1)/2,fy*j+(dj+1)/2,fz*k+(dk+1)/2) = Fields.at(i+di*0.25,j+dj*0.25,k+dk*0.25);
            FieldsDR(fx*i+(di+1)/2,fy*j+(dj+1)/2,fz*k+(dk+1)/2).finalize();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

NodePF PhaseField::Dot(const int i, const int j, const int k, const double dt) const
{
    NodePF value;
    if(Fields(i,j,k).flag)
    {
        NodePF OldFields = Fields(i,j,k);
        NodePF NewFields = Fields(i,j,k);
        NodePF tmp;

        for(auto psi = FieldsDot(i,j,k).cbegin(); psi != FieldsDot(i,j,k).cend(); ++psi)
        {
            if(psi->value1 != 0.0)
            {
                NewFields.add_value(psi->indexA,  psi->value1*dt);
                NewFields.add_value(psi->indexB, -psi->value1*dt);
                tmp.add_value(psi->indexA, 1);
                tmp.add_value(psi->indexB, 1);
            }
            if(psi->value2 != 0.0)
            {
                NewFields.add_value(psi->indexA,  psi->value2*dt);
                NewFields.add_value(psi->indexB, -psi->value2*dt);
                tmp.add_value(psi->indexA, 1);
                tmp.add_value(psi->indexB, 1);
            }
        }
        NewFields.finalize();

        for (auto alpha = tmp.cbegin(); alpha != tmp.end(); ++alpha)
        {
            value.add_value(alpha->index, (NewFields.get_value(alpha->index)-OldFields.get_value(alpha->index))/dt);
        }
    }
    return value;
}

NodePF PhaseField::Dot1(const int i, const int j, const int k, const double dt) const
{
    NodePF value;
    if(Fields(i,j,k).flag)
    {
        NodePF OldFields = Fields(i,j,k);
        NodePF NewFields = Fields(i,j,k);
        //NodePF tmp;

        for(auto psi = FieldsDot(i,j,k).cbegin(); psi != FieldsDot(i,j,k).cend(); ++psi)
        {
            if(psi->value1 != 0.0)
            {
                NewFields.add_value(psi->indexA,  psi->value1*dt);
                NewFields.add_value(psi->indexB, -psi->value1*dt);
                //tmp.add_value(psi->indexA, 1);
                //tmp.add_value(psi->indexB, 1);
            }
        }
        //NewFields.finalize();

        //for (auto alpha = tmp.cbegin(); alpha != tmp.end(); ++alpha)
        for(auto alpha = NewFields.cbegin(); alpha != NewFields.cend(); ++alpha)
        {
            value.add_value(alpha->index, (NewFields.get_value(alpha->index)-OldFields.get_value(alpha->index))/dt);
        }
    }
    return value;
}

NodePF PhaseField::Dot2(const int i, const int j, const int k, const double dt) const
{
    NodePF value;
    if(Fields(i,j,k).flag)
    {
        NodePF locDot1 = Dot1(i,j,k,dt);
        NodePF locDot  = Dot (i,j,k,dt);

        NodePF tmp;
        for (auto it : locDot1) if (it.value != 0.0) tmp.add_value(it.index, 1);
        for (auto it : locDot ) if (it.value != 0.0) tmp.add_value(it.index, 1);

        for (auto alpha = tmp.cbegin(); alpha != tmp.end(); ++alpha)
        {
            double locDot2Alpha = locDot.get_value(alpha->index)-locDot1.get_value(alpha->index);
            value.add_value(alpha->index,locDot2Alpha);
        }
    }
    return value;
}

void PhaseField::MergeIncrementsSR(const BoundaryConditions& BC,
                                           const double dt, const bool finalize)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        if(Fields(i,j,k).flag)
        {
            for(auto psi = FieldsDot(i,j,k).cbegin();
                     psi != FieldsDot(i,j,k).cend(); ++psi)
            {
                if(psi->value1 != 0.0)
                {
                    Fields(i,j,k).add_value(psi->indexA,  psi->value1 * dt);
                    Fields(i,j,k).add_value(psi->indexB, -psi->value1 * dt);
                }
                if(psi->value2 != 0.0)
                {
                    Fields(i,j,k).add_value(psi->indexA,  psi->value2 * dt);
                    Fields(i,j,k).add_value(psi->indexB, -psi->value2 * dt);
                }
            }

            FieldsDot(i,j,k).clear();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Finalize(BC, finalize);
}

void PhaseField::MergeIncrementsDR(const BoundaryConditions& BC,
                                           const double dt, const bool finalize)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,0,)
    {
        if(FieldsDR(i,j,k).flag)
        {
            for(auto psi = FieldsDotDR(i,j,k).cbegin();
                     psi != FieldsDotDR(i,j,k).cend(); ++psi)
            {
                if(psi->value1 != 0.0)
                {
                    FieldsDR(i,j,k).add_value(psi->indexA,  psi->value1 * dt);
                    FieldsDR(i,j,k).add_value(psi->indexB, -psi->value1 * dt);
                }
                if(psi->value2 != 0.0)
                {
                    FieldsDR(i,j,k).add_value(psi->indexA,  psi->value2 * dt);
                    FieldsDR(i,j,k).add_value(psi->indexB, -psi->value2 * dt);
                }
            }

            FieldsDotDR(i,j,k).clear();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Finalize(BC, finalize);
}

void PhaseField::NormalizeIncrementsSR(const BoundaryConditions& BC, const double dt)
{
    /** This function limits phase-field increments for all present phase-field
    pairs, so that the actual phase-field values are within their natural
    limits of 0.0 and 1.0.*/

    double precision = FLT_EPSILON;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    if (Fields(i,j,k).flag)
    {
        for(auto alpha = Fields(i,j,k).cbegin();
                 alpha != Fields(i,j,k).cend(); ++alpha)
        {
            double dPsiAlphaPos = 0.0;
            double dPsiAlphaNeg = 0.0;
            for(auto beta = Fields(i,j,k).cbegin();
                     beta != Fields(i,j,k).cend(); ++beta)
            if(alpha != beta)
            {
                double locdPsi = FieldsDot(i,j,k).get_asym1(alpha->index,
                                                           beta->index);
                if(locdPsi > 0.0) dPsiAlphaPos += locdPsi;
                if(locdPsi < 0.0) dPsiAlphaNeg += locdPsi;
            }
            /* Set out of bounds sets to zero */
            if((alpha->value == 0.0 and dPsiAlphaNeg < 0.0 and dPsiAlphaPos == 0.0)
            or (alpha->value == 1.0 and dPsiAlphaPos > 0.0 and dPsiAlphaNeg == 0.0))
            {
                for(auto beta = Fields(i,j,k).cbegin();
                         beta != Fields(i,j,k).cend(); ++beta)
                if(alpha != beta)
                {
                    FieldsDot(i,j,k).set_sym1(alpha->index, beta->index, 0.0);
                }
            }
        }
        for(auto alpha = FieldsDot(i,j,k).begin();
                 alpha != FieldsDot(i,j,k).end();)
        {
            /* Remove zero-sets from the storage */
            if(alpha->value1 == 0.0)
            {
                alpha = FieldsDot(i,j,k).erase(alpha);
            }
            else
            {
                ++alpha;
            }
        }
        /* Limit increments! This is done in a while loop to account for all
        existing pair-contributions.*/
        int number_of_iterations = 0;
        bool LimitingNeeded = true;
        while (LimitingNeeded)
        {
            number_of_iterations++;
            LimitingNeeded = false;
            NodeAB locIncrements;
            for(auto it = FieldsDot(i,j,k).cbegin();
                     it != FieldsDot(i,j,k).cend(); ++it)
            {
                /* Collect increments */
                if(it->value1 < 0.0)
                {
                    locIncrements.add_sym2(it->indexA,0, it->value1);
                    locIncrements.add_sym1(it->indexB,0, -it->value1);
                }
                if(it->value1 > 0.0)
                {
                    locIncrements.add_sym1(it->indexA,0, it->value1);
                    locIncrements.add_sym2(it->indexB,0, -it->value1);
                }
            }
            NodeAB locLimits;
            for(auto alpha = Fields(i,j,k).cbegin();
                     alpha != Fields(i,j,k).cend(); ++alpha)
            {
                /* Calculate limits */
                double posIncrement = locIncrements.get_sym1(alpha->index,0);
                double negIncrement = locIncrements.get_sym2(alpha->index,0);
                double newPFvalue = alpha->value+(posIncrement+negIncrement)*dt;
                locLimits.set_sym2(alpha->index,0, 1.0);
                locLimits.set_sym1(alpha->index,0, 1.0);
                if(newPFvalue < 0.0)
                {
                    double tmpLim = locLimits.get_sym2(alpha->index,0);
                    double tmpLim2 = min(tmpLim,-(alpha->value+posIncrement*dt)
                                                  /(negIncrement*dt));
                    locLimits.set_sym2(alpha->index,0, tmpLim2);
                }
                if(newPFvalue > 1.0)
                {
                    double tmpLim = locLimits.get_sym1(alpha->index,0);
                    double tmpLim2 = min(tmpLim,(1.0-(alpha->value+negIncrement
                                                 *dt))/(posIncrement*dt));
                    locLimits.set_sym1(alpha->index,0, tmpLim2);
                }
            }
            for(auto it = FieldsDot(i,j,k).begin();
                     it != FieldsDot(i,j,k).end(); ++it)
            {
                /* Limit increments */
                if(it->value1 < 0.0)
                {
                    double tmpLim = min(locLimits.get_sym2(it->indexA,0),
                                        locLimits.get_sym1(it->indexB,0));
                    it->value1 *= tmpLim;
                    if (tmpLim < 1.0)
                    LimitingNeeded = true;
                }
                if(it->value1 > 0.0)
                {
                    double tmpLim = min(locLimits.get_sym1(it->indexA,0),
                                        locLimits.get_sym2(it->indexB,0));
                    it->value1 *= tmpLim;
                    if (tmpLim < 1.0)
                    LimitingNeeded = true;
                }
            }
            /* Exit limiting loop, if no convergence after 24 iterations*/
            if (number_of_iterations > 24)
            LimitingNeeded = false;
        }
        for(auto alpha = FieldsDot(i,j,k).begin();
                 alpha != FieldsDot(i,j,k).end(); )
        {
            /* Clean up values which are too small */
            if(fabs(alpha->value1) < DBL_EPSILON)
            {
                alpha = FieldsDot(i,j,k).erase(alpha);
            }
            else
            {
                ++alpha;
            }
        }

        /* End plausibility check */
        NodePF locField = Fields(i,j,k);
        for(auto psi = FieldsDot(i,j,k).cbegin();
                 psi != FieldsDot(i,j,k).cend(); ++psi)
        if(psi->value1 != 0.0)
        {
            locField.add_value(psi->indexA,  psi->value1 * dt);
            locField.add_value(psi->indexB, -psi->value1 * dt);
        }
        for(auto it = locField.cbegin();
                 it != locField.cend(); ++it)
        if(it->value < -precision or it->value > 1.0 + precision)
        {
            double oldFields = Fields(i,j,k).get_value(it->index);
            string msg = "Normalizing of phase field increments failed in point ("
                       + to_string(i) + "," + to_string(j) + "," + to_string(k)
                       + "). " + to_string(Fields(i,j,k).size())
                       + " fields present. Grain "
                       + to_string(it->index) + " with a fields-value of "
                       + to_string(oldFields)
                       + " is incremented by "
                       + to_string(it->value-oldFields)
                       + ", which results in a fields-value of "
                       + to_string(it->value)
                       + ". This will result in undefined behavior!";
            Info::WriteWarning(msg,thisclassname,"NormalizeIncrements");
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    SetIncrementsBoundaryConditions(BC);

#ifdef DEBUG

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    if(Fields(i,j,k).flag)
    {
        double zero = 1.0;

        for(size_t n = 0; n < Nphases; n++)
        {
            double fraction = Fractions(i,j,k)({n});
            double newfraction = fraction;

            zero -= fraction;

            for(auto it = FieldsDot(i,j,k).begin();
                     it != FieldsDot(i,j,k).end(); ++it)
            {
                if((FieldsStatistics[it->indexA].Phase == n)
                and(FieldsStatistics[it->indexB].Phase != n))
                {
                    newfraction += it->value1*dt;
                }
                else if((FieldsStatistics[it->indexA].Phase != n)
                     and(FieldsStatistics[it->indexB].Phase == n))
                {
                    newfraction -= it->value1*dt;
                }
            }

            if((newfraction < -precision)
            or (newfraction >  precision + 1.0)
            or (fraction < -precision)
            or (fraction >  precision + 1.0))
            {
                string msg = "Redistribution of composition not possible. "
                             "Phase-field broke in point ("
                           + to_string(i) + "," + to_string(j) + ","
                           + to_string(k) + ") for phase " + to_string(n)
                           + ". Old fraction "
                           + to_string(fraction)
                           + " New fraction " + to_string(newfraction)
                           + ". This will break mass conservation!";

                Info::WriteWarning(msg,thisclassname,"NormalizeIncrements");
                //exit(1);
            }
        }

        if(zero < -precision or zero > precision)
        {
            string msg = "Redistribution of composition not possible. "
                         "Phase-field broke in point ("
                       + to_string(i) + "," + to_string(j) + ","
                       + to_string(k) + "). Phase fractions don't add up to "
                         "unity, difference is " + to_string(zero);

           Info::WriteWarning(msg,thisclassname,"NormalizeIncrements");
           exit(1);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
#endif
}

void PhaseField::NormalizeIncrementsDR(const BoundaryConditions& BC, const double dt)
{
    /** This function limits phase-field increments for all present phase-field
    pairs, so that the actual phase-field values are within their natural
    limits of 0.0 and 1.0.*/

    double precision = FLT_EPSILON;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,0,)
    if (FieldsDR(i,j,k).flag)
    {
        for(auto alpha = FieldsDR(i,j,k).cbegin();
                 alpha != FieldsDR(i,j,k).cend(); ++alpha)
        {
            double dPsiAlphaPos = 0.0;
            double dPsiAlphaNeg = 0.0;
            for(auto beta = FieldsDR(i,j,k).cbegin();
                     beta != FieldsDR(i,j,k).cend(); ++beta)
            if(alpha != beta)
            {
                double locdPsi = FieldsDotDR(i,j,k).get_asym1(alpha->index,
                                                             beta->index);
                if(locdPsi > 0.0) dPsiAlphaPos += locdPsi;
                if(locdPsi < 0.0) dPsiAlphaNeg += locdPsi;
            }
            /* Set out of bounds sets to zero */
            if((alpha->value == 0.0 and dPsiAlphaNeg < 0.0 and dPsiAlphaPos == 0.0)
            or (alpha->value == 1.0 and dPsiAlphaPos > 0.0 and dPsiAlphaNeg == 0.0))
            {
                for(auto beta = FieldsDR(i,j,k).cbegin();
                         beta != FieldsDR(i,j,k).cend(); ++beta)
                if(alpha != beta)
                {
                    FieldsDotDR(i,j,k).set_sym1(alpha->index, beta->index, 0.0);
                }
            }
        }
        for(auto alpha = FieldsDotDR(i,j,k).begin();
                 alpha != FieldsDotDR(i,j,k).end();)
        {
            /* Remove zero-sets from the storage */
            if(alpha->value1 == 0.0)
            {
                alpha = FieldsDotDR(i,j,k).erase(alpha);
            }
            else
            {
                ++alpha;
            }
        }
        /* Limit increments! This is done in a while loop, to acknowledge all
        existing pair-contributions.*/
        int number_of_iterations = 0;
        bool LimitingNeeded = true;
        while (LimitingNeeded)
        {
            number_of_iterations++;
            LimitingNeeded = false;
            NodeAB locIncrements;
            for(auto it = FieldsDotDR(i,j,k).cbegin();
                     it != FieldsDotDR(i,j,k).cend(); ++it)
            {
                /* Collect increments */
                if(it->value1 < 0.0)
                {
                    locIncrements.add_sym2(it->indexA,0, it->value1);
                    locIncrements.add_sym1(it->indexB,0, -it->value1);
                }
                if(it->value1 > 0.0)
                {
                    locIncrements.add_sym1(it->indexA,0, it->value1);
                    locIncrements.add_sym2(it->indexB,0, -it->value1);
                }
            }
            NodeAB locLimits;
            for(auto alpha = FieldsDR(i,j,k).cbegin();
                     alpha != FieldsDR(i,j,k).cend(); ++alpha)
            {
                /* Calculate limits */
                double posIncrement = locIncrements.get_sym1(alpha->index,0);
                double negIncrement = locIncrements.get_sym2(alpha->index,0);
                double newPFvalue = alpha->value+(posIncrement+negIncrement)*dt;
                locLimits.set_sym2(alpha->index,0,1.0);
                locLimits.set_sym1(alpha->index,0,1.0);
                if(newPFvalue < 0.0)
                {
                    double tmpLim = locLimits.get_sym2(alpha->index,0);
                    double tmpLim2 = min(tmpLim,-(alpha->value+posIncrement*dt)
                                                  /(negIncrement*dt));
                    locLimits.set_sym2(alpha->index,0, tmpLim2);
                }
                if(newPFvalue > 1.0)
                {
                    double tmpLim = locLimits.get_sym1(alpha->index,0);
                    double tmpLim2 = min(tmpLim,(1.0-(alpha->value+negIncrement
                                                 *dt))/(posIncrement*dt));
                    locLimits.set_sym1(alpha->index,0,tmpLim2);
                }
            }
            for(auto it = FieldsDotDR(i,j,k).begin();
                     it != FieldsDotDR(i,j,k).end(); ++it)
            {
                /* Limit increments */
                if(it->value1 < 0.0)
                {
                    double tmpLim= min(locLimits.get_sym2(it->indexA,0),
                                       locLimits.get_sym1(it->indexB,0));
                    it->value1 *= tmpLim;
                    if (tmpLim < 1.0)
                    LimitingNeeded = true;
                }
                if(it->value1 > 0.0)
                {
                    double tmpLim = min(locLimits.get_sym1(it->indexA,0),
                                        locLimits.get_sym2(it->indexB,0));
                    it->value1 *= tmpLim;
                    if (tmpLim < 1.0)
                    LimitingNeeded = true;
                }
            }
            /* Exit limiting loop, if no convergence after 24 iterations*/
            if (number_of_iterations > 24)
            LimitingNeeded = false;
        }
        for(auto alpha = FieldsDotDR(i,j,k).begin();
                 alpha != FieldsDotDR(i,j,k).end(); )
        {
            /* Clean up values which are too small */
            if(fabs(alpha->value1) < DBL_EPSILON)
            {
                alpha = FieldsDotDR(i,j,k).erase(alpha);
            }
            else
            {
                ++alpha;
            }
        }

        /* End plausibility check */
        NodePF locField = FieldsDR(i,j,k);
        for(auto psi = FieldsDotDR(i,j,k).cbegin();
                 psi != FieldsDotDR(i,j,k).cend(); ++psi)
        if(psi->value1 != 0.0)
        {
            locField.add_value(psi->indexA,  psi->value1 * dt);
            locField.add_value(psi->indexB, -psi->value1 * dt);
        }
        for(auto it = locField.cbegin();
                 it != locField.cend(); ++it)
        if(it->value < -precision or it->value > 1.0 + precision)
        {
            double oldFields = FieldsDR(i,j,k).get_value(it->index);
            string msg = "Normalizing of phase field increments failed in point ("
                       + to_string(i) + "," + to_string(j) + "," + to_string(k)
                       + "). " + to_string(FieldsDR(i,j,k).size())
                       + " fields present. Grain "
                       + to_string(it->index) + " with a fields-value of "
                       + to_string(oldFields)
                       + " is incremented by "
                       + to_string(it->value-oldFields)
                       + ", which results in a fields-value of "
                       + to_string(it->value)
                       + ". This will result in undefined behavior!";
            Info::WriteWarning(msg,thisclassname,"NormalizeIncrements");
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    SetIncrementsBoundaryConditionsDR(BC);
    CoarsenDot();
    SetIncrementsBoundaryConditions(BC);

#ifdef DEBUG

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    if(Fields(i,j,k).flag)
    {
        double zero = 1.0;

        for(size_t n = 0; n < Nphases; n++)
        {
            double fraction = Fractions(i,j,k)({n});
            double newfraction = fraction;

            zero -= fraction;

            for(auto it = FieldsDot(i,j,k).begin();
                     it != FieldsDot(i,j,k).end(); ++it)
            {
                if((FieldsStatistics[it->indexA].Phase == n)
                and(FieldsStatistics[it->indexB].Phase != n))
                {
                    newfraction += it->value1*dt;
                }
                else if((FieldsStatistics[it->indexA].Phase != n)
                     and(FieldsStatistics[it->indexB].Phase == n))
                {
                    newfraction -= it->value1*dt;
                }
            }

            if((newfraction < -precision*DBL_EPSILON)
            or (newfraction >  precision*DBL_EPSILON+1.0)
            or (fraction < -precision*DBL_EPSILON)
            or (fraction >  precision*DBL_EPSILON+1.0))
            {
                string msg = "Redistribution of composition not possible. "
                             "Phase-field broke in point ("
                           + to_string(i) + "," + to_string(j) + ","
                           + to_string(k) + ") for phase " + to_string(n)
                           + ". Old fraction "
                           + to_string(fraction)
                           + " New fraction " + to_string(newfraction)
                           + ". This will break mass conservation!";

                Info::WriteWarning(msg,thisclassname,"NormalizeIncrements");
                //exit(1);
            }
        }

        if(zero < -DBL_EPSILON*precision or zero > DBL_EPSILON*precision)
        {
            string msg = "Redistribution of composition not possible. "
                         "Phase-field broke in point ("
                       + to_string(i) + "," + to_string(j) + ","
                       + to_string(k) + "). Phase fractions don't add up to "
                         "unity, difference is " + to_string(zero);

           Info::WriteWarning(msg,thisclassname,"NormalizeIncrements");
           exit(1);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
#endif
}

void PhaseField::KeepPhaseFieldsVolume(void)
{
    Tensor<bool,2> AllowedTransitions({Nphases,Nphases});
    // Disallow all transitions
    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    {
        AllowedTransitions({n,m}) = false;
    }

    KeepPhaseVolume(AllowedTransitions);
}

void PhaseField::KeepPhaseVolume(Tensor<bool,2> AllowedTransitions)
{
    vector<vector<NodeAB>> PairwiseVolumeChange;

    size_t Nthreads = 1;
#ifdef _OPENMP
    Nthreads = omp_get_max_threads();
#endif
    PairwiseVolumeChange.resize(Nthreads);
    for(size_t th = 0; th < Nthreads; th++)
    {
        PairwiseVolumeChange[th].resize(FieldsStatistics.size());
    }

    // Collect pairwise increments and corresponding weights
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, 0, )
    {
        int thNum = 0;
#ifdef _OPENMP
        thNum = omp_get_thread_num();
#endif
        if(Fields(i,j,k).flag)
        for(auto alpha = Fields(i,j,k).cbegin();  alpha != Fields(i,j,k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1; beta != Fields(i,j,k).cend(); ++beta)
        {
            size_t pIndexA = FieldsStatistics[alpha->index].Phase;
            size_t pIndexB = FieldsStatistics[ beta->index].Phase;
            if(not AllowedTransitions({pIndexA, pIndexB}))
            {
                double locPhiDot = FieldsDot(i,j,k).get_asym1(alpha->index, beta->index);
                PairwiseVolumeChange[thNum][alpha->index].add_sym1(beta->index,0, locPhiDot);
                PairwiseVolumeChange[thNum][alpha->index].add_sym2(beta->index,0, sqrt(alpha->value*beta->value));
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    for(size_t th = 1; th < Nthreads; th++)
    for(size_t pf = 0; pf < PairwiseVolumeChange[0].size(); pf++)
    {
        PairwiseVolumeChange[0][pf].add_sym1(PairwiseVolumeChange[th][pf]);
        PairwiseVolumeChange[0][pf].add_sym2(PairwiseVolumeChange[th][pf]);
    }
    // Distribute pairwise increments using corresponding weights
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, 0, )
    {
        if(Fields(i,j,k).flag == 2)
        for(auto alpha = Fields(i,j,k).cbegin();  alpha != Fields(i,j,k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1; beta != Fields(i,j,k).cend(); ++beta)
        if(alpha->value*beta->value != 0.0)
        {
            size_t pIndexA = FieldsStatistics[alpha->index].Phase;
            size_t pIndexB = FieldsStatistics[ beta->index].Phase;
            if(not AllowedTransitions({pIndexA, pIndexB}))
            {
                if(PairwiseVolumeChange[0][alpha->index].get_sym2(beta->index,0) != 0.0)
                {
                    FieldsDot(i,j,k).add_asym1(alpha->index, beta->index,
                        -PairwiseVolumeChange[0][alpha->index].get_sym1(beta->index,0)*sqrt(alpha->value*beta->value)/
                         PairwiseVolumeChange[0][alpha->index].get_sym2(beta->index,0));
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::SetBoundaryConditions(const BoundaryConditions& BC)
{
    if(Resolution == Resolutions::Double)
    {
        SetBoundaryConditionsDR(BC);
    }
    if(dNx) BC.SetX(Fields);
    if(dNy) BC.SetY(Fields);
    if(dNz) BC.SetZ(Fields);
}

void PhaseField::SetBoundaryConditionsDR(const BoundaryConditions& BC)
{
    if(dNx) BC.SetX(FieldsDR);
    if(dNy) BC.SetY(FieldsDR);
    if(dNz) BC.SetZ(FieldsDR);
}

void PhaseField::SetIncrementsBoundaryConditions(const BoundaryConditions& BC)
{
    if(dNx) BC.SetX(FieldsDot);
    if(dNy) BC.SetY(FieldsDot);
    if(dNz) BC.SetZ(FieldsDot);
}

void PhaseField::SetIncrementsBoundaryConditionsDR(const BoundaryConditions& BC)
{
    if(dNx) BC.SetX(FieldsDotDR);
    if(dNy) BC.SetY(FieldsDotDR);
    if(dNz) BC.SetZ(FieldsDotDR);
}

void PhaseField::PrintPointStatistics(const int x, const int y, const int z) const
{
    cout << "Point: " << x << " " << y << " " << z << endl;
    cout << "Phase Fields: " << endl;
    cout << "Index: ";
    for (auto beta = Fields(x,y,z).cbegin();
              beta != Fields(x,y,z).cend();  ++beta)
    {
        cout << setw(12) << setfill(' ') << beta->index << " ";
    }
    cout << endl;
    cout << "Value: ";
    for (auto beta = Fields(x,y,z).cbegin();
              beta != Fields(x,y,z).cend();  ++beta)
    {
        cout.precision(6);
        cout << setw(12) << setfill(' ') << beta->value << " ";
    }
    cout << endl;
    cout << "Phase Fractions: " << endl;
    cout << "Index: ";
    for (size_t alpha = 0; alpha != Nphases; alpha++)
    {
        cout << setw(12) << setfill(' ') << alpha << " ";
    }
    cout << endl;
    cout << "Value: ";

    auto locFractions = Fractions(x,y,z);
    for (size_t alpha = 0; alpha != Nphases; alpha++)
    {
        cout.precision(6);
        cout << setw(12) << setfill(' ') << locFractions({alpha}) << " ";
    }
    cout << endl;
}

void PhaseField::WriteVTK(const int tStep, const Settings& locSettings, const bool CurvatureOutput,const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = UserInterface::MakeFileName(locSettings.VTKDir, thisclassname+"_", tStep, ".vts");
    switch (Resolution)
    {
        case Resolutions::Single:
        {
            ListOfFields.push_back((VTK::Field_t) {"Interfaces",  [this](int i,int j,int k){return Interfaces(i,j,k);}});
            ListOfFields.push_back((VTK::Field_t) {"Flags",       [this](int i,int j,int k){return Fields(i,j,k).flag;}});
            ListOfFields.push_back((VTK::Field_t) {"PhaseFields", [this](int i,int j,int k){return LocalIndex(i,j,k);}});
            for(size_t n = 0; n < Nphases; n++)
            {
                ListOfFields.push_back((VTK::Field_t) {"PhaseFraction_" + std::to_string(n), [n,this](int i,int j,int k){return Fractions(i,j,k)({n});}});
            }
            ListOfFields.push_back((VTK::Field_t) {"Junctions", [this](int i,int j,int k){return Junctions(i,j,k);}});
            ListOfFields.push_back((VTK::Field_t) {"Variants",  [this](int i,int j,int k){return Variants(i,j,k);}});
            if(CurvatureOutput)
            for (size_t n = 0; n < Nphases; n++)
            {
                ListOfFields.push_back((VTK::Field_t) {"Curvature_" + std::to_string(n), [n,this](int i,int j,int k){return Curvature(i,j,k,n);}});
                ListOfFields.push_back((VTK::Field_t) {"PrincipalCurvature_1_" + std::to_string(n), [n,this](int i,int j,int k){return PrincipalCurvatures(i,j,k,n)[0];}});
                ListOfFields.push_back((VTK::Field_t) {"PrincipalCurvature_2_" + std::to_string(n), [n,this](int i,int j,int k){return PrincipalCurvatures(i,j,k,n)[1];}});
            }
            VTK::Write(Filename, locSettings, ListOfFields, precision);
            break;
        }
        case Resolutions::Double:
        {
            ListOfFields.push_back((VTK::Field_t) {"Interfaces",  [this](int i,int j,int k){return InterfacesDR(i,j,k);}});
            ListOfFields.push_back((VTK::Field_t) {"Flags",       [this](int i,int j,int k){return FieldsDR(i,j,k).flag;}});
            ListOfFields.push_back((VTK::Field_t) {"PhaseFields", [this](int i,int j,int k){return LocalIndexDR(i,j,k);}});
            for(size_t n = 0; n < Nphases; n++)
            {
                ListOfFields.push_back((VTK::Field_t) {"PhaseFraction_" + std::to_string(n), [n,this](int i,int j,int k){
                double FractionsValue = 0.0;
                for (auto it = FieldsDR(i,j,k).cbegin();
                          it != FieldsDR(i,j,k).cend(); ++it)
                {
                    size_t pIndex = FieldsStatistics[it->index].Phase;
                    if(pIndex == n)
                    {
                        FractionsValue += it->value;
                    }
                }
                return FractionsValue;}});
            }
            ListOfFields.push_back((VTK::Field_t) {"Junctions", [this](int i,int j,int k){return JunctionsDR(i,j,k);}});
            ListOfFields.push_back((VTK::Field_t) {"Variants",  [this](int i,int j,int k){return VariantsDR(i,j,k);}});
            VTK::Write(Filename, locSettings, ListOfFields, precision, 2);
            break;
        }
    }
}

void PhaseField::WriteDistortedVTK(
        const int tStep,
        const Settings& locSettings,
        const ElasticProperties& EP,
        const bool CurvatureOutput,
        const int precision) const
{
    std::string Filename = UserInterface::MakeFileName(locSettings.VTKDir, thisclassname+"Distorted_", tStep, ".vts");
    std::vector<VTK::Field_t> ListOfFields;
    switch (Resolution)
    {
        case Resolutions::Single:
        {
            ListOfFields.push_back((VTK::Field_t) {"Interfaces",  [this](int i,int j,int k){return Interfaces(i,j,k);}});
            ListOfFields.push_back((VTK::Field_t) {"Flags",       [this](int i,int j,int k){return Fields(i,j,k).flag;}});
            ListOfFields.push_back((VTK::Field_t) {"PhaseFields", [this](int i,int j,int k){return LocalIndex(i,j,k);}});
            for(size_t n = 0; n < Nphases; n++)
            {
                ListOfFields.push_back((VTK::Field_t) {"PhaseFraction_" + std::to_string(n), [n,this](int i,int j,int k){return Fractions(i,j,k)({n});}});
            }
            ListOfFields.push_back((VTK::Field_t) {"Junctions", [this](int i,int j,int k){return Junctions(i,j,k);}});
            ListOfFields.push_back((VTK::Field_t) {"Variants",  [this](int i,int j,int k){return Variants(i,j,k);}});
            if(CurvatureOutput)
            for (size_t n = 0; n < Nphases; n++)
            {
                ListOfFields.push_back((VTK::Field_t) {"Curvature_" + std::to_string(n), [n,this](int i,int j,int k){return Curvature(i,j,k,n);}});
                ListOfFields.push_back((VTK::Field_t) {"PrincipalCurvature_1_" + std::to_string(n), [n,this](int i,int j,int k){return PrincipalCurvatures(i,j,k,n)[0];}});
                ListOfFields.push_back((VTK::Field_t) {"PrincipalCurvature_2_" + std::to_string(n), [n,this](int i,int j,int k){return PrincipalCurvatures(i,j,k,n)[1];}});
            }
            VTK::WriteDistorted(Filename, locSettings, EP, ListOfFields, precision);
            break;
        }
        case Resolutions::Double:
        {
            ListOfFields.push_back((VTK::Field_t) {"Interfaces",  [this](int i,int j,int k){return InterfacesDR(i,j,k);}});
            ListOfFields.push_back((VTK::Field_t) {"Flags",       [this](int i,int j,int k){return FieldsDR(i,j,k).flag;}});
            ListOfFields.push_back((VTK::Field_t) {"PhaseFields", [this](int i,int j,int k){return LocalIndexDR(i,j,k);}});
            for(size_t n = 0; n < Nphases; n++)
            {
                ListOfFields.push_back((VTK::Field_t) {"PhaseFraction_" + std::to_string(n), [n,this](int i,int j,int k){
                double FractionsValue = 0.0;
                for (auto it = FieldsDR(i,j,k).cbegin();
                          it != FieldsDR(i,j,k).cend(); ++it)
                {
                    size_t pIndex = FieldsStatistics[it->index].Phase;
                    if(pIndex == n)
                    {
                        FractionsValue += it->value;
                    }
                }
                return FractionsValue;}});
            }
            ListOfFields.push_back((VTK::Field_t) {"Junctions", [this](int i,int j,int k){return JunctionsDR(i,j,k);}});
            ListOfFields.push_back((VTK::Field_t) {"Variants",  [this](int i,int j,int k){return VariantsDR(i,j,k);}});
            VTK::WriteDistorted(Filename, locSettings, EP, ListOfFields, precision, 2);
            break;
        }
    }
}

void PhaseField::WriteLaplacianVTK(const int tStep, const Settings& locSettings,
                                   size_t PhiIndex, const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = UserInterface::MakeFileName(locSettings.VTKDir, "Laplacian_" + to_string(PhiIndex) + "_", tStep, ".vts");
    switch (Resolution)
    {
        case Resolutions::Single:
        {
            ListOfFields.push_back((VTK::Field_t) {"Laplacian_" + std::to_string(PhiIndex), [PhiIndex,this](int i,int j,int k){return Fields(i,j,k).get_laplacian(PhiIndex);}});
            VTK::Write(Filename, locSettings, ListOfFields, precision);
            break;
        }
        case Resolutions::Double:
        {
            ListOfFields.push_back((VTK::Field_t) {"Laplacian_" + std::to_string(PhiIndex), [PhiIndex,this](int i,int j,int k){return FieldsDR(i,j,k).get_laplacian(PhiIndex);}});
            VTK::Write(Filename, locSettings, ListOfFields, precision, 2);
            break;
        }
    }
}

void PhaseField::WriteIndividualPhaseFieldValuesVTK(const int tStep,
                            const Settings& locSettings,
                            const std::initializer_list<size_t> FieldIndices,
                            const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = UserInterface::MakeFileName(locSettings.VTKDir, "PhaseFieldValues_", tStep, ".vts");
    switch (Resolution)
    {
        case Resolutions::Single:
        {
            for (auto IteratorIndex = FieldIndices.begin();
                      IteratorIndex != FieldIndices.end(); IteratorIndex++)
            {
                ListOfFields.push_back((VTK::Field_t) {"FieldValue_" + std::to_string(*IteratorIndex), [IteratorIndex,this](int i,int j,int k){return Fields(i,j,k)[*IteratorIndex];}});
            }
            VTK::Write(Filename, locSettings, ListOfFields, precision);
            break;
        }
        case Resolutions::Double:
        {
            for (auto IteratorIndex = FieldIndices.begin();
                      IteratorIndex != FieldIndices.end(); IteratorIndex++)
            {
                ListOfFields.push_back((VTK::Field_t) {"FieldValue_" + std::to_string(*IteratorIndex), [IteratorIndex,this](int i,int j,int k){return FieldsDR(i,j,k)[*IteratorIndex];}});
            }
            VTK::Write(Filename, locSettings, ListOfFields, precision, 2);
            break;
        }
    }
}

void PhaseField::Write(const std::string& FileName) const
{
    fstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be opened",
                thisclassname, "Write()");
        exit(EXIT_FAILURE);
    };

    out.write(reinterpret_cast<const char*>(&Nx), sizeof(int));
    out.write(reinterpret_cast<const char*>(&Ny), sizeof(int));
    out.write(reinterpret_cast<const char*>(&Nz), sizeof(int));

    switch (Resolution)
    {
        case Resolutions::Single:
        {
            STORAGE_LOOP_BEGIN(i,j,k,Fields,0)
            {
                Fields(i,j,k).Write(out);
            }
            STORAGE_LOOP_END
            break;
        }
        case Resolutions::Double:
        {
            STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,0)
            {
                FieldsDR(i,j,k).Write(out);
            }
            STORAGE_LOOP_END
            break;
        }
    }
    out.close();
}

void PhaseField::Write(const int tStep) const
{
    #ifdef MPI_PARALLEL
    string FileName =
        UserInterface::MakeFileName(RawDataDir,thisclassname+"_"+ std::to_string(MPI_RANK) + "_", tStep, ".dat");
    #else
    string FileName =
        UserInterface::MakeFileName(RawDataDir,thisclassname+"_", tStep, ".dat");
    #endif
    Write(FileName);
    FieldsStatistics.Write(tStep);
}

void PhaseField::Write() const
{
    std::string FileName = RawDataDir+thisclassname+".dat";
    Write(FileName);
    FieldsStatistics.Write();
}

void PhaseField::WriteH5(const int tStep, H5Interface& H5)
{
    #ifdef H5OP
    std::vector<double> dbuffer;
    dbuffer.push_back(Nx);
    dbuffer.push_back(Ny);
    dbuffer.push_back(Nz);
    H5.WriteCheckPoint(tStep, "PFDomain", dbuffer);
    dbuffer.clear();
    switch (Resolution)
    {
        case Resolutions::Single:
        {
            dbuffer = Fields.pack();
            H5.WriteCheckPoint(tStep, "Fields", dbuffer);
            break;
        }
        case Resolutions::Double:
        {
            dbuffer = FieldsDR.pack();
            H5.WriteCheckPoint(tStep, "FieldsDR", dbuffer);
            break;
        }
    }
    FieldsStatistics.WriteH5(tStep,H5);
    #else
    std::cerr << "OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"" << std::endl;
    exit(5);
    #endif
}

bool PhaseField::ReadH5(const BoundaryConditions& BC, int tStep, H5Interface& H5)
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
    switch (Resolution)
    {
        case Resolutions::Single:
        {
            H5.ReadCheckPoint(tStep, "Fields", dbuffer);
            Fields.unpack(dbuffer);
            break;
        }
        case Resolutions::Double:
        {
            H5.ReadCheckPoint(tStep, "FieldsDR", dbuffer);
            FieldsDR.unpack(dbuffer);
            break;
        }
    }
    FieldsStatistics.ReadH5(tStep,H5);
    Finalize(BC);
    return true;
    #else
    std::cerr << "OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"" << std::endl;
    exit(5);
    #endif
}

bool PhaseField::Read(string FileName)
{
    ifstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteExit(FileName + " could not be opened",
                thisclassname, "Read()");
        return false;
    };

    int locNx = Nx;
    int locNy = Ny;
    int locNz = Nz;
    inp.read(reinterpret_cast<char*>(&locNx), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNy), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNz), sizeof(int));
    if(locNx != Nx or locNy != Ny or locNz != Nz)
    {
        stringstream message;
        message << "Inconsistent system dimensions!\n"
                << "Input data dimensions: (" << locNx
                << ", " << locNy << ", " << locNz << ") grid points.\n"
                << "Required data dimensions: (" << Nx
                << ", " << Ny << ", " << Nz << ") grid points.\n";
        Info::WriteExit(message.str(), thisclassname, "Read()");
        return false;
    }
    switch (Resolution)
    {
        case Resolutions::Single:
        {
            STORAGE_LOOP_BEGIN(i,j,k,Fields,0)
            {
                Fields(i,j,k).Read(inp);
            }
            STORAGE_LOOP_END
            break;
        }
        case Resolutions::Double:
        {
            STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,0)
            {
                FieldsDR(i,j,k).Read(inp);
            }
            STORAGE_LOOP_END
            break;
        }
    }
    inp.close();
    Info::WriteStandard(thisclassname, "Binary input loaded");
    return true;
}

bool PhaseField::Read(const BoundaryConditions& BC, int tStep, const bool finalize)
{
#ifdef MPI_PARALLEL
    string FileName =
        UserInterface::MakeFileName(RawDataDir,thisclassname+"_"+ std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName =
        UserInterface::MakeFileName(RawDataDir,thisclassname+"_", tStep, ".dat");
#endif

    bool phi_ok = Read(FileName);
    phi_ok = (phi_ok && FieldsStatistics.Read(tStep));
    if (phi_ok) Finalize(BC);
    return phi_ok;
}

bool PhaseField::Read(const BoundaryConditions& BC, const bool finalize)
{
    std::string FileName = RawDataDir+thisclassname + ".dat";

    bool phi_ok = Read(FileName);
    phi_ok = (phi_ok && FieldsStatistics.Read());
    if (phi_ok) Finalize(BC);
    return phi_ok;
}

void PhaseField::WriteAverageVolume(const int tStep, const size_t PhaseIndex) const
{
    stringstream converter;
    converter << PhaseIndex;

    string FileName = string("AverageVolumeOfPhase_")
        + converter.str() + string(".txt");

    if(!tStep)
    {
        fstream tout(FileName.c_str(), ios::out);
        tout << "Time\tAvgVolume" << endl;
        tout.close();
    }

    double avgVol = 0;
    size_t count = 0;
    for(size_t n = 0; n < FieldsStatistics.size(); n++)
    if(FieldsStatistics[n].Exist and FieldsStatistics[n].Phase == PhaseIndex)
    {
        avgVol += FieldsStatistics[n].Volume;
        count += 1.0;
    }
#ifdef MPI_PARALLEL
auto rVolume = avgVol;
MPI_Allreduce(&rVolume, &(avgVol), 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
size_t tmp = count;
MPI_Allreduce(&tmp, &(count), 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

    fstream out(FileName.c_str(), ios::out | ios::app);
    if(count)
    {
        out << tStep << "\t" << avgVol/count << endl;
    }
    else
    {
        out << tStep << "\t" << 0.0 << endl;
    }
    out.close();
}

void PhaseField::MoveFrame(const int dx, const int dy, const int dz,
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
        Fields(i, j, k) = Fields(i+dx,j+dy,k+dz);
    }
    Finalize(BC);
    Info::WriteStandard(thisclassname, "Frame moved");
}

void PhaseField::MoveFrameDR(const int dx, const int dy, const int dz,
                           const BoundaryConditions& BC)
{
    for(int n = 0; n <= 1; n++)
    {
        int xBeg = (dx >= 0) + (dx < 0)*(FieldsDR.sizeX()) - 1;
        int xEnd = (dx >= 0)*(FieldsDR.sizeX()) + (dx < 0) - 1;
        int xInc = 1 - 2*(dx < 0);

        int yBeg = (dy >= 0) + (dy < 0)*(FieldsDR.sizeY()) - 1;
        int yEnd = (dy >= 0)*(FieldsDR.sizeY()) + (dy < 0) - 1;
        int yInc = 1 - 2*(dy < 0);

        int zBeg = (dz >= 0) + (dz < 0)*(FieldsDR.sizeZ()) - 1;
        int zEnd = (dz >= 0)*(FieldsDR.sizeZ()) + (dz < 0) - 1;
        int zInc = 1 - 2*(dz < 0);

        for(int i = xBeg; ((dx >= 0) and (i <= xEnd)) or ((dx < 0) and (i >= xEnd)); i += xInc)
        for(int j = yBeg; ((dy >= 0) and (j <= yEnd)) or ((dy < 0) and (j >= yEnd)); j += yInc)
        for(int k = zBeg; ((dz >= 0) and (k <= zEnd)) or ((dz < 0) and (k >= zEnd)); k += zInc)
        {
            FieldsDR(i, j, k) = FieldsDR(i+dx,j+dy,k+dz);
        }
        FinalizeDR(BC);
    }
    Info::WriteStandard(thisclassname, "Frame moved");
}

void PhaseField::ConsumePlane(const int dx, const int dy, const int dz,
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
        Fields(i, j, k) = Fields(i+dx,j+dy,k+dz);
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
        Fields(i, j, k) = Fields(i-dx,j-dy,k-dz);
    }
    Finalize(BC);
    Info::WriteStandard(thisclassname, "Plane consumed");
}

void PhaseField::PrintPFVolumes() const
{
    Info::WriteLineInsert("Phase-field volumes","=");

    for(unsigned int idx = 0; idx < FieldsStatistics.size(); idx++)
    {
        if(FieldsStatistics[idx].Volume > 0.0)
        {

            Info::WriteStandardNarrow("PF", std::to_string(idx));
            Info::WriteStandardNarrow("Variant",
                    std::to_string(FieldsStatistics[idx].Variant));
            Info::WriteStandardNarrow("Stage",
                    std::to_string(FieldsStatistics[idx].Stage));
            Info::WriteStandardNarrow("Volume",
                    std::to_string(FieldsStatistics[idx].Volume));
            Info::WriteLine("-");
        }
    }
    Info::WriteLine("=");
}

void PhaseField::PrintVolumeFractions()
{
    Info::WriteSimple("Phase fractions:");
    for(size_t n = 0; n < Nphases; n++)
    {
        Info::WriteStandardNarrow("Phase " + std::to_string(n) + " (" + PhaseNames[n] + ")",
                std::to_string(100.0*FractionsTotal[n]) + " %" );
    }
}

void PhaseField::Remesh(int newNx, int newNy, int newNz,
        const BoundaryConditions& BC)
{
    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    Fields.Remesh(Nx, Ny, Nz);
    FieldsDot.Reallocate(Nx, Ny, Nz);
    Fractions.Reallocate(Nx, Ny, Nz);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        Fields(i,j,k).flag = 2*(Fields(i,j,k).size() > 1);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Finalize(BC);

    RefVolume = Pi;

    if(Nx < iWidth) RefVolume *= (Nx+1)/2.0;
    else RefVolume *= 1.1*iWidth;
    if(Ny < iWidth) RefVolume *= (Ny+1)/2.0;
    else RefVolume *= 1.1*iWidth;
    if(Nz < iWidth) RefVolume *= (Nz+1)/2.0;
    else RefVolume *= 1.1*iWidth;

    Info::WriteStandard(thisclassname, "Remeshed");
}

size_t PhaseField::PlantGrainNucleus(size_t PhaseIndex, int x, int y, int z)
{
    size_t locIndex = FieldsStatistics.add_nucleus(PhaseIndex);

    FieldsStatistics[locIndex].Rcm[0] = x;
    FieldsStatistics[locIndex].Rcm[1] = y;
    FieldsStatistics[locIndex].Rcm[2] = z;
    FieldsStatistics[locIndex].RefVolume = RefVolume;
    FieldsStatistics[locIndex].VolumeRatio = 0.0;
    FieldsStatistics[locIndex].State = PhaseAggregateStates[PhaseIndex];

    NucleationPresent = true;

    if (x - OffsetX >= 0 && x - OffsetX < Nx and
        y - OffsetY >= 0 && y - OffsetY < Ny and
        z - OffsetZ >= 0 && z - OffsetZ < Nz)
    {
        Fields(x - OffsetX, y - OffsetY, z - OffsetZ).set_value(locIndex, 0.0);

        long int fx = 1 + dNx;
        long int fy = 1 + dNy;
        long int fz = 1 + dNz;

        if(Resolution == Resolutions::Double)
        {
            for(int di = -dNx; di <= dNx; di+=2)
            for(int dj = -dNy; dj <= dNy; dj+=2)
            for(int dk = -dNz; dk <= dNz; dk+=2)
            {
                FieldsDR(fx*(x - OffsetX)+(di+1)/2,fy*(y - OffsetY)+(dj+1)/2,fz*(z - OffsetZ)+(dk+1)/2).set_value(locIndex, 0.0);
            }
        }
    }

    return locIndex;
}

size_t PhaseField::AddGrainInfo(size_t PhaseIndex)
{
    size_t locIndex = FieldsStatistics.add_grain(PhaseIndex);
    FieldsStatistics[locIndex].RefVolume = RefVolume;
    FieldsStatistics[locIndex].VolumeRatio = 0.0;
    FieldsStatistics[locIndex].State = PhaseAggregateStates[PhaseIndex];

    return locIndex;
}

PhaseField& PhaseField::operator= (const PhaseField& rhs)
{
    // protect against self-assignment and copy of uninitialized object
    if (this != &rhs and rhs.thisclassname == "PhaseField")
    {
        thisclassname = rhs.thisclassname;
        NucleationPresent = rhs.NucleationPresent;

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

        dx  = rhs.dx;
        Eta = rhs.Eta;
        Nphases = rhs.Nphases;
        iWidth = rhs.iWidth;
        RefVolume = rhs.RefVolume;
        PhaseNames = rhs.PhaseNames;
        PhaseAggregateStates = rhs.PhaseAggregateStates;

        FractionsTotal = rhs.FractionsTotal;

        LStencil = rhs.LStencil;
        Resolution = rhs.Resolution;

        GrainPairLimits = rhs.GrainPairLimits;

        VTKDir = rhs.VTKDir;
        RawDataDir = rhs.RawDataDir;

        if (FieldsStatistics.size() != rhs.FieldsStatistics.size())
        {
            FieldsStatistics.Allocate(rhs.FieldsStatistics.size());
        }

        if (Fields.IsNotAllocated())
        {
            FieldsStatistics.Allocate(rhs.FieldsStatistics.size());

            Fields.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, rhs.Fields.Bcells());
            FieldsDot.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, rhs.FieldsDot.Bcells());
            Fractions.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Nphases}, rhs.FieldsDot.Bcells());
        }
        else if (not Fields.IsSize(Nx, Ny, Nz))
        {
            Fields.Reallocate(Nx, Ny, Nz);
            FieldsDot.Reallocate(Nx, Ny, Nz);
            Fractions.Reallocate(Nx, Ny, Nz);
        }

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,Fields.Bcells(),)
        {
            Fields(i,j,k) = rhs.Fields(i,j,k);
            FieldsDot(i,j,k) = rhs.FieldsDot(i,j,k);
            Fractions(i,j,k) = rhs.Fractions(i,j,k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        FieldsStatistics = rhs.FieldsStatistics;

        if(Resolution == Resolutions::Double)
        {
            if (FieldsDR.IsNotAllocated())
            {
                FieldsDR.Allocate(2*Nx, 2*Ny, 2*Nz, dNx, dNy, dNz, rhs.Fields.Bcells());
                FieldsDotDR.Allocate(2*Nx, 2*Ny, 2*Nz, dNx, dNy, dNz, rhs.FieldsDot.Bcells());
            }
            else if (not FieldsDR.IsSize(2*Nx, 2*Ny, 2*Nz))
            {
                Fields.Reallocate(2*Nx, 2*Ny, 2*Nz);
                FieldsDot.Reallocate(2*Nx, 2*Ny, 2*Nz);
            }

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,FieldsDR.Bcells(),)
            {
                FieldsDR(i,j,k) = rhs.FieldsDR(i,j,k);
                FieldsDotDR(i,j,k) = rhs.FieldsDotDR(i,j,k);
            }
            OMP_PARALLEL_STORAGE_LOOP_END
        }
    }
    return *this;
}

std::vector<int> PhaseField::GetPresentPhaseFields() const
{
    std::vector<int> indeces;
    for(unsigned int n = 0; n < FieldsStatistics.size(); n++)
    if(FieldsStatistics[n].Exist)
    {
        indeces.push_back(n);
    }

    return indeces;
}

std::vector<int> PhaseField::ReturnVicinityPhaseFields(const int i,
                                           const int j, const int k) const
{
    // Analyze vicinity
    std::vector<int> tempPFindex;
    tempPFindex.push_back (Fields(i,j,k).front().index);
    for(int ii = -dNx; ii <= dNx; ii++)
    for(int jj = -dNy; jj <= dNy; jj++)
    for(int kk = -dNz; kk <= dNz; kk++)
    {
        for(auto alpha = Fields(i+ii,j+jj,k+kk).cbegin();
                 alpha != Fields(i+ii,j+jj,k+kk).cend(); ++alpha)
        {
            if(alpha->index != Fields(i,j,k).front().index)
            {
                tempPFindex.push_back (alpha->index);
            }
        }
    }

    // Erase double entries
    std::sort(tempPFindex.begin(), tempPFindex.end() );
    tempPFindex.erase( std::unique( tempPFindex.begin(),
                tempPFindex.end() ), tempPFindex.end() );

    return tempPFindex;
}

void PhaseField::CombinePhaseFields(const size_t PhaseIndex)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,Fields.Bcells(),)
    {
        double phasevalue = 0.0;
        bool ChangedPF = false;
        if (Interface(i,j,k))
        {
            for (auto it = Fields(i,j,k).begin();
                      it != Fields(i,j,k).end(); ++it)
            {
                if(FieldsStatistics[it->index].Phase and it->index != PhaseIndex)
                {
                    it->index = PhaseIndex;
                    phasevalue += it->value;
                    it->value = 0.0;
                    ChangedPF = true;
                }
            }
        }
        else
        {
            if(FieldsStatistics[Fields(i,j,k).front().index].Phase and
                    Fields(i,j,k).front().index != PhaseIndex)
            {
                Fields(i,j,k).begin()->index = PhaseIndex;
                phasevalue += Fields(i,j,k).begin()->value;
                Fields(i,j,k).begin()->value = 0.0;
                Fields(i,j,k).front().index = PhaseIndex;
                ChangedPF = true;
            }
        }

        if (ChangedPF == true) Fields(i,j,k).set_value(PhaseIndex,
                                Fields(i,j,k).get_value(PhaseIndex) + phasevalue);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::SelectiveCombinePhaseFieldsSR(BoundaryConditions& BC,
                                               const size_t TargetPFIndex,
                                               const size_t SourcePFIndex)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,Fields.Bcells(),)
    {
        double loc_value = 0.0;
        bool ChangedPF = false;

        for (auto it = Fields(i,j,k).begin();
                  it != Fields(i,j,k).end(); ++it)
        {
            if(it->index == SourcePFIndex)
            {
                loc_value += it->value;
                it->value = 0.0;
                ChangedPF = true;
            }
        }
        if(ChangedPF == true)
        {
            Fields(i,j,k).add_value(TargetPFIndex, loc_value);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Finalize(BC);
}

void PhaseField::SelectiveCombinePhaseFieldsDR(BoundaryConditions& BC,
                                               const size_t TargetPFIndex,
                                               const size_t SourcePFIndex)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,FieldsDR.Bcells(),)
    {
        double loc_value = 0.0;
        bool ChangedPF = false;

        for (auto it = FieldsDR(i,j,k).begin();
                  it != FieldsDR(i,j,k).end(); ++it)
        {
            if(it->index == SourcePFIndex)
            {
                loc_value += it->value;
                it->value = 0.0;
                ChangedPF = true;
            }
        }
        if(ChangedPF == true)
        {
            FieldsDR(i,j,k).add_value(TargetPFIndex, loc_value);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Finalize(BC);
}

pair<NodeA, NodeA> PhaseField::PrincipalCurvatures(const int i,
                                                   const int j,
                                                   const int k) const
{
    NodeA kappa1;
    NodeA kappa2;

    kappa1.clear();
    kappa2.clear();

    const double GWeights[3] = {-0.5/dx, 0.0, 0.5/dx};

    // Check if neighbour cell are also in the interface
    if (Interface(i,j,k))
    {
        // Calculate gradients of phase normal fields (Jacobian matrix)
        Tensor<NodeA, 2 > NormalGradient;
        NormalGradient.Allocate({3,3});
        NormalGradient.set_to_zero();

        // Calculate gradients of normals
        for (int ii = -1; ii <= +1; ii += 2)
        {
            const double GWeight = GWeights[ii+1];

            if (dNx)
            {
                NodeV3 locNormalsX = NormalsPhase(i+ii,j,k);
                for (auto it = locNormalsX.cbegin(); it < locNormalsX.cend(); ++it)
                {
                    NormalGradient({0,0}).add_value(it->index, GWeight*(it->vector3[0]));
                    NormalGradient({1,0}).add_value(it->index, GWeight*(it->vector3[1]));
                    NormalGradient({2,0}).add_value(it->index, GWeight*(it->vector3[2]));
                }
            }
            if (dNy)
            {
                NodeV3 locNormalsY = NormalsPhase(i,j+ii,k);
                for (auto it = locNormalsY.cbegin(); it < locNormalsY.cend(); ++it)
                {
                    NormalGradient({0,1}).add_value(it->index, GWeight*(it->vector3[0]));
                    NormalGradient({1,1}).add_value(it->index, GWeight*(it->vector3[1]));
                    NormalGradient({2,1}).add_value(it->index, GWeight*(it->vector3[2]));
                }
            }
            if (dNz)
            {
                NodeV3 locNormalsZ = NormalsPhase(i,j,k+ii);
                for (auto it = locNormalsZ.cbegin(); it < locNormalsZ.cend(); ++it)
                {
                    NormalGradient({0,2}).add_value(it->index, GWeight*(it->vector3[0]));
                    NormalGradient({1,2}).add_value(it->index, GWeight*(it->vector3[1]));
                    NormalGradient({2,2}).add_value(it->index, GWeight*(it->vector3[2]));
                }
            }
        }

        // Calculate basis of tangent space at (i,j,k)
        // the basis of a spherical coordinate system at (i,j,k) will be
        // used here. The phase normal vector is in this case equal to the
        // radial normal vector of the spherical coordinate and the
        // remaining basis vectors will be the basis of the tangent space

        NodeV3 locNormals = NormalsPhase(i,j,k);
        for (auto it = locNormals.cbegin(); it < locNormals.cend(); ++it)
        {
            double phi   = atan2(it->vector3[1],it->vector3[0]);
            double theta = acos(it->vector3[2]);

            // Calculate basis vectors of tangent space
            dVector3 e_n;
            e_n.set_to_zero();
            e_n[0] = it->vector3[0];
            e_n[1] = it->vector3[1];
            e_n[2] = it->vector3[2];

            dVector3 e_theta;
            e_theta.set_to_zero();
            e_theta[0] = cos(theta) * cos(phi);
            e_theta[1] = cos(theta) * sin(phi);
            e_theta[2] =            - sin(theta);

            dVector3 e_phi;
            e_phi.set_to_zero();
            e_phi[0] = - sin(phi);
            e_phi[1] =   cos(phi);

            // Calculate projection and inclusion matrices
            double Projection [2][3];
            double Inclusion  [3][2];
            for (size_t m = 0; m < 3; ++m)
            {
                Projection[1][m] = e_phi  [m];
                Projection[0][m] = e_theta[m];

                Inclusion[m][1] = Projection[1][m];
                Inclusion[m][0] = Projection[0][m];
            }

            // Calculate local Weingarten map W
            complex<double> W[2][2];
            for (size_t l = 0; l < 2; ++l)
            for (size_t m = 0; m < 2; ++m)
            {
                W[l][m] = 0.0;
                for (size_t n = 0; n < 3; ++n)
                for (size_t o = 0; o < 3; ++o)
                {
                    W[l][m] -= Projection[l][n] *
                        NormalGradient({n,o}).get_value(it->index)
                        * Inclusion[o][m];
                }
            }

            // Calculate eigenvalues of local Weingarten map
            const double part1 = real((W[0][0] + W[1][1])/2.0);
            const double W2    = real((W[0][0] - W[1][1]) * (W[0][0] - W[1][1]));
            const double part2 = real(0.5 * sqrt(W2) + 4.0*W[0][1] * W[1][0]);

            const double locKappa1 = real(part1 - part2);
            const double locKappa2 = real(part1 + part2);

            // Add eigenvalues to Node storage
            kappa1.add_value(it->index, locKappa1);
            kappa2.add_value(it->index, locKappa2);
        }
    }

    return make_pair(kappa1, kappa2);
}

std::array<double,2> PhaseField::PrincipalCurvatures(const int i, const int j, const int k, const size_t phase) const
{
    double phaseKappa1 = 0; // principle curvature phase
    double phaseKappa2 = 0; // principle curvature phase

    if (Interface(i,j,k))
    {
          NodeA locKappa1; // principle curvature of each phase field
          NodeA locKappa2; // principle curvature of each phase field
          // Calculate curvature of all phase fields
          tie(locKappa1, locKappa2) = PrincipalCurvatures(i,j,k);

          // Calculate curvature of desired phase
          for (auto it = locKappa1.cbegin();
                    it != locKappa1.cend(); it++)
          {
              size_t pIndex = FieldsStatistics[it->index].Phase;
              if (pIndex == phase) phaseKappa1 += it->value;
          }
          for (auto it = locKappa2.cbegin();
                    it != locKappa2.cend(); it++)
          {
              size_t pIndex = FieldsStatistics[it->index].Phase;
              if (pIndex == phase) phaseKappa2 += it->value;
          }
    }

    std::array<double,2> result({phaseKappa1,phaseKappa2});
    return result;
}

double PhaseField::Curvature(const int i, const int j, const int k, const size_t phase) const
{
    auto kappa = PrincipalCurvatures(i,j,k,phase);
    return 0.5*(kappa[0]+kappa[1]);
}

NodeA PhaseField::OrdinaryCurvatures(const int i, const int j,
        const int k) const
{
    NodeA locCurvature;

    const double GWeights[3] = {-0.5/dx, 0.0, 0.5/dx};

    // Check if neighbor cell are also in the interface
    if (Fields(i,j,k).flag == 2)
    {
        // Calculate gradients of phase fields
        NodeV3 locGradients = Fields(i,j,k).get_gradients();
        NodeA norms;
        for (auto it = locGradients.cbegin(); it != locGradients.cend(); ++it)
        {
            double loc_norm = it->vector3.abs();
            norms.set_value(it->index,loc_norm);
        }

        // Calculate normals divergence
        for (int ii = -1; ii <= +1; ii += 2)
        {
            const double GWeight = GWeights[ii+1];

            if(dNx)
            {
                NodeV3 locNormalsX = NormalsPhase(i+ii,j,k);
                for (auto it = locNormalsX.cbegin(); it != locNormalsX.cend(); ++it)
                {
                    locCurvature.add_value(it->index, GWeight*(it->vector3[0]));
                }
            }
            if(dNy)
            {
                NodeV3 locNormalsY = NormalsPhase(i,j+ii,k);
                for (auto it = locNormalsY.cbegin(); it < locNormalsY.cend(); ++it)
                {
                    locCurvature.add_value(it->index, GWeight*(it->vector3[1]));
                }
            }
            if(dNz)
            {
                NodeV3 locNormalsZ = NormalsPhase(i,j,k+ii);
                for (auto it = locNormalsZ.cbegin(); it < locNormalsZ.cend(); ++it)
                {
                    locCurvature.add_value(it->index, GWeight*(it->vector3[2]));
                }
            }
        }
        for (auto it = locCurvature.begin(); it != locCurvature.end(); ++it)
        {
            double loc_norm = norms[it->index];
            if(loc_norm > DBL_EPSILON)
            {
                it->value *= loc_norm;
            }
        }
    }
    return locCurvature;
}

NodeA PhaseField::OrdinaryCurvaturesDR(const int i, const int j, const int k) const
{
    NodeA locCurvature;

    const double GWeights[3] = {-0.5/dx, 0.0, 0.5/dx};

    // Check if neighbour cell are also in the interface
    if (Fields(i,j,k).flag == 2)
    {
        // Calculate gradients of phase fields
        NodeV3 locGradients = FieldsDR(i,j,k).get_gradients();
        NodeA norms;
        for (auto it = locGradients.cbegin(); it < locGradients.cend(); ++it)
        {
            double loc_norm = it->vector3.abs();
            norms.set_value(it->index,loc_norm);
        }

        // Calculate normals divergence
        for (int ii = -1; ii <= +1; ii += 2)
        {
            const double GWeight = GWeights[ii+1];

            if(dNx)
            {
                NodeV3 locNormalsX = NormalsPhaseDR(i+ii,j,k);
                for (auto it = locNormalsX.cbegin(); it < locNormalsX.cend(); ++it)
                {
                    locCurvature.add_value(it->index, GWeight*(it->vector3[0]));
                }
            }
            if(dNy)
            {
                NodeV3 locNormalsY = NormalsPhaseDR(i,j+ii,k);
                for (auto it = locNormalsY.cbegin(); it < locNormalsY.cend(); ++it)
                {
                    locCurvature.add_value(it->index, GWeight*(it->vector3[1]));
                }
            }
            if(dNz)
            {
                NodeV3 locNormalsZ = NormalsPhaseDR(i,j,k+ii);
                for (auto it = locNormalsZ.cbegin(); it < locNormalsZ.cend(); ++it)
                {
                    locCurvature.add_value(it->index, GWeight*(it->vector3[2]));
                }
            }
        }
        for (auto it = locCurvature.begin(); it < locCurvature.end(); ++it)
        {
            double loc_norm = norms[it->index];
            if(loc_norm > DBL_EPSILON)
            {
                it->value *= loc_norm;
            }
        }
    }
    return locCurvature;
}

bool PhaseField::PhaseFieldPresent(const int i, const int j, const int k,
                                                        const size_t Index) const
{
    if(!Interface(i,j,k))
    {
        if (Fields(i,j,k).front().index == Index) return true;
    }
    else
    {
        for (auto alpha = Fields(i,j,k).cbegin();
                  alpha != Fields(i,j,k).cend(); ++alpha)
        {
            if (alpha->index == Index) return true;
        }
    }
    return false;
}

bool PhaseField::ThermodynamicPhasePresent(size_t alpha)
{
    bool result = false;

    for(size_t beta = 0; beta < FieldsStatistics.size(); beta++)
    {
        if(FieldsStatistics[beta].Phase == alpha) result = true;
    }

    return result;
}

bool PhaseField::ThermodynamicPhasePairPresent(size_t alpha, size_t beta)
{
    bool result = false;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(x,y,z,Fields,0,)
    if (Interface(x,y,z))
    {
        bool alphapresent = false;
        bool betapresent  = false;
        for (auto it = Fields(x,y,z).cbegin();
                  it != Fields(x,y,z).cend(); ++it)
        {
            if(it->index == alpha)
            {
                alphapresent = true;
            }
            else if(it->index == beta)
            {
                betapresent = true;
            }
        }
        if((alphapresent) and (betapresent)) result = true;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    return result;
}

std::vector<int> PhaseField::GetMaxPhaseFieldOverlap(const size_t thPhase1,
        const size_t thPhase2)
{
    std::vector<int> maxOverlap;
    Matrix<int> Overlap;
    size_t numberOfGrains = FieldsStatistics.size();
    Overlap.Allocate(numberOfGrains, numberOfGrains);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        if (Interface(i,j,k))
        {
            for (auto it = Fields(i,j,k).cbegin();
                      it != Fields(i,j,k).cend(); ++it)
            {
                size_t locIndex1 = it->index;
                size_t thPhaseIndex1 = FieldsStatistics[locIndex1].Phase;
                if (thPhaseIndex1 == thPhase1)
                {
                    for (auto jt = it+1; jt != Fields(i,j,k).cend(); ++jt)
                    {
                        size_t locIndex2 = jt->index;
                        size_t thPhaseIndex2 = FieldsStatistics[locIndex2].Phase;
                        if (thPhaseIndex2 == thPhase2)
                        {
                            Overlap.add(locIndex1, locIndex2, 1);
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    int maxNx = -1;
    int maxNy = -1;
    for (size_t it = 0; it < numberOfGrains-1; it++)
    for (size_t jt = it+1; jt < numberOfGrains; jt++)
    {
        if (maxNx < Overlap.get(it, jt))
        {
            maxNx = it;
            maxNy = jt;
        }
    }

    if (maxNx == -1 or maxNy == -1)
    {
        maxOverlap.push_back(-1);
        maxOverlap.push_back(-1);
        maxOverlap.push_back(-1);
        return maxOverlap;
    }
    maxOverlap.push_back(maxNx);
    maxOverlap.push_back(maxNy);
    maxOverlap.push_back(Overlap.get(maxNx, maxNy));
    return maxOverlap;
}

Matrix<int> PhaseField::GetPhaseFieldOverlap(const size_t thPhase1,
        const size_t thPhase2)
{
    Matrix<int> Overlap;
    int numberOfGrains = FieldsStatistics.size();
    Overlap.Allocate(numberOfGrains, numberOfGrains);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        if (Interface(i,j,k))
        for (auto it = Fields(i,j,k).cbegin();
                  it != Fields(i,j,k).cend(); ++it)
        {
            size_t phaseIndex = it->index;
            size_t thPhaseIndex = FieldsStatistics[phaseIndex].Phase;
            if (thPhaseIndex == thPhase1)
            {
                for (auto jt = it+1;
                          jt < Fields(i,j,k).cend(); ++jt)
                {
                    size_t phaseIndex2 = jt->index;
                    size_t thPhaseIndex2 = FieldsStatistics[phaseIndex2].Phase;

                    #ifdef _OPENMP
                    #pragma omp critical
                    #endif
                    if (thPhaseIndex2 == thPhase2)
                    {
                        Overlap.add(phaseIndex, phaseIndex2, 1);
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    return Overlap;
}

void PhaseField::Advect(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi, const BoundaryConditions& BC, const double dt, const double tStep)
{
    Adv.AdvectPhaseField(Phi, Vel, BC, dx, dt, tStep, true);
}

void PhaseField::Advect(AdvectionHR& Adv, Velocities& Vel, BoundaryConditions& BC, FlowSolverLBM& LBM,
                                 double dt, double tStep)
{
    LBM.DetectObstacles(*this);
    Adv.AdvectPhaseField(*this, Vel, BC, dx, dt, tStep, true);
    LBM.DetectObstaclesAdvection(*this,Vel,BC);
}

Tensor<double,1> PhaseField::CalculateNewFractions(Tensor<double,1> oldFractions,
                                               NodeAB& FieldsDot, double locdt)
{
    Tensor<double,1> newFractions = oldFractions;
    for(auto it = FieldsDot.begin(); it != FieldsDot.end(); ++it)
    {
        size_t pIndexA = FieldsStatistics[it->indexA].Phase;                    // Index of thermodynamic phase alpha in the phase pair alpha-beta
        size_t pIndexB = FieldsStatistics[it->indexB].Phase;                    // Index of thermodynamic phase beta in the phase pair alpha-beta

        if(pIndexA != pIndexB)
        {
            newFractions({pIndexA}) += it->value1*locdt;
            newFractions({pIndexB}) -= it->value1*locdt;
        }
    }
    return newFractions;
}

Tensor<double,2> PhaseField::CalculatePsi(NodeAB& FieldsDot, double locdt)
{
    Tensor<double,2> Psi({Nphases,Nphases});
    Psi.set_to_zero();
    for(auto it = FieldsDot.begin(); it != FieldsDot.end(); ++it)
    {
        size_t pIndexA = FieldsStatistics[it->indexA].Phase;                    // Index of thermodynamic phase alpha in the phase pair alpha-beta
        size_t pIndexB = FieldsStatistics[it->indexB].Phase;                    // Index of thermodynamic phase beta in the phase pair alpha-beta

        if(pIndexA != pIndexB)
        {
            Psi({pIndexA,pIndexB}) += it->value1*locdt;
            Psi({pIndexB,pIndexA}) -= it->value1*locdt;
        }
    }
    return Psi;
}

Storage3D<double,1> PhaseField::NewFractions(double dt)
{
    Storage3D<double,1> result(Nx,Ny,Nz,dNx,dNy,dNz,{Nphases},0);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(x,y,z,result,0,)
    if(Fields(x,y,z).flag)
    {
        result(x,y,z) = CalculateNewFractions(Fractions(x,y,z),
                                              FieldsDot(x,y,z),dt);
    }
    else
    {
        result(x,y,z) = Fractions(x,y,z);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    return result;
}

double PhaseField::Interfaces(const int i, const int j, const int k) const
{
    double sum = 0.5;
    if(Interface(i,j,k))
    for (auto alpha = Fields(i,j,k).cbegin();
              alpha != Fields(i,j,k).cend() - 1; ++alpha)
    for (auto  beta = alpha + 1;
               beta != Fields(i,j,k).cend(); ++beta)
    {
        sum -= alpha->value*beta->value;
    }
    return 1.0/(sum*2.0);
};

double PhaseField::InterfacesDR(const int i, const int j, const int k) const
{
    double sum = 0.5;
    if(InterfaceDR(i,j,k))
    for (auto alpha = FieldsDR(i,j,k).cbegin();
              alpha != FieldsDR(i,j,k).cend() - 1; ++alpha)
    for (auto  beta = alpha + 1;
               beta != FieldsDR(i,j,k).cend(); ++beta)
    {
        sum -= alpha->value*beta->value;
    }
    return 1.0/(sum*2.0);
};

double PhaseField::LocalIndex(const int i, const int j, const int k) const
{
    int locIndex = 0;
    double locValue = 0.0;
    for (auto alpha  = Fields(i,j,k).cbegin();
              alpha != Fields(i,j,k).cend(); ++alpha)
    {
        if(alpha->value > locValue)
        {
            locValue = alpha->value;
            locIndex = alpha->index;
        }
    }
    return locIndex;
};

double PhaseField::LocalIndexDR(const int i, const int j, const int k) const
{
    int locIndex = 0;
    double locValue = 0.0;
    for (auto alpha  = FieldsDR(i,j,k).cbegin();
              alpha != FieldsDR(i,j,k).cend(); ++alpha)
    {
        if(alpha->value > locValue)
        {
            locValue = alpha->value;
            locIndex = alpha->index;
        }
    }
    return locIndex;
};

double PhaseField::Junctions(const int i, const int j, const int k) const
{
    double sum = 0.0;
    if(Interface(i,j,k))
    {
        for (auto alpha = Fields(i,j,k).cbegin();
                  alpha != Fields(i,j,k).cend(); ++alpha)
        if(alpha->value != 0.0)
        {
            sum += 1.0;
        }
    }
    else
    {
        sum = 1.0;
    }
    return sum;
}

double PhaseField::JunctionsDR(const int i, const int j, const int k) const
{
    double sum = 0.0;

    if(InterfaceDR(i,j,k))
    {
        for (auto alpha  = FieldsDR(i,j,k).cbegin();
                  alpha != FieldsDR(i,j,k).cend(); ++alpha)
        if(alpha->value != 0.0)
        {
            sum += 1.0;
        }
    }
    else
    {
        sum = 1.0;
    }
    return sum;
}

double PhaseField::Variants(const int i, const int j, const int k) const
{
    int locVariant = 0;
    double locValue = 0.0;
    for (auto alpha  = Fields(i,j,k).cbegin();
              alpha != Fields(i,j,k).cend(); ++alpha)
    {
        if(alpha->value > locValue)
        {
            locValue = alpha->value;
            locVariant = FieldsStatistics[alpha->index].Variant;
        }
    }
    return locVariant;
}

double PhaseField::VariantsDR(const int i, const int j, const int k) const
{
    int locVariant = 0;
    double locValue = 0.0;
    for (auto alpha  = FieldsDR(i,j,k).cbegin();
              alpha != FieldsDR(i,j,k).cend(); ++alpha)
    {
        if(alpha->value > locValue)
        {
            locValue = alpha->value;
            locVariant = FieldsStatistics[alpha->index].Variant;
        }
    }
    return locVariant;
}

double PhaseField::ColorScale(const int i, const int j, const int k) const
{
    int index = 0;
    double locValue = 0.0;
    for (auto beta = Fields(i,j,k).cbegin();
              beta != Fields(i,j,k).cend(); ++beta)
    {
        if(beta->value > locValue)
        {
            locValue = beta->value;
            index = beta->index;
        }
    }
    return index;
}

double PhaseField::ColorScaleDR(const int i, const int j, const int k) const
{
    int index = 0;
    double locValue = 0.0;
    for (auto beta = FieldsDR(i,j,k).cbegin();
              beta != FieldsDR(i,j,k).cend(); ++beta)
    {
        if(beta->value > locValue)
        {
            locValue = beta->value;
            index = beta->index;
        }
    }
    return index;
}


}// namespace openphase
