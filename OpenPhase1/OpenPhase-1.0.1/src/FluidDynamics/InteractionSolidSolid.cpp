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
 *   File created :   2010
 *   Main contributors :   Oleg Shchyglo; Dmitry Medvedev;
 *                         Marvin Tegeler; Raphael Schiedung
 *
 */

#include "Base/CommonFunctions.h"
#include "FluidDynamics/InteractionSolidSolid.h"
#include "FluidDynamics/FlowSolverLBM.h"

namespace openphase
{

void InteractionSolidSolid::Calculate(PhaseField& Phase,
        const FlowSolverLBM& LBM,
        const double dt, const int order, const int cutoff,
        const double lbStrength, const double elastic)
{
    //NOTE: lbStrength is assumed to be in lattice units!
    const double strength = lbStrength * LBM.dRho*LBM.dx*LBM.dx/(dt*dt);

    // NOTE this let the solid bodies come close but, as a result the simulation
    // is less stable
    if (LBM.Do_SolidSolid)
    {
        // Adjust strength if no lattice units are used
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
        {
            CalculateLocal(Phase,i,j,k,dt,order,cutoff,strength,elastic);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

void InteractionSolidSolid::Calculate(PhaseField& Phase,
        const double dt, const int order, const int cutoff,
        const double strength, const double elastic)
{
    //NOTE: strength has the physical dimension of an energy density!!
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        CalculateLocal(Phase,i,j,k,dt,order,cutoff,strength,elastic);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void InteractionSolidSolid::CalculateLocal(
         PhaseField& Phase, const int i, const int j, const int k,
         const double dt, const int order, const int cutoff,
         const double strength, const double elastic)
{
    //NOTE: strength has the physical dimension of an energy density!!
    const int Nx = Phase.Nx;
    const int Ny = Phase.Ny;
    const int Nz = Phase.Nz;

    const int dNx = Phase.dNx;
    const int dNy = Phase.dNy;
    const int dNz = Phase.dNz;

    const double dx = Phase.dx;
    const double dV = dx*dx*dx;

    for (auto it = Phase.Fields(i,j,k).cbegin();
              it != Phase.Fields(i,j,k).cend(); ++it)
    if (Phase.FieldsStatistics[it->index].State == AggregateStates::Solid)
    for (int ii = -cutoff*dNx; ii <= cutoff*dNx; ++ii)
    for (int jj = -cutoff*dNy; jj <= cutoff*dNy; ++jj)
    for (int kk = -cutoff*dNz; kk <= cutoff*dNz; ++kk)
    if (ii*ii+jj*jj+kk*kk > 0 and ii*ii+jj*jj+kk*kk < cutoff*cutoff)
    {
        const int xx  = ((i+ii)%Nx+Nx)%Nx;
        const int yy  = ((j+jj)%Ny+Ny)%Ny;
        const int zz  = ((k+kk)%Nz+Nz)%Nz;
        const int xx2 = (int(i+0.5*ii)%Nx+Nx)%Nx;
        const int yy2 = (int(j+0.5*jj)%Ny+Ny)%Ny;
        const int zz2 = (int(k+0.5*kk)%Nz+Nz)%Nz;
        for(auto nit = Phase.Fields(xx,yy,zz).cbegin();
                nit != Phase.Fields(xx,yy,zz).cend(); ++nit)
        if(Phase.FieldsStatistics[nit->index].State == AggregateStates::Solid &&
                it->index != nit->index )
        {
            dVector3 BounceForceDensity;
            double dist = sqrt(ii*ii+jj*jj+kk*kk);
            double Vol1 = 4.0/3.0*Pi*it->value*it->value*it->value;
            double Vol2 = 4.0/3.0*Pi*nit->value*nit->value*nit->value;
            double factor = 0.5*strength*(Vol1*Vol2)*order*std::pow((dist-cutoff)/cutoff,order-1)/(dist*cutoff*dx);
            dVector3 Pos1;
            Pos1[0] = xx2; Pos1[1] = yy2; Pos1[2] = zz2;
            dMatrix3x3 Wi;
            Wi(0,0) = 0.0;
            Wi(1,1) = 0.0;
            Wi(2,2) = 0.0;
            Wi(0,1) = -Phase.FieldsStatistics[it->index].aVel[2];
            Wi(0,2) =  Phase.FieldsStatistics[it->index].aVel[1];
            Wi(1,2) = -Phase.FieldsStatistics[it->index].aVel[0];
            Wi(1,0) =  Phase.FieldsStatistics[it->index].aVel[2];
            Wi(2,0) = -Phase.FieldsStatistics[it->index].aVel[1];
            Wi(2,1) =  Phase.FieldsStatistics[it->index].aVel[0];
            dVector3 Ri;
            dVector3 Veli;
            dVector3 distanceCMi;
            CommonFunctions::CalculateDistancePeriodic(Pos1,Phase.FieldsStatistics[it->index].Rcm,distanceCMi, Nx, Ny, Nz);
            Ri = distanceCMi*dx;
            Veli = Phase.FieldsStatistics[it->index].Vcm + Wi*Ri;
            dVector3 Pos2;
            Pos2[0] = xx2; Pos2[1] = yy2; Pos2[2] = zz2;
            dMatrix3x3 Wn;
            Wn(0,0) = 0.0;
            Wn(1,1) = 0.0;
            Wn(2,2) = 0.0;
            Wn(0,1) = -Phase.FieldsStatistics[nit->index].aVel[2];
            Wn(0,2) =  Phase.FieldsStatistics[nit->index].aVel[1];
            Wn(1,2) = -Phase.FieldsStatistics[nit->index].aVel[0];
            Wn(1,0) =  Phase.FieldsStatistics[nit->index].aVel[2];
            Wn(2,0) = -Phase.FieldsStatistics[nit->index].aVel[1];
            Wn(2,1) =  Phase.FieldsStatistics[nit->index].aVel[0];
            dVector3 Rn;
            dVector3 Veln;
            dVector3 distanceCMn;
            CommonFunctions::CalculateDistancePeriodic(Pos2,Phase.FieldsStatistics[nit->index].Rcm,distanceCMn, Nx, Ny, Nz);
            Rn = distanceCMn*dx;
            Veln = Phase.FieldsStatistics[nit->index].Vcm + Wn*Rn;
            auto normal = Phase.Normals(xx2,yy2,zz2).get_asym(it->index,nit->index);
            double v_delta = normal[0]*(Veli[0]-Veln[0])+normal[1]*(Veli[1]-Veln[1])+normal[2]*(Veli[2]-Veln[2]);
            double v_c = dx/dt; // Characteristic velocity to get physical dimensions right
            factor = elastic*factor+(1.0-elastic)*v_delta/v_c*factor;
            BounceForceDensity[0] = ii*factor;
            BounceForceDensity[1] = jj*factor;
            BounceForceDensity[2] = kk*factor;
            dVector3 locR;
            dVector3 pos;
            pos[0] = i; pos[1] = j; pos[2] = k;
            dVector3 distanceCM;
            CommonFunctions::CalculateDistancePeriodic(pos,Phase.FieldsStatistics[it->index].Rcm,distanceCM, Nx, Ny, Nz);
            locR = distanceCM*dx;
            dVector3 nlocR;
            pos[0] = xx; pos[1] = yy; pos[2] = zz;
            CommonFunctions::CalculateDistancePeriodic(pos,Phase.FieldsStatistics[it->index].Rcm,distanceCM, Nx, Ny, Nz);
            nlocR = distanceCM*dx;
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                Phase.FieldsStatistics[it->index].Force  += BounceForceDensity * dV;
                Phase.FieldsStatistics[it->index].Torque += locR.cross(BounceForceDensity) * dV;
                BounceForceDensity *= -1.0;
                Phase.FieldsStatistics[nit->index].Force  += BounceForceDensity * dV;
                Phase.FieldsStatistics[nit->index].Torque += nlocR.cross(BounceForceDensity) * dV;
            }
        }
    }
}

void InteractionSolidSolid::AdvectSolid(PhaseField& Phase,
        const BoundaryConditions& BC, const double dt)
{
    int Nx = Phase.Nx;
    int Ny = Phase.Ny;
    int Nz = Phase.Nz;
    double dx = Phase.dx;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Interface(i,j,k))
        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha < Phase.Fields(i,j,k).cend(); ++alpha)
        {
            if (Phase.FieldsStatistics[alpha->index].State == AggregateStates::Solid )
            {
                double dG_AB = 0.;
                dVector3 Velocity;
                dMatrix3x3 W;
                W(0,0) = 0.0;
                W(1,1) = 0.0;
                W(2,2) = 0.0;
                W(0,1) = -Phase.FieldsStatistics[alpha->index].aVel[2];
                W(0,2) =  Phase.FieldsStatistics[alpha->index].aVel[1];
                W(1,2) = -Phase.FieldsStatistics[alpha->index].aVel[0];
                W(1,0) =  Phase.FieldsStatistics[alpha->index].aVel[2];
                W(2,0) = -Phase.FieldsStatistics[alpha->index].aVel[1];
                W(2,1) =  Phase.FieldsStatistics[alpha->index].aVel[0];
                dVector3 R;
                dVector3 pos;
                pos[0] = i; pos[1] = j; pos[2] = k;
                dVector3 distanceCM;
                CommonFunctions::CalculateDistancePeriodic(pos,Phase.FieldsStatistics[alpha->index].Rcm,distanceCM, Nx, Ny, Nz);
                R = distanceCM*dx;
                Velocity = Phase.FieldsStatistics[alpha->index].Vcm + W*R;
                if (Velocity[0] > 0)
                    dG_AB  -= Velocity[0]*(Phase.Fields(i,j,k).get_value(alpha->index)-Phase.Fields(i-1,j,k).get_value(alpha->index))/Phase.dx;
                else
                    dG_AB  += Velocity[0]*(Phase.Fields(i,j,k).get_value(alpha->index)-Phase.Fields(i+1,j,k).get_value(alpha->index))/Phase.dx;
                if (Velocity[1] > 0)
                    dG_AB  -= Velocity[1]*(Phase.Fields(i,j,k).get_value(alpha->index)-Phase.Fields(i,j-1,k).get_value(alpha->index))/Phase.dx;
                else
                    dG_AB  += Velocity[1]*(Phase.Fields(i,j,k).get_value(alpha->index)-Phase.Fields(i,j+1,k).get_value(alpha->index))/Phase.dx;
                if (Velocity[2] > 0)
                    dG_AB  -= Velocity[2]*(Phase.Fields(i,j,k).get_value(alpha->index)-Phase.Fields(i,j,k-1).get_value(alpha->index))/Phase.dx;
                else
                    dG_AB  += Velocity[2]*(Phase.Fields(i,j,k).get_value(alpha->index)-Phase.Fields(i,j,k+1).get_value(alpha->index))/Phase.dx;
                Phase.FieldsDot(i,j,k).add_asym1(alpha->index,  0, dG_AB);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void InteractionSolidSolid::PreserveVolume(PhaseField& Phase,
            const std::vector<double>& RefVolume, const size_t _index,
            const double dt)
{
    std::vector<double> Area(Phase.FieldsStatistics.size(),0.0);
    std::vector<double> VolumeUpdate(Phase.FieldsStatistics.size(),0.0);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Interface(i,j,k))
        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        for(auto  beta = Phase.Fields(i,j,k).cbegin();
                  beta != Phase.Fields(i,j,k).cend(); ++beta)
        {
            if (alpha->index != beta->index && Phase.FieldsStatistics[alpha->index].Phase == _index)
            {
                {
                    #ifdef _OPENMP
                    #pragma omp atomic
                    #endif
                    Area[alpha->index] += 1.0;
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    for(size_t i = 0; i < Phase.FieldsStatistics.size(); ++i)
    {
        double VolumeDot = (Phase.FieldsStatistics[i].Volume - RefVolume[i])/dt;
        if (Area[i] > 0.0)
        {
            VolumeDot /= Area[i] + 1.0e-16;
        }
        else VolumeDot = 0.0;
        VolumeUpdate[i] = VolumeDot;
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Interface(i,j,k))
        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1;
                  beta != Phase.Fields(i,j,k).cend(); ++beta)
        {
            size_t index1 = Phase.FieldsStatistics[alpha->index].Phase;
            size_t index2 = Phase.FieldsStatistics[beta->index].Phase;
            double dG_AB = 0.0;
            if (index1 == _index)
                dG_AB += VolumeUpdate[alpha->index];
            if (index2 == _index)
                dG_AB -= VolumeUpdate[beta->index];
            dG_AB *= sqrt(std::max(0.,alpha->value*(1.0-alpha->value)))*8.0/Pi;
            Phase.FieldsDot(i,j,k).add_asym1(alpha->index,  beta->index, -dG_AB);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void InteractionSolidSolid::SetRefVolume(const PhaseField& Phase,
        std::vector<double>& RefVolume)
{
    RefVolume.resize(Phase.FieldsStatistics.size(),0.0);
    for(size_t i = 0; i < Phase.FieldsStatistics.size(); ++i)
    {
        RefVolume[i] = Phase.FieldsStatistics[i].Volume;
    }
}

double InteractionSolidSolid::VolumeError(const PhaseField& Phase,
        const std::vector<double>& RefVolume, const size_t index)
{
    double error = 0.0;

    for(size_t i = 0; i < Phase.FieldsStatistics.size(); ++i)
    {
        if ( Phase.FieldsStatistics[i].Phase == index)
        {
            error += std::fabs(RefVolume[i] - Phase.FieldsStatistics[i].Volume);
        }
    }
    return error;
}

} //namespace openphase
