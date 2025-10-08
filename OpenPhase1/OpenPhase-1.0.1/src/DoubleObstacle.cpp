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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev;
 *                         Johannes Goerler; Raphael Schiedung; Stephan Hubig
 *
 */

#include "Base/UserInterface.h"
#include "DoubleObstacle.h"
#include "Settings.h"
#include "DrivingForce.h"
#include "GrainInfo.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "VTK.h"

namespace openphase
{
using namespace std;

void DoubleObstacle::Initialize(Settings& locSettings)
{
    thisclassname = "DoubleObstacle";

    VTKDir = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void DoubleObstacle::CalculatePhaseFieldIncrementsSharp(PhaseField& Phase, InterfaceProperties& IP)
{
    switch(Phase.Resolution)
    {
        case Resolutions::Single:
        {
            CalculatePhaseFieldIncrementsSharpSR(Phase, IP);
            break;
        }
        case Resolutions::Double:
        {
            CalculatePhaseFieldIncrementsSharpDR(Phase, IP);
            break;
        }
    }
}

void DoubleObstacle::CalculatePhaseFieldIncrementsSharpSR(PhaseField& Phase,
                                                      InterfaceProperties& IP)
{
    dVector3 n_vector = {1.0, 0.0, 0.0};
    n_vector.normalize();

    double dx = Phase.dx;
    vector<PotentialCorrections> StencilDirections;

    for(auto ds = Phase.LStencil.cbegin(); ds != Phase.LStencil.cend(); ds++)
    {
        double d_x = ds->di;
        double d_y = ds->dj;
        double d_z = ds->dk;

        PotentialCorrections dirNew;

        dirNew.d_x = d_x;
        dirNew.d_y = d_y;
        dirNew.d_z = d_z;

        dirNew.stencil_weight = ds->weight;

        dirNew.scal_prod = (d_x*n_vector[0] + d_y*n_vector[1] + d_z*n_vector[2])*dx;

        dirNew.alpha     = cos(Pi*(dirNew.scal_prod)/Phase.Eta);
        dirNew.beta      = sin(Pi*(dirNew.scal_prod)/Phase.Eta);
        dirNew.phi_tl    = 0.5*cos(Pi*(1.0 - (dirNew.scal_prod)/Phase.Eta)) + 0.5;
        dirNew.phi_tu    = 0.5*cos(Pi*(dirNew.scal_prod)/Phase.Eta) + 0.5;

        StencilDirections.push_back(dirNew);
    }

    const double Prefactor = Pi*Pi/(Phase.Eta*Phase.Eta);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if (Phase.Fields(i,j,k).flag)
        {
            double norm_1 = 1.0/Phase.LocalNumberOfPhaseFieldsSR(i,j,k);

            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
            for(auto  beta = alpha + 1;
                      beta != Phase.Fields(i,j,k).cend(); ++beta)
            {
                double dPhi_dt = 0.0;

                double pot_term_a = 0.0;
                for(auto StDir = StencilDirections.cbegin(); StDir != StencilDirections.cend(); StDir++)
                {
                    double fik_term = 0;
                    if(StDir->scal_prod < 0 && alpha->value > StDir->phi_tu)
                    {
                        fik_term = 1.0;
                    }
                    else if(StDir->scal_prod > 0 && alpha->value < StDir->phi_tl)
                    {
                        fik_term = 0.0;
                    }
                    else
                    {
                        fik_term = 0.5 + StDir->alpha*(alpha->value - 0.5) - StDir->beta*sqrt(alpha->value*(1.0 - alpha->value));
                    }
                    pot_term_a += -StDir->stencil_weight*(fik_term - alpha->value);
                }

                double pot_term_b = 0.0;
                for(auto StDir = StencilDirections.cbegin(); StDir != StencilDirections.cend() ; StDir++)
                {
                    double fik_term = 0.0;
                    if(StDir->scal_prod < 0 && beta->value > StDir->phi_tu)
                    {
                        fik_term = 1.0;
                    }
                    else if(StDir->scal_prod > 0 && beta->value < StDir->phi_tl)
                    {
                        fik_term = 0.0;
                    }
                    else
                    {
                        fik_term = 0.5 + StDir->alpha*(beta->value - 0.5) - StDir->beta*sqrt(beta->value*(1.0 - beta->value));
                    }
                    pot_term_b += -StDir->stencil_weight*(fik_term - beta->value);
                }

                double scale = sqrt(Phase.FieldsStatistics[alpha->index].VolumeRatio*
                                    Phase.FieldsStatistics[ beta->index].VolumeRatio);

                dPhi_dt = IP.get_energy(i,j,k,alpha->index, beta->index)
                         *((alpha->laplacian - beta->laplacian) + 0.5*(1.0 + scale)*(pot_term_a - pot_term_b));

                if(Phase.Fields(i,j,k).size() > 2)
                for(auto gamma = Phase.Fields(i,j,k).cbegin();
                         gamma != Phase.Fields(i,j,k).cend(); ++gamma)
                if((gamma != alpha) && (gamma != beta))
                {
                    double pot_term_c = 0.0;
                    for(auto StDir = StencilDirections.cbegin(); StDir != StencilDirections.cend(); StDir++)
                    {
                        double fik_term = 0.0;
                        if(StDir->scal_prod < 0 && gamma->value > StDir->phi_tu)
                        {
                            fik_term = 1.0;
                        }
                        else if(StDir->scal_prod > 0 && gamma->value < StDir->phi_tl)
                        {
                            fik_term = 0.0;
                        }
                        else
                        {
                            fik_term = 0.5 + StDir->alpha*(gamma->value - 0.5) - StDir->beta*sqrt(gamma->value*(1.0 - gamma->value));
                        }
                        pot_term_c += -StDir->stencil_weight*(fik_term - gamma->value);
                    }

                    dPhi_dt += (IP.get_energy(i, j, k,  beta->index, gamma->index) -
                                IP.get_energy(i, j, k, alpha->index, gamma->index))
                              *(gamma->laplacian + pot_term_c + Prefactor*0.5);

                    if(IP.TripleJunctionFactor != 0.0 or scale < 1.0)
                    {
                        dPhi_dt += max(IP.TripleJunctionFactor, (1.0 - scale))*Prefactor
                             *(IP.get_energy(i,j,k,alpha->index, beta->index) +
                               IP.get_energy(i,j,k, beta->index,gamma->index) +
                               IP.get_energy(i,j,k,alpha->index,gamma->index))
                             *(alpha->value*gamma->value - beta->value*gamma->value);
                    }
                }
                dPhi_dt *= IP.get_mobility(i,j,k, alpha->index, beta->index)*norm_1*scale;
                Phase.FieldsDot(i,j,k).add_asym1(alpha->index, beta->index, dPhi_dt);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void DoubleObstacle::CalculatePhaseFieldIncrementsSharpDR(PhaseField& Phase,
                                                      InterfaceProperties& IP)
{
    dVector3 n_vector = {1.0, 0.0, 0.0};
    n_vector.normalize();

    double dx = Phase.dx;
    vector<PotentialCorrections> StencilDirections;

    for(auto ds = Phase.LStencil.cbegin(); ds != Phase.LStencil.cend(); ds++)
    {
        double d_x = ds->di;
        double d_y = ds->dj;
        double d_z = ds->dk;

        PotentialCorrections dirNew;

        dirNew.d_x = d_x;
        dirNew.d_y = d_y;
        dirNew.d_z = d_z;

        dirNew.stencil_weight = ds->weight;

        dirNew.scal_prod = (d_x*n_vector[0] + d_y*n_vector[1] + d_z*n_vector[2])*dx;

        dirNew.alpha     = cos(Pi*(dirNew.scal_prod)/Phase.Eta);
        dirNew.beta      = sin(Pi*(dirNew.scal_prod)/Phase.Eta);
        dirNew.phi_tl    = 0.5*cos(Pi*(1.0 - (dirNew.scal_prod)/Phase.Eta)) + 0.5;
        dirNew.phi_tu    = 0.5*cos(Pi*(dirNew.scal_prod)/Phase.Eta) + 0.5;

        StencilDirections.push_back(dirNew);
    }

    const double Prefactor = Pi*Pi/(Phase.Eta*Phase.Eta);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.FieldsDR,0,)
    {
        if (Phase.FieldsDR(i,j,k).flag)
        {
            double norm_1 = 1.0/Phase.LocalNumberOfPhaseFieldsDR(i,j,k);

            for(auto alpha = Phase.FieldsDR(i,j,k).cbegin();
                     alpha != Phase.FieldsDR(i,j,k).cend() - 1; ++alpha)
            for(auto  beta = alpha + 1;
                      beta != Phase.FieldsDR(i,j,k).cend(); ++beta)
            {
                double dPhi_dt = 0.0;

                double pot_term_a = 0.0;
                for(auto StDir = StencilDirections.cbegin(); StDir != StencilDirections.cend(); StDir++)
                {
                    double fik_term = 0;
                    if(StDir->scal_prod < 0 && alpha->value > StDir->phi_tu)
                    {
                        fik_term = 1.0;
                    }
                    else if(StDir->scal_prod > 0 && alpha->value < StDir->phi_tl)
                    {
                        fik_term = 0.0;
                    }
                    else
                    {
                        fik_term = 0.5 + StDir->alpha*(alpha->value - 0.5) - StDir->beta*sqrt(alpha->value*(1.0 - alpha->value));
                    }
                    pot_term_a += -StDir->stencil_weight*(fik_term - alpha->value);
                }

                double pot_term_b = 0.0;
                for(auto StDir = StencilDirections.cbegin(); StDir != StencilDirections.cend() ; StDir++)
                {
                    double fik_term = 0.0;
                    if(StDir->scal_prod < 0 && beta->value > StDir->phi_tu)
                    {
                        fik_term = 1.0;
                    }
                    else if(StDir->scal_prod > 0 && beta->value < StDir->phi_tl)
                    {
                        fik_term = 0.0;
                    }
                    else
                    {
                        fik_term = 0.5 + StDir->alpha*(beta->value - 0.5) - StDir->beta*sqrt(beta->value*(1.0 - beta->value));
                    }
                    pot_term_b += -StDir->stencil_weight*(fik_term - beta->value);
                }

                double scale = sqrt(Phase.FieldsStatistics[alpha->index].VolumeRatio*
                                    Phase.FieldsStatistics[ beta->index].VolumeRatio);

                dPhi_dt = 0.5*(1.0 + scale)*IP.get_energy_DR(i,j,k,alpha->index, beta->index)
                         *((alpha->laplacian - beta->laplacian) + (pot_term_a - pot_term_b));

                if(Phase.FieldsDR(i,j,k).size() > 2)
                for(auto gamma = Phase.FieldsDR(i,j,k).cbegin();
                         gamma != Phase.FieldsDR(i,j,k).cend(); ++gamma)
                if((gamma != alpha) && (gamma != beta))
                {
                    double pot_term_c = 0.0;
                    for(auto StDir = StencilDirections.cbegin(); StDir != StencilDirections.cend(); StDir++)
                    {
                        double fik_term = 0.0;
                        if(StDir->scal_prod < 0 && gamma->value > StDir->phi_tu)
                        {
                            fik_term = 1.0;
                        }
                        else if(StDir->scal_prod > 0 && gamma->value < StDir->phi_tl)
                        {
                            fik_term = 0.0;
                        }
                        else
                        {
                            fik_term = 0.5 + StDir->alpha*(gamma->value - 0.5) - StDir->beta*sqrt(gamma->value*(1.0 - gamma->value));
                        }
                        pot_term_c += -StDir->stencil_weight*(fik_term - gamma->value);
                    }

                    dPhi_dt += (IP.get_energy_DR(i, j, k,  beta->index, gamma->index) -
                                IP.get_energy_DR(i, j, k, alpha->index, gamma->index))
                              *(gamma->laplacian + pot_term_c + Prefactor*0.5);

                    if(IP.TripleJunctionFactor != 0.0 or scale < 1.0)
                    {
                        dPhi_dt += max(IP.TripleJunctionFactor, (1.0 - scale))*Prefactor
                             *(IP.get_energy_DR(i,j,k,alpha->index, beta->index) +
                               IP.get_energy_DR(i,j,k, beta->index,gamma->index) +
                               IP.get_energy_DR(i,j,k,alpha->index,gamma->index))
                             *(alpha->value*gamma->value - beta->value*gamma->value);
                    }
                }
                dPhi_dt *= IP.get_mobility_DR(i,j,k, alpha->index, beta->index)*norm_1*scale;
                Phase.FieldsDotDR(i,j,k).add_asym1(alpha->index, beta->index, dPhi_dt);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DoubleObstacle::CalculatePhaseFieldIncrements(PhaseField& Phase, InterfaceProperties& IP)
{
    switch(Phase.Resolution)
    {
        case Resolutions::Single:
        {
            CalculatePhaseFieldIncrementsSR(Phase, IP);
            break;
        }
        case Resolutions::Double:
        {
            CalculatePhaseFieldIncrementsDR(Phase, IP);
            break;
        }
    }
}

void DoubleObstacle::CalculatePhaseFieldIncrementsSR(PhaseField& Phase,
                                                     InterfaceProperties& IP)
{
    const double Prefactor = Pi*Pi/(Phase.Eta*Phase.Eta);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if (Phase.Fields(i,j,k).flag)
        {
            double norm_1 = 1.0/Phase.LocalNumberOfPhaseFieldsSR(i,j,k);

            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
            for(auto  beta = alpha + 1;
                      beta != Phase.Fields(i,j,k).cend(); ++beta)
            {
                double scale = sqrt(Phase.FieldsStatistics[alpha->index].VolumeRatio *
                                    Phase.FieldsStatistics[ beta->index].VolumeRatio);

                double dPhi_dt = IP.get_energy(i,j,k,alpha->index,beta->index)
                                 *((alpha->laplacian + 0.5*(1.0 + scale)*Prefactor*alpha->value) -
                                   ( beta->laplacian + 0.5*(1.0 + scale)*Prefactor* beta->value));

                if(Phase.Fields(i,j,k).size() > 2)
                for(auto gamma = Phase.Fields(i,j,k).cbegin();
                         gamma != Phase.Fields(i,j,k).cend(); ++gamma)
                if((gamma != alpha) && (gamma != beta))
                {
                    dPhi_dt += (IP.get_energy(i, j, k,  beta->index, gamma->index) -
                                IP.get_energy(i, j, k, alpha->index, gamma->index))
                              *(gamma->laplacian + Prefactor*gamma->value);

                    if((IP.TripleJunctionFactor != 0.0) or scale < 1.0)
                    {
                        dPhi_dt += max(IP.TripleJunctionFactor, (1.0 - scale))*Prefactor
                             *(IP.get_energy(i,j,k,alpha->index, beta->index) +
                               IP.get_energy(i,j,k, beta->index,gamma->index) +
                               IP.get_energy(i,j,k,alpha->index,gamma->index))
                             *(alpha->value*gamma->value - beta->value*gamma->value);
                    }
                }
                dPhi_dt *= IP.get_mobility(i,j,k, alpha->index, beta->index)*norm_1*scale;
                Phase.FieldsDot(i,j,k).add_asym1(alpha->index, beta->index, dPhi_dt);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DoubleObstacle::CalculatePhaseFieldIncrementsDR(PhaseField& Phase,
                                                     InterfaceProperties& IP)
{
    const double Prefactor = Pi*Pi/(Phase.Eta*Phase.Eta);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.FieldsDR,0,)
    {
        if(Phase.FieldsDR(i,j,k).flag)
        {
            double norm_1 = 1.0/Phase.LocalNumberOfPhaseFieldsDR(i,j,k);

            for(auto alpha = Phase.FieldsDR(i,j,k).cbegin();
                     alpha != Phase.FieldsDR(i,j,k).cend() - 1; ++alpha)
            for(auto  beta = alpha + 1;
                      beta != Phase.FieldsDR(i,j,k).cend(); ++beta)
            {
                double scale = sqrt(Phase.FieldsStatistics[alpha->index].VolumeRatio *
                                    Phase.FieldsStatistics[ beta->index].VolumeRatio);

                double dPhi_dt = 0.5*(1.0 + scale)*IP.get_energy_DR(i,j,k,alpha->index,beta->index)
                                 *((alpha->laplacian + Prefactor*alpha->value) -
                                   ( beta->laplacian + Prefactor* beta->value));

                if(Phase.FieldsDR(i,j,k).size() > 2)
                for(auto gamma = Phase.FieldsDR(i,j,k).cbegin();
                         gamma != Phase.FieldsDR(i,j,k).cend(); ++gamma)
                if((gamma != alpha) && (gamma != beta))
                {
                    dPhi_dt += (IP.get_energy_DR(i, j, k,  beta->index, gamma->index) -
                                IP.get_energy_DR(i, j, k, alpha->index, gamma->index))
                              *(gamma->laplacian + Prefactor*gamma->value);

                    if(IP.TripleJunctionFactor != 0.0 or scale < 1.0)
                    {
                        dPhi_dt += max(IP.TripleJunctionFactor, (1.0 - scale))*Prefactor
                             *(IP.get_energy_DR(i,j,k,alpha->index, beta->index) +
                               IP.get_energy_DR(i,j,k, beta->index,gamma->index) +
                               IP.get_energy_DR(i,j,k,alpha->index,gamma->index))
                             *(alpha->value*gamma->value - beta->value*gamma->value);
                    }
                }
                dPhi_dt *= IP.get_mobility_DR(i,j,k, alpha->index, beta->index)*norm_1*scale;
                Phase.FieldsDotDR(i,j,k).add_asym1(alpha->index, beta->index, dPhi_dt);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

double DoubleObstacle::PointEnergy(const PhaseField& Phase,
                                   const InterfaceProperties& IP,
                                   int i, int j, int k) const
{
    double energy = 0.0;
    double Prefactor = Phase.Eta*Phase.Eta/(Pi*Pi);

    if (Phase.Fields(i,j,k).flag)
    {
        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1;
                  beta != Phase.Fields(i,j,k).cend(); ++beta)
        {
            energy += 4.0*IP.get_energy(i,j,k,alpha->index, beta->index)/Phase.Eta*
            (alpha->value*beta->value - Prefactor*(alpha->gradient*beta->gradient));
        }
    }
    return energy;
}

double DoubleObstacle::PointEnergy(const PhaseField& Phase,
                                   const InterfaceProperties& IP,
                                   const ElasticProperties& EP,
                                   const int i, const int j, const int k) const
{
    double energy = 0.0;
    double Prefactor = Phase.Eta*Phase.Eta/(Pi*Pi);

    if (Phase.Fields(i,j,k).flag)
    {
        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1;
                  beta != Phase.Fields(i,j,k).cend(); ++beta)
        {
            energy += 4.0*IP.get_energy(EP,i,j,k,alpha->index, beta->index)/Phase.Eta*
            (alpha->value*beta->value - Prefactor*(alpha->gradient*beta->gradient));
        }
    }
    return energy;
}

double DoubleObstacle::Energy(const PhaseField& Phase,
                              const InterfaceProperties& IP) const
{
    double Energy = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(+:Energy))
    {
        Energy += PointEnergy(Phase, IP, i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    double loc_Energy = Energy;
    MPI_Allreduce(&loc_Energy, &Energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    const double dx = Phase.dx;
    return Energy*dx*dx*dx;
}

double DoubleObstacle::Energy(const PhaseField& Phase,
                              const InterfaceProperties& IP,
                              const ElasticProperties& EP) const
{
    double Energy = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(+:Energy))
    {
        Energy += PointEnergy(Phase, IP, EP, i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    double loc_Energy = Energy;
    MPI_Allreduce(&loc_Energy, &Energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    const double dx = Phase.dx;
    return Energy*dx*dx*dx;
}

double DoubleObstacle::AverageEnergyDensity(const PhaseField& Phase,
                                            const InterfaceProperties& IP) const
{
    double Energy = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(+:Energy))
    {
        Energy += PointEnergy(Phase, IP, i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    double loc_Energy = Energy;
    MPI_Allreduce(&loc_Energy, &Energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    const double Nx = Phase.TotalNx;
    const double Ny = Phase.Ny;
    const double Nz = Phase.Nz;
    return Energy/(Nx*Ny*Nz);
}

void DoubleObstacle::WriteEnergyVTK(const int tStep,
                                    const Settings& locSettings,
                                    const PhaseField& Phase,
                                    const InterfaceProperties& IP) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t){"Free Energy Density (Interface)", [this, &Phase, &IP](int i,int j,int k){return double(PointEnergy(Phase, IP, i,j,k));}});

    std::string Filename = UserInterface::MakeFileName(VTKDir, thisclassname, tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}

void DoubleObstacle::WriteEnergyVTK(const int tStep,
                                    const Settings& locSettings,
                                    const PhaseField& Phase,
                                    const InterfaceProperties& IP,
                                    const ElasticProperties& EP) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t){"Free Energy Density (Interface)", [this, &Phase, &IP, &EP](int i,int j,int k){return double(PointEnergy(Phase, IP, EP, i,j,k));}});

    std::string Filename = UserInterface::MakeFileName(VTKDir, thisclassname, tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}

}// namespace openphase
