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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Raphael Schiedung
 *
 */

#include "BoundaryConditions.h"
#include "InterfaceProperties.h"
#include "Mechanics/ElasticProperties.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Temperature.h"
#include "Info.h"
#include "Orientations.h"
#include "Nucleation.h"

namespace openphase
{
using namespace std;

void InterfaceProperties::Initialize(Settings& locSettings)
{
    thisclassname = "InterfaceProperties";

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;

    dNx = locSettings.dNx;
    dNy = locSettings.dNy;
    dNz = locSettings.dNz;

    dx = locSettings.dx;

    Nphases = locSettings.Nphases;

    maxSigma = 0.0;
    maxMu = 0.0;
    maxSigmaPhase = 0.0;

    R = PhysicalConstants::R;

    size_t Bcells = locSettings.Bcells;
    IntProperties.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);

    Resolution = locSettings.Resolution;
    if(Resolution == Resolutions::Double)
    {
        IntPropertiesDR.Allocate((1+dNx)*Nx, (1+dNy)*Ny, (1+dNz)*Nz, dNx, dNy, dNz, Bcells*2);
    }
    InterfaceEnergy.Allocate(Nphases, Nphases);
    InterfaceMobility.Allocate(Nphases, Nphases);
    RespectParentBoundaries.Allocate(Nphases, Nphases);

    VTKDir = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void InterfaceProperties::ReadInput(const string InputFileName)
{
    Info::WriteLineInsert("InterfaceProperties input");
    Info::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(1);
    }

    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);

    inp.close();
}

void InterfaceProperties::ReadInput(stringstream& inp)
{
    //  Read values of interface energy and mobility for pairs of phases

   int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    for (size_t alpha = 0; alpha < Nphases; alpha++)
    for (size_t beta  = alpha; beta < Nphases; beta++)
    {
        stringstream converter;
        converter << "_" << alpha << "_" << beta;
        string counter = converter.str();

        string locEnergyAnisotropyModel = UserInterface::ReadParameterK(inp, moduleLocation, "EnergyModel" + counter, false, "Iso");

        if (locEnergyAnisotropyModel == "EXT")
        {
            InterfaceEnergy(alpha,beta).Model = InterfaceEnergyModels::Ext;
        }
        if (locEnergyAnisotropyModel == "ISO")
        {
            InterfaceEnergy(alpha,beta).ReadInputIso(inp, moduleLocation, counter);
        }
        if (locEnergyAnisotropyModel == "CUBIC")
        {
            InterfaceEnergy(alpha,beta).ReadInputCubic(inp, moduleLocation, counter);
        }
        if (locEnergyAnisotropyModel == "HEXBOETTGER")
        {
            InterfaceEnergy(alpha,beta).ReadInputHexBoettger(inp, moduleLocation, counter);
        }
        if (locEnergyAnisotropyModel == "HEXYANG")
        {
            InterfaceEnergy(alpha,beta).ReadInputHexYang(inp, moduleLocation, counter);
        }
        if (locEnergyAnisotropyModel == "HEXSUN")
        {
            InterfaceEnergy(alpha,beta).ReadInputHexSun(inp, moduleLocation, counter);
        }
        if (locEnergyAnisotropyModel == "FACETED")
        {
            InterfaceEnergy(alpha,beta).ReadInputFaceted(inp, moduleLocation, counter);
        }
        InterfaceEnergy(beta,alpha) = InterfaceEnergy(alpha,beta);
        maxSigmaPhase = std::max(InterfaceEnergy(alpha,beta).MaxEnergy, maxSigmaPhase);

        string locMobilityAnisotropyModel = UserInterface::ReadParameterK(inp, moduleLocation, "MobilityModel" + counter, false, "Iso");

        if (locMobilityAnisotropyModel == "EXT")
        {
            InterfaceMobility(alpha,beta).Model = InterfaceMobilityModels::Ext;
        }
        if (locMobilityAnisotropyModel == "ISO")
        {
            InterfaceMobility(alpha,beta).ReadInputIso(inp, moduleLocation, counter);
        }
        if (locMobilityAnisotropyModel == "CUBIC")
        {
            InterfaceMobility(alpha,beta).ReadInputCubic(inp, moduleLocation, counter);
        }
        if (locMobilityAnisotropyModel == "HEXBOETTGER")
        {
            InterfaceMobility(alpha,beta).ReadInputHexBoettger(inp, moduleLocation, counter);
        }
        if (locMobilityAnisotropyModel == "HEXYANG")
        {
            InterfaceMobility(alpha,beta).ReadInputHexYang(inp, moduleLocation, counter);
        }
        if (locMobilityAnisotropyModel == "HEXSUN")
        {
            InterfaceMobility(alpha,beta).ReadInputHexSun(inp, moduleLocation, counter);
        }
        if (locMobilityAnisotropyModel == "FACETED")
        {
            InterfaceMobility(alpha,beta).ReadInputFaceted(inp, moduleLocation, counter);
        }
        InterfaceMobility(beta,alpha) = InterfaceMobility(alpha,beta);

        RespectParentBoundaries(alpha, beta) = UserInterface::ReadParameterB(inp, moduleLocation, string("RPB") + counter, false, false);
        RespectParentBoundaries(beta, alpha) = RespectParentBoundaries(alpha, beta);

    }
    TripleJunctionFactor = UserInterface::ReadParameterD(inp, moduleLocation, string("TripleJunctionFactor"), false, 0.0);

    Info::WriteLine();
    Info::WriteBlankLine();
}

void InterfaceProperties::ReInitialize(const PhaseField& Phase)
{
    Nx = Phase.Nx;
    Ny = Phase.Ny;
    Nz = Phase.Nz;

    IntProperties.Reallocate(Nx, Ny, Nz);
    if(Resolution == Resolutions::Double)
    {
        IntPropertiesDR.Reallocate((1+dNx)*Nx, (1+dNy)*Ny, (1+dNz)*Nz);
    }
    Info::WriteStandard(thisclassname, "Reinitialized");
}

void InterfaceProperties::Coarsen(const PhaseField& Phase)
{
    double norm = 1.0/pow(2.0,dNx+dNy+dNz);
    long int fx = 1 + dNx;
    long int fy = 1 + dNy;
    long int fz = 1 + dNz;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Fields(i,j,k).flag)
        {
            IntProperties(i,j,k).clear();

            for(int di = -dNx; di <= dNx; di+=2)
            for(int dj = -dNy; dj <= dNy; dj+=2)
            for(int dk = -dNz; dk <= dNz; dk+=2)
            {
                IntProperties(i,j,k).add_sym_pairs(IntPropertiesDR(fx*i+(di+1)/2,fy*j+(dj+1)/2,fz*k+(dk+1)/2));
            }
            IntProperties(i,j,k) *= norm;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void InterfaceProperties::Refine(const PhaseField& Phase)
{
    long int fx = 1 + dNx;
    long int fy = 1 + dNy;
    long int fz = 1 + dNz;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,IntProperties,0,)
    {
        if(Phase.Fields(i,j,k).flag)
        for(int di = -dNx; di <= dNx; di+=2)
        for(int dj = -dNy; dj <= dNy; dj+=2)
        for(int dk = -dNz; dk <= dNz; dk+=2)
        {
            IntPropertiesDR(fx*i+(di+1)/2,fy*j+(dj+1)/2,fz*k+(dk+1)/2) = IntProperties.at(i+di*0.25,j+dj*0.25,k+dk*0.25);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}


void InterfaceProperties::SetFacetOrientation(const PhaseField& Phase)
{
    for (size_t alpha = 0; alpha < Nphases; alpha++)
    for (size_t beta  = alpha; beta < Nphases; beta++)
    {
        int pIndexA = Phase.FieldsStatistics[alpha].Phase;
        int pIndexB = Phase.FieldsStatistics[beta].Phase;

        for(size_t i = 0; i < InterfaceEnergy(pIndexA, pIndexB).locFacets.size(); i++)
        {
            for(size_t j = 0 ; j < InterfaceEnergy(pIndexA, pIndexB).FacetVector[i].size() ; j++)
            {
                for(int unsigned n = 0; n < Phase.FieldsStatistics.size(); n++)
                {
                    InterfaceEnergy(pIndexA, pIndexB).FacetVector[i][j].rotate
                            (Phase.FieldsStatistics[n].Orientation.RotationMatrix);
                }
            }
        }
        for(size_t i = 0; i < InterfaceMobility(pIndexA, pIndexB).locFacets.size(); i++)
        {
            for(size_t j = 0 ; j < InterfaceMobility(pIndexA, pIndexB).FacetVector[i].size() ; j++)
            {
                for(int unsigned n = 0; n < Phase.FieldsStatistics.size(); n++)
                {
                    InterfaceMobility(pIndexA, pIndexB).FacetVector[i][j].rotate
                            (Phase.FieldsStatistics[n].Orientation.RotationMatrix);
                }
            }
        }
    }
}

void InterfaceProperties::SetSR(const PhaseField& Phase)
{
    double locMaxSigma = 0.0;
    double locMaxMu = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields, Bcells(), reduction(max:locMaxSigma) reduction(max:locMaxMu))
    {
        if (Phase.Fields(i,j,k).flag)
        {
            clear(i,j,k);

            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
            for(auto  beta = alpha + 1;
                      beta != Phase.Fields(i,j,k).cend(); ++beta)
            {
                int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

                double locEnergy = 0.0;

                switch(InterfaceEnergy(pIndexA,pIndexB).Model)
                {
                    case InterfaceEnergyModels::Cubic :
                    case InterfaceEnergyModels::HexBoettger :
                    case InterfaceEnergyModels::HexYang :
                    case InterfaceEnergyModels::HexSun :
                    {
                        dVector3 Norm = Phase.Normal(i,j,k, alpha->index, beta->index);

                        if(Phase.FieldsStatistics[alpha->index].State == AggregateStates::Solid and
                           Phase.FieldsStatistics[ beta->index].State != AggregateStates::Solid)
                        {
                            Norm = Phase.FieldsStatistics[alpha->index].Orientation.RotationMatrix*Norm;
                        }
                        if(Phase.FieldsStatistics[alpha->index].State != AggregateStates::Solid and
                           Phase.FieldsStatistics[ beta->index].State == AggregateStates::Solid)
                        {
                            Norm = Phase.FieldsStatistics[ beta->index].Orientation.RotationMatrix*Norm;
                        }

                        locEnergy = InterfaceEnergy(pIndexA, pIndexB).Calculate(Norm);
                        break;
                    }
                    case InterfaceEnergyModels::Faceted :
                    {
                        dVector3 Grad_alpha = Phase.Fields(i,j,k).get_gradient(alpha->index);
                        Grad_alpha.normalize();
                        dVector3 Grad_beta = Phase.Fields(i,j,k).get_gradient(beta->index);
                        Grad_beta.normalize();

                        if(Grad_alpha.length() > DBL_EPSILON and Grad_beta.length() > DBL_EPSILON)
                        {
                            double locEnergyAlpha = InterfaceEnergy(pIndexA, pIndexB).Calculate(Grad_alpha);
                            double locEnergyBeta = InterfaceEnergy(pIndexA, pIndexB).Calculate(Grad_beta);
                            locEnergy = 0.5*(locEnergyAlpha + locEnergyBeta);
                        }
                        else
                        {
                            locEnergy = maxSigmaPhase;
                        }
                        break;
                    }
                    case InterfaceEnergyModels::Iso :
                    {
                        locEnergy = InterfaceEnergy(pIndexA, pIndexB).MaxEnergy;
                        break;
                    }
                    case InterfaceEnergyModels::Ext :
                    {
                        break;
                    }
                }

                double locMobility = 0.0;
                switch(InterfaceMobility(pIndexA,pIndexB).Model)
                {
                    case InterfaceMobilityModels::Cubic :
                    case InterfaceMobilityModels::HexBoettger :
                    case InterfaceMobilityModels::HexYang :
                    case InterfaceMobilityModels::HexSun :
                    {
                        dVector3 Norm = Phase.Normal(i,j,k, alpha->index, beta->index);

                        if(Phase.FieldsStatistics[alpha->index].State == AggregateStates::Solid and
                           Phase.FieldsStatistics[ beta->index].State != AggregateStates::Solid)
                        {
                            Norm = Phase.FieldsStatistics[alpha->index].Orientation.RotationMatrix*Norm;
                        }
                        if(Phase.FieldsStatistics[alpha->index].State != AggregateStates::Solid and
                           Phase.FieldsStatistics[ beta->index].State == AggregateStates::Solid)
                        {
                            Norm = Phase.FieldsStatistics[ beta->index].Orientation.RotationMatrix*Norm;
                        }

                        locMobility = InterfaceMobility(pIndexA, pIndexB).Calculate(Norm);
                        break;
                    }
                    case InterfaceMobilityModels::Faceted :
                    {
                        dVector3 Grad_alpha = Phase.Fields(i,j,k).get_gradient(alpha->index);
                        Grad_alpha.normalize();
                        dVector3 Grad_beta = Phase.Fields(i,j,k).get_gradient(beta->index);
                        Grad_beta.normalize();

                        if(Grad_alpha.length() > DBL_EPSILON and Grad_beta.length() > DBL_EPSILON)
                        {
                            double locMobilityAlpha = InterfaceMobility(pIndexA, pIndexB).Calculate(Grad_alpha);
                            double locMobilityBeta = InterfaceMobility(pIndexA, pIndexB).Calculate(Grad_beta);
                            locMobility = 0.5*(locMobilityAlpha + locMobilityBeta);
                        }
                        else
                        {
                            locMobility = InterfaceMobility(pIndexA, pIndexB).MaxMobility;
                        }
                        break;
                    }
                    case InterfaceMobilityModels::Iso :
                    {
                        locMobility = InterfaceMobility(pIndexA, pIndexB).MaxMobility;
                        break;
                    }
                    case InterfaceMobilityModels::Ext :
                    {
                        break;
                    }
                }

                if(RespectParentBoundaries(pIndexA, pIndexB))
                {
                    if(pIndexA == pIndexB and
                       Phase.FieldsStatistics[alpha->index].Parent != Phase.FieldsStatistics[ beta->index].Parent)
                    {
                        locMobility = 0.0;
                    }

                    if(pIndexA != pIndexB and
                       Phase.FieldsStatistics[alpha->index].Parent !=  beta->index and
                       Phase.FieldsStatistics[ beta->index].Parent != alpha->index)
                    {
                        locMobility = 0.0;
                    }
                }

                if(locEnergy != 0.0 or locMobility != 0.0)
                {
                    set_energy_and_mobility(i,j,k,alpha->index, beta->index, locEnergy, locMobility);
                }

                locMaxSigma = max(locMaxSigma,locEnergy);
                locMaxMu = max(locMaxMu,locMobility);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxSigma = locMaxSigma;
    maxMu = locMaxMu;

#ifdef MPI_PARALLEL
    MPI_Allreduce(&locMaxSigma, &maxSigma, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&locMaxMu, &maxMu, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
}

void InterfaceProperties::SetDR(const PhaseField& Phase)
{
    double locMaxSigma = 0.0;
    double locMaxMu = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.FieldsDR, BcellsDR(), reduction(max:locMaxSigma) reduction(max:locMaxMu))
    {
        if (Phase.FieldsDR(i,j,k).flag)
        {
            clear_DR(i,j,k);

            for(auto alpha = Phase.FieldsDR(i,j,k).cbegin();
                     alpha != Phase.FieldsDR(i,j,k).cend() - 1; ++alpha)
            for(auto  beta = alpha + 1;
                      beta != Phase.FieldsDR(i,j,k).cend(); ++beta)
            {
                int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

                double locEnergy = 0.0;
                switch(InterfaceEnergy(pIndexA,pIndexB).Model)
                {
                    case InterfaceEnergyModels::Cubic :
                    case InterfaceEnergyModels::HexBoettger :
                    case InterfaceEnergyModels::HexYang :
                    case InterfaceEnergyModels::HexSun :
                    {
                        dVector3 Norm = Phase.NormalDR(i,j,k, alpha->index, beta->index);

                        if(Phase.FieldsStatistics[alpha->index].State == AggregateStates::Solid and
                           Phase.FieldsStatistics[ beta->index].State != AggregateStates::Solid)
                        {
                            Norm = Phase.FieldsStatistics[alpha->index].Orientation.RotationMatrix*Norm;
                        }
                        if(Phase.FieldsStatistics[alpha->index].State != AggregateStates::Solid and
                           Phase.FieldsStatistics[ beta->index].State == AggregateStates::Solid)
                        {
                            Norm = Phase.FieldsStatistics[ beta->index].Orientation.RotationMatrix*Norm;
                        }

                        locEnergy = InterfaceEnergy(pIndexA, pIndexB).Calculate(Norm);
                        break;
                    }
                    case InterfaceEnergyModels::Faceted :
                    {
                        dVector3 Grad_alpha = Phase.FieldsDR(i,j,k).get_gradient(alpha->index);
                        Grad_alpha.normalize();
                        dVector3 Grad_beta = Phase.FieldsDR(i,j,k).get_gradient(beta->index);
                        Grad_beta.normalize();

                        if(Grad_alpha.length() > DBL_EPSILON and Grad_beta.length() > DBL_EPSILON)
                        {
                            double locEnergyAlpha = InterfaceEnergy(pIndexA, pIndexB).Calculate(Grad_alpha);
                            double locEnergyBeta = InterfaceEnergy(pIndexA, pIndexB).Calculate(Grad_beta);
                            locEnergy = 0.5*(locEnergyAlpha + locEnergyBeta);
                        }
                        else
                        {
                            locEnergy = maxSigmaPhase;
                        }
                        break;
                    }
                    case InterfaceEnergyModels::Iso :
                    {
                        locEnergy = InterfaceEnergy(pIndexA, pIndexB).MaxEnergy;
                        break;
                    }
                    case InterfaceEnergyModels::Ext :
                    {
                        break;
                    }
                }

                double locMobility = 0.0;
                switch(InterfaceMobility(pIndexA,pIndexB).Model)
                {
                    case InterfaceMobilityModels::Cubic :
                    case InterfaceMobilityModels::HexBoettger :
                    case InterfaceMobilityModels::HexYang :
                    case InterfaceMobilityModels::HexSun :
                    {
                        dVector3 Norm = Phase.NormalDR(i,j,k, alpha->index, beta->index);

                        if(Phase.FieldsStatistics[alpha->index].State == AggregateStates::Solid and
                           Phase.FieldsStatistics[ beta->index].State != AggregateStates::Solid)
                        {
                            Norm = Phase.FieldsStatistics[alpha->index].Orientation.RotationMatrix*Norm;
                        }
                        if(Phase.FieldsStatistics[alpha->index].State != AggregateStates::Solid and
                           Phase.FieldsStatistics[ beta->index].State == AggregateStates::Solid)
                        {
                            Norm = Phase.FieldsStatistics[ beta->index].Orientation.RotationMatrix*Norm;
                        }

                        locMobility = InterfaceMobility(pIndexA, pIndexB).Calculate(Norm);
                        break;
                    }
                    case InterfaceMobilityModels::Faceted :
                    {
                        dVector3 Grad_alpha = Phase.Fields(i,j,k).get_gradient(alpha->index);
                        Grad_alpha.normalize();
                        dVector3 Grad_beta = Phase.Fields(i,j,k).get_gradient(beta->index);
                        Grad_beta.normalize();

                        if(Grad_alpha.length() > DBL_EPSILON and Grad_beta.length() > DBL_EPSILON)
                        {
                            double locMobilityAlpha = InterfaceMobility(pIndexA, pIndexB).Calculate(Grad_alpha);
                            double locMobilityBeta = InterfaceMobility(pIndexA, pIndexB).Calculate(Grad_beta);
                            locMobility = locMobilityAlpha + locMobilityBeta;
                        }
                        else
                        {
                            locMobility = InterfaceMobility(pIndexA, pIndexB).MaxMobility;
                        }
                        break;
                    }
                    case InterfaceMobilityModels::Iso :
                    {
                        locMobility = InterfaceMobility(pIndexA, pIndexB).MaxMobility;
                        break;
                    }
                    case InterfaceMobilityModels::Ext :
                    {
                        break;
                    }
                }

                if(RespectParentBoundaries(pIndexA, pIndexB))
                {
                    if(pIndexA == pIndexB and
                       Phase.FieldsStatistics[alpha->index].Parent != Phase.FieldsStatistics[ beta->index].Parent)
                    {
                        locMobility = 0.0;
                    }

                    if(pIndexA != pIndexB and
                       Phase.FieldsStatistics[alpha->index].Parent !=  beta->index and
                       Phase.FieldsStatistics[ beta->index].Parent != alpha->index)
                    {
                        locMobility = 0.0;
                    }
                }

                if(locEnergy != 0.0 or locMobility != 0.0)
                {
                    set_energy_and_mobility_DR(i,j,k,alpha->index, beta->index, locEnergy, locMobility);
                }
                locMaxSigma = max(locMaxSigma,locEnergy);
                locMaxMu = max(locMaxMu,locMobility);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxSigma = locMaxSigma;
    maxMu = locMaxMu;

#ifdef MPI_PARALLEL
    MPI_Allreduce(&locMaxSigma, &maxSigma, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&locMaxMu, &maxMu, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
}

void InterfaceProperties::SetwMD(const PhaseField& Phase, Orientations& OR)
{
    double locMaxSigma = 0.0;
    double locMaxMu = 0.0;
    double constant = 0.85;        // Shear Modulus*Burgers vector/4*Pi*(1-v)

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields, Bcells(), reduction(max:locMaxSigma) reduction(max:locMaxMu))
    {
        if (Phase.Fields(i,j,k).flag)
        {
            clear(i,j,k);

            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
            for(auto  beta = alpha + 1;
                      beta != Phase.Fields(i,j,k).cend(); ++beta)
            {
                int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

                double locEnergy = 0.0;
                double misorientation = 0;
                Quaternion OrientationA =  Phase.FieldsStatistics[alpha->index].Orientation;
                Quaternion OrientationB =  Phase.FieldsStatistics[beta->index].Orientation;
                misorientation = OR.getDisorientation(OrientationA, OrientationB);

                if (misorientation<=0.349066)
                    locEnergy = misorientation*constant*(0.8-log(misorientation));
                else
                    locEnergy = 0.349066*constant*(0.8-log(0.349066));

                double locMobility = 0.0;
                switch(InterfaceMobility(pIndexA,pIndexB).Model)
                {
                    case InterfaceMobilityModels::Cubic :
                    case InterfaceMobilityModels::HexBoettger :
                    case InterfaceMobilityModels::HexYang :
                    case InterfaceMobilityModels::HexSun :
                    {
                        dVector3 Norm = Phase.Normal(i,j,k, alpha->index, beta->index);

                        if(Phase.FieldsStatistics[alpha->index].State == AggregateStates::Solid and
                           Phase.FieldsStatistics[ beta->index].State != AggregateStates::Solid)
                        {
                            Norm = Phase.FieldsStatistics[alpha->index].Orientation.RotationMatrix*Norm;
                        }
                        if(Phase.FieldsStatistics[alpha->index].State != AggregateStates::Solid and
                           Phase.FieldsStatistics[ beta->index].State == AggregateStates::Solid)
                        {
                            Norm = Phase.FieldsStatistics[ beta->index].Orientation.RotationMatrix*Norm;
                        }

                        locMobility = InterfaceMobility(pIndexA, pIndexB).Calculate(Norm);
                        break;
                    }
                    case InterfaceMobilityModels::Faceted :
                    {
                        dVector3 Grad_alpha = Phase.Fields(i,j,k).get_gradient(alpha->index);
                        Grad_alpha.normalize();
                        dVector3 Grad_beta = Phase.Fields(i,j,k).get_gradient(beta->index);
                        Grad_beta.normalize();

                        if(Grad_alpha.length() > DBL_EPSILON and Grad_beta.length() > DBL_EPSILON)
                        {
                            double locMobilityAlpha = InterfaceMobility(pIndexA, pIndexB).Calculate(Grad_alpha);
                            double locMobilityBeta = InterfaceMobility(pIndexA, pIndexB).Calculate(Grad_beta);
                            locMobility = 0.5*(locMobilityAlpha + locMobilityBeta);
                        }
                        else
                        {
                            locMobility = InterfaceMobility(pIndexA, pIndexB).MaxMobility;
                        }
                        break;
                    }
                    case InterfaceMobilityModels::Iso :
                    {
                        locMobility = InterfaceMobility(pIndexA, pIndexB).MaxMobility;
                        break;
                    }
                    case InterfaceMobilityModels::Ext :
                    {
                        break;
                    }
                }

                if(RespectParentBoundaries(pIndexA, pIndexB))
                {
                    if(pIndexA == pIndexB and
                       Phase.FieldsStatistics[alpha->index].Parent != Phase.FieldsStatistics[ beta->index].Parent)
                    {
                        locMobility = 0.0;
                    }

                    if(pIndexA != pIndexB and
                       Phase.FieldsStatistics[alpha->index].Parent !=  beta->index and
                       Phase.FieldsStatistics[ beta->index].Parent != alpha->index)
                    {
                        locMobility = 0.0;
                    }
                }

                if(locEnergy != 0.0 or locMobility != 0.0)
                {
                    set_energy_and_mobility(i,j,k,alpha->index, beta->index, locEnergy, locMobility);
                }
                locMaxSigma = max(locMaxSigma,locEnergy);
                locMaxMu = max(locMaxMu,locMobility);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxSigma = locMaxSigma;
    maxMu = locMaxMu;

#ifdef MPI_PARALLEL
    MPI_Allreduce(&locMaxSigma, &maxSigma, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&locMaxMu, &maxMu, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
}

void InterfaceProperties::SetMobilityThermalEffectSR(const PhaseField& Phase, const Temperature& Tx)
{
    double locMaxMu = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields, Bcells(), reduction(max:locMaxMu))
    {
        if (Phase.Fields(i,j,k).flag)
        {
            double invT = 1.0/(R*Tx(i,j,k));
            for(auto it = IntProperties(i,j,k).begin();
                     it != IntProperties(i,j,k).end(); ++it)
            {
                int pIndexA = Phase.FieldsStatistics[it->indexA].Phase;
                int pIndexB = Phase.FieldsStatistics[it->indexB].Phase;
                it->value2 *= exp(-InterfaceMobility(pIndexA, pIndexB).ActivationEnergy*invT);
                locMaxMu = max(locMaxMu,it->value2);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxMu = locMaxMu;

#ifdef MPI_PARALLEL
    MPI_Allreduce(&locMaxMu, &maxMu, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
}

void InterfaceProperties::SetMobilityThermalEffectDR(const PhaseField& Phase, const Temperature& Tx)
{
    double locMaxMu = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.FieldsDR, BcellsDR(), reduction(max:locMaxMu))
    {
        if (Phase.FieldsDR(i,j,k).flag)
        {
            double invT = 1.0/(R*Tx.at(0.5*i - dNx*0.25, 0.5*j - dNy*0.25, 0.5*k - dNz*0.25));
            for(auto it = IntPropertiesDR(i,j,k).begin();
                     it != IntPropertiesDR(i,j,k).end(); ++it)
            {
                int pIndexA = Phase.FieldsStatistics[it->indexA].Phase;
                int pIndexB = Phase.FieldsStatistics[it->indexB].Phase;
                it->value2 *= exp(-InterfaceMobility(pIndexA, pIndexB).ActivationEnergy*invT);
                locMaxMu = max(locMaxMu,it->value2);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxMu = locMaxMu;

#ifdef MPI_PARALLEL
    MPI_Allreduce(&locMaxMu, &maxMu, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
}

void InterfaceProperties::ReduceMobilityForNucleation(PhaseField& Phi,
                                                    BoundaryConditions& BC,
                                                    Nucleation& Nuc)
{
    /** This function reduces the mobility of recently planted nuclei by a
    factor set in the Nucleation module. This can help avoid spreading of nuclei
    due to too high mobility during the nucleation stage. TODO: move to mobility
    calculation in InterfaceMobility::Set/Calculate! */

    for(size_t grain = 0; grain < Phi.FieldsStatistics.size(); grain++)
    if(Phi.FieldsStatistics[grain].IsNucleus())
    {
        double x = min(1.0,(double)fabs(Phi.FieldsStatistics[grain].Volume
                                       /Phi.RefVolume))*Pi;

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phi.Fields,Bcells(),)
        {
            if (Phi.Fields(i,j,k).flag)
            for(auto alpha = Phi.Fields(i,j,k).cbegin();
                     alpha != Phi.Fields(i,j,k).cend() - 1; ++alpha)
            for(auto  beta = alpha + 1;
                      beta != Phi.Fields(i,j,k).cend(); ++beta)
            {
                size_t indexA = alpha->index;
                size_t indexB = beta->index;

                if(indexA == grain or indexB == grain)
                {
                    int pIndexA = Phi.FieldsStatistics[indexA].Phase;
                    int pIndexB = Phi.FieldsStatistics[indexB].Phase;

                    double MobRed = std::min(Nuc.MobilityReduction(pIndexA,pIndexB),Nuc.MobilityReduction(pIndexB,pIndexA));
                    double corr = MobRed  + (1.0+(sin(x-Pi/2.)))*(1.0-MobRed);

                    double locIntMobility
                    = get_mobility(i,j,k,indexA,indexB)*corr;

                    set_mobility(i,j,k,indexA,indexB,locIntMobility);
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

void InterfaceProperties::WriteVTK(const PhaseField& Phase, const int tStep)
{
    int Lx = Nx;
    int Ly = Ny;
    int Lz = Nz;
    double a = 1.0;
    double b = 0.0;

    if(Resolution == Resolutions::Double)
    {
        Lx = (1+dNx)*Nx;
        Ly = (1+dNy)*Ny;
        Lz = (1+dNz)*Nz;
        a = 0.5;
        b = -0.25;
    }
    stringstream outbuffer;

    outbuffer << "# vtk DataFile Version 3.0\n";
    outbuffer << "InterfaceProperties\n";
    outbuffer << "ASCII\n";
    outbuffer << "DATASET STRUCTURED_GRID\n";
    outbuffer << "DIMENSIONS " << Lx << " " << Ly << " " << Lz << "\n";
    outbuffer << "POINTS " <<  Lx*Ly*Lz << " double\n";

    for(int k = 0; k < Lz; ++k)
    for(int j = 0; j < Ly; ++j)
    for(int i = 0; i < Lx; ++i)
    {
        outbuffer << i*a + b << " " << j*a+b << " " << k*a+b << "\n";
    }
    outbuffer << "\n";
    outbuffer << "POINT_DATA " << Lx*Ly*Lz << "\n";

    outbuffer << "SCALARS InterfaceEnergy double 1\n";
    outbuffer << "LOOKUP_TABLE default\n";

    for(int k = 0; k < Lz; ++k)
    for(int j = 0; j < Ly; ++j)
    for(int i = 0; i < Lx; ++i)
    {
        double sum = 0.0;
        double sum_weights = 0.0;
        if(Resolution == Resolutions::Single)
        {
            if(Phase.Interface(i,j,k))
            for (auto alpha = Phase.Fields(i,j,k).cbegin();
                    alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for (auto  beta = alpha + 1;
                        beta != Phase.Fields(i,j,k).cend(); ++beta)
                {
                    double weight = alpha->value + beta->value;
                    sum_weights += weight;
                    sum += weight*IntProperties(i,j,k).get_sym1(alpha->index, beta->index);
                }
        }
        else
        {
            if(Phase.InterfaceDR(i,j,k))
            for (auto alpha = Phase.FieldsDR(i,j,k).cbegin();
                    alpha != Phase.FieldsDR(i,j,k).cend() - 1; ++alpha)
                for (auto  beta = alpha + 1;
                        beta != Phase.FieldsDR(i,j,k).cend(); ++beta)
                {
                    double weight = alpha->value + beta->value;
                    sum_weights += weight;
                    sum += weight*IntPropertiesDR(i,j,k).get_sym1(alpha->index, beta->index);
                }
        }
        if(sum_weights > 0.0)
        {
            outbuffer << sum/(sum_weights) << "\n";
        }
        else
        {
            outbuffer << 0.0 << "\n";
        }
    }

    outbuffer << "SCALARS InterfaceMobility double 1\n";
    outbuffer << "LOOKUP_TABLE default\n";

    for(int k = 0; k < Lz; ++k)
    for(int j = 0; j < Ly; ++j)
    for(int i = 0; i < Lx; ++i)
    {
        double sum = 0.0;
        double sum_weights = 0.0;
        if(Resolution == Resolutions::Single)
        {
            if(Phase.Interface(i,j,k))
            for (auto alpha = Phase.Fields(i,j,k).cbegin();
                    alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for (auto  beta = alpha + 1;
                        beta != Phase.Fields(i,j,k).cend(); ++beta)
                {
                    double weight = alpha->value + beta->value;
                    sum_weights += weight;
                    sum += weight*IntProperties(i,j,k).get_sym2(alpha->index, beta->index);
                }
        }
        else
        {
            if(Phase.InterfaceDR(i,j,k))
            for (auto alpha = Phase.FieldsDR(i,j,k).cbegin();
                    alpha != Phase.FieldsDR(i,j,k).cend() - 1; ++alpha)
                for (auto  beta = alpha + 1;
                        beta != Phase.FieldsDR(i,j,k).cend(); ++beta)
                {
                    double weight = alpha->value + beta->value;
                    sum_weights += weight;
                    sum += weight*IntPropertiesDR(i,j,k).get_sym2(alpha->index, beta->index);
                }
        }
        if(sum_weights > 0.0)
        {
            outbuffer << sum/(sum_weights) << "\n";
        }
        else
        {
            outbuffer << 0.0 << "\n";
        }
    }
    stringstream tmp;
    tmp << "InterfaceProperties_";
    string FileName = UserInterface::MakeFileName(VTKDir, tmp.str(), tStep, ".vtk");

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbuffer.rdbuf();
    vtk_file.close();
}  //  WriteVTK

double InterfaceProperties::get_energy(const ElasticProperties& EP,
        const int i, const int j, const int k,
        const int alpha, const int beta) const
{
    double locVolumeChange = 0.0;
    for (int n = 0; n < 3; n++)
    {
        locVolumeChange += EP.TotalStrains(i,j,k).get_tensor(n,n);
        locVolumeChange -= EP.EigenStrains(i,j,k).get_tensor(n,n);
    }

    double locInterfaceEnergy = get_energy(i,j,k,alpha,beta);
    locInterfaceEnergy *= (1.0+locVolumeChange);
    return locInterfaceEnergy;
}

}// namespace openphase

