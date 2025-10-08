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
 *   File created :   2015
 *   Main contributors :   Philipp Engels; Marvin Tegeler; Raphael Schiedung
 *
 */

#include "AdvectionHR/AdvectionHR.h"
#include "Base/CommonFunctions.h"
#include "Base/Includes.h"
#include "Base/UserInterface.h"
#include "BoundaryConditions.h"
#include "Composition.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "Mechanics/ElasticProperties.h"
#include "Info.h"
#include "Orientations.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Temperature.h"
#include "Temperature.h"
#include "Velocities.h"

namespace openphase
{

AdvectionHR::AdvectionHR(Settings& locSettings, std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void AdvectionHR::Initialize(Settings& locSettings)
{
    thisclassname = "AdvectionHR";

    Scheme = AdvectionSchemes::Upwind;

    dNx = locSettings.dNx;
    dNy = locSettings.dNy;
    dNz = locSettings.dNz;

    if(locSettings.Bcells < 2)
    {
        std::string message  = "AdvectionHD requires at least 2 boundary cells for correct operation!\n";
                    message += "The actually selected number of boundary cells is ";
                    message += std::to_string(locSettings.Bcells);
                    message += ".\n";
                    message += "Adjust parameter $Bcells in @Settings section of the project input file accordingly!";
        Info::WriteExit(message, thisclassname, "Initialize()");
        exit(1);
    }
    Info::WriteStandard(thisclassname, "Initialized");
}

void AdvectionHR::ReadInput(const std::string InputFileName)
{
    Info::WriteLineInsert("AdvectionHR input");
    Info::WriteStandard("Source", InputFileName);

    std::fstream inpF(InputFileName.c_str(), std::ios::in | std::ios_base::binary);

    if (!inpF)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };
    std::stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();

	ReadInput(inp);	

    Info::WriteLine();
}

void AdvectionHR::ReadInput(std::stringstream& inp)
{
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
    std::string schemeString = UserInterface::ReadParameterK(inp, moduleLocation, "scheme", false, "UPWIND");

    if (schemeString == "MINMOD")
    {
        Scheme = AdvectionSchemes::Minmod;
    }
    else if (schemeString == "MC" or schemeString == "MONOTONISEDCENTRAL")
    {
        Scheme = AdvectionSchemes::MonotonizedCentral;
    }
    else if (schemeString == "SUPERBEE")
    {
        Scheme = AdvectionSchemes::Superbee;
    }
    /*else if (schemeString == "LAXWENDROFF")
    {
        Scheme = AdvectionSchemes::LaxWendroff;
    }*/
    else if (schemeString == "UPWIND")
    {
        Scheme = AdvectionSchemes::Upwind;
    }
    else
    {
        Info::WriteWarning("Wrong advection scheme specified!\nThe default \"Upwind\" scheme is used!", thisclassname, "ReadInput()");
    }
    Info::WriteLine();
    Info::WriteBlankLine();
}

void AdvectionHR::Advection(
        PhaseField& Phase,
        const BoundaryConditions& BC,
        const Velocities& Vel,
        const int direction,
        const double dt,
        const double dx,
        double (*Limiter)(const double,const double))
{
    Storage3D<NodePF, 0> PhaseFieldBackup(Phase.Fields);
    Storage3D<NodePF, 0> PhaseFieldUpdate(Phase.Fields);

    int CFL = 0;
    double locdt = dt;
    std::vector<int> dir(3);
    dir[0] = ((direction == 0)and(Phase.dNx))?1:0;
    dir[1] = ((direction == 1)and(Phase.dNy))?1:0;
    dir[2] = ((direction == 2)and(Phase.dNz))?1:0;
    const int offset = 0;
    int its = 1;
    do {
        CFL = 0;
        for (int it = 0; it < its; ++it)
        {
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,offset,reduction(+:CFL))
            if (Phase.Fields(i,j,k).flag)
            {
                PhaseFieldUpdate(i,j,k) = Phase.Fields(i,j,k);

                for (auto alpha  = Phase.Fields(i,j,k).begin();
                          alpha != Phase.Fields(i,j,k).end() - 1; alpha++)
                if(Phase.FieldsStatistics[alpha->index].State  == AggregateStates::Solid and
                   Phase.FieldsStatistics[alpha->index].Mobile == true)
                {
                    for (auto beta  = alpha + 1;
                              beta != Phase.Fields(i,j,k).end(); beta++)
                    if(Phase.FieldsStatistics[beta->index].State != AggregateStates::Solid)
                    {
                        const size_t Nx = Phase.Nx;
                        const size_t Ny = Phase.Ny;
                        const size_t Nz = Phase.Nz;

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
                        const dVector3 pos  = {double(i),double(j),double(k)};
                        const dVector3 posp = {double(i+dir[0]),double(j+dir[1]),double(k+dir[2])};
                        const dVector3 posm = {double(i-dir[0]),double(j-dir[1]),double(k-dir[2])};
                        dVector3 distanceCM;
                        dVector3 distanceCMp;
                        dVector3 distanceCMm;
                        CommonFunctions::CalculateDistancePeriodic(pos, Phase.FieldsStatistics[alpha->index].Rcm,distanceCM,  Nx, Ny, Nz);
                        CommonFunctions::CalculateDistancePeriodic(posp,Phase.FieldsStatistics[alpha->index].Rcm,distanceCMp, Nx, Ny, Nz);
                        CommonFunctions::CalculateDistancePeriodic(posm,Phase.FieldsStatistics[alpha->index].Rcm,distanceCMm, Nx, Ny, Nz);
                        const dVector3 R  = distanceCM *dx;
                        const dVector3 Rp = distanceCMp*dx;
                        const dVector3 Rm = distanceCMm*dx;
                        const dVector3 Vcm = Phase.FieldsStatistics[alpha->index].Vcm;
                        const dVector3 Velocity  = Vcm + W*R;
                        const dVector3 Velocityp = Vcm + W*Rp;
                        const dVector3 Velocitym = Vcm + W*Rm;

                        const double q   = alpha->value;
                        const double qp  = Phase.Fields(i+dir[0],j+dir[1],k+dir[2])[alpha->index];
                        const double qm  = Phase.Fields(i-dir[0],j-dir[1],k-dir[2])[alpha->index];
                        const double qpp = Phase.Fields(i+2*dir[0],j+2*dir[1],k+2*dir[2])[alpha->index];
                        const double qmm = Phase.Fields(i-2*dir[0],j-2*dir[1],k-2*dir[2])[alpha->index];
                        const double v   = Velocity[direction];
                        const double vp  = Velocityp[direction];
                        const double vm  = Velocitym[direction];
                        const double L0  = 0.5*Limiter(q-qm,qp-q);
                        const double Lp1 = 0.5*Limiter(qp-q,qpp-qp);
                        const double Lm1 = 0.5*Limiter(qm-qmm,q-qm);
                        const double ap  = 0.5*(vp+v);
                        const double am  = 0.5*(vm+v);
                        const double Fp  = std::max(ap,0.0)*(q+L0)
                                         + std::min(ap,0.0)*(qp-Lp1);
                        const double Fm  = std::max(am,0.0)*(qm+Lm1)
                                         + std::min(am,0.0)*(q-L0);
                        if (std::fabs(ap)*locdt > 0.49*dx || std::fabs(am)*locdt > 0.49*dx)
                        {
                            CFL = 1;
                        }
                        PhaseFieldUpdate(i,j,k).add_value(alpha->index, -(Fp-Fm)*locdt/dx);
                        PhaseFieldUpdate(i,j,k).add_value( beta->index,  (Fp-Fm)*locdt/dx);
                    }
                }
                else if (Phase.FieldsStatistics[alpha->index].State != AggregateStates::Solid)
                {
                    for (auto beta  = alpha + 1;
                              beta != Phase.Fields(i,j,k).end(); beta++)
                    if(Phase.FieldsStatistics[beta->index].State != AggregateStates::Solid)
                    {
                        const size_t pIndex = Phase.FieldsStatistics[alpha->index].Phase;
                        const double q   = Phase.Fields(i,j,k)[alpha->index];
                        const double qp  = Phase.Fields(i+dir[0],j+dir[1],k+dir[2])[alpha->index];
                        const double qm  = Phase.Fields(i-dir[0],j-dir[1],k-dir[2])[alpha->index];
                        const double qpp = Phase.Fields(i+2*dir[0],j+2*dir[1],k+2*dir[2])[alpha->index];
                        const double qmm = Phase.Fields(i-2*dir[0],j-2*dir[1],k-2*dir[2])[alpha->index];
                        const double v   = Vel.Phase(i,j,k)({pIndex})[direction];
                        const double vp  = Vel.Phase(i+dir[0],j+dir[1],k+dir[2])({pIndex})[direction];
                        const double vm  = Vel.Phase(i-dir[0],j-dir[1],k-dir[2])({pIndex})[direction];
                        const double L0  = 0.5*Limiter(q-qm,qp-q);
                        const double Lp1 = 0.5*Limiter(qp-q,qpp-qp);
                        const double Lm1 = 0.5*Limiter(qm-qmm,q-qm);
                        const double ap  = 0.5*(vp+v);
                        const double am  = 0.5*(vm+v);
                        const double Fp  = std::max(ap,0.0)*(q+L0)
                                         + std::min(ap,0.0)*(qp-Lp1);
                        const double Fm  = std::max(am,0.0)*(qm+Lm1)
                                         + std::min(am,0.0)*(q-L0);
                        if (std::fabs(ap)*locdt > 0.49*dx || std::fabs(am)*locdt > 0.49*dx)
                        {
                            CFL = 1;
                        }
                        PhaseFieldUpdate(i,j,k).add_value(alpha->index, -(Fp-Fm)*locdt/dx);
                        PhaseFieldUpdate(i,j,k).add_value( beta->index,  (Fp-Fm)*locdt/dx);
                    }
                    else if (Phase.FieldsStatistics[beta->index].State  == AggregateStates::Solid and
                             Phase.FieldsStatistics[beta->index].Mobile == true)
                    {
                        const size_t Nx = Phase.Nx;
                        const size_t Ny = Phase.Ny;
                        const size_t Nz = Phase.Nz;

                        dMatrix3x3 W;
                        W(0,0) = 0.0;
                        W(1,1) = 0.0;
                        W(2,2) = 0.0;
                        W(0,1) = -Phase.FieldsStatistics[beta->index].aVel[2];
                        W(0,2) =  Phase.FieldsStatistics[beta->index].aVel[1];
                        W(1,2) = -Phase.FieldsStatistics[beta->index].aVel[0];
                        W(1,0) =  Phase.FieldsStatistics[beta->index].aVel[2];
                        W(2,0) = -Phase.FieldsStatistics[beta->index].aVel[1];
                        W(2,1) =  Phase.FieldsStatistics[beta->index].aVel[0];
                        const dVector3 pos  = {double(i),double(j),double(k)};
                        const dVector3 posp = {double(i+dir[0]),double(j+dir[1]),double(k+dir[2])};
                        const dVector3 posm = {double(i-dir[0]),double(j-dir[1]),double(k-dir[2])};
                        dVector3 distanceCM;
                        dVector3 distanceCMp;
                        dVector3 distanceCMm;
                        CommonFunctions::CalculateDistancePeriodic(pos, Phase.FieldsStatistics[beta->index].Rcm,distanceCM,  Nx, Ny, Nz);
                        CommonFunctions::CalculateDistancePeriodic(posp,Phase.FieldsStatistics[beta->index].Rcm,distanceCMp, Nx, Ny, Nz);
                        CommonFunctions::CalculateDistancePeriodic(posm,Phase.FieldsStatistics[beta->index].Rcm,distanceCMm, Nx, Ny, Nz);
                        const dVector3 R  = distanceCM *dx;
                        const dVector3 Rp = distanceCMp*dx;
                        const dVector3 Rm = distanceCMm*dx;
                        const dVector3 Vcm = Phase.FieldsStatistics[beta->index].Vcm;
                        const dVector3 Velocity  = Vcm + W*R;
                        const dVector3 Velocityp = Vcm + W*Rp;
                        const dVector3 Velocitym = Vcm + W*Rm;

                        const double q   = beta->value;
                        const double qp  = Phase.Fields(i+dir[0],j+dir[1],k+dir[2])[beta->index];
                        const double qm  = Phase.Fields(i-dir[0],j-dir[1],k-dir[2])[beta->index];
                        const double qpp = Phase.Fields(i+2*dir[0],j+2*dir[1],k+2*dir[2])[beta->index];
                        const double qmm = Phase.Fields(i-2*dir[0],j-2*dir[1],k-2*dir[2])[beta->index];
                        const double v   = Velocity[direction];
                        const double vp  = Velocityp[direction];
                        const double vm  = Velocitym[direction];
                        const double L0  = 0.5*Limiter(q-qm,qp-q);
                        const double Lp1 = 0.5*Limiter(qp-q,qpp-qp);
                        const double Lm1 = 0.5*Limiter(qm-qmm,q-qm);
                        const double ap  = 0.5*(vp+v);
                        const double am  = 0.5*(vm+v);
                        const double Fp  = std::max(ap,0.0)*(q+L0)
                                         + std::min(ap,0.0)*(qp-Lp1);
                        const double Fm  = std::max(am,0.0)*(qm+Lm1)
                                         + std::min(am,0.0)*(q-L0);
                        if (std::fabs(ap)*locdt > 0.49*dx || std::fabs(am)*locdt > 0.49*dx)
                        {
                            CFL = 1;
                        }
                        PhaseFieldUpdate(i,j,k).add_value( beta->index, -(Fp-Fm)*locdt/dx);
                        PhaseFieldUpdate(i,j,k).add_value(alpha->index,  (Fp-Fm)*locdt/dx);
                    }
                }
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            if (CFL == 0)
            {
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,offset,)
                if (Phase.Fields(i,j,k).flag)
                {
                    Phase.Fields(i,j,k) = PhaseFieldUpdate(i,j,k);
                }
                OMP_PARALLEL_STORAGE_LOOP_END
                Phase.Finalize(BC);
                //BC.SetX(Phase.Fields);
                //BC.SetY(Phase.Fields);
                //BC.SetZ(Phase.Fields);
            }
            else
            {
                locdt *= 0.5;
                Phase.Fields = PhaseFieldBackup;
            }
        }
        its *= 2;
    }
    while (CFL > 0);
}

//void AdvectionHR::Advection(
//        PhaseField& Phase,
//        const BoundaryConditions& BC,
//        const Velocities& Vel,
//        const int direction,
//        const double dt,
//        const double dx,
//        double (*Limiter)(const double,const double))
//{
//    Storage3D<NodePF, 0> PhaseFieldBackup(Phase.Fields);
//    Storage3D<NodePF, 0> PhaseFieldUpdate;
//    PhaseFieldUpdate.Allocate(Phase.Fields);
//
//    int CFL = 0;
//    double locdt = dt;
//    std::vector<int> dir(3);
//    dir[0] = ((direction == 0)and(Phase.dNx))?1:0;
//    dir[1] = ((direction == 1)and(Phase.dNy))?1:0;
//    dir[2] = ((direction == 2)and(Phase.dNz))?1:0;
//    const int offset = 0;// (Phase.Fields.Bcells()-2);
//    int its = 1;
//    do {
//        CFL = 0;
//        for (int it = 0; it < its; ++it)
//        {
//            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,offset,reduction(+:CFL))
//            if (Phase.Fields(i,j,k).flag)
//            {
////                size_t solids = 0;
////                for(auto iit = Phase.Fields(i,j,k).cbegin();
////                        iit != Phase.Fields(i,j,k).cend(); iit++)
////                if (Phase.FieldsStatistics[iit->index].State == AggregateStates::Solid)
////                {
////                    solids++;
////                }
//
////                if (solids <= 1)
////                {
////                    int pInd = 0;
////                    for (auto alpha = Phase.Fields(i,j,k).cbegin(); alpha < Phase.Fields(i,j,k).cend(); alpha++)
////                    {
////                        if(Phase.FieldsStatistics[alpha->index].State == AggregateState::Solid)
////                        {
////                            pInd = alpha->index;
////                        }
////                    }
////
////                    pInd = Phase.FieldsStatistics[pInd].Phase;
////
////                    for (auto alpha = Phase.Fields(i,j,k).begin(); alpha < Phase.Fields(i,j,k).end(); alpha++)
////                    {
////                        int ind = alpha->index;
////                        double q   = Phase.Fields(i,j,k)[ind];
////                        double qp  = Phase.Fields(i+dir[0],j+dir[1],k+dir[2])[ind];
////                        double qm  = Phase.Fields(i-dir[0],j-dir[1],k-dir[2])[ind];
////                        double qpp = Phase.Fields(i+2*dir[0],j+2*dir[1],k+2*dir[2])[ind];
////                        double qmm = Phase.Fields(i-2*dir[0],j-2*dir[1],k-2*dir[2])[ind];
////                        double v  = Vel.Phase(i,j,k)({pInd})[direction];
////                        double vp = Vel.Phase(i+dir[0],j+dir[1],k+dir[2])({pInd})[direction];
////                        double vm = Vel.Phase(i-dir[0],j-dir[1],k-dir[2])({pInd})[direction];
////                        double L0  = 0.5*Limiter(q-qm,qp-q);
////                        double Lp1 = 0.5*Limiter(qp-q,qpp-qp);
////                        double Lm1 = 0.5*Limiter(qm-qmm,q-qm);
////                        double ap = 0.5*(vp+v);
////                        double am = 0.5*(vm+v);
////                        double Fp = std::max(ap,0.)*(q+L0)
////                                  + std::min(ap,0.)*(qp-Lp1);
////                        double Fm = std::max(am,0.)*(qm+Lm1)
////                                  + std::min(am,0.)*(q-L0);
////                        if (std::fabs(ap)*locdt > 0.49*dx || std::fabs(am)*locdt > 0.49*dx)
////                        {
////                            CFL = 1;
////                        }
////                        PhaseFieldUpdate(i,j,k).set(ind, Phase.Fields(i,j,k)[ind] - locdt/dx*(Fp-Fm));
////                    }
////                }
////                else
////                {
//                    for (auto alpha  = Phase.Fields(i,j,k).begin();
//                              alpha != Phase.Fields(i,j,k).end(); alpha++)
//                    if (Phase.FieldsStatistics[alpha->index].State == AggregateStates::Solid and
//                        Phase.FieldsStatistics[alpha->index].Mobile)
//                    {
//                        const size_t Nx = Phase.Nx;
//                        const size_t Ny = Phase.Ny;
//                        const size_t Nz = Phase.Nz;
//
//                        dMatrix3x3 W;
//                        W(0,0) = 0.0;
//                        W(1,1) = 0.0;
//                        W(2,2) = 0.0;
//                        W(0,1) = -Phase.FieldsStatistics[alpha->index].aVel[2];
//                        W(0,2) =  Phase.FieldsStatistics[alpha->index].aVel[1];
//                        W(1,2) = -Phase.FieldsStatistics[alpha->index].aVel[0];
//                        W(1,0) =  Phase.FieldsStatistics[alpha->index].aVel[2];
//                        W(2,0) = -Phase.FieldsStatistics[alpha->index].aVel[1];
//                        W(2,1) =  Phase.FieldsStatistics[alpha->index].aVel[0];
//                        const dVector3 pos  = {double(i),double(j),double(k)};
//                        const dVector3 posp = {double(i+dir[0]),double(j+dir[1]),double(k+dir[2])};
//                        const dVector3 posm = {double(i-dir[0]),double(j-dir[1]),double(k-dir[2])};
//                        dVector3 distanceCM;
//                        dVector3 distanceCMp;
//                        dVector3 distanceCMm;
//                        CommonFunctions::CalculateDistancePeriodic(pos, Phase.FieldsStatistics[alpha->index].Rcm,distanceCM,  Nx, Ny, Nz);
//                        CommonFunctions::CalculateDistancePeriodic(posp,Phase.FieldsStatistics[alpha->index].Rcm,distanceCMp, Nx, Ny, Nz);
//                        CommonFunctions::CalculateDistancePeriodic(posm,Phase.FieldsStatistics[alpha->index].Rcm,distanceCMm, Nx, Ny, Nz);
//                        const dVector3 R  = distanceCM *dx;
//                        const dVector3 Rp = distanceCMp*dx;
//                        const dVector3 Rm = distanceCMm*dx;
//                        const dVector3 Vcm = Phase.FieldsStatistics[alpha->index].Vcm;
//                        const dVector3 Velocity  = Vcm + W*R;
//                        const dVector3 Velocityp = Vcm + W*Rp;
//                        const dVector3 Velocitym = Vcm + W*Rm;
//
//                        const size_t ind = alpha->index;
//                        const double q   = Phase.Fields(i,j,k)[ind];
//                        const double qp  = Phase.Fields(i+dir[0],j+dir[1],k+dir[2])[ind];
//                        const double qm  = Phase.Fields(i-dir[0],j-dir[1],k-dir[2])[ind];
//                        const double qpp = Phase.Fields(i+2*dir[0],j+2*dir[1],k+2*dir[2])[ind];
//                        const double qmm = Phase.Fields(i-2*dir[0],j-2*dir[1],k-2*dir[2])[ind];
//                        const double v   = Velocity[direction];
//                        const double vp  = Velocityp[direction];
//                        const double vm  = Velocitym[direction];
//                        const double L0  = 0.5*Limiter(q-qm,qp-q);
//                        const double Lp1 = 0.5*Limiter(qp-q,qpp-qp);
//                        const double Lm1 = 0.5*Limiter(qm-qmm,q-qm);
//                        const double ap  = 0.5*(vp+v);
//                        const double am  = 0.5*(vm+v);
//                        const double Fp  = std::max(ap,0.0)*(q+L0)
//                                         + std::min(ap,0.0)*(qp-Lp1);
//                        const double Fm  = std::max(am,0.0)*(qm+Lm1)
//                                         + std::min(am,0.0)*(q-L0);
//                        if (std::fabs(ap)*locdt > 0.49*dx || std::fabs(am)*locdt > 0.49*dx)
//                        {
//                            CFL = 1;
//                        }
//                        PhaseFieldUpdate(i,j,k).set_value(ind, Phase.Fields(i,j,k)[ind] - (Fp-Fm)*locdt/dx);
//                    }
//                    else
//                    {
//                        const size_t ind = alpha->index;
//                        const double q   = Phase.Fields(i,j,k)[ind];
//                        const double qp  = Phase.Fields(i+dir[0],j+dir[1],k+dir[2])[ind];
//                        const double qm  = Phase.Fields(i-dir[0],j-dir[1],k-dir[2])[ind];
//                        const double qpp = Phase.Fields(i+2*dir[0],j+2*dir[1],k+2*dir[2])[ind];
//                        const double qmm = Phase.Fields(i-2*dir[0],j-2*dir[1],k-2*dir[2])[ind];
//                        const double v   = Vel.Phase(i,j,k)({ind})[direction];
//                        const double vp  = Vel.Phase(i+dir[0],j+dir[1],k+dir[2])({ind})[direction];
//                        const double vm  = Vel.Phase(i-dir[0],j-dir[1],k-dir[2])({ind})[direction];
//                        const double L0  = 0.5*Limiter(q-qm,qp-q);
//                        const double Lp1 = 0.5*Limiter(qp-q,qpp-qp);
//                        const double Lm1 = 0.5*Limiter(qm-qmm,q-qm);
//                        const double ap  = 0.5*(vp+v);
//                        const double am  = 0.5*(vm+v);
//                        const double Fp  = std::max(ap,0.0)*(q+L0)
//                                         + std::min(ap,0.0)*(qp-Lp1);
//                        const double Fm  = std::max(am,0.0)*(qm+Lm1)
//                                         + std::min(am,0.0)*(q-L0);
//                        if (std::fabs(ap)*locdt > 0.49*dx || std::fabs(am)*locdt > 0.49*dx)
//                        {
//                            CFL = 1;
//                        }
//                        PhaseFieldUpdate(i,j,k).set_value(ind, Phase.Fields(i,j,k)[ind] - (Fp-Fm)*locdt/dx);
//                    }
//                //}
//            }
//            OMP_PARALLEL_STORAGE_LOOP_END
//            if (CFL == 0)
//            {
//                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,offset,)
//                if (Phase.Fields(i,j,k).flag)
//                {
//                    Phase.Fields(i,j,k) = PhaseFieldUpdate(i,j,k);
//                }
//                /*else
//                {
//                    Phase.Fields(i,j,k).clear();
//                    Phase.Fields(i,j,k).add_value(Phase.Fields(i,j,k).front().index, 1.0);
//                }*/
//                OMP_PARALLEL_STORAGE_LOOP_END
//                BC.SetX(Phase.Fields);
//                BC.SetY(Phase.Fields);
//                BC.SetZ(Phase.Fields);
//            }
//            else
//            {
//                locdt *= 0.5;
//                Phase.Fields = PhaseFieldBackup;
//            }
//        }
//        its *= 2;
//    }
//    while (CFL > 0);
//}

void AdvectionHR::AdvectPhaseField(
        PhaseField& Phase,
        const Velocities& Vel,
        const BoundaryConditions& BC,
        const double dx,
        const double dt,
        const int tStep,
        const bool finalize)
{
    if (Phase.Fields.Bcells() < 2)
    {
        Info::WriteExit("Number of Bcells for storage PhaseField::Fields needs to be 2 or higher.",
                        "AdvectionHR", "AdvectPhaseField()");
        raise(SIGABRT);
    }

    BC.SetX(Phase.Fields);
    BC.SetY(Phase.Fields);
    BC.SetZ(Phase.Fields);
    double (*Limiter)(const double,const double);
    switch(Scheme)
    {
        case AdvectionSchemes::Minmod:
        {
            Limiter = &Slope_Limiter_Minmod;
            break;
        }
        case AdvectionSchemes::MonotonizedCentral:
        {
            Limiter = &Slope_Limiter_MC;
            break;
        }
        case AdvectionSchemes::Superbee:
        {
            Limiter = &Slope_Limiter_Superbee;
            break;
        }
        /*case AdvectionSchemes::LaxWendroff:
        {
            Limiter = &Slope_Limiter_LaxWendroff;
            break;
        }*/
        case AdvectionSchemes::Upwind:
        default:
        {
            Limiter = &Slope_Limiter_Upwind;
            break;
        }
    }
    if (tStep % 2 == 0)
    {
        Advection(Phase, BC, Vel, 0, dt, dx, Limiter);
        Advection(Phase, BC, Vel, 1, dt, dx, Limiter);
        Advection(Phase, BC, Vel, 2, dt, dx, Limiter);
    }
    else
    {
        Advection(Phase, BC, Vel, 2, dt, dx, Limiter);
        Advection(Phase, BC, Vel, 1, dt, dx, Limiter);
        Advection(Phase, BC, Vel, 0, dt, dx, Limiter);
    }

    // Do special kind of finalization
    const int offset = Phase.Fields.Bcells();
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,offset,)
    {
        // restrict the phase-field values to the range [0,1]
        if(Phase.Fields(i,j,k).flag)
        {
            for(auto it = Phase.Fields(i,j,k).begin();
                     it != Phase.Fields(i,j,k).end();)
            {
                bool erase = false;
                if (it->value >= 1.0 - DBL_EPSILON) it->value = 1.0;
                if (it->value <= 0.0 + DBL_EPSILON)
                {
                    it->value = 0.0;
                    erase = true;
                }

                if (erase)
                {
                    it = Phase.Fields(i,j,k).erase(it);
                }
                else
                {
                    ++it;
                }
            }

            double total = 0.0;
            for(auto it = Phase.Fields(i,j,k).begin();
                     it != Phase.Fields(i,j,k).end(); it++)
            {
                total += it->value;
            }

            // Determine number solid and fluid phase-fields and the
            // respective solid and fluid fractions
            size_t FluidPhases = 0;
            //size_t SolidPhases = 0;
            double SolidFraction = 0;
            //double FluidFraction = 0;
            for(auto it = Phase.Fields(i,j,k).cbegin();
                     it != Phase.Fields(i,j,k).cend(); it++)
            if (Phase.FieldsStatistics[it->index].State != AggregateStates::Solid)
            {
                //FluidFraction += it->value;
                FluidPhases++;
            }
            else
            {
                SolidFraction += it->value;
                //SolidPhases++;
            }

            if (SolidFraction > 1)
            {
                for(auto it = Phase.Fields(i,j,k).begin();
                        it != Phase.Fields(i,j,k).end();)
                if (Phase.FieldsStatistics[it->index].State != AggregateStates::Solid)
                {
                    it = Phase.Fields(i,j,k).erase(it);
                }
                else
                {
                    ++it;
                }
            }
            else if (FluidPhases > 0)
            {
                //if there are solid and fluid phase-fields present
                // normalize only the fluid phases
                const double norm = (1.0 - total)/FluidPhases;
                for(auto it = Phase.Fields(i,j,k).begin();
                            it != Phase.Fields(i,j,k).end(); it++)
                if (Phase.FieldsStatistics[it->index].State != AggregateStates::Solid)
                {
                    it->value += norm;
                    it->laplacian = 0.0;
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    BC.SetX(Phase.Fields);
    BC.SetY(Phase.Fields);
    BC.SetZ(Phase.Fields);
    Phase.Finalize(BC,finalize);
}
}// namespace openphase
