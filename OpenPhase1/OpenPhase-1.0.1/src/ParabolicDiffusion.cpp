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
 *   File created :   2022
 *   Main contributors :   Raphael Schiedung
 *
 */

#include "Base/UserInterface.h"
#include "BoundaryConditions.h"
#include "Composition.h"
#include "DrivingForce.h"
#include "ParabolicDiffusion.h"
#include "GrainInfo.h"
#include "Info.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "PhysicalConstants.h"
#include "Settings.h"
#include "Temperature.h"
#include "VTK.h"
#include <cassert>

namespace openphase
{

ParabolicDiffusion::ParabolicDiffusion(Settings& locSettings, std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}
void ParabolicDiffusion::Initialize(Settings& locSettings)
{
    TotalNx = locSettings.TotalNx;

    Nx  = locSettings.Nx;
    Ny  = locSettings.Ny;
    Nz  = locSettings.Nz;

    dNx = locSettings.dNx;
    dNy = locSettings.dNy;
    dNz = locSettings.dNz;

    ElementNames = locSettings.ElementNames;
    Ncomp        = locSettings.Ncomp;
    Nphases      = locSettings.Nphases;
    PhaseNames   = locSettings.PhaseNames;
    RawDataDir   = locSettings.RawDataDir;
    TextDir      = locSettings.TextDir;
    VTKDir       = locSettings.VTKDir;
    dx           = locSettings.dx;

    C0                       .Allocate({Nphases, Ncomp});
    ConstantIsotropicMobility.Allocate({Ncomp});
    EPS0                     .Allocate({Nphases, Ncomp});
    EPS1                     .Allocate({Nphases, Ncomp});
    EPS2                     .Allocate({Nphases, Ncomp});
    PhaseDiffusionCoefficient.Allocate({Nphases, Ncomp});
    PhaseMobility            .Allocate({Nphases, Ncomp});

    C0                       .set_to_value(0.0);
    EPS0                     .set_to_value(0.0);
    EPS1                     .set_to_value(0.0);
    EPS2                     .set_to_value(0.0);
    PhaseDiffusionCoefficient.set_to_value(0.0);
    PhaseMobility            .set_to_value(0.0);

    switch(locSettings.ActiveDimensions())
    {
        case 1:
        {
            LStencil.Set(LaplacianStencil1D_3, dx,dNx,dNy,dNz);
            GStencil.Set(GradientStencil1D,    dx,dNx,dNy,dNz);
            break;
        }
        case 2:
        {
            switch(locSettings.DiffusionStencil)
            {
                case LaplacianStencils::Simple:
                {
                    LStencil.Set(LaplacianStencil2D_5, dx,dNx,dNy,dNz);
                    GStencil.Set(GradientStencil1D,    dx,dNx,dNy,dNz);
                    break;
                }
                case LaplacianStencils::Isotropic:
                {
                    LStencil.Set(LaplacianStencil2D_9, dx,dNx,dNy,dNz);
                    GStencil.Set(GradientStencil2D,    dx,dNx,dNy,dNz);
                    break;
                }
                case LaplacianStencils::LB:
                {
                    LStencil.Set(LaplacianStencil2D_LB, dx,dNx,dNy,dNz);
                    GStencil.Set(GradientStencil2D_LB,  dx,dNx,dNy,dNz);
                    break;
                }
            }
            break;
        }
        case 3:
        {
            switch(locSettings.DiffusionStencil)
            {
                case LaplacianStencils::Simple:
                {
                    LStencil.Set(LaplacianStencil3D_7, dx);
                    GStencil.Set(GradientStencil1D,    dx);
                    break;
                }
                case LaplacianStencils::Isotropic:
                {
                    LStencil.Set(LaplacianStencil3D_27a, dx);
                    GStencil.Set(GradientStencil3D,      dx);
                    break;
                }
                case LaplacianStencils::LB:
                {
                    LStencil.Set(LaplacianStencil3D_LB, dx);
                    GStencil.Set(GradientStencil3D_LB,  dx);
                    break;
                }
            }
            break;
        }
    }

    Bcells = locSettings.Bcells;
    if (Bcells < 2)
    {
        Info::WriteExit("ParabolicDiffusion::Initialize: Too few boundary cell! Two are required!");
        exit(EXIT_FAILURE);
    }

    ChemicalPotentialOld.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Ncomp}, Bcells);
    ChemicalPotential   .Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Ncomp}, Bcells);
    Flux                .Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Ncomp}, Bcells);
    Mobility            .Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Ncomp}, Bcells);

    Info::WriteStandard("ParabolicDiffusion", "Initialized");
}
void ParabolicDiffusion::ReadInput(const std::string InputFileName)
{
    // Check if there is reference component for which the diffusion should not be solved
    Info::WriteLineInsert("Parabolic Diffusion input");
    Info::WriteStandard("Source", InputFileName);

    std::fstream inpF(InputFileName.c_str(), std::ios::in | std::ios_base::binary);
    if (!inpF)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        std::exit(1);
    };
    std::stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();

    int moduleLocation = UserInterface::FindModuleLocation(inp, "ParabolicDiffusion");

    UseLaplacian = UserInterface::ReadParameterB(inp, moduleLocation, "UseLaplacian");

    std::string RefCompName = UserInterface::ReadParameterK(inp, moduleLocation, "RefElement");
    for(size_t n = 0; n < Ncomp; n++)
    if(ElementNames[n] == RefCompName) RefComp = n;

    for(size_t PhaseIdx = 0; PhaseIdx < Nphases; PhaseIdx++)
    for(size_t comp = 0; comp < Ncomp; comp++)
    if (comp != RefComp)
    {
        std::string converter = ""+PhaseNames[PhaseIdx]+"_"+ElementNames[comp];
        EPS0({PhaseIdx, comp}) = UserInterface::ReadParameterD(inp, moduleLocation, "EPS0_"+converter);
        EPS1({PhaseIdx, comp}) = UserInterface::ReadParameterD(inp, moduleLocation, "EPS1_"+converter);
        EPS2({PhaseIdx, comp}) = UserInterface::ReadParameterD(inp, moduleLocation, "EPS2_"+converter);
    }
    for(size_t PhaseIdx = 0; PhaseIdx < Nphases; PhaseIdx++)
    for(size_t comp = 0; comp < Ncomp; comp++)
    if (comp != RefComp)
    {
        std::string converter = ""+PhaseNames[PhaseIdx]+"_"+ElementNames[comp];
        C0({PhaseIdx, comp}) = UserInterface::ReadParameterD(inp, moduleLocation, "C0_"+converter);
    }
    for(size_t PhaseIdx = 0; PhaseIdx < Nphases; PhaseIdx++)
    for(size_t comp     = 0; comp     < Ncomp;       comp++)
    if (comp != RefComp)
    {
        PhaseDiffusionCoefficient({PhaseIdx, comp}) = UserInterface::ReadParameterD(inp, moduleLocation, "D0_"+PhaseNames[PhaseIdx]+"_"+ElementNames[comp]);
        PhaseMobility ({PhaseIdx, comp}) = PhaseDiffusionCoefficient({PhaseIdx, comp})/EPS2({PhaseIdx, comp});
    }
    for(size_t comp = 0; comp < Ncomp; comp++)
    if (comp != RefComp)
    {
        auto GetConstantIsotropicMobility = [&]()
        {
            for(size_t alpha = 0;         alpha < Nphases; alpha++)
            for(size_t beta  = alpha + 1; beta  < Nphases; beta++)
            {
                if (PhaseMobility({alpha, comp}) != PhaseMobility({beta, comp}))
                {
                    return 0.0;
                }
            }
            return PhaseMobility({0, comp});
        };
        ConstantIsotropicMobility({comp}) = GetConstantIsotropicMobility();
    }

    Info::WriteLine();
    Info::WriteBlankLine();
}
void ParabolicDiffusion::CalculateMobility(PhaseField& Phase, size_t comp)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Mobility,Mobility.Bcells(),)
    {
        Mobility(i,j,k)({comp}) = 0;
        for(auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        {
            size_t indexA = Phase.FieldsStatistics[alpha->index].Phase;
            Mobility(i,j,k)({comp}) += alpha->value*PhaseMobility ({indexA, comp});
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void ParabolicDiffusion::GetDrivingForce(PhaseField& Phase, Composition& Elements, Temperature& Tx, DrivingForce& dGab)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,)
    {
        if(Phase.Interface(i,j,k))
        {
            for(size_t comp = 0; comp < Ncomp; comp++)
            if (comp != RefComp)
            for(auto alpha  = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
            for(auto  beta  = alpha + 1; beta < Phase.Fields(i,j,k).cend(); ++beta)
            {
                size_t indexA = Phase.FieldsStatistics[alpha->index].Phase;
                size_t indexB = Phase.FieldsStatistics[ beta->index].Phase;
                if(indexA != indexB)
                {
                    double dCA = Elements.MoleFractionsTotal(i,j,k)({comp}) - C0({indexA, comp});
                    double dCB = Elements.MoleFractionsTotal(i,j,k)({comp}) - C0({indexB, comp});
                    double f_A = 0.5 * EPS2({indexA, comp}) * dCA *dCA
                                     + EPS1({indexA, comp}) * dCA
                                     + EPS0({indexA, comp});
                    double f_B = 0.5 * EPS2({indexB, comp}) * dCB *dCB
                                     + EPS1({indexB, comp}) * dCB
                                     + EPS0({indexB, comp});
                    double dG_AB = f_B - f_A;
                    dGab.Raw(i,j,k).add_asym1(alpha->index,  beta->index, dG_AB);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void ParabolicDiffusion::CalculateChemicalPotential(PhaseField& Phase, Composition& Elements, Temperature& Tx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,Bcells,)
    {
        for(size_t comp = 0; comp < Ncomp; comp++)
        if (comp != RefComp)
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            {
                size_t indexA = Phase.FieldsStatistics[alpha->index].Phase;
                double dC     = Elements.MoleFractionsTotal(i,j,k)({comp}) - C0({indexA, comp});
                ChemicalPotential(i,j,k)({comp}) += alpha->value*(EPS2({indexA, comp})*dC + EPS1({indexA, comp}));
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void ParabolicDiffusion::CaluculateDiffusionFlux(size_t comp)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,Bcells-1,)
    {
        Flux(i,j,k)({comp})[0] = 0.0;
        Flux(i,j,k)({comp})[1] = 0.0;
        Flux(i,j,k)({comp})[2] = 0.0;
        for (auto gs = GStencil.cbegin(); gs != GStencil.cend(); ++gs)
        {
            const int ii = i + gs->di;
            const int jj = j + gs->dj;
            const int kk = k + gs->dk;
            Flux(i,j,k)({comp})[0] -= gs->weightX * Mobility(i,j,k)({comp}) * ChemicalPotential(ii,jj,kk)({comp});
            Flux(i,j,k)({comp})[1] -= gs->weightY * Mobility(i,j,k)({comp}) * ChemicalPotential(ii,jj,kk)({comp});
            Flux(i,j,k)({comp})[2] -= gs->weightZ * Mobility(i,j,k)({comp}) * ChemicalPotential(ii,jj,kk)({comp});
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void ParabolicDiffusion::ApplyIncrements(PhaseField& Phase, Composition& Elements, double dt, size_t comp)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,)
    {
        if (dNx > 0) Elements.MoleFractionsTotal(i,j,k)({comp}) -= (Flux(i+1,j,k)({comp})[0] - Flux(i-1,j,k)({comp})[0])/dx/2*dt;
        if (dNy > 0) Elements.MoleFractionsTotal(i,j,k)({comp}) -= (Flux(i,j+1,k)({comp})[1] - Flux(i,j-1,k)({comp})[1])/dx/2*dt;
        if (dNz > 0) Elements.MoleFractionsTotal(i,j,k)({comp}) -= (Flux(i,j,k+1)({comp})[2] - Flux(i,j,k-1)({comp})[2])/dx/2*dt;

        assert(Elements.MoleFractionsTotal(i,j,k)({comp}) <=  1.001);
        assert(Elements.MoleFractionsTotal(i,j,k)({comp}) >= -0.001);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void ParabolicDiffusion::ApplyIncrementsConstantIsotropicMobility(PhaseField& Phase, Composition& Elements, double dt, size_t comp)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,)
    {
        double DeltaMFT = 0.0;
        for(auto ds = LStencil.begin(); ds != LStencil.end(); ds++)
        {
            const int ii = i + ds->di;
            const int jj = j + ds->dj;
            const int kk = k + ds->dk;
            DeltaMFT += ds->weight*ChemicalPotential(ii,jj,kk)({comp});
        }
        DeltaMFT *= ConstantIsotropicMobility({comp})*dt;
        Elements.MoleFractionsTotal(i,j,k)({comp}) += DeltaMFT;

        assert(Elements.MoleFractionsTotal(i,j,k)({comp}) <=  1.001);
        assert(Elements.MoleFractionsTotal(i,j,k)({comp}) >= -0.001);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void ParabolicDiffusion::CalculateReferenceElementMoleFractions(Composition& Elements)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,)
    {
        Elements.MoleFractionsTotal(i,j,k)({RefComp}) = 1.0;
        for(size_t comp = 0; comp < Ncomp; comp++)
        if (comp != RefComp)
        {
            Elements.MoleFractionsTotal(i,j,k)({RefComp}) -= Elements.MoleFractionsTotal(i,j,k)({comp});
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void ParabolicDiffusion::SetPhaseMolarFractions(Composition& Elements)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,)
    {
        for(size_t comp = 0; comp < Ncomp; comp++)
        for(size_t phase = 0; phase < Nphases; phase++)
        {
            //NOTE this model has no phase molar fractions! They set to unsure compatibility with the rest of openphase
            Elements.MoleFractions(i,j,k)({phase,comp}) = Elements.MoleFractionsTotal(i,j,k)({comp});
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ParabolicDiffusion::ResetChemicalPotential()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,Bcells,)
    {
        for(size_t comp = 0; comp < Ncomp; comp++)
        if (comp != RefComp)
        {
            ChemicalPotentialOld(i,j,k)({comp}) = ChemicalPotential(i,j,k)({comp});
            ChemicalPotential   (i,j,k)({comp}) = 0.0;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void ParabolicDiffusion::Solve(PhaseField& Phase, Composition& Elements, Temperature& Tx, BoundaryConditions& BC, double dt)
{
    Elements.SetBoundaryConditions(BC);
    CalculateChemicalPotential(Phase, Elements, Tx);
    for(size_t comp = 0; comp < Ncomp; comp++)
    if (comp != RefComp)
    {
         if (UseLaplacian and ConstantIsotropicMobility({comp}) != 0.0)
         {
             ApplyIncrementsConstantIsotropicMobility(Phase, Elements, dt, comp);
         }
         else
         {
             CalculateMobility(Phase, comp);
             CaluculateDiffusionFlux(comp);
             ApplyIncrements(Phase, Elements, dt, comp);
         }
    }
    ResetChemicalPotential();
    CalculateReferenceElementMoleFractions(Elements);
    SetPhaseMolarFractions(Elements);
}
void ParabolicDiffusion::WriteVTK(const int tStep, const Settings& locSettings, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        ListOfFields.push_back((VTK::Field_t) {"Chemical_Potential_" + ElementNames[comp], [this, comp](int i,int j,int k){return ChemicalPotential(i,j,k)({comp});}});
        ListOfFields.push_back((VTK::Field_t) {"Flux_"               + ElementNames[comp], [this, comp](int i,int j,int k){return              Flux(i,j,k)({comp});}});
    }
    std::string Filename = UserInterface::MakeFileName(VTKDir, "ParabolicDiffusion_", tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}
}// namespace openphase
