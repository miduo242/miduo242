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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */


#include "EquilibriumPartitionDiffusionBinary.h"
#include "Base/UserInterface.h"
#include "BoundaryConditions.h"
#include "Composition.h"
#include "DrivingForce.h"
#include "GrainInfo.h"
#include "Info.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "PhysicalConstants.h"
#include "Settings.h"
#include "Temperature.h"
#include "VTK.h"

namespace openphase
{

using namespace std;
/*************************************************************************/
void EquilibriumPartitionDiffusionBinary::Initialize(Settings& locSettings)
{
    thisclassname = "EquilibriumPartitionDiffusion";

    Nphases = locSettings.Nphases;
    Ncomp   = locSettings.Ncomp;
    dx      = locSettings.dx;
    Eta     = locSettings.Eta;
    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;

    dNx = locSettings.dNx;
    dNy = locSettings.dNy;
    dNz = locSettings.dNz;

    ElementNames = locSettings.ElementNames;

    Comp = 0;
    R = PhysicalConstants::R;
    TotalMass = 0.0;
    maxDC = 0.0;

    Precision = DBL_EPSILON;
    EnableAntiTrapping = true;

    size_t Bcells = locSettings.Bcells;
    dMu.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Nphases}, Bcells);
    DC.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Nphases}, Bcells);

    Pt.Allocate({Nphases, Nphases});
    Ts.Allocate({Nphases, Nphases});
    Cs.Allocate({Nphases, Nphases});
    mL.Allocate({Nphases, Nphases});
    mL_1.Allocate({Nphases, Nphases});

    DC0.resize(Nphases);
    AE.resize(Nphases);
    Stoichiometric.resize(Nphases);
    S.resize(Nphases);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,dMu, Bcells,)
    {
        for (size_t n = 0; n < Nphases; ++n)
        {
            dMu(i,j,k)({n}) = 0.0;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    LimitingNeeded = false;
    switch(locSettings.ActiveDimensions())
    {
        case 1:
        {
            DStencil.SetNoCenter(LaplacianStencil1D_3, dx,dNx,dNy,dNz);
            break;
        }
        case 2:
        {
            if(locSettings.DiffusionStencil == LaplacianStencils::Simple)
            {
                DStencil.SetNoCenter(LaplacianStencil2D_5, dx,dNx,dNy,dNz);
            }
            else
            {
                DStencil.SetNoCenter(LaplacianStencil2D_9, dx,dNx,dNy,dNz);
            }
            break;
        }
        case 3:
        {
            if(locSettings.DiffusionStencil == LaplacianStencils::Simple)
            {
                DStencil.SetNoCenter(LaplacianStencil3D_7, dx);
            }
            else
            {
                DStencil.SetNoCenter(LaplacianStencil3D_27a, dx);
            }
            break;
        }
        default:
        {

        }
    }
    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void EquilibriumPartitionDiffusionBinary::ReadInput(const string InputFileName)
{
    Info::WriteLineInsert("EquilibriumPartitionDiffusionBinary input");
    Info::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };

    Info::WriteLineInsert("EquilibriumPartitionDiffusion properties");
    Info::WriteStandard("Source", InputFileName);

    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);

    inp.close();
}

void EquilibriumPartitionDiffusionBinary::ReadInput(stringstream& inp)
{
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    string RefCompName = UserInterface::ReadParameterK(inp, moduleLocation, string("RefElement"), true, "NN");
    EnableAntiTrapping = UserInterface::ReadParameterB(inp, moduleLocation, string("EnableAT"), false, true);

    bool comp_set = false;
    for(size_t n = 0; n < Ncomp; n++)
    if(ElementNames[n] == RefCompName)
    {
        RefComp = n;
        comp_set = true;
    }
    else
    {
        Comp = n;
    }
    if(!comp_set)
    {
        Info::WriteExit("Nonexistent reference element is selected",thisclassname, "ReadInput()");
        exit(1);
    }
    for(size_t n = 0; n < Nphases; ++n)
    {
        stringstream converter;
        converter << n;
        string counter = converter.str();
        Stoichiometric[n] = UserInterface::ReadParameterB(inp, moduleLocation, string("Flag_") + counter);
        DC0[n]            = UserInterface::ReadParameterD(inp, moduleLocation, string("DC_") + counter);
        AE[n]             = UserInterface::ReadParameterD(inp, moduleLocation, string("AE_") + counter);
        S[n]              = UserInterface::ReadParameterD(inp, moduleLocation, string("EF_") + counter);
    }
    ThereAreStoichiometricPhases = false;
    for(size_t n = 0; n < Nphases; ++n)
    {
        if (Stoichiometric[n])
        {
            ThereAreStoichiometricPhases = true;
            break;
        }
    }

    for(size_t n =   0; n < Nphases-1; ++n)
    for(size_t m = n+1; m < Nphases  ; ++m)
    {
        stringstream converter;
        converter << n << "_" << m;
        string counter = converter.str();
        Cs({n, m})    = Cs({m, n})   = UserInterface::ReadParameterD(inp, moduleLocation, string("Cs_") + counter);
        Ts({n, m})    = Ts({m, n})   = UserInterface::ReadParameterD(inp, moduleLocation, string("Ts_") + counter);
        mL({n, m})    = UserInterface::ReadParameterD(inp, moduleLocation, string("ML_") + counter);
        stringstream converter2;
        converter2 << m << "_" << n;
        string counter2 = converter2.str();
        mL({m, n})    = UserInterface::ReadParameterD(inp, moduleLocation, string("ML_") + counter2);
    }

    for(size_t n = 0; n < Nphases; ++n)
    {
        Pt({n, n}) = 1.0;
        mL({n, n}) = 1.0;
        mL_1({n, n}) = 1.0;
        Ts({n, n}) = 0.0;
        Cs({n, n}) = 0.0;
    };

    for(size_t n = 0; n < Nphases; ++n)
    for(size_t m = 0; m < Nphases; ++m)
    if(n != m)
    {
        if(mL({n, m}) == 0.0)
        {
            Pt({m, n}) = 0.0;
            mL_1({n, m}) = 0.0;
        }
        else
        {
            mL_1({n, m}) = 1.0/mL({n, m});
        }
        Pt({n, m}) = fabs(mL({m, n}) * mL_1({n, m}));
    }

    Info::WriteLine();
    Info::WriteBlankLine();
}

inline void EquilibriumPartitionDiffusionBinary::CalculateLocalPhaseConcentrations(PhaseField& Phase,
                                                             Temperature& Tx, Composition& Elements,
                                                             int i, int j, int k)
{
    size_t locNphases = 0;
    size_t locPindex = 0;
    for(size_t n = 0; n < Nphases; n++)
    if(Phase.Fractions(i,j,k)({n}) > 0.0)
    {
        locNphases += 1;
        locPindex = n;
    }

    if(locNphases == 1)
    //if(!Phase.Interface(i,j,k))
    {
        size_t pIndex = locPindex;
        Elements.MoleFractions(i,j,k)({pIndex, Comp}) = Elements.MoleFractionsTotal(i,j,k)({Comp});
        if(!Stoichiometric[pIndex])
        {
            for(size_t n = 0; n < Nphases; ++n)
            if(n != pIndex)
            {
                Elements.MoleFractions(i,j,k)({n, Comp}) = Cs({n, pIndex}) +
                         (Elements.MoleFractions(i,j,k)({pIndex, Comp}) - Cs({n, pIndex})) *
                         mL({pIndex, n})*mL_1({n, pIndex});
            }
        }
        else
        {
            for(size_t n = 0; n < Nphases; ++n)
            if(n != pIndex)
            {
                Elements.MoleFractions(i,j,k)({n, Comp}) = Cs({n, pIndex}) +
                                    (Tx(i,j,k) - Ts({n, pIndex}))*mL_1({n, pIndex});
            }
        }
    }
    else
    {
        Tensor <double, 2> eqCx({Nphases, Nphases});

        for(size_t n = 0; n < Nphases; ++n)
        for(size_t m = 0; m < Nphases; ++m)
        if(n != m)
        {
            eqCx({n, m}) = Cs({n, m}) + (Tx(i,j,k) - Ts({n, m}))*mL_1({n, m});
        }

        double dCx = Elements.MoleFractionsTotal(i,j,k)({Comp});
        for(size_t n = 0; n < Nphases; ++n)
        {
            Elements.MoleFractions(i,j,k)({n, Comp}) = 0.0;
            eqCx({n,n}) = 0.0;
            double div = 0.0;
            for(size_t m = 0; m < Nphases; ++m)
            if(n != m)
            {
                eqCx({n,n}) += eqCx({n, m})*Phase.Fractions(i,j,k)({m});
                div += Phase.Fractions(i,j,k)({m});
            }

            if(div > Precision)
            {
                eqCx({n,n}) /= div;
            }
            else
            {
                eqCx({n,n}) = Elements.MoleFractionsTotal(i,j,k)({Comp});
            }
            dCx -= eqCx({n,n})*Phase.Fractions(i,j,k)({n});
        }

        for(size_t n = 0; n < Nphases; ++n)
        if (!Stoichiometric[n])
        {
            double denom = 0.0;
            double num = 1.0;
            for(size_t m = 0; m < Nphases; ++m)
            {
                if(!Stoichiometric[m])
                {
                    denom += Pt({m, n})*Phase.Fractions(i,j,k)({m});
                }
                else
                {
                    num -= Phase.Fractions(i,j,k)({m});
                }
            }

            Elements.MoleFractions(i,j,k)({n, Comp}) = eqCx({n,n});

            if(fabs(denom) > 0.0)
            {
                Elements.MoleFractions(i,j,k)({n, Comp}) += num*dCx/denom;
            }
        }
        else
        {
            Elements.MoleFractions(i,j,k)({n, Comp}) = eqCx({n,n});
        }
    }
}

void EquilibriumPartitionDiffusionBinary::CalculatePhaseConcentrations(PhaseField& Phase,
                                    Temperature& Tx, Composition& Elements)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,Elements.MoleFractionsTotal.Bcells(),)
    {
        if(Phase.Fields(i,j,k).flag)
        {
            CalculateLocalPhaseConcentrations(Phase, Tx, Elements, i, j, k);
        }
        else
        {
            size_t pIndex = Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase;
            Elements.MoleFractions(i,j,k)({pIndex, Comp}) = Elements.MoleFractionsTotal(i,j,k)({Comp});
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::SetDiffusionCoefficients(PhaseField& Phase,
                                                                Temperature& Tx)
{
    maxDC = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DC,DC.Bcells(),reduction(max:maxDC))
    {
        double invRT  = 1.0/(R*Tx(i,j,k));
        for (size_t n = 0; n < Nphases; ++n)
        {
            DC(i, j, k)({n}) = DC0[n] * exp(-AE[n] * invRT);
            maxDC = max(DC(i, j, k)({n}), maxDC);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::GetDrivingForce(PhaseField& Phase,
                                                    Composition& Elements,
                                                    Temperature& Tx,
                                                    DrivingForce& dGab)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,)
    {
        if(Phase.Interface(i,j,k))
        {
            CalculateLocalPhaseConcentrations(Phase, Tx, Elements, i, j, k);

            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
            for(auto  beta = alpha + 1;
                      beta != Phase.Fields(i,j,k).cend(); ++beta)
            {
                size_t index1 = Phase.FieldsStatistics[alpha->index].Phase;
                size_t index2 = Phase.FieldsStatistics[ beta->index].Phase;

                if(index1 != index2)
                {
                    double dG_AB = 0.5*((Ts({index1, index2}) + mL({index1, index2})*
                                       (Elements.MoleFractions(i,j,k)({index1, Comp}) - Cs({index1, index2})) - Tx(i,j,k))*
                                       (1 - Stoichiometric[index1])*(1 + Stoichiometric[index2])
                                        +
                                       (Ts({index2, index1}) + mL({index2, index1})*
                                       (Elements.MoleFractions(i,j,k)({index2, Comp}) - Cs({index2, index1})) - Tx(i,j,k))*
                                       (1 - Stoichiometric[index2])*(1 + Stoichiometric[index1]))*
                                       (S[index2] - S[index1]);

                    dGab.Raw(i,j,k).add_asym1(alpha->index,  beta->index,  dG_AB);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

double EquilibriumPartitionDiffusionBinary::GetDrivingForceAlpha(PhaseField& Phase,
                                                         Composition& Elements,
                                                         Temperature& Tx, size_t alpha,
                                                         int i, int j, int k)
{
    double dG = 0.0;
    CalculateLocalPhaseConcentrations(Phase, Tx, Elements, i, j, k);

    for(auto  beta = Phase.Fields(i,j,k).cbegin();
              beta < Phase.Fields(i,j,k).cend(); ++beta)
    {
        size_t index1 = alpha;
        size_t index2 = Phase.FieldsStatistics[beta->index].Phase;
        if(index1 != index2 and beta->value != 0.0)
        {
            double dG_AB = 0.5*((Ts({index1, index2}) + mL({index1, index2})*
                               (Elements.MoleFractions(i,j,k)({index1, Comp}) - Cs({index1, index2})) - Tx(i,j,k))*
                               (1 - Stoichiometric[index1])*(1 + Stoichiometric[index2])
                                +
                               (Ts({index2, index1}) + mL({index2, index1})*
                               (Elements.MoleFractions(i,j,k)({index2, Comp}) - Cs({index2, index1})) - Tx(i,j,k))*
                               (1 - Stoichiometric[index2])*(1 + Stoichiometric[index1]))*
                               (S[index2] - S[index1]);

            dG += dG_AB;
        }
    }
    return dG;
}

void EquilibriumPartitionDiffusionBinary::RestoreStoichiometric(PhaseField& Phase,
                                                          Temperature& Tx,
                                                          Composition& Elements)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,)
    {
        if(Phase.Fields(i,j,k).flag == 1)
        {
            size_t alpha = Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase;
            if (Stoichiometric[alpha] &&
                Elements.MoleFractionsTotal(i,j,k)({Comp}) != Elements.Initial({alpha, Comp}))
            {
                double mass = 0.0;
                double dC = Elements.MoleFractionsTotal(i,j,k)({Comp}) - Elements.Initial({alpha, Comp});
                for(int x = -dNx; x <= dNx; x++)
                for(int y = -dNy; y <= dNy; y++)
                for(int z = -dNz; z <= dNz; z++)
                if(Phase.Interface(i + x, j + y, k + z) &&
                   (int)i+x >= 0 && (int)i+x < Nx &&
                   (int)j+y >= 0 && (int)j+y < Ny &&
                   (int)k+z >= 0 && (int)k+z < Nz)
                {
                    for(size_t beta = 0; beta < Nphases; ++beta)
                    if(!Stoichiometric[beta])
                    {
                        mass += Phase.Fractions(i + x, j + y, k + z)({beta});
                    }
                }

                if(mass > 0.0)
                {
                    mass = 1.0/(mass);
                    Elements.MoleFractionsTotal(i, j, k)({Comp}) -= dC;

                    for(int x = -dNx; x <= dNx; x++)
                    for(int y = -dNy; y <= dNy; y++)
                    for(int z = -dNz; z <= dNz; z++)
                    if(Phase.Interface(i + x, j + y, k + z) &&
                       (int)i+x >= 0 && (int)i+x < Nx &&
                       (int)j+y >= 0 && (int)j+y < Ny &&
                       (int)k+z >= 0 && (int)k+z < Nz)
                    {
                        for(size_t beta = 0; beta < Nphases; ++beta)
                        if(!Stoichiometric[beta])
                        {
                            Elements.MoleFractionsTotal(i + x, j + y, k + z)({Comp}) += mass *
                                dC * Phase.Fractions(i + x, j + y, k + z)({beta});
                        }
                        CalculateLocalPhaseConcentrations(Phase, Tx, Elements,i,j,k);
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::RestoreStoichiometricThreadSafe(PhaseField& Phase,
                                                          Temperature& Tx,
                                                          Composition& Elements)
{
    int Bcells = Elements.MoleFractionsTotal.Bcells();
    #ifdef _OPENMP
    int locksize = (Nx+2*Bcells*dNx)*(Ny+2*Bcells*dNy)*(Nz+2*Bcells*dNz);
    std::vector<omp_lock_t> writelock(locksize);
    for (int i = 0; i < locksize; ++i)
    {
        omp_init_lock(&writelock[i]);
    }
    #endif
    Storage3D< double, 0 > dC;
    dC.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Elements.MoleFractionsTotal.Bcells());

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,)
    {
        if(Phase.Fields(i,j,k).flag == 1)
        {
            size_t alpha = Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase;
            if (Stoichiometric[alpha] &&
                Elements.MoleFractionsTotal(i,j,k)({Comp}) != Elements.Initial({alpha, Comp}))
            {
                dC(i,j,k) = Elements.MoleFractionsTotal(i,j,k)({Comp}) - Elements.Initial({alpha, Comp});
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,)
    {
        if(Phase.Fields(i,j,k).flag == 1)
        {
            size_t alpha = Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase;
            if (Stoichiometric[alpha] &&
                Elements.MoleFractionsTotal(i,j,k)({Comp}) != Elements.Initial({alpha, Comp}))
            {
                double mass = 0.0;
                for(int x = -dNx; x <= dNx; x++)
                for(int y = -dNy; y <= dNy; y++)
                for(int z = -dNz; z <= dNz; z++)
                if(Phase.Interface(i + x, j + y, k + z) &&
                   (int)i+x >= -Bcells*dNx && (int)i+x < Nx + Bcells*dNx &&
                   (int)j+y >= -Bcells*dNy && (int)j+y < Ny + Bcells*dNy &&
                   (int)k+z >= -Bcells*dNz && (int)k+z < Nz + Bcells*dNz)
                {
                    for(size_t beta = 0; beta < Nphases; ++beta)
                    if(!Stoichiometric[beta])
                    {
                        mass += Phase.Fractions(i + x, j + y, k + z)({beta});
                    }
                }
                if(mass > 0.0)
                {
                    mass = 1.0/(mass);
                    #ifdef _OPENMP
                    omp_set_lock(&writelock[((Ny + 2*Bcells*dNy)*(i + Bcells*dNx) + j + Bcells*dNy)*(Nz + 2*Bcells*dNz) + k + Bcells*dNz]);
                    #endif
                    Elements.MoleFractionsTotal(i, j, k)({Comp}) -= dC(i,j,k);
                    CalculateLocalPhaseConcentrations(Phase, Tx, Elements,i, j, k);
                    #ifdef _OPENMP
                    omp_unset_lock(&writelock[((Ny + 2*Bcells*dNy)*(i + Bcells*dNx) + j + Bcells*dNy)*(Nz + 2*Bcells*dNz) + k + Bcells*dNz]);
                    #endif
                    for(int x = -dNx; x <= dNx; x++)
                    for(int y = -dNy; y <= dNy; y++)
                    for(int z = -dNz; z <= dNz; z++)
                    if(Phase.Interface(i + x, j + y, k + z) &&
                        (int)i+x >= -Bcells*dNx && (int)i+x < Nx + Bcells*dNx &&
                        (int)j+y >= -Bcells*dNy && (int)j+y < Ny + Bcells*dNy &&
                        (int)k+z >= -Bcells*dNz && (int)k+z < Nz + Bcells*dNz)
                    {
                        #ifdef _OPENMP
                        omp_set_lock(&writelock[((Ny + 2*Bcells*dNy)*(i + Bcells*dNx) + j + Bcells*dNy)*(Nz + 2*Bcells*dNz) + k + Bcells*dNz]);
                        #endif
                        for(size_t beta = 0; beta < Nphases; ++beta)
                        if(!Stoichiometric[beta])
                        {
                            Elements.MoleFractionsTotal(i + x, j + y, k + z)({Comp}) += mass *
                                dC(i,j,k) * Phase.Fractions(i + x, j + y, k + z)({beta});
                        }
                        CalculateLocalPhaseConcentrations(Phase, Tx, Elements,i + x, j + y, k + z);
                        #ifdef _OPENMP
                        omp_unset_lock(&writelock[((Ny + 2*Bcells*dNy)*(i + Bcells*dNx) + j + Bcells*dNy)*(Nz + 2*Bcells*dNz) + k + Bcells*dNz]);
                        #endif
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)
{
    DC.Reallocate(newNx, newNy, newNz);
    dMu.Reallocate(newNx, newNy, newNz);
    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    Info::WriteStandard(thisclassname, "Remeshed!");
}

void EquilibriumPartitionDiffusionBinary::CalculateInterfaceMobility(PhaseField& Phase,
        Composition& Elements, Temperature& Tx, BoundaryConditions& BC, InterfaceProperties& IP,
        bool InternalDiffusivities)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,Elements.MoleFractionsTotal.Bcells(),)
    {
        if(Phase.Fields(i,j,k).flag)
        {
            if(Phase.Fields(i,j,k).flag)
            for(auto alpha  = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
            for(auto beta   = alpha + 1;
                     beta  != Phase.Fields(i,j,k).cend(); ++beta)
            {
                size_t index1 = Phase.FieldsStatistics[alpha->index].Phase;
                size_t index2 = Phase.FieldsStatistics[ beta->index].Phase;
                if(IP.InterfaceMobility(index1,index2).Model == InterfaceMobilityModels::Ext)
                {
                    double mob = 0.0;
                    if(!Stoichiometric[index1] or !Stoichiometric[index2])
                    {
                        double dS = (S[index1] - S[index2]);
                        double dC = (Cs({index1, index2}) + (Tx(i,j,k) - Ts({index1, index2}))*mL_1({index1, index2}) -
                                    (Cs({index2, index1}) + (Tx(i,j,k) - Ts({index2, index1}))*mL_1({index2, index1})));

                        double mLtmp = 0.0;
                        if(!Stoichiometric[index1])
                        {
                            mLtmp = min(fabs(mL({index1, index2})), fabs(mL({index2, index1})));
                        }
                        else
                        {
                            mLtmp = mL({index2, index1});
                        }

                        double denominator = mLtmp*Eta*dS*dC;

                        if(denominator != 0.0)
                        {
                            mob = fabs(8.0*(DC(i,j,k)({index1}) + DC(i,j,k)({index2}))/denominator);
                        }
                    }
                    IP.set_mobility(i, j, k, alpha->index, beta->index, mob);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

void EquilibriumPartitionDiffusionBinary::CalculateAntitrappingIncrements(PhaseField& Phase,
                                                                Composition& Elements)
{
    const double koef = Eta/Pi;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,)
    {
        if(Phase.Fields(i,j,k).flag)
        for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
        {
            int di = ds->di;
            int dj = ds->dj;
            int dk = ds->dk;
            if(Phase.Fields(i+di, j+dj, k+dk).flag)
            {
                dVector3 StencilDir{double(di),double(dj),double(dk)};
                StencilDir.normalize();

                NodePF locPhaseFields = Phase.Fields(i,j,k)*0.5;
                locPhaseFields.add_existing_values(Phase.Fields(i+di, j+dj, k+dk)*0.5);

                NodeV3 locNormals = (Phase.Normals(i,j,k) + Phase.Normals(i+di,j+dj,k+dk));
                if(locPhaseFields.size() > 1)
                for(auto alpha = locPhaseFields.cbegin();
                         alpha != locPhaseFields.cend(); ++alpha)
                {
                    size_t pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                    if(!Stoichiometric[pIndexA])
                    {
                        for(auto beta = locPhaseFields.cbegin();
                                 beta != locPhaseFields.cend(); ++beta)
                        if(alpha != beta)
                        {
                            size_t pIndexB = Phase.FieldsStatistics[beta->index].Phase;
                            if (pIndexA != pIndexB and (DC0[pIndexA] + DC0[pIndexB]) != 0.0)
                            {
                                double scale = 1.0;

                                if(Phase.FieldsStatistics[alpha->index].Stage == 1)
                                {
                                    scale *= Phase.FieldsStatistics[alpha->index].VolumeRatio;
                                }
                                if(Phase.FieldsStatistics[ beta->index].Stage == 1)
                                {
                                    scale *= Phase.FieldsStatistics[ beta->index].VolumeRatio;
                                }

                                dVector3 locNormal = locNormals.get_asym(alpha->index, beta->index).normalized();
                                double Projection = locNormal*StencilDir;

                                double locPsi = 0.5*(Phase.FieldsDot(i   , j   , k   ).get_asym1(alpha->index, beta->index) +
                                                     Phase.FieldsDot(i+di, j+dj, k+dk).get_asym1(alpha->index, beta->index));

                                double locDCalpha = 0.5*(DC(i,j,k)({pIndexA}) + DC(i+di, j+dj, k+dk)({pIndexA}));
                                double locDCbeta  = 0.5*(DC(i,j,k)({pIndexB}) + DC(i+di, j+dj, k+dk)({pIndexB}));

                                double localDelta = -dx*ds->weight*
                                                    sqrt(alpha->value*beta->value)*0.5*
                                                    ((Elements.MoleFractions(i   , j   , k   )({pIndexA, Comp}) +
                                                      Elements.MoleFractions(i+di, j+dj, k+dk)({pIndexA, Comp})) -
                                                     (Elements.MoleFractions(i   , j   , k   )({pIndexB, Comp}) +
                                                      Elements.MoleFractions(i+di, j+dj, k+dk)({pIndexB, Comp}))) *
                                                    locPsi * Projection * koef * scale *
                                                        (locDCalpha - locDCbeta)/
                                                        (locDCalpha + locDCbeta);

                                if (localDelta < -Precision)
                                {
                                    Elements.MoleFractionsDot(i, j, k)({pIndexA, Comp}).out += localDelta;
                                }
                                if (localDelta > Precision)
                                {
                                    Elements.MoleFractionsDot(i, j, k)({pIndexA, Comp}).in += localDelta;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::LimitAntitrappingIncrements(PhaseField& Phase,
                                                          Composition& Elements,
                                                          Temperature& Tx)
{
    const double koef = Eta/Pi;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,)
    {
        if(Elements.Limiting(i, j, k))
        {
            if(Phase.Fields(i,j,k).flag)
            {
                for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
                {
                    int di = ds->di;
                    int dj = ds->dj;
                    int dk = ds->dk;
                    if(Phase.Fields(i+di, j+dj, k+dk).flag)
                    {
                        dVector3 StencilDir{double(di),double(dj),double(dk)};
                        StencilDir.normalize();

                        NodePF locPhaseFields = Phase.Fields(i,j,k)*0.5;
                        locPhaseFields.add_existing_values(Phase.Fields(i+di, j+dj, k+dk)*0.5);

                        NodeV3 locNormals = (Phase.Normals(i,j,k) + Phase.Normals(i+di,j+dj,k+dk));

                        if(locPhaseFields.size() > 1)
                        for(auto alpha = locPhaseFields.cbegin();
                                 alpha != locPhaseFields.cend(); ++alpha)
                        {
                            size_t pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                            if(!Stoichiometric[pIndexA])
                            for(auto beta = locPhaseFields.cbegin();
                                     beta < locPhaseFields.cend(); ++beta)
                            if(alpha != beta)
                            {
                                size_t pIndexB = Phase.FieldsStatistics[beta->index].Phase;
                                if (pIndexA != pIndexB and (DC0[pIndexA] + DC0[pIndexB]) != 0.0)

                                {
                                    double scale = 1.0;

                                    if(Phase.FieldsStatistics[alpha->index].Stage == 1)
                                    {
                                        scale *= (Phase.FieldsStatistics[alpha->index].Volume/Phase.RefVolume);
                                    }
                                    if(Phase.FieldsStatistics[ beta->index].Stage == 1)
                                    {
                                        scale *= (Phase.FieldsStatistics[beta->index].Volume/Phase.RefVolume);
                                    }

                                    dVector3 locNormal = locNormals.get_asym(alpha->index, beta->index).normalized();
                                    double Projection = locNormal*StencilDir;

                                    double locPsi = 0.5*(Phase.FieldsDot(i   , j   , k   ).get_asym1(alpha->index, beta->index) +
                                                         Phase.FieldsDot(i+di, j+dj, k+dk).get_asym1(alpha->index, beta->index));

                                    double locDCalpha = 0.5*(DC(i,j,k)({pIndexA}) + DC(i+di, j+dj, k+dk)({pIndexA}));
                                    double locDCbeta  = 0.5*(DC(i,j,k)({pIndexB}) + DC(i+di, j+dj, k+dk)({pIndexB}));

                                    double localDelta = -dx*ds->weight*
                                                        sqrt(alpha->value*beta->value)*0.5*
                                                        ((Elements.MoleFractions(i   , j   , k   )({pIndexA, Comp}) +
                                                          Elements.MoleFractions(i+di, j+dj, k+dk)({pIndexA, Comp})) -
                                                         (Elements.MoleFractions(i   , j   , k   )({pIndexB, Comp}) +
                                                          Elements.MoleFractions(i+di, j+dj, k+dk)({pIndexB, Comp}))) *
                                                        locPsi * Projection * koef * scale *
                                                            (locDCalpha - locDCbeta)/
                                                            (locDCalpha + locDCbeta);

                                    if (Elements.Norm(i, j, k)({pIndexA, Comp}).in <=
                                        Elements.Norm(i+di, j+dj, k+dk)({pIndexA, Comp}).out and
                                        localDelta > Precision)
                                    {
                                        double localNorm = Elements.Norm(i, j, k)({pIndexA, Comp}).in;
                                        double localIncrement = localDelta*(1.0 - localNorm);
                                        if(localIncrement > Precision)
                                        {
                                            Elements.MoleFractionsDot(   i,    j,    k)({pIndexA, Comp}).in  -= localIncrement;
                                            Elements.MoleFractionsDot(i+di, j+dj, k+dk)({pIndexA, Comp}).out += localIncrement;
                                        }
                                    }
                                    if (Elements.Norm(i, j, k)({pIndexA, Comp}).out <
                                        Elements.Norm(i+di, j+dj, k+dk)({pIndexA, Comp}).in and
                                        localDelta < -Precision)
                                    {
                                        double localNorm = Elements.Norm(i, j, k)({pIndexA, Comp}).out;
                                        double localIncrement = localDelta*(1.0 - localNorm);
                                        if(localIncrement < -Precision)
                                        {
                                            Elements.MoleFractionsDot(   i,    j,    k)({pIndexA, Comp}).out -= localIncrement;
                                            Elements.MoleFractionsDot(i+di, j+dj, k+dk)({pIndexA, Comp}).in  += localIncrement;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::LimitAntitrappingIncrementsThreadSafe(PhaseField& Phase,
                                                          Composition& Elements,
                                                          Temperature& Tx)
{
    #ifdef _OPENMP
    int Bcells = Elements.MoleFractionsTotal.Bcells();
    int locksize = (Nx+2*Bcells*dNx)*(Ny+2*Bcells*dNy)*(Nz+2*Bcells*dNz);
    std::vector<omp_lock_t> writelock(locksize);
    for (int i = 0; i < locksize; ++i)
    {
        omp_init_lock(&writelock[i]);
    }
    #endif

    const double koef = Eta/Pi;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,)
    {
        if(Elements.Limiting(i, j, k))
        {
            if(Phase.Fields(i,j,k).flag)
            {
                for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
                {
                    int di = ds->di;
                    int dj = ds->dj;
                    int dk = ds->dk;
                    if(Phase.Fields(i+di, j+dj, k+dk).flag)
                    {
                        dVector3 StencilDir{double(di),double(dj),double(dk)};
                        StencilDir.normalize();

                        NodePF locPhaseFields = Phase.Fields(i,j,k)*0.5;
                        locPhaseFields.add_existing_values(Phase.Fields(i+di, j+dj, k+dk)*0.5);

                        NodeV3 locNormals = (Phase.Normals(i,j,k) + Phase.Normals(i+di, j+dj, k+dk));
                        if(locPhaseFields.size() > 1)
                        for(auto alpha = locPhaseFields.cbegin();
                                 alpha != locPhaseFields.cend(); ++alpha)
                        {
                            size_t pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                            if(!Stoichiometric[pIndexA])
                            for(auto beta = locPhaseFields.cbegin();
                                     beta < locPhaseFields.cend(); ++beta)
                            if(alpha != beta)
                            {
                                size_t pIndexB = Phase.FieldsStatistics[beta->index].Phase;
                                if (pIndexA != pIndexB and (DC0[pIndexA] + DC0[pIndexB]) != 0.0)
                                {
                                    double scale = 1.0;

                                    if(Phase.FieldsStatistics[alpha->index].Stage == 1)
                                    {
                                        scale *= (Phase.FieldsStatistics[alpha->index].Volume/Phase.RefVolume);
                                    }
                                    if(Phase.FieldsStatistics[ beta->index].Stage == 1)
                                    {
                                        scale *= (Phase.FieldsStatistics[beta->index].Volume/Phase.RefVolume);
                                    }

                                    dVector3 locNormal = locNormals.get_asym(alpha->index, beta->index).normalized();
                                    double Projection = locNormal*StencilDir;

                                    double locPsi = 0.5*(Phase.FieldsDot(i   , j   , k   ).get_asym1(alpha->index, beta->index) +
                                                         Phase.FieldsDot(i+di, j+dj, k+dk).get_asym1(alpha->index, beta->index));

                                    double locDCalpha = 0.5*(DC(i,j,k)({pIndexA}) + DC(i+di, j+dj, k+dk)({pIndexA}));
                                    double locDCbeta  = 0.5*(DC(i,j,k)({pIndexB}) + DC(i+di, j+dj, k+dk)({pIndexB}));

                                    double localDelta = -dx*ds->weight*
                                                        sqrt(alpha->value*beta->value)*0.5*
                                                        ((Elements.MoleFractions(i   , j   , k   )({pIndexA, Comp}) +
                                                          Elements.MoleFractions(i+di, j+dj, k+dk)({pIndexA, Comp})) -
                                                         (Elements.MoleFractions(i   , j   , k   )({pIndexB, Comp}) +
                                                          Elements.MoleFractions(i+di, j+dj, k+dk)({pIndexB, Comp}))) *
                                                        locPsi * Projection * koef * scale *
                                                            (locDCalpha - locDCbeta)/
                                                            (locDCalpha + locDCbeta);

                                    if (Elements.Norm(i, j, k)({pIndexA, Comp}).in <=
                                        Elements.Norm(i+di, j+dj, k+dk)({pIndexA, Comp}).out and
                                        localDelta > Precision)
                                    {
                                        double localNorm = Elements.Norm(i, j, k)({pIndexA, Comp}).in;
                                        double localIncrement = localDelta*(1.0 - localNorm);
                                        if(localIncrement > Precision)
                                        {
                                            Elements.MoleFractionsDot(  i,   j,   k)({pIndexA, Comp}).in -= localIncrement;
                                            #ifdef _OPENMP
                                            omp_set_lock(&writelock[((Ny + 2*Bcells*dNy)*(i + Bcells*dNx) + j + Bcells*dNy)*(Nz + 2*Bcells*dNz) + k + Bcells*dNz]);
                                            #endif
                                            Elements.MoleFractionsDot(i+di, j+dj, k+dk)({pIndexA, Comp}).out += localIncrement;
                                            #ifdef _OPENMP
                                            omp_unset_lock(&writelock[((Ny + 2*Bcells*dNy)*(i + Bcells*dNx) + j + Bcells*dNy)*(Nz + 2*Bcells*dNz) + k + Bcells*dNz]);
                                            #endif
                                        }
                                    }
                                    if (Elements.Norm(i, j, k)({pIndexA, Comp}).out <
                                        Elements.Norm(i+di, j+dj, k+dk)({pIndexA, Comp}).in and
                                        localDelta < -Precision)
                                    {
                                        double localNorm = Elements.Norm(i, j, k)({pIndexA, Comp}).out;
                                        double localIncrement = localDelta*(1.0 - localNorm);
                                        if(localIncrement < -Precision)
                                        {
                                            Elements.MoleFractionsDot(  i,   j,  k)({pIndexA, Comp}).out -= localIncrement;
                                            #ifdef _OPENMP
                                            omp_set_lock(&writelock[((Ny + 2*Bcells*dNy)*(i + Bcells*dNx) + j + Bcells*dNy)*(Nz + 2*Bcells*dNz) + k + Bcells*dNz]);
                                            #endif
                                            Elements.MoleFractionsDot(i+di, j+dj, k+dk)({pIndexA, Comp}).in += localIncrement;
                                            #ifdef _OPENMP
                                            omp_unset_lock(&writelock[((Ny + 2*Bcells*dNy)*(i + Bcells*dNx) + j + Bcells*dNy)*(Nz + 2*Bcells*dNz) + k + Bcells*dNz]);
                                            #endif
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::CalculateDiffusionIncrements(PhaseField& Phase,
                                                          Composition& Elements,
                                                          Temperature& Tx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,)
    {
        for(size_t alpha = 0; alpha != Nphases; ++alpha)
        {
            Elements.MoleFractionsDot(i, j, k)({alpha, Comp}).in = 0.0;
            Elements.MoleFractionsDot(i, j, k)({alpha, Comp}).out = 0.0;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,)
    {
        if(Phase.Fields(i,j,k).flag)
        {
            for(size_t alpha = 0; alpha != Nphases; ++alpha)
            if(!Stoichiometric[alpha] and Phase.Fractions(i, j, k)({alpha}) > 0.0)
            {
                double invRT  = Elements.TotalMolarVolume/(R*Tx(i,j,k));
                for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
                {
                    int di = ds->di;
                    int dj = ds->dj;
                    int dk = ds->dk;
                    if (Phase.Fractions(i+di, j+dj, k+dk)({alpha}) > 0.0)
                    {
                        double invRTxyz  = Elements.TotalMolarVolume/(R*Tx(i+di, j+dj, k+dk));
                        double localDelta = -ds->weight*
                                sqrt(Phase.Fractions(i, j, k)({alpha}) *
                                     Phase.Fractions(i+di, j+dj, k+dk)({alpha}))*

                               0.5*((DC(i,j,k)({alpha}) + DC(i+di, j+dj, k+dk)({alpha}))*    // D * Grad ( C {alpha})
                                    (Elements.MoleFractions(   i,    j,    k)({alpha, Comp}) -
                                     Elements.MoleFractions(i+di, j+dj, k+dk)({alpha, Comp})) +

                                    (invRT*DC(i,j,k)({alpha}) +
                                     invRTxyz*DC(i+di, j+dj, k+dk)({alpha}))*
                                    (dMu(   i,    j,    k)({alpha}) -
                                     dMu(i+di, j+dj, k+dk)({alpha})));          // dMu = 0;

                        if (localDelta < -Precision)
                        {
                            Elements.MoleFractionsDot(i, j, k)({alpha, Comp}).out += localDelta;
                        }
                        if (localDelta > Precision)
                        {
                            Elements.MoleFractionsDot(i, j, k)({alpha, Comp}).in += localDelta;
                        }
                    }
                }
            }
        }
        else
        {
            size_t alpha = Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase;
            if(!Stoichiometric[alpha])
            {
                double invRT  = Elements.TotalMolarVolume/(R*Tx(i,j,k));
                for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
                {
                    int di = ds->di;
                    int dj = ds->dj;
                    int dk = ds->dk;
                    double invRTxyz  = Elements.TotalMolarVolume/(R*Tx(i+di, j+dj, k+dk));
                    double localDelta = - ds->weight*
                           0.5*((DC(i,j,k)({alpha}) + DC(i+di, j+dj, k+dk)({alpha}))*
                                (Elements.MoleFractions(  i,   j,   k)({alpha, Comp}) -
                                 Elements.MoleFractions(i+di, j+dj, k+dk)({alpha, Comp})) +

                                (invRT*DC(i,j,k)({alpha}) +
                                 invRTxyz*DC(i+di, j+dj, k+dk)({alpha}))*
                                (dMu(   i,    j,    k)({alpha}) -
                                 dMu(i+di, j+dj, k+dk)({alpha})));

                    if (localDelta < -Precision)
                    {
                        Elements.MoleFractionsDot(i, j, k)({alpha, Comp}).out += localDelta;
                    }
                    if (localDelta > Precision)
                    {
                        Elements.MoleFractionsDot(i, j, k)({alpha, Comp}).in += localDelta;
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::CalculateLimits(PhaseField& Phase,
                                              Composition& Elements, double dt)
{
    bool ompLimitingNeeded = false;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Norm,Elements.Norm.Bcells(),)
    {
        for(size_t alpha = 0; alpha != Nphases; ++alpha)
        {
            Elements.Norm(i,j,k)({alpha, Comp}).in = 1.0;
            Elements.Norm(i,j,k)({alpha, Comp}).out = 1.0;
        }
        Elements.Limiting(i,j,k) = 0;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,reduction(||:ompLimitingNeeded))
    {
        if(Phase.Fields(i,j,k).flag)
        {
            for(size_t alpha = 0; alpha != Nphases; ++alpha)
            if(!Stoichiometric[alpha] and Phase.Fractions(i,j,k)({alpha}) > 0.0)
            {
                double dPhaseIn  = Elements.MoleFractionsDot(i,j,k)({alpha, Comp}).in*dt;
                double dPhaseOut = Elements.MoleFractionsDot(i,j,k)({alpha, Comp}).out*dt;

                if(dPhaseIn <= Precision) dPhaseIn = 0.0;
                if(dPhaseOut >= -Precision) dPhaseOut = 0.0;

                double PhaseOld = Elements.MoleFractions(i,j,k)({alpha, Comp});

                if(PhaseOld + dPhaseIn > Elements.MoleFractionsMAX({alpha, Comp}) && dPhaseIn != 0.0)
                {
                    Elements.Norm(i,j,k)({alpha, Comp}).in =
                                          (Elements.MoleFractionsMAX({alpha, Comp}) - PhaseOld)/dPhaseIn;
                    ompLimitingNeeded = true;

                    for (int x = -dNx; x <= dNx; ++x)
                    for (int y = -dNy; y <= dNy; ++y)
                    for (int z = -dNz; z <= dNz; ++z)
                    {
                        Elements.Limiting(i+x, j+y, k+z) = 1;
                    }
                }
                if(PhaseOld + dPhaseOut < Elements.MoleFractionsMIN({alpha, Comp}) && dPhaseOut != 0.0)
                {
                    Elements.Norm(i,j,k)({alpha, Comp}).out =
                                          (Elements.MoleFractionsMIN({alpha, Comp}) - PhaseOld)/dPhaseOut;
                    ompLimitingNeeded = true;
                    for (int x = -dNx; x <= dNx; ++x)
                    for (int y = -dNy; y <= dNy; ++y)
                    for (int z = -dNz; z <= dNz; ++z)
                    {
                        Elements.Limiting(i+x, j+y, k+z) = 1;
                    }
                }
            }
        }
        else
        {
            size_t alpha = Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase;
            if(!Stoichiometric[alpha])
            {
                double dPhaseIn  = Elements.MoleFractionsDot(i,j,k)({alpha, Comp}).in*dt;
                double dPhaseOut = Elements.MoleFractionsDot(i,j,k)({alpha, Comp}).out*dt;

                if(dPhaseIn <= Precision) dPhaseIn = 0.0;
                if(dPhaseOut >= -Precision) dPhaseOut = 0.0;

                double PhaseOld = Elements.MoleFractions(i,j,k)({alpha, Comp});

                if(PhaseOld + dPhaseIn > Elements.MoleFractionsMAX({alpha, Comp}) && dPhaseIn != 0.0)
                {
                    Elements.Norm(i,j,k)({alpha, Comp}).in =
                                          (Elements.MoleFractionsMAX({alpha, Comp}) - PhaseOld)/dPhaseIn;
                    ompLimitingNeeded = true;
                    for (int x = -dNx; x <= dNx; ++x)
                    for (int y = -dNy; y <= dNy; ++y)
                    for (int z = -dNz; z <= dNz; ++z)
                    {
                        Elements.Limiting(i+x, j+y, k+z) = 1;
                    }
                }
                if(PhaseOld + dPhaseOut < Elements.MoleFractionsMIN({alpha, Comp}) && dPhaseOut != 0.0)
                {
                    Elements.Norm(i,j,k)({alpha, Comp}).out =
                                          (Elements.MoleFractionsMIN({alpha, Comp}) - PhaseOld)/dPhaseOut;
                    ompLimitingNeeded = true;
                    for (int x = -dNx; x <= dNx; ++x)
                    for (int y = -dNy; y <= dNy; ++y)
                    for (int z = -dNz; z <= dNz; ++z)
                    {
                        Elements.Limiting(i+x, j+y, k+z) = 1;
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    LimitingNeeded = ompLimitingNeeded;
#ifdef MPI_PARALLEL
    bool rLimitingNeeded = LimitingNeeded;
    MPI_Allreduce(&rLimitingNeeded, &LimitingNeeded, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
#endif
}

void EquilibriumPartitionDiffusionBinary::LimitDiffusionIncrements(PhaseField& Phase,
                                                          Composition& Elements,
                                                          Temperature& Tx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,)
    {
        if(Elements.Limiting(i, j, k))
        {
            if(Phase.Fields(i,j,k).flag)
            {
                for(size_t alpha = 0; alpha != Nphases; ++alpha)
                if(!Stoichiometric[alpha] and Phase.Fractions(i, j, k)({alpha}) > 0.0)
                {
                    double invRT = Elements.TotalMolarVolume/(R*Tx(i, j, k));
                    for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
                    {
                        int di = ds->di;
                        int dj = ds->dj;
                        int dk = ds->dk;
                        if (Phase.Fractions(i+di, j+dj, k+dk)({alpha}) > 0.0)
                        {
                            double invRTxyz  = Elements.TotalMolarVolume/(R*Tx(i+di, j+dj, k+dk));
                            double localDelta = -ds->weight*
                                    sqrt(Phase.Fractions(i, j, k)({alpha}) *
                                         Phase.Fractions(i+di, j+dj, k+dk)({alpha}))*

                                   0.5*((DC(i,j,k)({alpha}) + DC(i+di, j+dj, k+dk)({alpha}))*
                                        (Elements.MoleFractions(  i,   j,   k)({alpha, Comp}) -
                                         Elements.MoleFractions(i+di, j+dj, k+dk)({alpha, Comp})) +

                                        (invRT*DC(i,j,k)({alpha}) +
                                         invRTxyz*DC(i+di, j+dj, k+dk)({alpha}))*
                                        (dMu(  i,   j,   k)({alpha}) -
                                         dMu(i+di, j+dj, k+dk)({alpha})));

                            if (Elements.Norm(i, j, k)({alpha, Comp}).in <=
                                Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).out)
                            {
                                double localNorm = Elements.Norm(i, j, k)({alpha, Comp}).in;
                                double localIncrement = localDelta*(1.0 - localNorm);
                                if(localIncrement > Precision)
                                {
                                    Elements.MoleFractionsDot(  i,   j,   k)({alpha, Comp}).in -=
                                                                      localIncrement;
                                }
                            }
                            if (Elements.Norm(i, j, k)({alpha, Comp}).in >
                                Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).out)
                            {
                                double localNorm = Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).out;
                                double localIncrement = localDelta*(1.0 - localNorm);
                                if(localIncrement > Precision)
                                {
                                    Elements.MoleFractionsDot(  i,   j,   k)({alpha, Comp}).in -=
                                                                      localIncrement;
                                }
                            }

                            if (Elements.Norm(i, j, k)({alpha, Comp}).out <
                                Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).in)
                            {
                                double localNorm = Elements.Norm(i, j, k)({alpha, Comp}).out;
                                double localIncrement = localDelta*(1.0 - localNorm);
                                if(localIncrement < -Precision)
                                {
                                    Elements.MoleFractionsDot(  i,   j,   k)({alpha, Comp}).out -=
                                                                      localIncrement;
                                }
                            }
                            if (Elements.Norm(i, j, k)({alpha, Comp}).out >=
                                Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).in)
                            {
                                double localNorm = Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).in;
                                double localIncrement = localDelta*(1.0 - localNorm);
                                if(localIncrement < -Precision)
                                {
                                    Elements.MoleFractionsDot(  i,   j,   k)({alpha, Comp}).out -=
                                                                      localIncrement;
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                size_t alpha = Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase;
                if(!Stoichiometric[alpha])
                {
                    double invRT  = Elements.TotalMolarVolume/(R*Tx(i, j, k));
                    for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
                    {
                        int di = ds->di;
                        int dj = ds->dj;
                        int dk = ds->dk;
                        double invRTxyz  = Elements.TotalMolarVolume/(R*Tx(i+di, j+dj, k+dk));
                        double localDelta = - ds->weight*

                               0.5*((DC(i,j,k)({alpha}) + DC(i+di, j+dj, k+dk)({alpha}))*
                                    (Elements.MoleFractions(  i,   j,   k)({alpha, Comp}) -
                                     Elements.MoleFractions(i+di, j+dj, k+dk)({alpha, Comp})) +

                                    (invRT*DC(i,j,k)({alpha}) +
                                     invRTxyz*DC(i+di, j+dj, k+dk)({alpha}))*
                                    (dMu(  i,   j,   k)({alpha}) -
                                     dMu(i+di, j+dj, k+dk)({alpha})));

                        if (Elements.Norm(i, j, k)({alpha, Comp}).in <=
                            Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).out)
                        {
                            double localNorm = Elements.Norm(i, j, k)({alpha, Comp}).in;
                            double localIncrement = localDelta*(1.0 - localNorm);
                            if(localIncrement > Precision)
                            {
                                Elements.MoleFractionsDot(  i,   j,   k)({alpha, Comp}).in -=
                                                                  localIncrement;
                            }
                        }
                        if (Elements.Norm(i, j, k)({alpha, Comp}).in >
                            Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).out)
                        {
                            double localNorm = Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).out;
                            double localIncrement = localDelta*(1.0 - localNorm);
                            if(localIncrement > Precision)
                            {
                                Elements.MoleFractionsDot(  i,   j,   k)({alpha, Comp}).in -=
                                                                  localIncrement;
                            }
                        }

                        if (Elements.Norm(i, j, k)({alpha, Comp}).out <
                            Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).in)
                        {
                            double localNorm = Elements.Norm(i, j, k)({alpha, Comp}).out;
                            double localIncrement = localDelta*(1.0 - localNorm);
                            if(localIncrement < -Precision)
                            {
                                Elements.MoleFractionsDot(  i,   j,   k)({alpha, Comp}).out -=
                                                                  localIncrement;
                            }
                        }
                        if (Elements.Norm(i, j, k)({alpha, Comp}).out >=
                            Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).in)
                        {
                            double localNorm = Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).in;
                            double localIncrement = localDelta*(1.0 - localNorm);
                            if(localIncrement < -Precision)
                            {
                                Elements.MoleFractionsDot(  i,   j,   k)({alpha, Comp}).out -=
                                                                  localIncrement;
                            }
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::ApplyIncrements(PhaseField& Phase,
                                               Composition& Elements, double dt)
{
    double ompTotalMass = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,0,reduction(+:ompTotalMass))
    {
        for(size_t alpha = 0; alpha < Nphases; ++alpha)
        if(!Stoichiometric[alpha])
        {
            double dPhase = (Elements.MoleFractionsDot(i,j,k)({alpha, Comp}).in +
                             Elements.MoleFractionsDot(i,j,k)({alpha, Comp}).out)*dt;
            if(fabs(dPhase) <= Precision) dPhase = 0.0;
            Elements.MoleFractionsTotal(i,j,k)({Comp}) += dPhase;
        }
        ompTotalMass += Elements.MoleFractionsTotal(i,j,k)({Comp});
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    TotalMass = ompTotalMass;
}

void EquilibriumPartitionDiffusionBinary::CalculateReferenceElementMoleFractions(Composition& Elements)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.MoleFractionsTotal,Elements.MoleFractionsTotal.Bcells(),)
    {
        for(size_t alpha = 0; alpha < Nphases; ++alpha)
        {
            Elements.MoleFractions(i,j,k)({alpha,RefComp}) = 1.0 - Elements.MoleFractions(i,j,k)({alpha,Comp});
        }
        Elements.MoleFractionsTotal(i,j,k)({RefComp}) = 1.0 - Elements.MoleFractionsTotal(i,j,k)({Comp});
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::Solve(PhaseField& Phase,
                    Composition& Elements, Temperature& Tx,
                    BoundaryConditions& BC, double dt,
                    bool InternalDiffusivities)
{
    Elements.SetBoundaryConditions(BC);

    CalculatePhaseConcentrations(Phase, Tx, Elements);

    if(InternalDiffusivities)
    {
        SetDiffusionCoefficients(Phase, Tx);
    }
    CalculateDiffusionIncrements(Phase, Elements, Tx);
    if(EnableAntiTrapping)
    {
        CalculateAntitrappingIncrements(Phase, Elements);
    }
    CalculateLimits(Phase, Elements, dt);

    if(LimitingNeeded)
    {
        Elements.SetLimitsBoundaryConditions(BC);
        LimitDiffusionIncrements(Phase, Elements, Tx);
        if(EnableAntiTrapping)
        {
            LimitAntitrappingIncrements(Phase, Elements, Tx);
        }
    }
    ApplyIncrements(Phase, Elements, dt);

    Elements.SetBoundaryConditions(BC);

    CalculatePhaseConcentrations(Phase, Tx, Elements);
    if(ThereAreStoichiometricPhases)
    {
        RestoreStoichiometric(Phase, Tx, Elements);
        Elements.SetBoundaryConditions(BC);
    }
    CalculateReferenceElementMoleFractions(Elements);
    Elements.CollectStatistics(Phase);
}

void EquilibriumPartitionDiffusionBinary::PrintPointStatistics(int x, int y, int z)
{
    cout << "Diffusion Coefficients:" << endl;
    for (size_t n = 0; n < Nphases; n++)
    {
        cout << "D[" << n << "] = " << DC(x, y, z)({n}) << endl;
    }
    cout << endl;
}

}// namespace openphase
