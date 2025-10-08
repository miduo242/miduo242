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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Philipp Engels;
 *                         Raphael Schiedung
 *
 */

#include "Base/UserInterface.h"
#include "Composition.h"
#include "DrivingForce.h"
#include "EquilibriumPartitionDiffusionBinary.h"
#include "Info.h"
#include "Mechanics/ElasticProperties.h"
#include "Mechanics/ElasticityModels/ElasticityKhachaturyan.h"
#include "Orientations.h"
#include "ParabolicDiffusion.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Temperature.h"
#include "VTK.h"

namespace openphase
{

using namespace std;

vStrain ElasticityKhachaturyan::CalculateLocEigenStrainDifference(
        const ElasticProperties& EP,
        const dMatrix3x3 locStretchesAlpha,
        const dMatrix3x3 locStretchesBeta,
        const int i, const int j, const int k)
{
    vStrain locEigenStrainDifference;
    locEigenStrainDifference = EP.StrainSmall(locStretchesBeta)
                             - EP.StrainSmall(locStretchesAlpha);
    return locEigenStrainDifference;
}

dMatrix6x6 ElasticityKhachaturyan::CalculateLocElasticConstants(
            const PhaseField& Phase,
            const ElasticProperties& EP, const Composition& Cx,
            const int i, const int j, const int k,
            const size_t idx)
{
    const size_t pIndex = Phase.FieldsStatistics[idx].Phase;
    dMatrix6x6 locElasticConstants = EP.ElasticConstants[idx];
    for(int ii = 0; ii < 6; ii++)
    for(int jj = 0; jj < 6; jj++)
    {
        double delta = 1.0;
        for(size_t comp = 0; comp < Cx.Ncomp; comp ++)
        {
            delta += EP.Kappa({idx, comp})(ii,jj)*
                (Cx.MoleFractions(i,j,k)({pIndex, comp}) - EP.Cref({pIndex, comp}));
        }
        locElasticConstants(ii,jj) *= delta;
    }
    return locElasticConstants;
}

void ElasticityKhachaturyan::CalculateDrivingForce(const PhaseField& Phase,
        const ElasticProperties& EP, DrivingForce& dGab)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.DeformationGradientsTotal,0,)
    if(Phase.Interface(i,j,k))
    {
        dMatrix6x6 locCompliance = EP.EffectiveElasticConstants(i,j,k).inverted();
        vStrain ElasticStrains = locCompliance*EP.Stresses(i, j, k);

        for(auto alpha = Phase.Fields(i, j, k).cbegin();
                 alpha != Phase.Fields(i, j, k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1;
                  beta != Phase.Fields(i, j, k).cend();  ++beta)
        {
            size_t pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
            size_t pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

            vStrain ElasticStrainsCorr = EP.CalculateNeuberCorrection(ElasticStrains,locCompliance, i, j, k, pIndexA, pIndexB);

            dMatrix3x3 locStretchesAlpha = EP.TransformationStretches[alpha->index];
            dMatrix3x3 locStretchesBeta  = EP.TransformationStretches[ beta->index];

            vStrain locEigenStrainDifference =
                CalculateLocEigenStrainDifference(EP, locStretchesAlpha, locStretchesBeta, i,j,k);

            double dG_AB = 0.5*(ElasticStrainsCorr*
                               ((EP.ElasticConstants[beta->index] -
                                 EP.ElasticConstants[alpha->index])*
                                ElasticStrainsCorr))

                    - locEigenStrainDifference*(EP.EffectiveElasticConstants(i,j,k)*ElasticStrainsCorr);

            dGab.Raw(i,j,k).add_asym1(alpha->index, beta->index, dG_AB);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticityKhachaturyan::CalculateDrivingForce(const PhaseField& Phase,
        const ElasticProperties& EP, const Composition& Cx, DrivingForce& dGab)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.DeformationGradientsTotal,0,)
    if(Phase.Interface(i,j,k))
    {
        dMatrix6x6 locCompliance = EP.EffectiveElasticConstants(i,j,k).inverted();
        vStrain ElasticStrains = locCompliance*EP.Stresses(i, j, k);

        for(auto alpha = Phase.Fields(i, j, k).cbegin();
                 alpha != Phase.Fields(i, j, k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1;
                  beta != Phase.Fields(i, j, k).cend();  ++beta)
        {
            size_t pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
            size_t pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

            vStrain ElasticStrainsCorr = EP.CalculateNeuberCorrection(ElasticStrains,locCompliance, i, j, k, pIndexA, pIndexB);

            dMatrix3x3 locStretchesAlpha = EP.TransformationStretches[alpha->index];
            dMatrix3x3 locStretchesBeta  = EP.TransformationStretches[ beta->index];

            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
            {
                locStretchesAlpha += EP.Lambda({alpha->index, comp})*(Cx.MoleFractions(i,j,k)({pIndexA, comp}) - EP.Cref({pIndexA, comp}));
                locStretchesBeta  += EP.Lambda({ beta->index, comp})*(Cx.MoleFractions(i,j,k)({pIndexB, comp}) - EP.Cref({pIndexB, comp}));
            }

            vStrain locEigenStrainDifference =
                    CalculateLocEigenStrainDifference(EP, locStretchesAlpha, locStretchesBeta, i,j,k);

            const dMatrix6x6 locElasticConstantsAlpha =
                CalculateLocElasticConstants(Phase,EP,Cx,i,j,k,alpha->index);
            const dMatrix6x6 locElasticConstantsBeta  =
                CalculateLocElasticConstants(Phase,EP,Cx,i,j,k,beta->index );

            const double dG_AB = (ElasticStrainsCorr*
                    ((locElasticConstantsAlpha - locElasticConstantsBeta)*
                            ElasticStrainsCorr))*0.5
                    - locEigenStrainDifference*(EP.EffectiveElasticConstants(i,j,k)*ElasticStrainsCorr);

            dGab.Raw(i,j,k).add_asym1(alpha->index, beta->index, dG_AB);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticityKhachaturyan::CalculateDrivingForce(const PhaseField& Phase,
        const ElasticProperties& EP, const Temperature& Tx, DrivingForce& dGab)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.DeformationGradientsTotal,0,)
    if(Phase.Interface(i,j,k))
    {
        dMatrix6x6 locCompliance = EP.EffectiveElasticConstants(i,j,k).inverted();
        vStrain ElasticStrains = locCompliance*EP.Stresses(i, j, k);

        for(auto alpha = Phase.Fields(i, j, k).cbegin();
                 alpha != Phase.Fields(i, j, k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1;
                  beta != Phase.Fields(i, j, k).cend();  ++beta)
        {
            size_t pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
            size_t pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

            double dG_AB = 0.0;

            vStrain ElasticStrainsCorr = EP.CalculateNeuberCorrection(ElasticStrains,locCompliance, i, j, k, pIndexA, pIndexB);

            dMatrix3x3 locStretchesAlpha = EP.TransformationStretches[alpha->index];
            dMatrix3x3 locStretchesBeta  = EP.TransformationStretches[ beta->index];
            for(int ii = 0; ii < 3; ii++)
            for(int jj = 0; jj < 3; jj++)
            {
                locStretchesAlpha(ii,jj) *= (1.0 + EP.Alpha[alpha->index](ii,jj)*(Tx(i,j,k) - EP.Tref[pIndexA]));
                locStretchesBeta(ii,jj)  *= (1.0 + EP.Alpha[ beta->index](ii,jj)*(Tx(i,j,k) - EP.Tref[pIndexB]));
            }
            vStrain locEigenStrainDifference =
                    CalculateLocEigenStrainDifference(EP, locStretchesAlpha, locStretchesBeta, i,j,k);

            dMatrix6x6 C_alpha = EP.ElasticConstants[alpha->index];
            dMatrix6x6 C_beta  = EP.ElasticConstants[beta->index ];

            for(int ii = 0; ii < 6; ii++)
            for(int jj = 0; jj < 6; jj++)
            {
                double delta_alpha = 1.0;
                double delta_beta  = 1.0;
                delta_alpha += EP.Gamma[alpha->index](ii,jj)*
                    (Tx(i,j,k) - EP.Tref[pIndexA]);
                delta_beta  += EP.Gamma[ beta->index](ii,jj)*
                    (Tx(i,j,k) - EP.Tref[pIndexB]);
                C_alpha(ii,jj) *= delta_alpha;
                C_beta(ii,jj)  *= delta_beta;
            }

            dG_AB += (ElasticStrainsCorr*((C_beta - C_alpha)*ElasticStrainsCorr))*0.5
                    - locEigenStrainDifference*(EP.EffectiveElasticConstants(i,j,k)*ElasticStrainsCorr);//eigenstrain difference

            dGab.Raw(i,j,k).add_asym1(alpha->index, beta->index, dG_AB);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}


void ElasticityKhachaturyan::CalculateChemicalPotentialContribution(
        const PhaseField& Phase, const ElasticProperties& EP, EquilibriumPartitionDiffusionBinary& DF)
{
    /** This function calculates the partial derivative of the mechanical
    energy with respect to the composition. This is later used in the EQPD-model
    to calculate an additional diffusion flux. */

    const size_t comp = DF.Comp;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DF.dMu, DF.dMu.Bcells(),)
    {
        vStrain locStrain = (EP.TotalStrains(i,j,k) - EP.EigenStrains(i,j,k));

        for (size_t n = 0; n < EP.Nphases; n++)
        {
            DF.dMu(i,j,k)({n}) = 0.0;
        }

        for(auto alpha = Phase.Fields(i, j, k).cbegin();
                 alpha != Phase.Fields(i, j, k).cend(); ++alpha)
        {
            double locdMuEl = 0.0;
            vStrain locLambda = VoigtStrain(EP.Lambda({alpha->index, comp}));

            for(int ii = 0; ii < 6; ii++)
            {
                locdMuEl -= locLambda[ii]*alpha->value*EP.Stresses(i, j, k)[ii];
                for(int jj = 0; jj < 6; jj++)
                {
                    locdMuEl += 0.5 * EP.Kappa({alpha->index, comp})(ii,jj)*alpha->value* locStrain[ii]*EP.ElasticConstants[alpha->index](ii,jj)*locStrain[jj];
                }
            }
            size_t pIndex = Phase.FieldsStatistics[alpha->index].Phase;
            DF.dMu(i,j,k)({pIndex}) += locdMuEl;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticityKhachaturyan::CalculateChemicalPotentialContribution(
        const PhaseField& Phase, const ElasticProperties& EP, ParabolicDiffusion& DF)
{
    /** This function calculates the partial derivative of the mechanical
    energy with respect to the composition. This is later used in the EQP-model
    to calculate an additional diffusion flux. */

    const size_t Ncomp = EP.Ncomp;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DF.ChemicalPotential, DF.ChemicalPotential.Bcells(),)
    {
        vStrain locStrain = (EP.TotalStrains(i,j,k) - EP.EigenStrains(i,j,k));

        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            for(auto alpha  = Phase.Fields(i, j, k).cbegin();
                     alpha != Phase.Fields(i, j, k).cend(); ++alpha)
            {
                double locdMuEl = 0.0;
                vStrain locLambda = VoigtStrain(EP.Lambda({alpha->index, comp}));

                for(int ii = 0; ii < 6; ii++)
                {
                    locdMuEl -= locLambda[ii]*alpha->value*EP.Stresses(i, j, k)[ii];
                    for(int jj = 0; jj < 6; jj++)
                    {
                        locdMuEl += 0.5 * EP.Kappa({alpha->index, comp})(ii,jj)*alpha->value* locStrain[ii]*EP.ElasticConstants[alpha->index](ii,jj)*locStrain[jj];
                    }
                }

                DF.ChemicalPotential(i,j,k)({comp}) += locdMuEl;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticityKhachaturyan::SetEffectiveElasticConstants(
        const PhaseField& Phase, ElasticProperties& EP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveElasticConstants, EP.EffectiveElasticConstants.Bcells(),)
    {
        EP.EffectiveElasticConstants(i,j,k).set_to_zero();
        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        {
            EP.EffectiveElasticConstants(i,j,k) +=
                         EP.ElasticConstants[alpha->index]*alpha->value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticityKhachaturyan::SetEffectiveElasticConstants(
        const PhaseField& Phase, ElasticProperties& EP, const Orientations& OR)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveElasticConstants,  EP.EffectiveElasticConstants.Bcells(),)
    {
        EP.EffectiveElasticConstants(i,j,k).set_to_zero();
        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        {
            EP.EffectiveElasticConstants(i,j,k) +=
                         EP.ElasticConstants[alpha->index].rotated(
                         OR.Quaternions(i,j,k).RotationMatrix)*alpha->value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticityKhachaturyan::SetEffectiveElasticConstants(
        const PhaseField& Phase, ElasticProperties& EP, const Composition& Cx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveElasticConstants, EP.EffectiveElasticConstants.Bcells(),)
    {
        EP.EffectiveElasticConstants(i,j,k).set_to_zero();
        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        {
            EP.EffectiveElasticConstants(i,j,k) +=
                CalculateLocElasticConstants(Phase,EP,Cx,i,j,k,alpha->index)*alpha->value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticityKhachaturyan::SetEffectiveElasticConstants(const PhaseField& Phase,
                                         ElasticProperties& EP, const Temperature& Tx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveElasticConstants, EP.EffectiveElasticConstants.Bcells(),)
    {
        EP.EffectiveElasticConstants(i,j,k).set_to_zero();

        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        {
            size_t pIndex = Phase.FieldsStatistics[alpha->index].Phase;

            for(int ii = 0; ii < 6; ii++)
            for(int jj = 0; jj < 6; jj++)
            {
                double deltaCij = EP.Gamma[alpha->index](ii,jj)*(Tx(i,j,k) - EP.Tref[pIndex])
                                   *EP.ElasticConstants[alpha->index](ii,jj);

                EP.EffectiveElasticConstants(i,j,k)(ii,jj) += alpha->value *
                        (EP.ElasticConstants[alpha->index](ii,jj) + deltaCij);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
}// namespace openphase
