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
#include "Info.h"
#include "Mechanics/ElasticityModels/ElasticitySteinbach.h"
#include "Mechanics/PlasticFlow/PlasticFlowNeuberMethods.h"
#include "Mechanics/ElasticProperties.h"
#include "Orientations.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Temperature.h"
#include "EquilibriumPartitionDiffusionBinary.h"
#include "VTK.h"

namespace openphase
{

using namespace std;

vStrain ElasticitySteinbach::CalculateLocEigenStrainDifference(
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

void ElasticitySteinbach::CalculateNeuberCorrection(
        vStress& ElasticStresses, const vStrain& ElasticStrains,
        const ElasticProperties& EP, const dMatrix6x6 locCompliance,
        const int i, const int j, const int k,
        const int pIndexA, const int pIndexB)
{
    if (EP.NeuberCorrection[pIndexA].Active or EP.NeuberCorrection[pIndexB].Active)
    {
        vStress deviatoricStresses;
        vStrain deviatoricStrains;
        vStrain dummyStrain;
        double StressTrace = EP.Stresses(i, j, k).trace();
        double StrainTrace = ElasticStrains.trace();
        for(int n = 0; n < 3; n++)
        {
            deviatoricStresses[n] = EP.Stresses(i, j, k)[n] - 0.3*StressTrace;
            deviatoricStrains[n] = ElasticStrains[n] - 0.3*StrainTrace;
        }
        if(EP.NeuberCorrection[pIndexA].Active and !EP.NeuberCorrection[pIndexB].Active)
        {
            PlasticFlowNeuberMethods::getNeuberDataRO(
                    deviatoricStresses,
                    deviatoricStrains,
                    ElasticStresses,
                    dummyStrain,
                    EP.NeuberCorrection[pIndexA].YoungsModulus,
                    EP.NeuberCorrection[pIndexA].YieldStrength,
                    EP.NeuberCorrection[pIndexA].Hardening);
        }
        if(!EP.NeuberCorrection[pIndexA].Active and EP.NeuberCorrection[pIndexB].Active)
        {
            PlasticFlowNeuberMethods::getNeuberDataRO(
                    deviatoricStresses,
                    deviatoricStrains,
                    ElasticStresses,
                    dummyStrain,
                    EP.NeuberCorrection[pIndexB].YoungsModulus,
                    EP.NeuberCorrection[pIndexB].YieldStrength,
                    EP.NeuberCorrection[pIndexB].Hardening);
        }
        if(EP.NeuberCorrection[pIndexA].Active and EP.NeuberCorrection[pIndexB].Active)
        {
            if(EP.NeuberCorrection[pIndexA].YieldStrength <= EP.NeuberCorrection[pIndexB].YieldStrength)
            {
                PlasticFlowNeuberMethods::getNeuberDataRO(
                        deviatoricStresses,
                        deviatoricStrains,
                        ElasticStresses,
                        dummyStrain,
                        EP.NeuberCorrection[pIndexA].YoungsModulus,
                        EP.NeuberCorrection[pIndexA].YieldStrength,
                        EP.NeuberCorrection[pIndexA].Hardening);
            }
            else
            {
                PlasticFlowNeuberMethods::getNeuberDataRO(
                        deviatoricStresses,
                        deviatoricStrains,
                        ElasticStresses,
                        dummyStrain,
                        EP.NeuberCorrection[pIndexB].YoungsModulus,
                        EP.NeuberCorrection[pIndexB].YieldStrength,
                        EP.NeuberCorrection[pIndexB].Hardening);
            }
        }
        for(int n = 0; n < 3; n++)
        {
            ElasticStresses[n] += 0.3*StressTrace;
        }
    }
    else
    {
        ElasticStresses = EP.Stresses(i, j, k);
    }
}

dMatrix6x6 ElasticitySteinbach::CalculateLocElasticConstants(
        const PhaseField& Phase, const ElasticProperties& EP,
        const Composition& Cx, const int i, const int j, const int k,
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

void ElasticitySteinbach::CalculateDrivingForce(const PhaseField& Phase,
        const ElasticProperties& EP, DrivingForce& dGab)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.Stresses, 0,)
    if(Phase.Interface(i,j,k))
    {
        dMatrix6x6 locCompliance = EP.EffectiveElasticConstants(i,j,k).inverted();
        vStrain ElasticStrains = locCompliance*EP.Stresses(i, j, k);
        vStress ElasticStresses = EP.Stresses(i, j, k);

        for(auto alpha = Phase.Fields(i, j, k).cbegin();
                 alpha != Phase.Fields(i, j, k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1;
                  beta != Phase.Fields(i, j, k).cend();  ++beta)
        {
            const size_t pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
            const size_t pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

            CalculateNeuberCorrection(ElasticStresses, ElasticStrains,
                    EP, locCompliance, i, j, k, pIndexA, pIndexB);

            dMatrix3x3 locStretchesAlpha = EP.TransformationStretches[alpha->index];
            dMatrix3x3 locStretchesBeta  = EP.TransformationStretches[ beta->index];

            vStrain locEigenStrainDifference =
                CalculateLocEigenStrainDifference(EP, locStretchesAlpha, locStretchesBeta, i,j,k);

            double dG_AB = 0.0;
            for(int ii = 0; ii < 6; ii++)
            {
                dG_AB -= ElasticStresses[ii]*locEigenStrainDifference[ii];

                for(int jj = 0; jj < 6; jj++)
                {
                    dG_AB += 0.5*ElasticStresses[ii]*
                                (EP.Compliances[alpha->index](ii,jj) -
                                 EP.Compliances[beta->index](ii,jj))*
                                 ElasticStresses[jj];
                }
            }
            dGab.Raw(i,j,k).add_asym1(alpha->index, beta->index, dG_AB);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySteinbach::CalculateDrivingForce(const PhaseField& Phase,
        const ElasticProperties& EP, const Composition& Cx, DrivingForce& dGab)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.Stresses, 0,)
    if(Phase.Interface(i,j,k))
    {
        dMatrix6x6 locCompliance = EP.EffectiveElasticConstants(i,j,k).inverted();
        vStrain ElasticStrains = locCompliance*EP.Stresses(i, j, k);
        vStress ElasticStresses = EP.Stresses(i, j, k);

        for(auto alpha = Phase.Fields(i, j, k).cbegin();
                 alpha != Phase.Fields(i, j, k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1;
                  beta != Phase.Fields(i, j, k).cend();  ++beta)
        {
            const size_t pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
            const size_t pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

            CalculateNeuberCorrection(ElasticStresses, ElasticStrains,
                    EP, locCompliance, i, j, k, pIndexA, pIndexB);

            dMatrix3x3 locStretchesAlpha = EP.TransformationStretches[alpha->index];
            dMatrix3x3 locStretchesBeta  = EP.TransformationStretches[ beta->index];

            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
            {
                locStretchesAlpha += EP.Lambda({alpha->index, comp})*(Cx.MoleFractions(i,j,k)({pIndexA, comp}) - EP.Cref({pIndexA, comp}));
                locStretchesBeta  += EP.Lambda({ beta->index, comp})*(Cx.MoleFractions(i,j,k)({pIndexB, comp}) - EP.Cref({pIndexB, comp}));
            }

            vStrain locEigenStrainDifference =
                    CalculateLocEigenStrainDifference(EP, locStretchesAlpha, locStretchesBeta, i,j,k);

            const dMatrix6x6 Compliences_alpha = CalculateLocElasticConstants(Phase,EP,Cx,i,j,k,alpha->index).inverted();
            const dMatrix6x6 Compliences_beta  = CalculateLocElasticConstants(Phase,EP,Cx,i,j,k,beta->index ).inverted();

            double dG_AB = 0.0;
            for(int ii = 0; ii < 6; ii++)
            {
                dG_AB -= ElasticStresses[ii]*locEigenStrainDifference[ii];

                for(int jj = 0; jj < 6; jj++)
                {
                    dG_AB += 0.5*ElasticStresses[ii]*
                                (Compliences_alpha(ii,jj) -
                                 Compliences_beta(ii,jj))*
                                 ElasticStresses[jj];
                }
            }
            dGab.Raw(i,j,k).add_asym1(alpha->index, beta->index, dG_AB);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySteinbach::CalculateDrivingForce(const PhaseField& Phase,
        const ElasticProperties& EP, const Temperature& Tx, DrivingForce& dGab)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.Stresses, 0,)
    if(Phase.Interface(i,j,k))
    {
        dMatrix6x6 locCompliance = EP.EffectiveElasticConstants(i,j,k).inverted();
        vStrain ElasticStrains = locCompliance*EP.Stresses(i, j, k);
        vStress ElasticStresses = EP.Stresses(i, j, k);

        for(auto alpha = Phase.Fields(i, j, k).cbegin();
                 alpha != Phase.Fields(i, j, k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1;
                  beta != Phase.Fields(i, j, k).cend();  ++beta)
        {
            const size_t pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
            const size_t pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

            CalculateNeuberCorrection(ElasticStresses, ElasticStrains,
                    EP, locCompliance, i, j, k, pIndexA, pIndexB);

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

            dMatrix6x6 Compliences_alpha = C_alpha.inverted();
            dMatrix6x6 Compliences_beta  = C_beta.inverted();

            double dG_AB = 0.0;
            for(int ii = 0; ii < 6; ii++)
            {
                dG_AB -= ElasticStresses[ii]*locEigenStrainDifference[ii];

                for(int jj = 0; jj < 6; jj++)
                {
                    dG_AB += 0.5*ElasticStresses[ii]*
                                (Compliences_alpha(ii,jj) -
                                 Compliences_beta(ii,jj))*
                                 ElasticStresses[jj];
                }
            }
            dGab.Raw(i,j,k).add_asym1(alpha->index, beta->index, dG_AB);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySteinbach::CalculateChemicalPotentialContribution(
        const PhaseField& Phase, const ElasticProperties& EP,
        EquilibriumPartitionDiffusionBinary& DF)
{
    /** This function calculates the partial derivative of the mechanical
    energy density with respect to the composition.*/

    const size_t comp = DF.Comp;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.DeformationGradientsTotal, EP.DeformationGradientsTotal.Bcells(),)
    {
        for(size_t n = 0; n < EP.Nphases; n++)
        {
            DF.dMu(i,j,k)({n}) = 0.0;
        }

        for(auto alpha = Phase.Fields(i, j, k).cbegin();
                 alpha != Phase.Fields(i, j, k).cend(); ++alpha)
        {
            double locdMuEl = 0.0;
            size_t index = alpha->index;
            size_t pIndex = Phase.FieldsStatistics[index].Phase;
            vStrain locLambda = VoigtStrain(EP.Lambda({index, comp}));
            for(int ii = 0; ii < 6; ii++)
            {
                locdMuEl -= locLambda[ii]*alpha->value*EP.Stresses(i, j, k)[ii];
                for(int jj = 0; jj < 6; jj++)
                {
                    locdMuEl += 0.5 * EP.Kappa({index, comp})(ii,jj)*alpha->value*
                                     (EP.TotalStrains(i,j,k)[ii] -
                                      EP.EigenStrains(i,j,k)[ii]) *
                                      EP.ElasticConstants[index](ii,jj) *
                                     (EP.TotalStrains(i,j,k)[jj] -
                                      EP.EigenStrains(i,j,k)[jj]);
                }
            }
            DF.dMu(i,j,k)({pIndex}) += locdMuEl;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void ElasticitySteinbach::CalculateChemicalPotentialContribution(
        const PhaseField& Phase, const ElasticProperties& EP,
        ParabolicDiffusion& DF)
{
    //TODO
}
/*void ElasticitySteinbach::CalculateChemicalPotentialContrib(PhaseField& Phase,
                            ElasticProperties& EP, ThermodynamicProperties& TP)
{
    const size_t Ncomp = EP.Ncomp;
    #pragma omp parallel for collapse(OMP_COLLAPSE_LOOPS) schedule(dynamic, OMP_DYNAMIC_CHUNKSIZE)
    //OP_STORAGE_LOOP_ENTIRE(i,j,k,EP.Strains)
    OP_STORAGE_LOOP_INTERIOR(i,j,k,EP.Strains)
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        if(Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                double locdMuEl = 0.0;
                size_t index = alpha->index;
                size_t pIndex = Phase.FieldsStatistics[index].Phase;

                for(int ii = 0; ii < 6; ii++)
                for(int jj = 0; jj < 6; jj++)
                {
                    locdMuEl += 0.5 * EP.Kappa({index, comp})(ii,jj)*alpha->value*
                                     (EP.Strains(i,j,k)[ii] -
                                      EP.EffectiveEigenStrains(i,j,k)[ii]) *
                                     (EP.Strains(i,j,k)[jj] -
                                      EP.EffectiveEigenStrains(i,j,k)[jj]) +

                                      EP.Lambda({index, comp})[ii]*alpha->value*
                                      EP.EffectiveElasticConstants(i,j,k)(ii,jj)*
                                      (EP.Strains(i,j,k)[jj] -
                                       EP.EffectiveEigenStrains(i,j,k)[jj]);
                }
                TP.ChemicalPotential(i,j,k)({pIndex, comp}) += locdMuEl;
            }
        }
        else
        {
            double locdMuEl = 0.0;
            size_t index = Phase.Fields(i,j,k).front().index;
            size_t pIndex = Phase.FieldsStatistics[index].Phase;

            for(int ii = 0; ii < 6; ii++)
            for(int jj = 0; jj < 6; jj++)
            {
                locdMuEl += 0.5 * EP.Kappa({index,comp})(ii,jj)*
                                 (EP.Strains(i,j,k)[ii] -
                                  EP.EffectiveEigenStrains(i,j,k)[ii]) *
                                 (EP.Strains(i,j,k)[jj] -
                                  EP.EffectiveEigenStrains(i,j,k)[jj]) +

                                  EP.Lambda({index, comp})[ii]*
                                  EP.EffectiveElasticConstants(i,j,k)(ii,jj)*
                                  (EP.Strains(i,j,k)[jj] -
                                   EP.EffectiveEigenStrains(i,j,k)[jj]);
            }
            TP.ChemicalPotential(i,j,k)({pIndex, comp}) += locdMuEl;
        }
    }
}*/


void ElasticitySteinbach::SetEffectiveElasticConstants(const PhaseField& Phase,
                                                       ElasticProperties& EP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveElasticConstants, EP.EffectiveElasticConstants.Bcells(),)
    {
        EP.EffectiveElasticConstants(i,j,k).set_to_zero();
        if(Phase.Interface(i,j,k))
        {
            dMatrix6x6 TempCompliences;
            TempCompliences.set_to_zero();

            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            {
                TempCompliences += (EP.Compliances[alpha->index]*alpha->value);
            }
            EP.EffectiveElasticConstants(i,j,k) = TempCompliences.inverted();
        }
        else
        {
            EP.EffectiveElasticConstants(i,j,k) =
              EP.ElasticConstants[Phase.Fields(i,j,k).front().index];
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySteinbach::SetEffectiveElasticConstants(const PhaseField& Phase,
                                        ElasticProperties& EP, const Orientations& OR)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveElasticConstants, EP.EffectiveElasticConstants.Bcells(),)
    {
        EP.EffectiveElasticConstants(i,j,k).set_to_zero();
        if(Phase.Interface(i,j,k))
        {
            dMatrix6x6 TempCompliences;
            TempCompliences.set_to_zero();

            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            {
                dMatrix6x6 tempStiff =
                   EP.ElasticConstants[alpha->index].rotated(OR.Quaternions(i,j,k).RotationMatrix);
                TempCompliences += (tempStiff.inverted()*alpha->value);
            }
            EP.EffectiveElasticConstants(i,j,k) = TempCompliences.inverted();
        }
        else
        {
            EP.EffectiveElasticConstants(i,j,k) =
               EP.ElasticConstants[Phase.Fields(i,j,k).front().index].rotated(OR.Quaternions(i,j,k).RotationMatrix);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySteinbach::SetEffectiveElasticConstants(
        const PhaseField& Phase, ElasticProperties& EP, const Composition& Cx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveElasticConstants, EP.EffectiveElasticConstants.Bcells(),)
    {
        EP.EffectiveElasticConstants(i,j,k).set_to_zero();
        if(Phase.Interface(i,j,k))
        {
            dMatrix6x6 locEffectiveCompliances;
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            {
                locEffectiveCompliances +=
                    CalculateLocElasticConstants(Phase,EP,Cx,i,j,k,alpha->index).inverted()*alpha->value;
            }
            EP.EffectiveElasticConstants(i,j,k) = locEffectiveCompliances.inverted();
        }
        else
        {
            const size_t index = Phase.Fields(i,j,k).front().index;
            EP.EffectiveElasticConstants(i,j,k) =
                CalculateLocElasticConstants(Phase,EP,Cx,i,j,k,index);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySteinbach::SetEffectiveElasticConstants(const PhaseField& Phase,
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
