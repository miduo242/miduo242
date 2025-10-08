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

#ifndef ELASTICITYSTEINBACH_H
#define ELASTICITYSTEINBACH_H

#include "Base/Includes.h"

namespace openphase
{
class Composition;
class DrivingForce;
class ElasticProperties;
class EquilibriumPartitionDiffusionBinary;
class Orientations;
class ParabolicDiffusion;
class PhaseField;
class Plasticity;
class Settings;
class Temperature;
class dMatrix6x6;
class vStrain;
class vStress;


class ElasticitySteinbach                                                       ///< Steinbach's approximation for elasticity parameters
{
 public:

    /// Calculates and sets effective elastic constants in ElasticProperties
    static void SetEffectiveElasticConstants(const PhaseField& Phase, ElasticProperties& EP);
    static void SetEffectiveElasticConstants(const PhaseField& Phase, ElasticProperties& EP, const Orientations& OR);
    static void SetEffectiveElasticConstants(const PhaseField& Phase, ElasticProperties& EP, const Composition& Cx);
    static void SetEffectiveElasticConstants(const PhaseField& Phase, ElasticProperties& EP, const Temperature& Tx);

    static void CalculateDrivingForce(const PhaseField& Phase, const ElasticProperties& EP, DrivingForce& dGab);
    static void CalculateDrivingForce(const PhaseField& Phase, const ElasticProperties& EP, const Composition& Cx, DrivingForce& dGab);
    static void CalculateDrivingForce(const PhaseField& Phase, const ElasticProperties& EP, const Temperature& Tx, DrivingForce& dGab);

    static void CalculateChemicalPotentialContribution(const PhaseField& Phase, const ElasticProperties& EP, EquilibriumPartitionDiffusionBinary& DF);
    static void CalculateChemicalPotentialContribution(const PhaseField& Phase, const ElasticProperties& EP, ParabolicDiffusion& DF);

 protected:
 private:

    static vStrain CalculateLocEigenStrainDifference(
            const ElasticProperties& EP,
            const dMatrix3x3 locStretchesAlpha,
            const dMatrix3x3 locStretchesBeta,
            const int i, const int j, const int k);
    static void CalculateNeuberCorrection(
            vStress& ElasticStresses, const vStrain& ElasticStrains,
            const ElasticProperties& EP, const dMatrix6x6 locCompliance,
            const int i, const int j, const int k,
            const int pIndexA, const int pIndexB);

    static dMatrix6x6 CalculateLocElasticConstants(
            const PhaseField& Phase,
            const ElasticProperties& EP, const Composition& Cx,
            const int i, const int j, const int k,
            const size_t idx);
};
}// namespace openphase
#endif
