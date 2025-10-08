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
 *   File created :   2014
 *   Main contributors :   Oleg Shchyglo; Raphael Schiedung
 *
 */

#ifndef INTERACTIONSOLIDFLUID_H
#define INTERACTIONSOLIDFLUID_H

#include "Base/Includes.h"
namespace openphase
{

class BoundaryConditions;
class PhaseField;
class Settings;
class Velocities;

class InteractionSolidFluid                                                     ///< Processes solid liquid interaction
{
 public:
    static size_t CalculateSolidVelocities(PhaseField& Phase, Velocities& Vel,
            const BoundaryConditions& BC, const Settings& locSettings,
            const double dt, double* VLimit = nullptr,
            bool EnforceZeroTotalMomentum = false,
            bool EnforceZeroTotalAngularMomentum = false);                      ///< Calculates velocities of solid particles

    static void CollectGrainsStatistics(PhaseField& Phase,
            const BoundaryConditions& BC, const Settings& locSettings);         ///<  Collects center of mass position and moments of inertia for each phase field.
    static void CalculateCenterOfMass(PhaseField& Phase);                       ///<  Substep of CollectGrainStatistics
    static void CollectGrainsStatisticsStepTwo(PhaseField& Phase,
            const BoundaryConditions& BC, const Settings& locSettings);         ///<  Substep of CollectGrainStatistics
    static void CalculateCenterOfMassWithPeriodicBoundaryConditions(PhaseField& Phase);

 protected:

 private:
};

} //namespace openphase
#endif
