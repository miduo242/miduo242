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
 *   Main contributors :   Oleg Shchyglo; Dmitry Medvedev; Raphael Schiedung
 *
 */

#ifndef INTERACTIONSOLIDSOLID_H
#define INTERACTIONSOLIDSOLID_H

#include "PhaseField.h"

namespace openphase
{

class D3Q27;
class FlowSolverLBM;
class Settings;
class Velocities;

class InteractionSolidSolid : public OPObject                                   ///<  Handles the interaction between solid particles in the fluid environment
{
 public:

    void Initialize(Settings& locSettings) override;                            ///< Allocates memory and initiates global settings, structure described in OPSettings.h file
    void ReadInput(std::string InputFileName) override;                         ///< Reads input parameters from a file

    static void Calculate(PhaseField& Phase,
            const FlowSolverLBM& LBM,
            const double dt, const int order = 4, const int cutoff = 3,
            const double lbStrength = 100.0, const double elastic = 0.0);       ///< Calculates the solid-solid interaction the entire simulation domain, but only near LMB.Obstacle!

    static void Calculate(PhaseField& Phase,
            const double dt, const int order = 4, const int cutoff = 3,
            const double strength = 1.0, const double elastic = 0.0);           ///< Calculates the solid-solid interaction the entire simulation domain

    static void AdvectSolid(PhaseField& Phase, const BoundaryConditions& BC,    ///< NOTE: deprecated use AdvectionHR!
            const double dt);

    static void PreserveVolume(PhaseField& Phase,
            const std::vector<double>& RefVolume, const size_t _index,
            const double dt);                                                   ///< NOTE: deprecated use AdvectionHR!

    static void SetRefVolume(const PhaseField& Phase,
            std::vector<double>& RefVolume);                                    ///< NOTE: deprecated use AdvectionHR!

    static double VolumeError(const PhaseField& Phase,
            const std::vector<double>& RefVolume, const size_t index);          ///< NOTE: deprecated use AdvectionHR!
 protected:
 private:
    static void CalculateLocal(PhaseField& Phase,
            const int i, const int j, const int k,
            const double dt, const int order = 4, const int cutoff = 3,
            const double strength = 100.0, const double elastic = 0.0);         ///< Calculates the solid-solid interaction at the point (i,j,k)
};

} //namespace openphase
#endif
