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
 *   File created :   2021
 *   Main contributors :   Raphael Schiedung
 *
 */

#ifndef GRANDPOTENTIALDENSITY_H
#define GRANDPOTENTIALDENSITY_H

#include <cstddef>
#include <string>
#include <vector>

#include "PhaseField.h"
#include "Settings.h"
#include "GrandPotential/PhaseDensities/PhaseDensity.h"

namespace openphase
{
class Settings;
namespace GrandPotential
{
class Density                                                                   ///< Base class of grand potential densities / equations of state
{
 public:
    Density () = default;                                                       ///< Explicitly defaulted  automatically generated constructor
    Density (Settings& locSettings, std::string InputFileName);                 ///< Allocates and initialises memory
    void InitializeAndReadInput (Settings& locSettings, std::string FileName);
    PhaseDensity& operator() (size_t PhaseIdx) const
    {
        assert(PhaseIdx < Nphases);
        return *storage[PhaseIdx];
    };
    double operator() (int i, int j, int k, const PhaseField& Phase) const     ///< Computes grand potential density at (i,j,k)
    {
        double omega = 0.0;
        for(auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        {
            size_t PhaseIdx = Phase.FieldsStatistics[alpha->index].Phase;
            omega += alpha->value*(*storage[PhaseIdx])(i,j,k);
        }
        return omega;
    };

    //NOTE: The set-methods set references to external storage classes.
    // A set-method needs to be called explicitly only once but will it be
    // implicitly invoked by the compiler on each function call when the
    // external external memory is accessed during that function call!
    void Set(const double& Temperature, const Storage3D<double,1>& ChemicalPotential)
    {
         for(auto omega : storage) omega->Set(Temperature,ChemicalPotential);
    };
    void Set(const double& Temperature, const Storage3D<double,1>& ChemicalPotential, const int& GravityDirection)
    {
         for(auto omega : storage) omega->Set(Temperature,ChemicalPotential,GravityDirection);
    };

    void WriteVTK (long tStep, const Settings& locSettings, const PhaseField& Phase, long precision=16) const; ///< Writes mole fraction in VTK format (.vts file)

    static constexpr auto thisclassname = "GrandPotentialDensity";

 protected:
    size_t Nphases;                                                             ///< Number of thermodynamic phases
    std::vector<std::string> PhaseNames;                                        ///< Names of corresponding thermodynamic phases
    std::vector<PhaseDensity*> storage;                                         ///< Storage of grand potential densities for each phase
};
}// namespace openphase::GrandPotential
}// namespace openphase
#endif
