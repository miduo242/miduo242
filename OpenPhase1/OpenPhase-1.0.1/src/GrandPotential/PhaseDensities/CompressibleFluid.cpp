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


#include "Base/UserInterface.h"
#include "GrandPotential/PhaseDensities/CompressibleFluid.h"
#include "Info.h"
#include "PhysicalConstants.h"
#include "Settings.h"
#include <cassert>

namespace openphase::GrandPotential
{
void CompressibleFluid::Initialize(Settings& locSettings)
{
    PhaseDensity::Initialize(locSettings);

    if (Ncomp != 1)
    {
        Info::WriteExit("Model is only suitable for single component simulation!", thisclassname, "Initialize");
        std::terminate();
    }

    BulkModulus               .assign(Ncomp, 0.0);
    ReferenceChemicalPotential.assign(Ncomp, 0.0);
    ReferenceMolarVolume      .assign(Ncomp, 0.0);
    ReferencePressure         .assign(Ncomp, 0.0);
    MolarMass                 .assign(Ncomp, 0.0);
}
void CompressibleFluid::ReadInput(std::stringstream& InputFile, int moduleLocation)
{
    Info::Write("");
    Info::Write("Compressible non-gas grand potential density properties of "+PhaseNames[PhaseIdx]);
    std::string converter = "_"+PhaseNames[PhaseIdx]+"_"+ElementNames[0];
    BulkModulus                [0] = UserInterface::ReadParameterD(InputFile, moduleLocation, "K"  +converter);
    ReferenceChemicalPotential [0] = UserInterface::ReadParameterD(InputFile, moduleLocation, "MU0"+converter);
    ReferencePressure          [0] = UserInterface::ReadParameterD(InputFile, moduleLocation, "P0" +converter);
    ReferenceMolarVolume       [0] = UserInterface::ReadParameterD(InputFile, moduleLocation, "VM0"+converter);
    MolarMass                  [0] = UserInterface::ReadParameterD(InputFile, moduleLocation, "M_" +ElementNames[0]);
    GravityAcceleration            = UserInterface::ReadParameterD(InputFile, moduleLocation, "g");
    Info::Write("");
}
}// namespace openphase::GrandPotential
