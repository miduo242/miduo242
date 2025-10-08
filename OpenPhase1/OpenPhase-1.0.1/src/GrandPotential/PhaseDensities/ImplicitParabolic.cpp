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
#include "GrandPotential/PhaseDensities/ImplicitParabolic.h"
#include "Info.h"
#include "Settings.h"
#include <cassert>
namespace openphase::GrandPotential
{
void ImplicitParabolic::Initialize(Settings& locSettings)
{
    Parabolic::Initialize(locSettings);
    Zeta.assign(Ncomp, 0.0);
}
void ImplicitParabolic::ReadInput(std::stringstream& InputFile, int moduleLocation)
{
    Parabolic::ReadInput(InputFile,moduleLocation);
    Info::Write("");
    Info::Write("Implicit parabolic grand potential density material properties for "+PhaseNames[PhaseIdx]);
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        std::string converter = ""+PhaseNames[PhaseIdx]+"_"+ElementNames[comp];
        Zeta[comp] = UserInterface::ReadParameterD(InputFile, moduleLocation, "Zeta_"+converter);
        if (Zeta[comp] <= 0) Info::WriteWarning("Zeta needs to positive!",thisclassname,"Read Input");
    }
    std::string converter = ""+PhaseNames[PhaseIdx];
    precision     = UserInterface::ReadParameterD(InputFile, moduleLocation, "PREC_"+converter);
    MaxIterations = UserInterface::ReadParameterI(InputFile, moduleLocation, "MAXI_"+converter);
    Info::Write("");
}
}// namespace openphase::GrandPotential

