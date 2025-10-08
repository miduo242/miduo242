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
#include "GrandPotential/PhaseDensities/Parabolic.h"
#include "Info.h"
#include "Settings.h"
#include <cassert>

namespace openphase::GrandPotential
{
void Parabolic::Initialize(Settings& locSettings)
{
    PhaseDensity::Initialize(locSettings);

    EPS.assign(Ncomp, 0.0);
    C0 .assign(Ncomp, 0.0);
}
void Parabolic::ReadInput (std::stringstream& InputFile, int moduleLocation)
{
    Info::Write("");
    Info::Write("Parabolic grad potential density material properties for "+PhaseNames[PhaseIdx]);
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        std::string converter = ""+PhaseNames[PhaseIdx]+"_"+ElementNames[comp];
        C0[comp] = UserInterface::ReadParameterD(InputFile, moduleLocation, "C0_"+converter);
        if (C0[comp] < 0) Info::WriteWarning("C0 needs to positive!",thisclassname,"Read Input");
    }
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        std::string converter = ""+PhaseNames[PhaseIdx]+"_"+ElementNames[comp];
        EPS[comp] = UserInterface::ReadParameterD(InputFile, moduleLocation, "EPS_"+converter);
        if (EPS[comp] <= 0) Info::WriteWarning("EPS needs to positive!",thisclassname,"Read Input");
    }
    Info::Write("");
}
}// namespace openphase::GrandPotential
