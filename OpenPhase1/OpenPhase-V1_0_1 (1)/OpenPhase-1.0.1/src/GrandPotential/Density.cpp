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


#ifndef GRANDPOTENTIALDENSITIES_H
#define GRANDPOTENTIALDENSITIES_H

#include "Settings.h"
#include "GrandPotential/Density.h"
#include "GrandPotential/PhaseDensities/CompressibleFluid.h"
#include "GrandPotential/PhaseDensities/IdealGas.h"
#include "GrandPotential/PhaseDensities/ImplicitParabolic.h"
#include "GrandPotential/PhaseDensities/Parabolic.h"
#include "Base/UserInterface.h"
#include "VTK.h"
#include "Info.h"

namespace openphase::GrandPotential
{
Density::Density(Settings& locSettings, std::string filename)
{
    InitializeAndReadInput(locSettings, filename);
}
void Density::InitializeAndReadInput(Settings& locSettings, std::string filename)
{
    PhaseNames = locSettings.PhaseNames;
    Nphases    = locSettings.Nphases;

    std::fstream inp(filename.c_str(), std::ios::in);
    if (!inp)
    {
        Info::WriteExit("File \"" + filename + "\" could not be opened", thisclassname, "ReadInput");
        std::exit(EXIT_FAILURE);
    };
    std::stringstream inp_data;
    inp_data << inp.rdbuf();
    inp.close();

    Info::WriteLine();
    Info::WriteLineInsert("Grand Potential Density");
    Info::WriteStandard("Source", filename);

    int moduleLocation = UserInterface::FindModuleLocation(inp_data, "GrandPotentialDensity");

    for(size_t PhaseIdx = 0; PhaseIdx < Nphases; PhaseIdx++)
    {

        std::string converter = ""+PhaseNames[PhaseIdx];
        std::string str_omega = UserInterface::ReadParameterK(inp_data, moduleLocation, "OMEGA_"+converter);
        if      (str_omega == "PARABOLIC" )         storage.push_back(new Parabolic         (PhaseIdx));
        else if (str_omega == "IMPLICITPARABOLIC")  storage.push_back(new ImplicitParabolic (PhaseIdx));
        else if (str_omega == "IDEALGAS")           storage.push_back(new IdealGas          (PhaseIdx));
        else if (str_omega == "COMPRESSIBLEFLUID")  storage.push_back(new CompressibleFluid (PhaseIdx));
        else
        {
            std::cerr << "Phase Model "+str_omega+" is not implemented!\n";
            std::terminate();
        }
    }
    for (auto omega : storage) omega->Initialize(locSettings);
    for (auto omega : storage) omega->ReadInput(inp_data, moduleLocation);
}
}// namespace openphase::GrandPotential
#endif
