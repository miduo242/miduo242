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
 *   Main contributors :   Marvin Tegeler
 *
 */

#include "Termination.h"
#include "Base/UserInterface.h"
#include "Settings.h"
#include "Info.h"
#include "PhaseField.h"
#include "Temperature.h"

namespace openphase
{

using namespace std;

void Termination::Initialize(Settings& locSettings)
{
    thisclassname = "Termination";
    thisobjectname = thisclassname;
    criterion = Criterion::None;
    Info::WriteStandard(thisclassname, "Initialized");
}

void Termination::ReadInput(const string InputFileName)
{
    if(thisclassname == "")
    {
        Info::WriteExit("Attempt to use uninitialized object", "Settings", "ReadInput()");
        exit(1);
    }

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        Info::WriteExit("File " + InputFileName + " could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };

    Info::WriteLine();
    Info::WriteLine();
    Info::WriteStandard("Run time control parameters source", InputFileName);

    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);

    inp.close();
}

void Termination::ReadInput(stringstream& inp)
{
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    std::string crit = UserInterface::ReadParameterF(inp, moduleLocation, string("Criterion"));

    if (crit == "VolumeFraction")
    {
        criterion = Criterion::VolumeFraction;
    }
    if (crit == "AverageTemperature")
    {
        criterion = Criterion::AverageTemperature;
    }

    if (criterion == Criterion::VolumeFraction)
    {
        SignificantPhase = UserInterface::ReadParameterI(inp, moduleLocation, string("SignificantPhase"));
    }

    EndValue = UserInterface::ReadParameterD(inp, moduleLocation, string("EndValue"));

    std::string comp = UserInterface::ReadParameterF(inp, moduleLocation, string("Comparison"));

    if (comp == "Greater")
    {
        comparison = Comparison::Greater;
    }
    if (comp == "Smaller")
    {
        comparison = Comparison::Smaller;
    }

    Info::WriteLine();
}

bool Termination::ReachedTerminationCondition(PhaseField& Phi, Temperature& Tx)
{
    bool end = false;
    if (criterion == Criterion::VolumeFraction)
    {
        if (comparison == Comparison::Greater)
        {
            if (Phi.FractionsTotal[SignificantPhase] > EndValue) end = true;
        }
        if (comparison == Comparison::Smaller)
        {
            if (Phi.FractionsTotal[SignificantPhase] < EndValue) end = true;
        }
    }
    if (criterion == Criterion::AverageTemperature)
    {
        double aveT = Tx.Average();
        if (comparison == Comparison::Greater)
        {
            if (aveT > EndValue) end = true;
        }
        if (comparison == Comparison::Smaller)
        {
            if (aveT < EndValue) end = true;
        }
    }
    return end;
}

} //namespace openphase
