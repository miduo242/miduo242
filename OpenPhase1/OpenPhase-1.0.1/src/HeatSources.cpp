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
 *   Main contributors :   Oleg Shchyglo; Marvin Tegeler; Matthias Stratmann
 *
 */

#include "HeatSources.h"
#include "Info.h"
#include "Settings.h"
#include "PhaseField.h"
#include "Temperature.h"
#include "HeatDiffusion.h"
#include "RunTimeControl.h"

namespace openphase
{
using namespace std;

void HeatSources::Initialize(Settings& locSettings)
{
    thisclassname = "HeatSources";

    Nphases = locSettings.Nphases;
    ConsiderHeatSources = false;

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void HeatSources::ReadInput(const std::string InputFileName)
{
    Info::WriteLineInsert("HeatSources input");
    Info::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };

    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);
    inp.close();

    Info::WriteLine();
}

void HeatSources::ReadInput(std::stringstream& inp)
{


    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    bool ThereAreSources = false;
    if(moduleLocation != -1)
    {
        ThereAreSources = true;
    }

    size_t idx = 0;
    while(ThereAreSources)
    {
        stringstream converter;
        converter << idx;
        string counter = converter.str();
        ThereAreSources = UserInterface::ReadParameterB(inp, moduleLocation, string("Source_") + counter, false, false);
        if(not ThereAreSources)
        {
            break;
        }
        HeatSourceStructure locHeatSource;
        locHeatSource.Active = false;
        locHeatSource.Value  = UserInterface::ReadParameterD(inp, moduleLocation, string("Value_") + counter, true, 0.0);

        bool SourceTypeSet = false;
        string SourceType = UserInterface::ReadParameterK(inp, moduleLocation, string("Type_") + counter, true, "NONE");
        if(SourceType == "BCX0"       ) {locHeatSource.Type = HeatSourceTypes::BCX0;        SourceTypeSet = true;}
        if(SourceType == "BCXN"       ) {locHeatSource.Type = HeatSourceTypes::BCXN;        SourceTypeSet = true;}
        if(SourceType == "BCY0"       ) {locHeatSource.Type = HeatSourceTypes::BCY0;        SourceTypeSet = true;}
        if(SourceType == "BCYN"       ) {locHeatSource.Type = HeatSourceTypes::BCYN;        SourceTypeSet = true;}
        if(SourceType == "BCZ0"       ) {locHeatSource.Type = HeatSourceTypes::BCZ0;        SourceTypeSet = true;}
        if(SourceType == "BCZN"       ) {locHeatSource.Type = HeatSourceTypes::BCZN;        SourceTypeSet = true;}
        if(SourceType == "Phase"      ) {locHeatSource.Type = HeatSourceTypes::Phase;       SourceTypeSet = true;}
        if(SourceType == "Ellipsoidal") {locHeatSource.Type = HeatSourceTypes::Ellipsoidal; SourceTypeSet = true;}
        if(SourceType == "Rectangular") {locHeatSource.Type = HeatSourceTypes::Rectangular; SourceTypeSet = true;}

        if(not SourceTypeSet)
        {
            string message = "No valid type is provided for the heat source " + counter + " --> " + SourceType;
            Info::WriteExit(message,thisclassname,"ReadInput()");
            exit(1);
        }

        switch(locHeatSource.Type)
        {
            case HeatSourceTypes::Phase:
            {
                locHeatSource.PhaseIndex = UserInterface::ReadParameterI(inp, moduleLocation, string("PhaseIndex_") + counter, true, 0);
                break;
            }
            case HeatSourceTypes::Ellipsoidal:
            {
                locHeatSource.Position = UserInterface::ReadParameterV3(inp, moduleLocation, string("Position_") + counter, true, {0.0,0.0,0.0});
                locHeatSource.Size     = UserInterface::ReadParameterV3(inp, moduleLocation, string("Radii_"   ) + counter, true, {0.0,0.0,0.0});
                break;
            }
            case HeatSourceTypes::Rectangular:
            {
                locHeatSource.Position = UserInterface::ReadParameterV3(inp, moduleLocation, string("Position_") + counter, true, {0.0,0.0,0.0});
                locHeatSource.Size     = UserInterface::ReadParameterV3(inp, moduleLocation, string("Sizes_"   ) + counter, true, {0.0,0.0,0.0});
                break;
            }
            default:
            {
                break;
            }
        }

        bool SourceONtriggerSet = false;
        string SourceONtrigger = UserInterface::ReadParameterK(inp, moduleLocation, string("TriggerON_") + counter, true, "USER");
        if(SourceONtrigger == "USER"            ) {locHeatSource.TriggerON = EventTriggers::User;             SourceONtriggerSet = true;}
        if(SourceONtrigger == "TMAX"            ) {locHeatSource.TriggerON = EventTriggers::Tmax;             SourceONtriggerSet = true;}
        if(SourceONtrigger == "TMIN"            ) {locHeatSource.TriggerON = EventTriggers::Tmin;             SourceONtriggerSet = true;}
        if(SourceONtrigger == "TIME"            ) {locHeatSource.TriggerON = EventTriggers::Time;             SourceONtriggerSet = true;}
        if(SourceONtrigger == "TIMESTEP"        ) {locHeatSource.TriggerON = EventTriggers::TimeStep;         SourceONtriggerSet = true;}
        if(SourceONtrigger == "PHASEFRACTIONMAX") {locHeatSource.TriggerON = EventTriggers::PhaseFractionMax; SourceONtriggerSet = true;}
        if(SourceONtrigger == "PHASEFRACTIONMIN") {locHeatSource.TriggerON = EventTriggers::PhaseFractionMin; SourceONtriggerSet = true;}

        if(not SourceONtriggerSet)
        {
            string message = "No valid ON trigger is provided for the heat source " + counter + " --> " + SourceONtrigger;
            Info::WriteExit(message,thisclassname,"ReadInput()");
            exit(1);
        }

        if(locHeatSource.TriggerON != EventTriggers::User)
        {
            locHeatSource.ONtriggerValue = UserInterface::ReadParameterD(inp, moduleLocation, string("TriggerONvalue_") + counter, true, 0);
        }
        if(locHeatSource.TriggerON == EventTriggers::PhaseFractionMax or locHeatSource.TriggerON == EventTriggers::PhaseFractionMin)
        {
            locHeatSource.PhaseIndexON = UserInterface::ReadParameterI(inp, moduleLocation, string("PhaseIndexON_") + counter, true, 0);
        }

        bool SourceOFFtriggerSet = false;
        string SourceOFFtrigger = UserInterface::ReadParameterK(inp, moduleLocation, string("TriggerOFF_") + counter, false, SourceONtrigger);
        if(SourceOFFtrigger == "USER"            ) {locHeatSource.TriggerOFF = EventTriggers::User;             SourceOFFtriggerSet = true;}
        if(SourceOFFtrigger == "TMAX"            ) {locHeatSource.TriggerOFF = EventTriggers::Tmax;             SourceOFFtriggerSet = true;}
        if(SourceOFFtrigger == "TMIN"            ) {locHeatSource.TriggerOFF = EventTriggers::Tmin;             SourceOFFtriggerSet = true;}
        if(SourceOFFtrigger == "TIME"            ) {locHeatSource.TriggerOFF = EventTriggers::Time;             SourceOFFtriggerSet = true;}
        if(SourceOFFtrigger == "TIMESTEP"        ) {locHeatSource.TriggerOFF = EventTriggers::TimeStep;         SourceOFFtriggerSet = true;}
        if(SourceOFFtrigger == "PHASEFRACTIONMAX") {locHeatSource.TriggerOFF = EventTriggers::PhaseFractionMax; SourceOFFtriggerSet = true;}
        if(SourceOFFtrigger == "PHASEFRACTIONMIN") {locHeatSource.TriggerOFF = EventTriggers::PhaseFractionMin; SourceOFFtriggerSet = true;}

        if(not SourceOFFtriggerSet)
        {
            string message = "No valid OFF trigger is provided for the heat source " + counter + " --> " + SourceOFFtrigger;
            Info::WriteExit(message,thisclassname,"ReadInput()");
            exit(1);
        }

        if(locHeatSource.TriggerOFF != EventTriggers::User)
        {
            locHeatSource.OFFtriggerValue = UserInterface::ReadParameterD(inp, moduleLocation, string("TriggerOFFvalue_") + counter, true, 0);
        }

        if(locHeatSource.TriggerOFF == EventTriggers::PhaseFractionMax or locHeatSource.TriggerOFF == EventTriggers::PhaseFractionMin)
        {
            locHeatSource.PhaseIndexOFF = UserInterface::ReadParameterI(inp, moduleLocation, string("PhaseIndexOFF_") + counter, true, 0);
        }

        locHeatSource.Repeat = UserInterface::ReadParameterI(inp, moduleLocation, string("Repeat_") + counter, false, 1);
        if(locHeatSource.Repeat <= 0)
        {
            string message = "Zero or negative \"Repeat\" parameter for heat source " + counter + " --> " + to_string(locHeatSource.Repeat);
            Info::WriteExit(message,thisclassname,"ReadInput()");
            exit(1);
        }

        Sources.push_back(locHeatSource);
        idx++;
    }

    Info::WriteLine();
    Info::WriteBlankLine();
}

void HeatSources::Apply(PhaseField& Phase, Temperature& Tx, HeatDiffusion& HD)
{
    for(auto source = Sources.begin(); source != Sources.end(); source++)
    if(source->Active)
    {
        switch(source->Type)
        {
            case HeatSourceTypes::BCX0:
            {
                if(Tx.ExtensionX0.isActive())
                {
                    Tx.ExtensionX0.Qdot = source->Value;
                }
                else if(Tx.OffsetX == 0)
                {
                    int x = 0;
                    for(int y = 0; y < HD.Ny; y++)
                    for(int z = 0; z < HD.Nz; z++)
                    {
                        HD.Qdot(x,y,z) += source->Value/HD.dx;
                    }
                }
                break;
            }
            case HeatSourceTypes::BCXN:
            {
                if(Tx.ExtensionXN.isActive())
                {
                    Tx.ExtensionXN.Qdot = source->Value/HD.dx;
                }
                else if(Tx.OffsetX + Tx.Nx == Tx.TotalNx)
                {
                    int x = HD.Nx-1;
                    for(int y = 0; y < HD.Ny; y++)
                    for(int z = 0; z < HD.Nz; z++)
                    {
                        HD.Qdot(x,y,z) += source->Value/HD.dx;
                    }
                }
                break;
            }
            case HeatSourceTypes::BCY0:
            {
                if(Tx.ExtensionY0.isActive())
                {
                    Tx.ExtensionY0.Qdot = source->Value/HD.dx;
                }
                else if(Tx.OffsetY == 0)
                {
                    int y = 0;
                    for(int x = 0; x < HD.Nx; x++)
                    for(int z = 0; z < HD.Nz; z++)
                    {
                        HD.Qdot(x,y,z) += source->Value/HD.dx;
                    }
                }
                break;
            }
            case HeatSourceTypes::BCYN:
            {
                if(Tx.ExtensionYN.isActive())
                {
                    Tx.ExtensionYN.Qdot = source->Value/HD.dx;
                }
                else if(Tx.OffsetY + Tx.Ny == Tx.TotalNy)
                {
                    int y = HD.Ny-1;
                    for(int x = 0; x < HD.Nx; x++)
                    for(int z = 0; z < HD.Nz; z++)
                    {
                        HD.Qdot(x,y,z) += source->Value/HD.dx;
                    }
                }
                break;
            }
            case HeatSourceTypes::BCZ0:
            {
                if(Tx.ExtensionZ0.isActive())
                {
                    Tx.ExtensionZ0.Qdot = source->Value/HD.dx;
                }
                else if(Tx.OffsetZ == 0)
                {
                    int z = 0;
                    for(int x = 0; x < HD.Nx; x++)
                    for(int y = 0; y < HD.Ny; y++)
                    {
                        HD.Qdot(x,y,z) += source->Value/HD.dx;
                    }
                }
                break;
            }
            case HeatSourceTypes::BCZN:
            {
                if(Tx.ExtensionZN.isActive())
                {
                    Tx.ExtensionZN.Qdot = source->Value/HD.dx;
                }
                else if(Tx.OffsetZ + Tx.Nz == Tx.TotalNz)
                {
                    int z = HD.Nz-1;
                    for(int x = 0; x < HD.Nx; x++)
                    for(int y = 0; y < HD.Ny; y++)
                    {
                        HD.Qdot(x,y,z) += source->Value/HD.dx;
                    }
                }
                break;
            }
            default:
            {
                break;
            }
        }
    }
}

void HeatSources::Activate(PhaseField& Phase, Temperature& Tx, RunTimeControl& RTC)
{
    int counter = 0;
    for(auto source = Sources.begin(); source != Sources.end(); source++)
    {
        if(source->Active)
        {
            switch(source->TriggerOFF)
            {
                case EventTriggers::Tmin:
                {
                    if(Tx.Tmin < source->OFFtriggerValue)
                    {
                        source->Active = false;
                    }
                    break;
                }
                case EventTriggers::Tmax:
                {
                    if(Tx.Tmax > source->OFFtriggerValue)
                    {
                        source->Active = false;
                    }
                    break;
                }

                case EventTriggers::TimeStep:
                {
                    if(RTC.tStep == source->OFFtriggerValue)
                    {
                        source->Active = false;
                    }
                    break;
                }
                case EventTriggers::Time:
                {
                    if(RTC.SimulationTime > source->OFFtriggerValue - RTC.dt and
                       RTC.SimulationTime < source->OFFtriggerValue + RTC.dt)
                    {
                        source->Active = false;
                    }
                    break;
                }
                case EventTriggers::PhaseFractionMax:
                {
                    if(Phase.FractionsTotal[source->PhaseIndexOFF] > source->OFFtriggerValue)
                    {
                        source->Active = false;
                    }
                    break;
                }
                case EventTriggers::PhaseFractionMin:
                {
                    if(Phase.FractionsTotal[source->PhaseIndexOFF] < source->OFFtriggerValue)
                    {
                        source->Active = false;
                    }
                    break;
                }
                case EventTriggers::User:
                default:
                {
                    break;
                }
            }
            if(!source->Active)
            {
                std::stringstream message;
                message << " Activate Heat Sources- heat source " << counter << " is deactivated.";
                Info::WriteStandard(thisclassname, message.str());
            }
        }
        else if (source->Repeat != 0)
        {
            switch(source->TriggerON)
            {
                case EventTriggers::Tmin:
                {
                    if(Tx.Tmin < source->ONtriggerValue)
                    {
                        source->Active = true;
                        source->Repeat -= 1;
                    }
                    break;
                }
                case EventTriggers::Tmax:
                {
                    if(Tx.Tmax > source->ONtriggerValue)
                    {
                        source->Active = true;
                        source->Repeat -= 1;
                    }
                    break;
                }

                case EventTriggers::TimeStep:
                {
                    if(RTC.tStep == source->ONtriggerValue)
                    {
                        source->Active = true;
                        source->Repeat -= 1;
                    }
                    break;
                }
                case EventTriggers::Time:
                {
                    if(RTC.SimulationTime > source->ONtriggerValue - RTC.dt and
                       RTC.SimulationTime < source->ONtriggerValue + RTC.dt)
                    {
                        source->Active = true;
                        source->Repeat -= 1;
                    }
                    break;
                }
                case EventTriggers::PhaseFractionMax:
                {
                    if(Phase.FractionsTotal[source->PhaseIndexON] > source->ONtriggerValue)
                    {
                        source->Active = true;
                        source->Repeat -= 1;
                    }
                    break;
                }
                case EventTriggers::PhaseFractionMin:
                {
                    if(Phase.FractionsTotal[source->PhaseIndexON] < source->ONtriggerValue)
                    {
                        source->Active = true;
                        source->Repeat -= 1;
                    }
                    break;
                }
                case EventTriggers::User:
                default:
                {
                    break;
                }
            }
            if(source->Active)
            {
                std::stringstream message;
                message << " Activate Heat Sources- heat source " << counter << " is activated.";
                Info::WriteStandard(thisclassname, message.str());
            }
        }
        counter++;
    }
}

}// namespace openphase
