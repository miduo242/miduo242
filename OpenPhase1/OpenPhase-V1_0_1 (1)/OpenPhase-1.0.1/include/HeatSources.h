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

#ifndef HEATSOURCES_H
#define HEATSOURCES_H

#include "Base/Includes.h"

namespace openphase
{
class Settings;
class HeatDiffusion;
class Temperature;
class PhaseField;
class RunTimeControl;

enum class HeatSourceTypes                                                      ///< Types of available heat sources
{
    BCX0,                                                                       ///< Heat source as boundary condition for lower X boundary
    BCXN,                                                                       ///< Heat source as boundary condition for upper X boundary
    BCY0,                                                                       ///< Heat source as boundary condition for lower Y boundary
    BCYN,                                                                       ///< Heat source as boundary condition for upper Y boundary
    BCZ0,                                                                       ///< Heat source as boundary condition for lower Z boundary
    BCZN,                                                                       ///< Heat source as boundary condition for upper Z boundary
    Phase,                                                                      ///< Selected phase will act as a heat source
    Ellipsoidal,                                                                ///< Ellipsoidal shape heat source. Position and size specified by the user
    Rectangular,                                                                ///< Rectangular shape heat source. Position and size specified by the user
};

struct HeatSourceStructure
{
    HeatSourceTypes Type;                                                       ///< Heat source type
    EventTriggers TriggerON;                                                    ///< Heat source ON trigger condition selector
    EventTriggers TriggerOFF;                                                   ///< Heat source OFF trigger condition selector
    bool Active;                                                                ///< Indicate if the current heat source is active (true) or inactive (false)

    double Value;                                                               ///< Heat flux value [W/m^2]
    double ONtriggerValue;                                                      ///< Trigger parameter value to turn the heat source ON
    double OFFtriggerValue;                                                     ///< Trigger parameter value to turn the heat source OFF
    size_t PhaseIndex;                                                          ///< Phase index for the phase sensitive trigger conditions
    size_t PhaseIndexON;                                                        ///< Phase index for the phase fraction evaluation for ON trigger
    size_t PhaseIndexOFF;                                                       ///< Phase index for the phase fraction evaluation for OFF trigger
    dVector3 Position;                                                          ///< Position of the center of the heat source for ellipsoidal or rectangular shapes [grid coordinates]
    dVector3 Size;                                                              ///< Heat source size for ellipsoidal (three orthogonal radii) and rectangular (three orthogonal side lengths) shape sources [grid coordinates]
    int Repeat;                                                                 ///< How many times to repeat the ON/OFF sequence (-1 for infinite)
};

class OP_EXPORTS HeatSources : public OPObject                                             ///< Heat sources class
{
 public:
    HeatSources(){};                                                            ///< Constructor
    HeatSources(Settings& locSettings,
                std::string InputFileName = DefaultInputFileName)               ///< Initializes storages, sets internal variables.
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }
    void Initialize(Settings& locSettings) override;                            ///< Allocates memory, initializes the settings);
    void ReadInput(const std::string FileName) override;                        ///< Reads input parameters
	void ReadInput(std::stringstream& FileName);                        ///< Reads input parameters

    void Apply(PhaseField& Phase, Temperature& Tx, HeatDiffusion& HD);          ///< Applies active heat sources
    void Activate(PhaseField& Phase, Temperature& Tx, RunTimeControl& RTC);     ///< Activates heat sources based on selected trigger conditions

    size_t Add(HeatSourceStructure& myHeatSource)                               ///< Adds custom heat source to the list of sources
    {
        Sources.push_back(myHeatSource);
        return Sources.size() - 1;
    }
    std::vector<HeatSourceStructure> Sources;                                   ///< List of active heat sources
    bool ConsiderHeatSources;                                                   ///< Indicated if there are active heat sources
    size_t Nphases;                                                             ///< Number of phases

 protected:
 private:
};

} // namespace openphase
#endif
