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
 *   File created :   2013
 *   Main contributors :   Oleg Shchyglo
 *
 */

#ifndef TIMEINFO_H
#define TIMEINFO_H

#include "Base/Includes.h"
#include "Info.h"
namespace openphase
{

class Settings;

class TimeInfo                                                                  ///< Execution time statistics collection and output
{
 public:

    ~TimeInfo();
    TimeInfo(){};
    TimeInfo(const Settings& locSettings, const std::string Name = "Execution Time Statistics");
    void Initialize(const Settings& locSettings, const std::string Name);
    void SetStart(void);
    void SetTimeStamp(const std::string Message);
    void SkipToHere(void);
    void Reset(void);
    void PrintWallClockSummary(void) const;
    void PrintCPUClockSummary(void) const
    {
        Info::WriteWarning("TimeInfo::PrintCPUClockSummary() not yet implemented.",
                "TimeInfo", "PrintCPUClockSummary()");
    };
    void PrintFullSummary(void) const
    {
        Info::WriteWarning("TimeInfo::PrintFullSummary() not yet implemented.",
                "TimeInfo", "PrintFullSummary()");
    };

 protected:
 private:
    std::string thisclassname;
    std::vector<clock_t> TimeCollector;
    std::vector<std::string> Messages;

    std::string TimerName;
};
}// namespace openphase
#endif
