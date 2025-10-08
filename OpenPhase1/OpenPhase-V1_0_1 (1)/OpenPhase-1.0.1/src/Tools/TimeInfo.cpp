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
 *   Main contributors :
 *
 */

#include "Tools/TimeInfo.h"
#include "Settings.h"
#include "Info.h"

namespace openphase
{

using namespace std;

TimeInfo::TimeInfo(const Settings& locSettings, const std::string Name)
{
    Initialize(locSettings, Name);
}

TimeInfo::~TimeInfo(void)
{
    Info::WriteStandard(thisclassname, "Exited normally");
}

void TimeInfo::Initialize(const Settings& locSettings, const string Name)
{
    thisclassname = "TimeInfo";
    TimerName = Name;
}
void TimeInfo::Reset(void)
{
    TimeCollector.clear();
    Messages.clear();
}

void TimeInfo::SetStart(void)
{
    Reset();
    TimeCollector.push_back(clock());
    Messages.push_back("Start");
}

void TimeInfo::SetTimeStamp(const string Message)
{
    TimeCollector.push_back(clock());
    Messages.push_back(Message);
}

void TimeInfo::SkipToHere(void)
{
    TimeCollector.push_back(clock());
    Messages.push_back("#");
}

void TimeInfo::PrintWallClockSummary(void) const
{
    Info::WriteLine("=");
    Info::WriteSimple(TimerName);
    Info::WriteLine("-");

    clock_t TotalTime = 0;
    for(unsigned int it = 1; it < TimeCollector.size(); it++)
    {
        if(Messages[it]!="#")
        {
            TotalTime += (TimeCollector[it] - TimeCollector[it-1]);
        }
    }
    double TotalConsumedTime = double(TotalTime) / double(CLOCKS_PER_SEC*omp_get_max_threads());

    for(unsigned int it = 1; it < TimeCollector.size(); it++)
    {
        if(Messages[it]!="#")
        {
            stringstream showtime;
            double SectionConsumedTime = double(TimeCollector[it] - TimeCollector[it-1]) / double(CLOCKS_PER_SEC*omp_get_max_threads());
            showtime << std::fixed << std::setprecision(2) << std::scientific << SectionConsumedTime << "  "
                     << std::fixed << std::setprecision(2) << ((SectionConsumedTime/TotalConsumedTime)*100.0) << " %";
            Info::WriteStandard(Messages[it], showtime.str());
        }
    }
    Info::WriteLine("-");
    Info::WriteStandard("Total", to_string(TotalConsumedTime));
    Info::WriteLine("=");
}

}// namespace openphase
