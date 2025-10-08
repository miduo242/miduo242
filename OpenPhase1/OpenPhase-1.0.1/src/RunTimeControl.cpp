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
 *   File created :   2018
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */

#include "RunTimeControl.h"
#include "Base/UserInterface.h"
#include "Settings.h"
#include "Info.h"
#ifdef _WIN32
#include <filesystem>
#endif


namespace openphase
{

using namespace std;

void RunTimeControl::Initialize(Settings& locSettings)
{
    thisclassname = "RunTimeControl";
    thisobjectname = thisclassname;

    nSteps          = 0;
    dt              = 0.0;
    SimulationTime  = 0.0;
    Restart         = false;
    tStart          = 0;
    tStep           = 0;
    tRestartWrite   = 1;
    tFileWrite      = 1;
    tScreenWrite    = 1;
    nOMP            = 1;

    VTKDir = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;
    TextDir = locSettings.TextDir;

    CheckTime = mygettime();
    bStop = false;

    Info::WriteStandard(thisclassname, "Initialized");
}

void RunTimeControl::ReadInput(const string InputFileName)
{
    Info::WriteLineInsert("RunTimeControl input");
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
}

void RunTimeControl::ReadInput(stringstream& inp)
{
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    SimulationTitle    = UserInterface::ReadParameterF(inp, moduleLocation, string("SimTtl"));

    UnitsOfLength      = UserInterface::ReadParameterS(inp, moduleLocation, string("LUnits"), false, "m");
    UnitsOfMass        = UserInterface::ReadParameterS(inp, moduleLocation, string("MUnits"), false, "kg");
    UnitsOfTime        = UserInterface::ReadParameterS(inp, moduleLocation, string("TUnits"), false, "s");
    UnitsOfEnergy      = UserInterface::ReadParameterS(inp, moduleLocation, string("EUnits"), false, "J");

#ifndef SERIAL
    nOMP               = UserInterface::ReadParameterI(inp, moduleLocation, string("nOMP"));
    omp_set_num_threads(nOMP);
#endif

    nSteps             = UserInterface::ReadParameterI(inp, moduleLocation, string("nSteps"));
    dt                 = UserInterface::ReadParameterD(inp, moduleLocation, string("dt"));

    Restart            = UserInterface::ReadParameterB(inp, moduleLocation, string("Restrt"), false, false);
    if(Restart)
    {
        tStart         = UserInterface::ReadParameterI(inp, moduleLocation, string("tStart"));
    }
    else
    {
        tStart = 0;
    }
    tRestartWrite      = UserInterface::ReadParameterI(inp, moduleLocation, string("tRstrt"));
    tFileWrite         = UserInterface::ReadParameterI(inp, moduleLocation, string("FTime"));
    tScreenWrite       = UserInterface::ReadParameterI(inp, moduleLocation, string("STime"));
    StopTime           = UserInterface::ReadParameterD(inp, moduleLocation, string("StopTime"), false, 0.0);

    LogModeScreen      = UserInterface::ReadParameterB(inp, moduleLocation, string("LogScreen"),  false, false);
    LogModeVTK         = UserInterface::ReadParameterB(inp, moduleLocation, string("LogVTK"),     false, false);
    LogScreenFactor    = UserInterface::ReadParameterI(inp, moduleLocation, string("LogScreenF"), false, 1);
    LogVTKFactor       = UserInterface::ReadParameterI(inp, moduleLocation, string("LogVTKF"),    false, 1);


    Info::WriteLine();
    Info::WriteLine();

    if (tFileWrite < 1)
    {
        Info::WriteSimple(string("Illegal value for FTime => ") + to_string(tFileWrite) + string(" ! FTime is set to 1"));
        tFileWrite = 1;
    }
    if (tScreenWrite < 1)
    {
        Info::WriteSimple(string("Illegal value for STime => ") + to_string(tScreenWrite) + string(" ! STime is set to 1"));
        tScreenWrite = 1;
    }

    struct stat st;

    if(stat(VTKDir.c_str(),&st) != 0)
    {
#ifdef _WIN32
 //       ignore_result( system(string("mkdir \\" + VTKDir).c_str())+"\\");
        std::filesystem::create_directories(VTKDir);
#else
        ignore_result( system(string("mkdir -p " + VTKDir).c_str()));
#endif
        Info::WriteStandard("Created directory", VTKDir);
    }

    struct stat st2;
    if(stat(RawDataDir.c_str(),&st2) != 0)
    {
#ifdef _WIN32
        //ignore_result( system(string("mkdir \\" + RawDataDir).c_str())+"\\");
        std::filesystem::create_directories(RawDataDir);
#else
        ignore_result( system(string("mkdir -p " + RawDataDir).c_str()));
#endif
        Info::WriteStandard("Created directory", RawDataDir);
    }

    struct stat st3;
    if(stat(TextDir.c_str(),&st3) != 0)
    {
#ifdef _WIN32
        //ignore_result( system(string("mkdir \\" + TextDir).c_str())+"\\");
        std::filesystem::create_directories(TextDir);
#else
        ignore_result( system(string("mkdir -p " + TextDir).c_str()));
#endif
        Info::WriteStandard("Created directory", TextDir);
    }

    Info::WriteLine();
    Info::WriteBlankLine();
}

RunTimeControl& RunTimeControl::operator= (const RunTimeControl& rhs)
{
    if (this != &rhs) // protect against self-assignment
    {
        thisclassname = rhs.thisclassname;

        SimulationTitle = rhs.SimulationTitle;

        UnitsOfLength = rhs.UnitsOfLength;
        UnitsOfMass = rhs.UnitsOfMass;
        UnitsOfTime = rhs.UnitsOfTime;
        UnitsOfEnergy = rhs.UnitsOfEnergy;

        dt = rhs.dt;
        tStep = rhs.tStep;
        SimulationTime = rhs.SimulationTime;
        nSteps = rhs.nSteps;
        tFileWrite = rhs.tFileWrite;
        tScreenWrite = rhs.tScreenWrite;
        tStart = rhs.tStart;
        Restart = rhs.Restart;
        tRestartWrite = rhs.tRestartWrite;
        nOMP = rhs.nOMP;

        VTKDir = rhs.VTKDir;
        RawDataDir = rhs.RawDataDir;
        TextDir = rhs.TextDir;

        omp_set_num_threads(nOMP);
    }
    return *this;
}
void RunTimeControl::Write()
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    {
        std::ofstream outp(RawDataDir + thisobjectname + ".dat");
        outp << tStep << std::flush;
    }
}
bool RunTimeControl::RestartPossible()
{
    std::ifstream inp(RawDataDir + thisobjectname + ".dat");
    if (inp)
    {
        inp >> tStart;
        Restart = true;
    }
    else
    {
        Restart = false;
    }
    return Restart;
}

} //namespace openphase
