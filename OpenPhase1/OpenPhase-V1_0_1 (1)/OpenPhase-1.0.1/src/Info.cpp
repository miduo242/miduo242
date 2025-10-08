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
 *   Main contributors :   Philipp Engels; Raphael Schiedung
 *
 */

#include "Info.h"
#include "RunTimeControl.h"
#include "Base/dVector6.h"
#include "Base/dMatrix6x6.h"

namespace openphase
{

using namespace std;

VerbosityLevels Info::OutputVerbosity = VerbosityLevels::Normal;

void Info::WriteTimeStep(const RunTimeControl& RTC, const string Message)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        time_t rawtime;
        time(&rawtime);

        string sTime;
        string thetime(ctime( &rawtime ));
        sTime = thetime.substr(0, thetime.length() -1);

        WriteLine("=");
        cout << setfill(' ') << setw(ColumnWidth)  << left <<
                "Time step" << ": " << to_string(RTC.tStep) + "/" + to_string(RTC.nSteps) << "\n";
        cout << setfill(' ') << setw(ColumnWidth)  << left <<
                "Simulation time" << ": " << RTC.SimulationTime << "\n";
        cout << setfill(' ') << setw(ColumnWidth)  << left <<
                "Wall clock time" << ": " << sTime << "\n";
        if (Message != "")
        {
            WriteLine("-");
            cout << Message;
        }
        WriteLine("=");
    }
}

void Info::WriteTimeStep(const int tStep, const int nSteps, const string Message)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        time_t rawtime;
        time(&rawtime);

        string sTime;
        string thetime(ctime( &rawtime ));
        sTime = thetime.substr(0, thetime.length() -1);

        WriteLine("_");
        cout << setfill(' ') << setw(ColumnWidth)  << left <<
                "TimeStep" << " " << to_string(tStep) + "/" + to_string(nSteps) << "\n";
        cout << setfill(' ') << setw(ColumnWidth)  << left <<
                "Time" << " " << sTime << "\n";
        if (Message != "")
        {
            cout << Message << "\n";
        }
        WriteLine("_");
    }
}

void Info::WriteTimeStep(const int tScreenWrite, const int tStep,
        const int nSteps, const string Message)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        if (!(tStep%tScreenWrite))
        {
            WriteTimeStep(tStep, nSteps, Message);
        }
    }
}

void Info::WriteLine(const string LineType)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        cout << setfill(LineType[0]) << setw(LineLength) << "" << "\n";
    }
}

void Info::WriteLineInsert(const string Insert, const string LineType)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    {
        cout << LineType[0] << LineType[0] << "< " << Insert << " >" << setfill(LineType[0]) <<
                setw(std::max<int>(0,LineLength-Insert.size()-6)) << "" << "\n";
    }
}

void Info::WriteBlankLine(void)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        cout << "\n";
    }
}

void Info::WriteSimple(const string Message)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        cout << Message << "\n";
    }
}

void Info::WriteStandard(const string Left, const bool Right)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0 )
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        std::cout << setfill(' ') << setw(ColumnWidth)  << left << Left << ":  ";
        std::cout << std::boolalpha;
        std::cout << Right << "\n";
    }
}

void Info::WriteStandard(const string Left, const string Right, const int precision, floatfield notation)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        std::cout << setfill(' ') << setw(ColumnWidth)  << left << Left << ":  ";
        switch(notation)
        {
            case floatfield::Scientific:
            {
                std::cout << std::scientific;
                break;
            }
            case floatfield::Fixed:
            {
                std::cout << std::fixed;
                break;
            }
        }
        std::cout << Right << "\n";
    }
}

void Info::WriteStandard(const std::string Left, const dVector3 vec, const int precision, floatfield notation)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        std::cout << setfill(' ') << setw(ColumnWidth)  << left << Left << ": ";
        std::cout << setprecision(precision);
        switch(notation)
        {
            case floatfield::Scientific:
            {
                std::cout << std::scientific;
                break;
            }
            case floatfield::Fixed:
            {
                std::cout << std::fixed;
                break;
            }
        }

        std::cout << "|"
                  << sgn(vec[0]) << vec[0] << " "
                  << sgn(vec[1]) << vec[1] << " "
                  << sgn(vec[2]) << vec[2] << " |\n";
    }
}

void Info::WriteStandard(const std::string Left, const dVector6 vec, const int precision, floatfield notation)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        std::cout << setfill(' ') << setw(ColumnWidth)  << left << Left << ": ";

        std::cout << setprecision(precision);
        switch(notation)
        {
            case floatfield::Scientific:
            {
                std::cout << std::scientific;
                break;
            }
            case floatfield::Fixed:
            {
                std::cout << std::fixed;
                break;
            }
        }

        std::cout << "|"
                  << sgn(vec[0]) << vec[0] << " "
                  << sgn(vec[1]) << vec[1] << " "
                  << sgn(vec[2]) << vec[2] << " "
                  << sgn(vec[3]) << vec[3] << " "
                  << sgn(vec[4]) << vec[4] << " "
                  << sgn(vec[5]) << vec[5] << " |\n";
    }
}

void Info::WriteStandard(const std::string Left, const dMatrix3x3 mat, const int precision, floatfield notation)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        std::cout << setfill(' ') << setw(ColumnWidth)  << left << Left << ":\n";
        //std::cout << setprecision(precision);
        //std::cout << std::scientific;

        std::cout << mat.print() << "\n";
    }
}

void Info::WriteStandard(const std::string Left, const dMatrix6x6 mat, const int precision, floatfield notation)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        std::cout << setfill(' ') << setw(ColumnWidth)  << left << Left << ":\n";
        //std::cout << setprecision(precision);
        //std::cout << std::scientific;

        std::cout << mat.print() << "\n";
    }
}

void Info::WriteStandard(const std::string Left, const vStress Stresses, const int precision, floatfield notation)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        std::cout << setfill(' ') << setw(ColumnWidth)  << left << Left << ": ";

        std::cout << setprecision(precision);
        switch(notation)
        {
            case floatfield::Scientific:
            {
                std::cout << std::scientific;
                break;
            }
            case floatfield::Fixed:
            {
                std::cout << std::fixed;
                break;
            }
        }

        std::cout << "\n|"
                << sgn(Stresses[0]) << Stresses[0] << " "
                << sgn(Stresses[1]) << Stresses[1] << " "
                << sgn(Stresses[2]) << Stresses[2] << " "
                << sgn(Stresses[3]) << Stresses[3] << " "
                << sgn(Stresses[4]) << Stresses[4] << " "
                << sgn(Stresses[5]) << Stresses[5] << " |\n";
    }
}

void Info::WriteStandard(const std::string Left, const vStrain Strain, const int precision, floatfield notation)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        std::cout << setfill(' ') << setw(ColumnWidth)  << left << Left << ": ";

        std::cout << setprecision(precision);
        switch(notation)
        {
            case floatfield::Scientific:
            {
                std::cout << std::scientific;
                break;
            }
            case floatfield::Fixed:
            {
                std::cout << std::fixed;
                break;
            }
        }

        std::cout << "\n|"
                << sgn(Strain[0]) << Strain[0] << " "
                << sgn(Strain[1]) << Strain[1] << " "
                << sgn(Strain[2]) << Strain[2] << " "
                << sgn(Strain[3]) << Strain[3] << " "
                << sgn(Strain[4]) << Strain[4] << " "
                << sgn(Strain[5]) << Strain[5] << " |\n";
    }
}

void Info::WriteStandardNarrow(const string Left, const string Right)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        cout << setfill(' ') << setw(20)  << left << Left << ": " << Right << "\n";
    }
}

string Info::GetStandard(const string Left, const string Right)
{
    stringstream ss;
    ss << setfill(' ') << setw(ColumnWidth)  << left << Left << ": " << setprecision(8) << Right << "\n";
    string returnStandard = ss.str();

    return returnStandard;
}

string Info::GetStandardNarrow(const string Left, const string Right)
{
    stringstream ss;
    ss << setfill(' ') << setw(20)  << left << Left << ": " << setprecision(8) << Right << "\n";
    string returnStandard = ss.str();

    return returnStandard;
}

void Info::WriteCoordinate(const int x, const int y, const int z, const double dx)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        cout << "Point:      (" << x << ", " << y << ", " << z << ")" << "\n";
        cout << "Coordinate: (" << x*dx << ", " << y*dx << ", " << z*dx << ")" << "\n";
    }
}

void Info::WriteWarning(const string Message, const string Instance,
                        const string Method)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        string thisInstance = Instance;
        WriteLine("~");
        if (Method != "")
        {
            thisInstance += "::" + Method;
        }
        cout << setfill(' ') << setw(10)  << left << "Warning: "  << thisInstance << "\n"
                                                  << "          " << Message  << "\n";
        WriteLine("~");
    }
}

void Info::WriteExit(const string Message, const string Instance, const string Method)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    {
        string thisInstance = Instance;
        time_t rawtime;
        time(&rawtime);

        string sTime;
        string thetime(ctime( &rawtime ));
        sTime = thetime.substr(0, thetime.length() -1);

        cerr << "\n";
        WriteLine("*");
        if (Method != "")
        {
            thisInstance += "::" + Method;
        }
        cerr << setfill(' ') << setw(10)  << left << "Calculation terminated!" << "\n"
                             << setw(10)  << left << "Instance:" << thisInstance << "\n"
                             << setw(10)  << left << "Reason: " << Message << "\n"
                             << setw(10)  << left << "Time: " << sTime << "\n";
        WriteLine("*");
        cerr << "\n";
    }
}

void Info::WriteStartScreen(void)
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    if(OutputVerbosity != VerbosityLevels::Silent)
    {
        cout << "\n";
        WriteLine(">");
        cout << "  OpenPhase software library 1.0.0 \n\n"
             << "  Copyright (c) Ruhr-Universitaet Bochum, Universitaetsstrasse 150, D-44801 Bochum, Germany\n"
             << "            AND OpenPhase Solutions GmbH, Universitaetsstrasse 136, D-44799 Bochum, Germany.\n"
             << "\n"
             << "  This program is free software: you can redistribute it and/or modify\n"
             << "  it under the terms of the GNU General Public License as published by\n"
             << "  the Free Software Foundation, either version 3 of the License, or\n"
             << "  (at your option) any later version.\n"
             << "\n"
             << "  This program is distributed in the hope that it will be useful,\n"
             << "  but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
             << "  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
             << "  GNU General Public License for more details.\n"
             << "\n"
             << "  You should have received a copy of the GNU General Public License\n"
             << "  along with this program.  If not, see <http://www.gnu.org/licenses/>.\n";
;
        WriteLine("<");
        cout << "\n";
    }
}

void Info::PressEnterToContinue(void)
{
    cout << "Press ENTER to continue... " << flush;
    cin.ignore( numeric_limits <streamsize> ::max(), '\n' );
}
}// namespace openphase
