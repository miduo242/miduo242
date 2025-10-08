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

#ifndef INFO_H
#define INFO_H

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>

//#include "Base/Includes.h"
#include "Base/Definitions.h"
#include "Base/vStrain.h"
#include "Base/vStress.h"

namespace openphase
{

class RunTimeControl;
class dMatrix3x3;
class dMatrix6x6;
class dVector3;
class dVector6;

enum class VerbosityLevels                                                      ///< Terminal output verbosity levels
{
    Silent,                                                                     ///< No terminal output except for exits.
    Normal,                                                                     ///< All outputs are enabled
    Debug                                                                       ///< All outputs plus debug outputs
};

enum class floatfield                                                           ///< Terminal output format for floating point numbers
{
    Fixed,                                                                      ///< Use fixed floating-point notation for output
    Scientific                                                                  ///< Use scientific floating-point notation for output
};

/// A Collection of method to format output to screen in a standardized format
class OP_EXPORTS Info
{
 public:

    static VerbosityLevels OutputVerbosity;                                     ///< Sets the level of simulation output verbosity

    static const size_t ColumnWidth = 40;
    static const size_t LineLength  = 80;

    /// Return white space for positive values to align positive ans negative
    /// numbers vertically
    template <typename T>
    static inline std::string sgn(const T value)
    {
         return (value >= 0)?" ":"";
    }

    /// Writes a Message, Blank Line, or Line of char to screen
    static void Write(const std::string message)
    {
        if (message.size() > 1)
        {
            WriteSimple(message);
        }
        else if (message.size() == 1)
        {
            WriteLine(message);
        }
        else
        {
            WriteBlankLine();
        }
    }

    /// Writes name and value to screen
    ///
    /// @param name of the value printed to screen
    /// @param value printed to screen
    /// @param precision of value printed to screen
    template <typename T>
    static void Write(const std::string name, const T value,
            const int precision=4)
    {
        WriteStandard(name, value, precision);
    }

    /// Writes name and value to screen and to line buffer for later file write
    ///
    /// @param line is a buffer which can be written to log file
    /// @param tStep current time step
    /// @param name of the value printed to screen
    /// @param value printed to screen
    /// @param precision of value printed to screen
    /// @param sep is separator char used in the log file
    template <typename T>
    static void WriteWithLog(std::array<std::stringstream,2> &line,
            const int tStep, const std::string name, const  T value,
            const int precision=4, const char sep = ',')
    {
        if (tStep == 0)
        {
            line[0] << name << sep;
        }
        line[1] << value << sep;
        WriteStandard(name, value, precision);
    }

    /// Writes line buffer to log file
    ///
    /// @param log is the log file where line will be written to
    /// @param line is a buffer which can be written to log file
    /// @param tStep current time step
    /// @param precision of value printed to screen
    /// @param sep is separator char used in the log file
    static void WriteLineToLogfile(std::fstream &log,
            const std::array<std::stringstream,2> &line,
            const int tStep, const char sep = ',')
    {
        if (tStep == 0)
        {
            log << tStep << sep << line[0].str() << std::endl;
        }
        log << tStep << sep << line[1].str() << std::endl;
    }

    // Standard output writing methods
    template <typename T>
    static void WriteStandard(const std::string Left, const T Right,
            [[maybe_unused]] const int precision = 6,
            [[maybe_unused]] floatfield notation = floatfield::Scientific)
    {
#ifdef MPI_PARALLEL
        if (MPI_RANK == 0)
#endif
        {
            if constexpr (std::is_floating_point<T>::value)
            {
                std::cout << std::setfill(' ') << std::setw(ColumnWidth)  << std::left << Left << ": ";
                std::cout << std::setprecision(precision);
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
                std::cout << sgn(Right) << Right << std::endl;
            }
            else if constexpr (std::is_signed<T>::value)
            {
                std::cout << std::setfill(' ') << std::setw(ColumnWidth)  << std::left << Left << ": ";
                std::cout << sgn(Right) << Right << std::endl;
            }
            else
            {
                std::cout << std::setfill(' ') << std::setw(ColumnWidth)  << std::left << Left << ": ";
                std::cout << " " << Right << std::endl;
            }
        }
    }
    static void WriteStandard(const std::string Left, const bool Right);
    static void WriteStandard(const std::string Left, const std::string Right,
            const int precision = 4,
            floatfield notation = floatfield::Scientific);
    static void WriteStandard(const std::string Left, const dVector3 vec,
            const int precision = 3,
            floatfield notation = floatfield::Scientific);
    static void WriteStandard(const std::string Left, const dVector6 vec,
            const int precision = 3,
            floatfield notation = floatfield::Scientific);
    static void WriteStandard(const std::string Left, const dMatrix3x3 mat,
            const int precision = 3,
            floatfield notation = floatfield::Scientific);
    static void WriteStandard(const std::string Left, const dMatrix6x6 mat,
            const int precision = 3,
            floatfield notation = floatfield::Scientific);
    static void WriteStandard(const std::string Left, const vStress Stresses,
            const int precision = 3,
            floatfield notation = floatfield::Scientific);
    static void WriteStandard(const std::string Left, const vStrain Strains,
            const int precision = 3,
            floatfield notation = floatfield::Scientific);

    static void WriteStresses(const std::string Left, const vStress Stresses,
            const int precision = 3,
            floatfield notation = floatfield::Scientific)
    {
        WriteStandard(Left, Stresses, precision, notation);
    }
    static void WriteStrains(const std::string Left, const vStrain Strains,
            const int precision = 3,
            floatfield notation = floatfield::Scientific)
    {
        WriteStandard(Left, Strains, precision, notation);
    }

    // Static console output methods
    static void WriteLine(const std::string LineType = "-");
    static void WriteLineInsert(const std::string Insert, const std::string LineType = "-");
    static void WriteBlankLine(void);
    static void WriteTimeStep(const RunTimeControl& RTC, const std::string Message = "");
    static void WriteTimeStep(const int tStep, const int nSteps, const std::string Message = "");
    static void WriteTimeStep(const int tScreenWrite, const int tStep, const int nSteps, const std::string Message = "");
    static void WriteSimple(const std::string Message);
    static void WriteStandardNarrow(const std::string Left, const std::string Right);
    static void WriteCoordinate(const int x, const int y, const int z, const double dx);
    static void WriteWarning(const std::string Message, const std::string Instance = "", const std::string Method = "");
    static void WriteExit(const std::string Message, const std::string Instance = "", const std::string Method = "");
    static void WriteStartScreen(void);
    static void PressEnterToContinue(void);
    static std::string GetStandard(const std::string Left, const std::string Right);
    static std::string GetStandardNarrow(const std::string Left, const std::string Right);

    template <typename T>
    static std::string to_string_with_precision(const T a_value, const int n = 6)
    {
        std::ostringstream out;
        out << std::setprecision(n) << a_value;
        return out.str();
    }
};
}// namespace openphase
#endif
