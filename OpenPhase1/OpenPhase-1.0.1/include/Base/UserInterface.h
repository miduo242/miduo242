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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev; Philipp Engels;
 *   Johannes Goerler
 *
 */

#ifndef USERINTERFACE_H
#define USERINTERFACE_H

#include "Base/Includes.h"

namespace openphase
{
class OP_EXPORTS UserInterface
{
 public:
    static std::string MakeFileName(std::string Directory,
                                      std::string NameBase,
                                      int Index,
                                      std::string FileExtension);               /// Creates full filename string using the location directory, name base, running index and file extension

    // Methods to read input parameters from OpenPhase input files.
    // They search for $KEY in the entire file using the following syntax:
    // $KEY    commment    :   value

    static bool ParameterPresent(std::stringstream& Inp,                             /// Checks if Key is present in the specified module, returns true if parameter is found and false otherwise
                                 const int location, std::string Key);

    static int FindParameter(std::stringstream& Inp,                                 /// Checks if Key is present in the specified module and returns location right before the key
                             const int location, std::string Key);

    static int FindParameterLocation(std::stringstream& Inp,                         /// Checks if Key is present in the specified module and returns location right after the key
                             const int location, std::string Key);

    static int FindModuleLocation(std::stringstream& Inp,                            /// Returns location of the module
                                const std::string module);

    static double ReadParameterD(std::stringstream& Inp,                             /// Read double precision floating point parameter value
        int currentLocation,
        const std::string Key,
        const bool mandatory = true,
        const double defaultval = 0.0);

    static dVector3 ReadParameterV3(std::stringstream& Inp,                          /// Reads dVector3 vector
        int currentLocation,
        std::string Key,
        const bool mandatory = true,
        const dVector3 defaultval = {0,0,0});

    static dMatrix3x3 ReadParameterM3x3(std::stringstream& Inp,                      /// Reads dMatrix3x3 tensor
            int currentLocation,
            const std::string Key,
            const bool mandatory = true,
            const dMatrix3x3 defaultval = dMatrix3x3::ZeroTensor());

    static dMatrix6x6 ReadParameterT6(std::stringstream& Inp,                        /// Read dMatrix6x6 tensor
        int currentLocation,
        const std::string Key,
        const bool mandatory = true,
        const dMatrix6x6 defaultval = dMatrix6x6::UnitTensor());

    static dMatrix6x6 ReadParameterM6x6(std::stringstream& Inp,                      /// Reads dMatrix6x6 tensor
                int currentLocation,
                const std::string Key,
                const bool mandatory = true,
                const dMatrix6x6 defaultval = dMatrix6x6::ZeroTensor());

    static int ReadParameterI(std::stringstream& Inp,                                /// Reads integer parameter value
        int currentLocation,
        std::string Key,
        const bool mandatory = true,
        int const defaultval = 0);

    static std::string ReadParameterS(std::stringstream& Inp,                        /// Reads string parameter value
        int currentLocation,
        const std::string Key,
        const bool mandatory = true,
        const std::string defaultval = "NN");

    static std::string ReadParameterK(std::stringstream& Inp,                        /// Reads keyword (removes spaces and converts to upper case)
            int currentLocation,
            const std::string Key,
            const bool mandatory = true,
            const std::string defaultval = "NN");

    static std::vector<std::string> ReadParameterVS(std::stringstream& Inp,          /// Reads list of strings as a vector separated by any punctuation marks
        int currentLocation,
        std::string Key,
        const bool mandatory = true,
        const std::vector<std::string> defaultval = {""});

    static bool ReadParameterB(std::stringstream& Inp,                               /// Read boolean parameter value
        int currentLocation,
        std::string Key,
        const bool mandatory = true,
        const bool defaultval = false);

    static std::string ReadParameterF(std::stringstream& Inp,                        /// Reads filename string (assumes no spaces in the filename)
        int currentLocation,
        const std::string Key,
        const bool mandatory = true,
        const std::string defaultval = "NN");

    static std::string ReadParameterFW(std::stringstream& Inp,                        /// Reads filename string (assumes no spaces in the filename)
        int currentLocation,
        const std::string Key,
        const bool mandatory = true,
        const std::string defaultval = "NN");

    static char ReadParameterC(std::stringstream& Inp,                               /// Reads char parameter value
        int currentLocation,
        const std::string Key,
        const bool mandatory = true,
        const char defaultval = 'X');
};
}// namespace openphase
#endif
