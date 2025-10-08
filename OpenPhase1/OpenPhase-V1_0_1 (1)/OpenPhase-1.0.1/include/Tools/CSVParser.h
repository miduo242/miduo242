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
 *   File created :   2015
 *   Main contributors :   Rajendran Mohan Kumar, Philipp Engels
 *
 */

#ifndef CSV_H
#define CSV_H

#include "Info.h"
#include "Base/Includes.h"

namespace openphase
{

class CSVParser
{
 public:
    static void WriteHeader(const std::string fileName,
        const std::vector<std::string> headerArray, const std::string seperator = ",");
    static void ClearContent(const std::string fileName);
    static void WriteData(const std::string fileName,
        const std::vector<int> dataArray, const std::string seperator = ",");
    static void WriteData(const std::string fileName,
        const std::vector<double> dataArray, const std::string seperator = ",");
    static void readFile(std::string& fileName,
        std::vector<std::vector<double>>& dataArray,
        std::vector<std::string>& headerArray, const std::string seperator = ",");
    static void readFile(std::string& fileName,
        std::vector<std::vector<std::string>>& dataArray,
        std::vector<std::string>& headerArray, std::string seperator = ",");

 protected:
 private:
};

} // namespace openphase
#endif //CSV_H
