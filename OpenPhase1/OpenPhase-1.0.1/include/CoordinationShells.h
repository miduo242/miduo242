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
 *   Main contributors :   Oleg Shchyglo;
 *
 */

#ifndef COORDINATIONSHELLS_H
#define COORDINATIONSHELLS_H

#include "Base/Includes.h"

namespace openphase
{

class Settings;

class CoordinationShells                                                        /// Stores the coordination shells for grains intialization algorithm.
{
 public:

    CoordinationShells(){};
    CoordinationShells(const Settings& locSettings)
    {
        Initialize(locSettings);
    }
    void Initialize(const Settings& locSettings);                               /// Initializes storages, sets internal variables.
    std::vector<iVector3> operator[](size_t n);                                 /// Returns a list of points belonging to a given coordination shell
    size_t Nshells()                                                            /// Returns the maximum number of generated coordination shells
    {
        return CShells.size();
    }

    static std::vector<dVector3> findPermutations(std::vector<dVector3>& locFacets, size_t n);

 protected:
 private:
    std::vector<iVector3> CShells;
};

} // namespace openphase

#endif
