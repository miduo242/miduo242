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
 *   Main contributors :   Johannes Goerler; Oleg Shchyglo
 *
 */

#ifndef PLASTICFLOWNEUBERMETHODS_H
#define PLASTICFLOWNEUBERMETHODS_H

#include "Base/Includes.h"

namespace openphase
{
class ElasticProperties;

class PlasticFlowNeuberMethods
{
 public:
    static void getNeuberDataRO(vStress& StressIn, vStrain& StrainIn,
                                vStress& StressOut, vStrain& StrainOut,
                                double YoungModulus, double YieldStress,
                                double alpha);                                  ///< Neuber correction method based on Ramberg-Osgood equation with exponent n=3

    static vStrain getNeuberStrains(vStrain StrainIn);
    static vStress getNeuberStresses(vStress StressIn);

    static void writeNeuberStressVTK(const Settings& locSettings, ElasticProperties& EP, const int tStep);
    static void writeNeuberStrainVTK(const Settings& locSettings, ElasticProperties& EP, const int tStep);

 protected:
    static void WriteStressesVTKData(ElasticProperties& EP, std::stringstream& buffer);
    static void WriteStrainsVTKData(ElasticProperties& EP, std::stringstream& buffer);

 private:
};
}// namespace openphase

#endif

