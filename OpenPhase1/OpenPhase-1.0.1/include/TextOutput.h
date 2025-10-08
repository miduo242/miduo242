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
 *   Main contributors :   Matthias Stratmann; Johannes Goerler
 *
 */

#ifndef TEXTOUTPUT_H
#define TEXTOUTPUT_H

#include "Base/Includes.h"

namespace openphase
{
class BoundaryConditions;
class ElasticProperties;
class GrainInfo;
class InterfaceField;
class PhaseField;
class PhaseField;
class Settings;
class Temperature;
class Composition;
class Nucleation;

class TextOutput : public OPObject                                              ///< Output of tabulated simulation data
{
 public:
    static void PhasePercent(Composition& Cx, PhaseField& Phi,
                               Settings& OPSettings, std::string filename,
                               double time);                                    ///< Write volume percent of thermodynamic phases
    static void GrainVolumes(PhaseField& Phi, std::string filename, double time);    ///< Write volume of each grain
    static void WriteValue(double Value,std::string filename, double time);
    static void WriteMultipleValues(std::vector<std::string> Names, std::vector<double> value, std::string filename, double time);
    static void LineConcentration(Composition& Cx, PhaseField& Phi,
                                  std::string filename, double timestep,
                                  std::string type, std::string axis, int x, int y,int z);///< Write total composition over a straight line in separate files
    static void AverageStress(ElasticProperties& EP, std::string filename,
                              double timeOrStrain);                             ///< Write average stress over time/strain
    static void AverageStrain(ElasticProperties& EP, std::string filename,
                              double timeOrStrain);                             ///< Write average strain over time
    static void AverageDMatrix3x3(Storage3D<dMatrix3x3, 0>& Matrix, Settings& OP,
                                  std::string filename, double timeOrStrain);
    static void AverageDouble(Storage3D<double,0>& value, Settings& OP,
                              std::string filename, double timeOrStrain);
    static void maxElasticRotation(ElasticProperties& EP, PhaseField& Phase,
                                   Orientations& OR, int tStep,
                                   std::string Filename);                       ///< Write max elastic rotations over time
    static void AverageTemp(Temperature& Tx, std::string filename, double time);///< Write average system temperature
    static bool FileExists(std::string filename);                               ///< Check availability of given filename

 private:
    static void LocalPhaseComposition(Composition& Cx,
                                      PhaseField& Phi, std::string filename,
                                      double time, int x, int y, int z);         ///<
};

} // namespace openphase
#endif
