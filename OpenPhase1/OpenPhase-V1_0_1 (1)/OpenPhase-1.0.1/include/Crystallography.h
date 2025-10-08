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
 *   Main contributors :   Philipp Engels; Hesham Salama
 *
 */

#ifndef CRYSTALLOGRAPHY_H
#define CRYSTALLOGRAPHY_H

#include "Base/Includes.h"

namespace openphase
{

enum class CrystalSymmetry
{
    CUBIC,
    HEXAGONAL,
    TETRAGONAL,
    TRIGONAL,
    ORTHORHOMBIC,
    MONOCLINIC,
    TRICLINIC
};


class Settings;

class OP_EXPORTS Crystallography : public OPObject
{
 public:

    CrystalSymmetry CrystalSym;

    static const std::vector<double> SlipSystem_Cubic;
    static const std::vector<double> SlipSystem_BCC;
    static const std::vector<double> SlipSystem_FCC_XY;
    static const std::vector<double> SlipSystem_FCC_XZ;
    static const std::vector<double> SlipSystem_FCC_YZ;
    static const std::vector<double> SlipSystem_FCC;
    static const std::vector<double> SlipSystem_PartialFCC;
    static const std::vector<double> SlipSystem_HCPBasal;
    static const std::vector<double> SlipSystem_HCPPrismatic;
    static const std::vector<double> SlipSystem_HCPPyramidal;
    static const std::vector<double> SlipSystem_HCPPyramidal_1st;
    static const std::vector<double> SlipSystem_HCPPyramidal_2nd;
    static const std::vector<double> TwinSystem_FCC;
    static const std::vector<double> TwinSystem_HCPTensileTwin1;
    static const std::vector<double> TwinSystem_HCPTensileTwin2;
    static const std::vector<double> TwinSystem_HCPCompressionTwin1;
    static const std::vector<double> TwinSystem_HCPCompressionTwin2;
    Crystallography(){};
    Crystallography(Settings& locSettings, std::string InputFileName);
    void ReadInput(std::string InputFileName) override;                         /// Reads input parameters
    void Initialize(Settings& locSettings) override;                            /// Initializes storages, sets internal variables.

    Storage<dMatrix3x3> SymmetriesCubic;                                        /// Stores cubic symmetry operations
    Storage<dMatrix3x3> CrystalSymmetries;

    size_t nsym = 1;

    double Phi = 0.0;
    double Phi_min = 0.0;
    double Phi_max = 0.0;
    double Phi_max2 = 0.0;
    double Theta = 0.0;
    double Theta_min = 0.0;
    double Theta_max = 0.0;

    const double a = std::sqrt(3.0) / 2.0;

 protected:
 private:
};

} // namespace openphase

#endif
