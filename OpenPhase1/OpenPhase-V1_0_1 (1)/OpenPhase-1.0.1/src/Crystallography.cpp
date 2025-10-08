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

#include "Crystallography.h"
#include "Info.h"
#include "Settings.h"
#include "Base/UserInterface.h"

namespace openphase
{
using namespace std;

Crystallography::Crystallography(Settings& locSettings, const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void Crystallography::ReadInput(string InputFileName)
{
    Info::WriteLineInsert("Crystal Symmetry input");
    Info::WriteStandard("Source", InputFileName);

    fstream inpF(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inpF)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };
    stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();

    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    string Symmetry = UserInterface::ReadParameterK(inp, moduleLocation, "CrystalSymmetry", false, "CUBIC");

    // Symmetry matrices for crystal lattice. Taken from Kocks, Tome and Wenk
    // 'Texture and Anisotropy', Cambridge University Press 1998
    if (Symmetry == "CUBIC")
    {
        // Load cubic parameters (class 432)
        CrystalSym =  CrystalSymmetry::CUBIC;

        nsym = 24;
        Theta_min = 0 * (M_PI / 180);
        Theta_max = 45 * (M_PI / 180);
        Phi_min = 0 * (M_PI);
        Phi_max = acos(sqrt(1.0 / (2.0 + (pow(tan(Theta_max),2.0)))));

        CrystalSymmetries.Allocate(nsym);

        CrystalSymmetries[0].set(1, 0, 0, 0, 1, 0, 0, 0, 1);
        CrystalSymmetries[1].set(0, 0, 1, 1, 0, 0, 0, 1, 0);
        CrystalSymmetries[2].set(0, 1, 0, 0, 0, 1, 1, 0, 0);
        CrystalSymmetries[3].set(0,-1, 0, 0, 0, 1,-1, 0, 0);
        CrystalSymmetries[4].set(0,-1, 0, 0, 0,-1, 1, 0, 0);
        CrystalSymmetries[5].set(0, 1, 0, 0, 0,-1,-1, 0, 0);

        CrystalSymmetries[6].set( 0, 0,-1, 1, 0, 0, 0,-1, 0);
        CrystalSymmetries[7].set( 0, 0,-1,-1, 0, 0, 0, 1, 0);
        CrystalSymmetries[8].set( 0, 0, 1,-1, 0, 0, 0,-1, 0);
        CrystalSymmetries[9].set(-1, 0, 0, 0, 1, 0, 0, 0,-1);
        CrystalSymmetries[10].set(-1, 0, 0, 0,-1, 0, 0, 0, 1);
        CrystalSymmetries[11].set( 1, 0, 0, 0,-1, 0, 0, 0,-1);

        CrystalSymmetries[12].set( 0, 0,-1, 0,-1, 0,-1, 0, 0);
        CrystalSymmetries[13].set( 0, 0, 1, 0,-1, 0, 1, 0, 0);
        CrystalSymmetries[14].set( 0, 0, 1, 0, 1, 0,-1, 0, 0);
        CrystalSymmetries[15].set( 0, 0,-1, 0, 1, 0, 1, 0, 0);
        CrystalSymmetries[16].set(-1, 0, 0, 0, 0,-1, 0,-1, 0);
        CrystalSymmetries[17].set( 1, 0, 0, 0, 0,-1, 0, 1, 0);

        CrystalSymmetries[18].set( 1, 0, 0, 0, 0, 1, 0,-1, 0);
        CrystalSymmetries[19].set(-1, 0, 0, 0, 0, 1, 0, 1, 0);
        CrystalSymmetries[20].set( 0,-1, 0,-1, 0, 0, 0, 0,-1);
        CrystalSymmetries[21].set( 0, 1, 0,-1, 0, 0, 0, 0,-1);
        CrystalSymmetries[22].set( 0, 1, 0, 1, 0, 0, 0, 0,-1);
        CrystalSymmetries[23].set( 0,-1, 0, 1, 0, 0, 0, 0, 1);
    }
    else if (Symmetry == "HEXAGONAL")
    {
        //  Load hexagonal parameters (class 622)
        CrystalSym = CrystalSymmetry::HEXAGONAL;

        nsym = 12;
        Theta_min = 0 * (M_PI / 180);
        Theta_max = 30 * (M_PI / 180);
        Phi_min = 0 * (M_PI / 180);
        Phi_max = M_PI / 2;

        CrystalSymmetries.Allocate(nsym);

        CrystalSymmetries[0].set( 1, 0, 0, 0, 1, 0, 0, 0, 1);
        CrystalSymmetries[1].set(-0.5, a, 0, -a, -0.5, 0, 0, 0, 1);
        CrystalSymmetries[2].set(-0.5,-a, 0,  a, -0.5, 0, 0, 0, 1);
        CrystalSymmetries[3].set( 0.5, a, 0, -a,  0.5, 0, 0, 0, 1);
        CrystalSymmetries[4].set(-1, 0, 0, 0,-1, 0, 0, 0, 1);
        CrystalSymmetries[5].set( 0.5,-a, 0,  a,  0.5, 0, 0, 0, 1);

        CrystalSymmetries[6].set(-0.5, -a, 0, -a, 0.5, 0, 0, 0, -1);
        CrystalSymmetries[7].set( 1, 0, 0, 0, -1, 0, 0, 0, -1);
        CrystalSymmetries[8].set(-0.5,  a, 0,  a, 0.5, 0, 0, 0, -1);
        CrystalSymmetries[9].set( 0.5,  a, 0,  a,-0.5, 0, 0, 0, -1);
        CrystalSymmetries[10].set(-1, 0, 0, 0, 1, 0, 0, 0, -1);
        CrystalSymmetries[11].set(0.5, -a, 0, -a,-0.5, 0, 0, 0, -1);
    }
    else if (Symmetry == "TETRAGONAL")
    {
        //  Load tetragonal parameters (class 422)
        CrystalSym = CrystalSymmetry::TETRAGONAL;

        nsym = 8;
        Theta_min = 0 * (M_PI / 180);
        Theta_max = 45 * (M_PI / 180);
        Phi_min = 0 * (M_PI / 180);
        Phi_max = M_PI / 2;

        CrystalSymmetries.Allocate(nsym);

        CrystalSymmetries[0].set( 1, 0, 0, 0, 1, 0, 0, 0, 1);
        CrystalSymmetries[1].set(-1, 0, 0, 0, 1, 0, 0, 0,-1);
        CrystalSymmetries[2].set( 1, 0, 0, 0,-1, 0, 0, 0,-1);
        CrystalSymmetries[3].set(-1, 0, 0, 0,-1, 0, 0, 0, 1);
        CrystalSymmetries[4].set( 0, 1, 0,-1, 0, 0, 0, 0, 1);
        CrystalSymmetries[5].set( 0,-1, 0, 1, 0, 0, 0, 0, 1);
        CrystalSymmetries[6].set( 0, 1, 0, 1, 0, 0, 0, 0,-1);
        CrystalSymmetries[7].set( 0,-1, 0,-1, 0, 0, 0, 0,-1);
    }
    else if (Symmetry == "TRIGONAL")
    {
        //  Load trigonal parameters (class 32)
        CrystalSym = CrystalSymmetry::TRIGONAL;

        nsym = 6;
        Theta_min = 0 * (M_PI / 180);
        Theta_max = 60 * (M_PI / 180);
        Phi_min = 0 * (M_PI / 180);
        Phi_max = M_PI / 2;

        CrystalSymmetries.Allocate(nsym);

        CrystalSymmetries[0].set( 1, 0, 0, 0, 1, 0, 0, 0, 1);
        CrystalSymmetries[1].set(-0.5, a, 0,-a, -0.5, 0, 0, 0, 1);
        CrystalSymmetries[2].set(-0.5,-a, 0, a, -0.5, 0, 0, 0, 1);
        CrystalSymmetries[3].set( 0.5, a, 0, a, -0.5, 0, 0, 0,-1);
        CrystalSymmetries[4].set(-1, 0, 0, 0, 1, 0, 0, 0, 1);
        CrystalSymmetries[5].set( 0.5,-a, 0,-a, -0.5, 0, 0, 0,-1);
    }
    else if (Symmetry == "ORTHORHOMBIC")
    {
        //  Load orthorhombic parameters (class 22)
        CrystalSym = CrystalSymmetry::ORTHORHOMBIC;

        nsym = 4;
        Theta_min = 0 * (M_PI / 180);
        Theta_max = 90 * (M_PI / 180);
        Phi_min = 0 * (M_PI / 180);
        Phi_max = M_PI / 2;

        CrystalSymmetries.Allocate(nsym);

        CrystalSymmetries[0].set( 1, 0, 0, 0, 1, 0, 0, 0, 1);
        CrystalSymmetries[1].set(-1, 0, 0, 0, 1, 0, 0, 0,-1);
        CrystalSymmetries[2].set( 1, 0, 0, 0,-1, 0, 0, 0, 1);
        CrystalSymmetries[3].set(-1, 0, 0, 0,-1, 0, 0, 0, 1);
    }
    else if (Symmetry == "MONOCLINIC")
    {
        //  Load monoclinic parameters (class 2)
        CrystalSym = CrystalSymmetry::MONOCLINIC;

        nsym = 2;
        Theta_min = 0 * (M_PI / 180);
        Theta_max = 180 * (M_PI / 180);
        Phi_min = 0 * (M_PI / 180);
        Phi_max = M_PI / 2;

        CrystalSymmetries.Allocate(nsym);

        CrystalSymmetries[0].set( 1, 0, 0, 0, 1, 0, 0, 0, 1);
        CrystalSymmetries[1].set(-1, 0, 0, 0, 1, 0, 0, 0,-1);

    }
    else if (Symmetry == "TRICLINIC")
    {
        //  Load triclinic parameters (class 1)
        CrystalSym = CrystalSymmetry::TRICLINIC;

        nsym = 1;
        Theta_min = 0 * (M_PI / 180);
        Theta_max = 360 * (M_PI / 180);
        Phi_min = 0 * (M_PI / 180);
        Phi_max = M_PI / 2;

        CrystalSymmetries.Allocate(nsym);

        CrystalSymmetries[0].set(1, 0, 0, 0, 1, 0, 0, 0, 1);
    }
    else
    {
        nsym = 0;
        Info::WriteWarning("No or wrong CrystalSymmetry specified!\nThe default \"Cubic\" model is used!", thisclassname, "ReadInput()");
    }

    Info::WriteLine();
    Info::WriteBlankLine();
}

void Crystallography::Initialize(Settings& locSettings)
{
    thisclassname = "Crystallography";

    // Symmetry matrices for cubic lattice. Taken from Kocks, Tome and Wenk
    // 'Texture and Anisotropy', Cambridge University Press 1998

    dMatrix3x3 SCubic1  { 1, 0, 0,   0, 1, 0,   0, 0, 1 };
    dMatrix3x3 SCubic2  { 0, 0,-1,   0,-1, 0,  -1, 0, 0 };
    dMatrix3x3 SCubic3  { 0, 0, 1,   1, 0, 0,   0, 1, 0 };
    dMatrix3x3 SCubic4  { 0, 0, 1,   0,-1, 0,   1, 0, 0 };
    dMatrix3x3 SCubic5  { 0, 1, 0,   0, 0, 1,   1, 0, 0 };
    dMatrix3x3 SCubic6  { 0, 0, 1,   0, 1, 0,  -1, 0, 0 };

    dMatrix3x3 SCubic7  { 0,-1, 0,   0, 0, 1,  -1, 0, 0 };
    dMatrix3x3 SCubic8  { 0, 0,-1,   0, 1, 0,   1, 0, 0 };
    dMatrix3x3 SCubic9  { 0,-1, 0,   0, 0,-1,   1, 0, 0 };
    dMatrix3x3 SCubic10 {-1, 0, 0,   0, 0,-1,   0,-1, 0 };
    dMatrix3x3 SCubic11 { 0, 1, 0,   0, 0,-1,  -1, 0, 0 };
    dMatrix3x3 SCubic12 { 1, 0, 0,   0, 0,-1,   0, 1, 0 };

    dMatrix3x3 SCubic13 { 0, 0,-1,   1, 0, 0,   0,-1, 0 };
    dMatrix3x3 SCubic14 { 1, 0, 0,   0, 0, 1,   0,-1, 0 };
    dMatrix3x3 SCubic15 { 0, 0,-1,  -1, 0, 0,   0, 1, 0 };
    dMatrix3x3 SCubic16 {-1, 0, 0,   0, 0, 1,   0, 1, 0 };
    dMatrix3x3 SCubic17 { 0, 0, 1,  -1, 0, 0,   0,-1, 0 };
    dMatrix3x3 SCubic18 { 0,-1, 0,  -1, 0, 0,   0, 0,-1 };

    dMatrix3x3 SCubic19 {-1, 0, 0,   0, 1, 0,   0, 0,-1 };
    dMatrix3x3 SCubic20 { 0, 1, 0,  -1, 0, 0,   0, 0, 1 };
    dMatrix3x3 SCubic21 {-1, 0, 0,   0,-1, 0,   0, 0, 1 };
    dMatrix3x3 SCubic22 { 0, 1, 0,   1, 0, 0,   0, 0,-1 };
    dMatrix3x3 SCubic23 { 1, 0, 0,   0,-1, 0,   0, 0,-1 };
    dMatrix3x3 SCubic24 { 0,-1, 0,   1, 0, 0,   0, 0, 1 };

    SymmetriesCubic.Allocate(24);

    SymmetriesCubic[0] = SCubic1;
    SymmetriesCubic[1] = SCubic2;
    SymmetriesCubic[2] = SCubic3;
    SymmetriesCubic[3] = SCubic4;
    SymmetriesCubic[4] = SCubic5;
    SymmetriesCubic[5] = SCubic6;
    SymmetriesCubic[6] = SCubic7;
    SymmetriesCubic[7] = SCubic8;
    SymmetriesCubic[8] = SCubic9;
    SymmetriesCubic[9] = SCubic10;
    SymmetriesCubic[10] = SCubic11;
    SymmetriesCubic[11] = SCubic12;
    SymmetriesCubic[12] = SCubic13;
    SymmetriesCubic[13] = SCubic14;
    SymmetriesCubic[14] = SCubic15;
    SymmetriesCubic[15] = SCubic16;
    SymmetriesCubic[16] = SCubic17;
    SymmetriesCubic[17] = SCubic18;
    SymmetriesCubic[18] = SCubic19;
    SymmetriesCubic[19] = SCubic20;
    SymmetriesCubic[20] = SCubic21;
    SymmetriesCubic[21] = SCubic22;
    SymmetriesCubic[22] = SCubic23;
    SymmetriesCubic[23] = SCubic24;

    Info::WriteStandard(thisclassname, "Initialized");
}

const std::vector<double> Crystallography::SlipSystem_Cubic
{ //  {  n1,  n2,  n3,  d1,  d2,  d3, }
    0.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 1.0,-1.0, 1.0, 0.0, 0.0,

    1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    -1.0, 0.0, 1.0, 0.0, 1.0, 0.0,

    1.0, 1.0, 0.0, 0.0, 0.0, 1.0,
    -1.0, 1.0, 0.0, 0.0, 0.0, 1.0
};
const std::vector<double> Crystallography::SlipSystem_BCC
{ //  {  n1,  n2,  n3,  d1,  d2,  d3 }
    0.0, 1.0,-1.0, 1.0, 1.0, 1.0,
    1.0, 0.0,-1.0, 1.0, 1.0, 1.0,
    1.0,-1.0, 0.0, 1.0, 1.0, 1.0,

    0.0, 1.0,-1.0,-1.0, 1.0, 1.0,
    1.0, 0.0, 1.0,-1.0, 1.0, 1.0,
    1.0, 1.0, 0.0,-1.0,-1.0, 1.0,

    0.0, 1.0, 1.0, 1.0, 1.0,-1.0,
    1.0, 0.0, 1.0, 1.0, 1.0, 1.0,
    1.0,-1.0, 0.0, 1.0, 1.0,-1.0,

    0.0, 1.0, 1.0, 1.0,-1.0, 1.0,
    1.0, 0.0,-1.0, 1.0,-1.0, 1.0,
    1.0, 1.0, 0.0, 1.0,-1.0, 1.0
};
const std::vector<double> Crystallography::SlipSystem_FCC_XY
{
    1.0, 1.0, 0.0, 1.0,-1.0, 0.0,
    1.0, 1.0, 0.0,-1.0, 1.0, 0.0,
    -1.0, 1.0, 0.0, 1.0, 1.0, 0.0,
    -1.0, 1.0, 0.0,-1.0,-1.0, 0.0,
};
const std::vector<double> Crystallography::SlipSystem_FCC_XZ
{
    1.0, 0.0, 1.0, 1.0, 0.0,-1.0,
    1.0, 0.0, 1.0,-1.0, 0.0, 1.0,
    -1.0, 0.0, 1.0, 1.0, 0.0, 1.0,
    -1.0, 0.0, 1.0,-1.0, 0.0,-1.0,
};
const std::vector<double> Crystallography::SlipSystem_FCC_YZ
{
    0.0, 1.0, 1.0, 0.0, 1.0,-1.0,
    0.0, 1.0, 1.0, 0.0,-1.0, 1.0,
    0.0,-1.0, 1.0, 0.0, 1.0, 1.0,
    0.0,-1.0, 1.0, 0.0,-1.0,-1.0
};
const std::vector<double> Crystallography::SlipSystem_FCC
{
    1.0, 1.0, 1.0, 1.0,-1.0, 0.0,
    -1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
    1.0, 1.0,-1.0, 1.0,-1.0, 0.0,
    1.0,-1.0, 1.0, 1.0, 1.0, 0.0,

    1.0,-1.0, 1.0, 1.0, 0.0,-1.0,
    -1.0, 1.0, 1.0, 1.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.0,-1.0,
    1.0, 1.0,-1.0, 1.0, 0.0, 1.0,

    -1.0, 1.0, 1.0, 0.0, 1.0,-1.0,
    1.0, 1.0,-1.0, 0.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 1.0,-1.0,
    1.0,-1.0, 1.0, 0.0, 1.0, 1.0
};
const std::vector<double> Crystallography::SlipSystem_PartialFCC =
{ //  {  n1,  n2,  n3,  d1,  d2,  d3 }
    -2.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0,-2.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0,-2.0, 1.0, 1.0, 1.0,

    2.0, 1.0, 1.0,-1.0, 1.0, 1.0,
    -1.0,-2.0, 1.0,-1.0, 1.0, 1.0,
    -1.0, 1.0,-2.0,-1.0, 1.0, 1.0,

    -2.0,-1.0, 1.0, 1.0,-1.0, 1.0,
    1.0, 2.0, 1.0, 1.0,-1.0, 1.0,
    1.0,-1.0,-2.0, 1.0,-1.0, 1.0,

    -2.0, 1.0,-1.0, 1.0, 1.0,-1.0,
    1.0,-2.0,-1.0, 1.0, 1.0,-1.0,
    1.0, 1.0, 2.0, 1.0, 1.0,-1.0
};
const std::vector<double> Crystallography::SlipSystem_HCPBasal
{ //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
    //Basal systems {0001}<11-20>
    0.0, 0.0, 0.0, 1.0, 2.0,-1.0,-1.0, 0.0,
    0.0, 0.0, 0.0, 1.0,-1.0, 2.0,-1.0, 0.0,
    0.0, 0.0, 0.0, 1.0,-1.0,-1.0, 2.0, 0.0
};
const std::vector<double> Crystallography::SlipSystem_HCPPrismatic
{ //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
    //Prismatic systems {10-10}<11-20>
    0.0, 1.0,-1.0, 0.0, 2.0,-1.0,-1.0, 0.0,
    -1.0, 0.0, 1.0, 0.0,-1.0, 2.0,-1.0, 0.0,
    1.0,-1.0, 0.0, 0.0,-1.0,-1.0, 2.0, 0.0
};
const std::vector<double> Crystallography::SlipSystem_HCPPyramidal
{ //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
    //Pyramidal systems <a> {10-11}<11-20>
    0.0, 1.0,-1.0, 1.0, 2.0,-1.0,-1.0, 0.0,
    -1.0, 0.0, 1.0, 1.0,-1.0, 2.0,-1.0, 0.0,
    1.0,-1.0, 0.0, 1.0,-1.0,-1.0, 2.0, 0.0,
    -1.0, 1.0, 0.0, 1.0, 1.0, 1.0,-2.0, 0.0,
    0.0,-1.0, 1.0, 1.0,-2.0, 1.0, 1.0, 0.0,
    1.0, 0.0,-1.0, 1.0, 1.0,-2.0, 1.0, 0.0
};
const std::vector<double> Crystallography::SlipSystem_HCPPyramidal_1st
{ //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
    //Pyramidal systems (1st) <c+a> {10-11}<11-23>
    -1.0, 1.0, 0.0, 1.0, 2.0, -1.0,-1.0, 3.0,
    -1.0, 1.0, 0.0, 1.0, 1.0, -2.0, 1.0, 3.0,
    1.0, 0.0,-1.0, 1.0,-1.0, -1.0, 2.0, 3.0,

    1.0, 0.0,-1.0, 1.0,-2.0, 1.0, 1.0, 3.0,
    0.0,-1.0, 1.0, 1.0,-1.0, 2.0,-1.0, 3.0,
    0.0,-1.0, 1.0, 1.0, 1.0, 1.0,-2.0, 3.0,

    1.0,-1.0, 0.0, 1.0,-2.0, 1.0, 1.0, 3.0,
    1.0,-1.0,-0.0, 1.0,-1.0, 2.0,-1.0, 3.0,
    -1.0, 0.0, 1.0, 1.0, 1.0, 1.0,-2.0, 3.0,

    -1.0, 0.0, 1.0, 1.0, 2.0,-1.0,-1.0, 3.0,
    0.0, 1.0,-1.0, 1.0, 1.0,-2.0, 1.0, 3.0,
    0.0, 1.0,-1.0, 1.0,-1.0,-1.0, 2.0, 3.0
};
const std::vector<double> Crystallography::SlipSystem_HCPPyramidal_2nd
{ //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
    //Pyramidal systems (2nd) <c+a> {11-22}<11-23>
    -2.0, 1.0, 1.0, 2.0, 2.0,-1.0,-1.0, 3.0,
    1.0,-2.0, 1.0, 2.0,-1.0, 2.0,-1.0, 3.0,
    1.0, 1.0,-2.0, 2.0,-1.0,-1.0, 2.0, 3.0,
    2.0,-1.0,-1.0, 2.0,-2.0, 1.0, 1.0, 3.0,
    -1.0, 2.0,-1.0, 2.0, 1.0,-2.0, 1.0, 3.0,
    -1.0,-1.0, 2.0, 2.0, 1.0, 1.0,-2.0, 3.0
};
const std::vector<double> Crystallography::TwinSystem_FCC
{
    -2.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0,-2.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0,-2.0, 1.0, 1.0, 1.0,

    2.0, 1.0, 1.0,-1.0, 1.0, 1.0,
    -1.0,-2.0, 1.0,-1.0, 1.0, 1.0,
    -1.0, 1.0,-2.0,-1.0, 1.0, 1.0,

    -2.0,-1.0, 1.0, 1.0,-1.0, 1.0,
    1.0, 2.0, 1.0, 1.0,-1.0, 1.0,
    1.0,-1.0,-2.0, 1.0,-1.0, 1.0,

    -2.0, 1.0,-1.0, 1.0, 1.0,-1.0,
    1.0,-2.0,-1.0, 1.0, 1.0,-1.0,
    1.0, 1.0, 2.0, 1.0, 1.0,-1.0
};
const std::vector<double> Crystallography::TwinSystem_HCPTensileTwin1
{ //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
    //Tensile twinning in Co, Mg, Zr, Ti, and Be; compressive twinning in Cd and Zn (T1)
    -1.0, 1.0, 0.0, 2.0, 1.0,-1.0, 0.0, 1.0,
    1.0, 0.0,-1.0, 2.0,-1.0, 0.0, 1.0, 1.0,
    0.0,-1.0, 1.0, 2.0, 0.0, 1.0,-1.0, 1.0,
    1.0,-1.0, 0.0, 2.0,-1.0, 1.0, 0.0, 1.0,
    -1.0, 0.0, 1.0, 2.0, 1.0, 0.0,-1.0, 1.0,
    0.0, 1.0,-1.0, 2.0, 0.0,-1.0, 1.0, 1.0
};
const std::vector<double> Crystallography::TwinSystem_HCPTensileTwin2
{ //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
    //Tensile twinning in Co, Re, and Zr (T2)
    -2.0, 1.0, 1.0, 1.0, 2.0,-1.0,-1.0, 6.0,
    1.0,-2.0, 1.0, 1.0,-1.0, 2.0,-1.0, 6.0,
    1.0, 1.0,-2.0, 1.0,-1.0,-1.0, 2.0, 6.0,
    2.0,-1.0,-1.0, 1.0,-2.0, 1.0, 1.0, 6.0,
    -1.0, 2.0,-1.0, 1.0, 1.0,-2.0, 1.0, 6.0,
    -1.0,-1.0, 2.0, 1.0, 1.0, 1.0,-2.0, 6.0
};
const std::vector<double> Crystallography::TwinSystem_HCPCompressionTwin1
{ //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
    //Compressive twinning in Ti and Zr (C1)
    2.0,-1.0,-1.0, 2.0, 2.0,-1.0,-1.0,-3.0,
    -1.0, 2.0,-1.0, 2.0,-1.0, 2.0,-1.0,-3.0,
    -1.0,-1.0, 2.0, 2.0,-1.0,-1.0, 2.0,-3.0,
    -2.0, 1.0, 1.0, 2.0,-2.0, 1.0, 1.0,-3.0,
    1.0,-2.0, 1.0, 2.0, 1.0,-2.0, 1.0,-3.0,
    1.0, 1.0,-2.0, 2.0, 1.0, 1.0,-2.0,-3.0
};
const std::vector<double> Crystallography::TwinSystem_HCPCompressionTwin2
{ //  {  n1,  n2,  n3, n4,  d1,  d2,  d3,  d4}
    //Compressive twinning in Mg (C1)
    -1.0, 1.0, 0.0,-2.0,-1.0, 1.0, 0.0, 1.0,
    1.0, 0.0,-1.0,-2.0, 1.0, 0.0,-1.0, 1.0,
    0.0,-1.0, 1.0,-2.0, 0.0,-1.0, 1.0, 1.0,
    1.0,-1.0, 0.0,-2.0, 1.0,-1.0, 0.0, 1.0,
    -1.0, 0.0, 1.0,-2.0,-1.0, 0.0, 1.0, 1.0,
    0.0, 1.0,-1.0,-2.0, 0.0, 1.0,-1.0, 1.0
};

}// namespace openphase
