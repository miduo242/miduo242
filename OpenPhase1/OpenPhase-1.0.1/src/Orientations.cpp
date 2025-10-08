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
 *   Main contributors :   Efim Borukhovich; Philipp Engels; Oleg Shchyglo
 *
 */

#include "Orientations.h"
#include "Info.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"
#include "Settings.h"
#include "Base/UserInterface.h"
#include "Base/Quaternion.h"
#include "Tools.h"
#include "Crystallography.h"
#include "Velocities.h"
#include "VTK.h"
#include "AdvectionHR/AdvectionHR.h"

namespace openphase
{
using namespace std;

Orientations::Orientations(Settings& locSettings)
{
    Initialize(locSettings);
}

double Orientations::getMisorientation(const dMatrix3x3 RotMatA, const dMatrix3x3 RotMatB)
{
    // Taken from Kocks, Tome and Wenk
    // 'Texture and Anisotropy', Cambride University Press 1998, p. 69, eq. 9b
    return acos(((RotMatB*RotMatA.transposed()).trace()-1.0)/2.0);
}

double Orientations::getMisorientationCubic(const dMatrix3x3 RotMatA, const dMatrix3x3 RotMatB, const Crystallography& CR)
{
    // Taken from Kocks, Tome and Wenk
    // 'Texture and Anisotropy', Cambride University Press 1998, p. 69, eq. 9b

    // Only for cubic system. Needs to be generalized!

    double misorientation = 10.0;
    for (int i = 0; i < 24; i++)
    {
        double misorientationloc = acos(((CR.SymmetriesCubic[i]*RotMatB*
                RotMatA.transposed()).trace()-1.0)/2.0);
        if (abs(misorientationloc) <= abs(misorientation))
        {
            misorientation = misorientationloc;
        }
    }
    return misorientation;
}

double Orientations::getDisorientation(const Quaternion OrientationA, const Quaternion OrientationB)
{
    // Taken from Grimmer 1974
    // Only for cubic system. Needs to be generalized!

    double Angle[24];
    double misorientation = 0.0;
    Quaternion deltaQ;
    deltaQ  = OrientationA* OrientationB.inverted();

    Angle[0]  = deltaQ[0];
    Angle[1]  = deltaQ[1];
    Angle[2]  = deltaQ[2];
    Angle[3]  = deltaQ[3];

    Angle[4]  = (deltaQ[0] + deltaQ[1])*0.707;
    Angle[5]  = (deltaQ[0] - deltaQ[1])*0.707;
    Angle[6]  = (deltaQ[2] + deltaQ[3])*0.707;
    Angle[7]  = (deltaQ[2] - deltaQ[3])*0.707;

    Angle[8]  = (deltaQ[0] + deltaQ[2])*0.707;
    Angle[9]  = (deltaQ[0] - deltaQ[2])*0.707;
    Angle[10] = (deltaQ[1] + deltaQ[3])*0.707;
    Angle[11] = (deltaQ[1] - deltaQ[3])*0.707;

    Angle[12] = (deltaQ[0] + deltaQ[3])*0.707;
    Angle[13] = (deltaQ[0] - deltaQ[3])*0.707;
    Angle[14] = (deltaQ[1] + deltaQ[2])*0.707;
    Angle[15] = (deltaQ[1] - deltaQ[2])*0.707;

    Angle[16] = (deltaQ[0] + deltaQ[1] + deltaQ[2] + deltaQ[3])*0.5;
    Angle[17] = (deltaQ[0] + deltaQ[1] - deltaQ[2] - deltaQ[3])*0.5;
    Angle[18] = (deltaQ[0] - deltaQ[1] + deltaQ[2] - deltaQ[3])*0.5;
    Angle[19] = (deltaQ[0] - deltaQ[1] - deltaQ[2] + deltaQ[3])*0.5;

    Angle[20] = (deltaQ[0] + deltaQ[1] + deltaQ[2] - deltaQ[3])*0.5;
    Angle[21] = (deltaQ[0] + deltaQ[1] - deltaQ[2] + deltaQ[3])*0.5;
    Angle[22] = (deltaQ[0] - deltaQ[1] + deltaQ[2] + deltaQ[3])*0.5;
    Angle[23] = (deltaQ[0] - deltaQ[1] - deltaQ[2] - deltaQ[3])*0.5;

    for (int a=0; a<24; a++)
    {
        Angle[a] = abs(Angle[a]);
    }
    double DisAngle = 0;
    for (int a=0; a<24; a++)
    {
        if (DisAngle < Angle[a])
            DisAngle =  Angle[a];
    }
    misorientation = 2*acos(DisAngle);
    return misorientation;
}

dVector3 Orientations::MillerConversion(const std::vector<double>& hkil, dVector3 hkl, bool PlaneNormal)
{
    if(PlaneNormal)
    {
        hkl[0] = hkil[0];
        hkl[1] = hkil[1];
        hkl[2] = hkil[3];
    }
    else
    {
        hkl[0] = hkil[0]-hkil[2];
        hkl[1] = hkil[1]-hkil[2];
        hkl[2] = hkil[3];
    }

    // Titanium c/a ratio ( to be generalized depends on the element )
    double a = 0.2950;
    double c = 0.4683;

    dMatrix3x3 locMat {a, -a/2.0, 0,
                       0, a*sqrt(3)/2.0, 0,
                       0, 0, c};
    hkl = locMat * hkl;

    hkl.normalize();
    return hkl;
}

void Orientations::Initialize(Settings& locSettings)
{
    thisclassname = "Orientations";

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;

    dNx = locSettings.dNx;
    dNy = locSettings.dNy;
    dNz = locSettings.dNz;

    dx = locSettings.dx;

    Nphases = locSettings.Nphases;

    size_t Bcells = locSettings.Bcells;
    Quaternions.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Quaternions, Bcells,)
    {
        Quaternions(i,j,k).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    VTKDir = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;

    initialized = true;
    Info::WriteStandard("Orientations", "Initialized");
}

void Orientations::SetRandomGrainOrientations(PhaseField& Phase, const int seed)
{
    // uniform sampling of quaternion parameters using the method presented in (Shoemake, 1992)
    // taken from :: http://planning.cs.uiuc.edu/node198.html#eqn:shoemake

    default_random_engine generator(seed);

    uniform_real_distribution <double> Q1Distribution(0, 1.0);
    uniform_real_distribution <double> Q2Distribution(0, 1.0);
    uniform_real_distribution <double> Q3Distribution(0, 1.0);

    uniform_real_distribution <double> A1Distribution(0.0, 2.0*Pi);
    uniform_real_distribution <double> A2Distribution(0.0, 2.0*Pi);
    uniform_real_distribution <double> A3Distribution(0.0, 2.0*Pi);

    for(size_t alpha = 0; alpha < Phase.FieldsStatistics.size(); alpha++)
    if(Phase.FieldsStatistics[alpha].Exist)
    {
        Quaternion tempQuat;
        switch(dNx + dNy + dNz)
        {
            case 1: // 1D
            {
                Phase.FieldsStatistics[alpha].Orientation.set(1.0, 0.0, 0.0, 0.0);
                break;
            }
            case 2: // 2D
            {
                double a1 = 0.0;
                double a2 = 0.0;
                double a3 = 0.0;
                if(Phase.dNx == 0) a1 = A1Distribution(generator);
                if(Phase.dNy == 0) a2 = A2Distribution(generator);
                if(Phase.dNz == 0) a3 = A3Distribution(generator);
                EulerAngles ph1({a1,a2,a3},XYZ);
                Phase.FieldsStatistics[alpha].Orientation = ph1.getQuaternion().normalized();
                break;
            }
            case 3: // 3D Full Rotation
            {
                double u1 = Q1Distribution(generator);
                double u2 = Q2Distribution(generator);
                double u3 = Q3Distribution(generator);

                tempQuat.set(sqrt(1.0-u1)*sin(2.0*Pi*u2),
                             sqrt(1.0-u1)*cos(2.0*Pi*u2),
                             sqrt(u1)*sin(2.0*Pi*u3),
                             sqrt(u1)*cos(2.0*Pi*u3));

                Phase.FieldsStatistics[alpha].Orientation = tempQuat.normalized();
                break;
            }
        }
    }
}

void Orientations::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(Quaternions);
    BC.SetY(Quaternions);
    BC.SetZ(Quaternions);
}

void Orientations::Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)
{
    Quaternions.Remesh(newNx, newNy, newNz);

    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    SetBoundaryConditions(BC);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Quaternions, Quaternions.Bcells(),)
    {
        Quaternions(i, j, k).setRotationMatrix();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Info::WriteStandard(thisclassname, "Remeshed");
}

Orientations& Orientations::operator= (const Orientations& rhs)
{
    // protect against invalid self-assignment and copy of unitialized object
    if (this != &rhs and rhs.thisclassname == "Orientations")
    {
        thisclassname = rhs.thisclassname;
        //DefaultInputFileName = rhs.DefaultInputFileName;

        Nx = rhs.Nx;
        Ny = rhs.Ny;
        Nz = rhs.Nz;

        dNx = rhs.dNx;
        dNy = rhs.dNy;
        dNz = rhs.dNz;

        dx = rhs.dx;

        Nphases = rhs.Nphases;

        Quaternions = rhs.Quaternions;
    }
    return *this;
}

void Orientations::Write(const int tStep) const
{
    string FileName = UserInterface::MakeFileName(RawDataDir,"Rotations_", tStep, ".dat");

    fstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be created", thisclassname, "WriteRotations");
        exit(1);
    };

    for(long int i = 0; i < Quaternions.sizeX(); ++i)
    for(long int j = 0; j < Quaternions.sizeY(); ++j)
    for(long int k = 0; k < Quaternions.sizeZ(); ++k)
    {
        for(int n = 0; n < 4; n++)
        {
            double tmp = Quaternions(i,j,k)[n];
            out.write(reinterpret_cast<char*>(&tmp), sizeof(double));
        }
    }
}

void Orientations::Read(const int tStep)
{
    string FileName = UserInterface::MakeFileName(RawDataDir,"Rotations_",
                                                  tStep, ".dat");

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be created",
                        thisclassname, "ReadRotations");
        exit(1);
    };

    for(long int i = 0; i < Quaternions.sizeX(); ++i)
    for(long int j = 0; j < Quaternions.sizeY(); ++j)
    for(long int k = 0; k < Quaternions.sizeZ(); ++k)
    {
        double tmp[4];

        for(int n = 0; n < 4; n++)
        {
            inp.read(reinterpret_cast<char*>(&tmp[n]), sizeof(double));
        }
        Quaternions(i,j,k).set(tmp[0],tmp[1],tmp[2],tmp[3]);
    }
}

void Orientations::WriteVTK(const int tStep, const Settings& locSettings,
                            const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"RQ_0", [this](int i,int j,int k){return Quaternions(i,j,k)[0];}});
    ListOfFields.push_back((VTK::Field_t) {"RQ_1", [this](int i,int j,int k){return Quaternions(i,j,k)[1];}});
    ListOfFields.push_back((VTK::Field_t) {"RQ_2", [this](int i,int j,int k){return Quaternions(i,j,k)[2];}});
    ListOfFields.push_back((VTK::Field_t) {"RQ_3", [this](int i,int j,int k){return Quaternions(i,j,k)[3];}});

    std::string Filename = UserInterface::MakeFileName(VTKDir, "Quaternions_", tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}

void Orientations::WriteTotalVTK(const PhaseField& Phase, const int tStep,
                                 const Settings& locSettings,
                                 const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"R", [this,Phase](int i,int j,int k){return getTotalRotation(Phase, i, j, k);}});

    std::string Filename = UserInterface::MakeFileName(VTKDir, "TotalRotations_", tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}

void Orientations::WriteMisorientationsVTK(const int tStep,
                                           const Settings& locSettings,
                                           const int precision,
                                           const std::string measure) const
{
    double mult = 1.0;
    if (!measure.compare("deg")) mult*=180.0/Pi;

    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"Omega", [this,mult](int i,int j,int k){return getMisorientation(dMatrix3x3::UnitTensor(), Quaternions(i,j,k).RotationMatrix)*mult;}});

    std::string Filename = UserInterface::MakeFileName(VTKDir, "Misorientations_", tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}

EulerAngles Orientations::RotationToEuler(const dMatrix3x3& Rot, const EulerConvention EConvention )
{
    EulerAngles Euler;
    Euler.Convention = EConvention;
    dMatrix3x3 RT = Rot.transposed();

    if(EConvention == ZXZ)
    {
        // Euler ZXZ passive (following formulation by Martin Boeff)
        double squvw = sqrt(RT(0,0)*RT(0,0)
            + RT(1,0)*RT(1,0) + RT(2,0)*RT(2,0));
        double sqhk = sqrt(RT(0,2)*RT(0,2)
            + RT(1,2)*RT(1,2));
        double sqhkl = sqrt(RT(0,2)*RT(0,2)
            + RT(1,2)*RT(1,2) + RT(2,2)*RT(2,2));
        double tempval = RT(2,2)/sqhkl;

        if(tempval >  1.0) {tempval =  1.0;}
        if(tempval < -1.0) {tempval = -1.0;}
        Euler.Q[1] = acos(tempval);

        if(Euler.Q[1] < 1.0e-8)
        {
            // calculate phi2
            Euler.Q[2] = 0.0;
            // calculate phi1
            tempval = RT(0,0)/squvw;
            if(tempval >  1.0) {tempval =  1.0;}
            if(tempval < -1.0) {tempval = -1.0;}

            Euler.Q[0] = acos(tempval);
            if(RT(1,0) > 0.0) {Euler.Q[0] = 2.0*Pi - Euler.Q[0];}
        }
        else
        {
            // calculate phi2
            tempval = RT(1,2)/sqhk;
            if(tempval >  1.0) {tempval =  1.0;}
            if(tempval < -1.0) {tempval = -1.0;}

            Euler.Q[2] = acos(tempval);
            if(RT(0,2) < 0.0) {Euler.Q[2] = 2.0*Pi - Euler.Q[2];}
            // calculate phi1
            tempval = - RT(2,1)/Euler.SinQ[1];
            if(tempval >  1.0) {tempval = 1.0;}
            if(tempval < -1.0) {tempval = -1.0;}

            Euler.Q[0] = acos(tempval);
            if(RT(2,0) < 0.0) {Euler.Q[0] = 2.0*Pi - Euler.Q[0];}
        }
        Euler.setTrigonometricFunctions();
    }
    else
    {
        std::stringstream message;
        message<< "Wrong/Unknown/None Euler convention used: " << EConvention ;

        Info::WriteExit(message.str(), "Orientations", "RotationToEulerAngles()");
        exit(13);
    }
    return Euler;
}

void Orientations::WriteGrainEBSDDataQuaternions(const PhaseField& Phase,
                                                        const int tStep)
{
    stringstream outbuffer;
    outbuffer << "Index\t"
              << "Phase\t"
              << "X\t"
              << "Y\t"
              << "Z\t"
              << "Quat real\t"
              << "Quat i\t"
              << "Quat j\t"
              << "Quat k\t" << endl;

    int index = 1;
    STORAGE_LOOP_BEGIN(i,j,k,Quaternions,0)
    {
        EulerAngles tempAng;
        dMatrix3x3  tempMatrix;
        Quaternion  tempQuat;

        outbuffer << index << "\t";
        if (Phase.Interface(i,j,k))
        {
            // linear interpolation in interface via quaternions
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha < Phase.Fields(i, j, k).cend(); ++alpha)
            {
                tempQuat += (Quaternions(i,j,k) + Phase.FieldsStatistics[alpha->index].Orientation)*alpha->value;
            }
            tempQuat.normalize();
            // Set phase index to 0 in interface
            outbuffer << 0 << "\t";
        }
        else
        {
            tempQuat = Quaternions(i,j,k) + Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Orientation;
            outbuffer << Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase + 1 << "\t";
        }

        outbuffer << i*Phase.dx << "\t"
                  << j*Phase.dx << "\t"
                  << k*Phase.dx << "\t"
                  << tempQuat[0] << "\t"
                  << tempQuat[1] << "\t"
                  << tempQuat[2] << "\t"
                  << tempQuat[3] << endl;
        index++;
    }
    STORAGE_LOOP_END

    string FileName = UserInterface::MakeFileName(RawDataDir, "EBSD_", tStep, ".dat");

    ofstream ebsd_file(FileName.c_str());
    ebsd_file << outbuffer.rdbuf();
    ebsd_file.close();

    Info::WriteStandard("EBSD file", FileName);
}

void Orientations::WriteRotated100Vector(const int tStep)
{
    stringstream outbufer;
    dVector3 X;
    X[0] = 1.0;
    X[1] = 0.0;
    X[2] = 0.0;

    outbufer << "# vtk DataFile Version 3.0\n";
    outbufer << "Rotated100Vector" << "\n";
    outbufer << "ASCII\n";
    outbufer << "DATASET RECTILINEAR_GRID\n";
    outbufer << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
    outbufer << "X_COORDINATES " << Nx << " double\n";
    for (int i = 0; i < Nx; i++) outbufer << i << " ";
    outbufer << "\n";
    outbufer << "Y_COORDINATES " << Ny << " double\n";
    for (int j = 0; j < Ny; j++) outbufer << j << " ";
    outbufer << "\n";
    outbufer << "Z_COORDINATES " << Nz << " double\n";
    for (int k = 0; k < Nz; k++) outbufer << k << " ";
    outbufer << "\n";
    outbufer << "POINT_DATA " << Nx*Ny*Nz << "\n";

    for (int dir = 0; dir < 3; ++dir)
    {
        outbufer << "SCALARS X_" << dir << " double 1\n";
        outbufer << "LOOKUP_TABLE default\n";

        for (int k = 0; k < Nz; k++)
        for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++)
        {
            outbufer << (Quaternions(i,j,k).RotationMatrix*X)[dir] <<"\n";
        }
    }
    outbufer << "VECTORS X double\n";

    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        outbufer << (Quaternions(i,j,k).RotationMatrix*X)[0] << " "
                 << (Quaternions(i,j,k).RotationMatrix*X)[1] << " "
                 << (Quaternions(i,j,k).RotationMatrix*X)[2] << "\n";
    }

    string FileName = UserInterface::MakeFileName(VTKDir,"Rotated100Vector_",
                                                  tStep, ".vtk");

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbufer.rdbuf();
    vtk_file.close();
}

void Orientations::Advect(AdvectionHR& Adv,
                          const Velocities& Vel,
                          PhaseField& Phi,
                          const BoundaryConditions& BC,
                          const double dt, const double tStep)
{
    Adv.AdvectField(Quaternions, Vel, BC, dx, dt, tStep);
}


void Orientations::PrintPointStatistics(int x, int y, int z)
{
    cout << "Point:      (" << x << ", " << y << ", " << z << ")" << endl;
    cout << "Quaternions:\n" << Quaternions(x,y,z).print() << endl;

    cout << endl;
}

}// namespace openphase
