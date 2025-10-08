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

#ifndef ORIENTATIONS_H
#define ORIENTATIONS_H

#include "Base/Includes.h"
#include "PhaseField.h"

namespace openphase
{

class Settings;
class BoundaryConditions;
class PhaseField;
class Quaternion;
class Crystallography;
class Velocities;
class AdvectionHR;

class OP_EXPORTS Orientations : public OPObject                                            ///< Module which stores and handles orientation data like rotation matrices, Euler angles etc.
{

 public:

    Orientations(){};
    Orientations(Settings& locSettings);

    void Initialize(Settings& locSettings) override;                            ///< Initializes the module, allocate the storage, assign internal variables
    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC) override;///< Remesh and reallocate orientations

    void SetBoundaryConditions(const BoundaryConditions& BC);                   ///< Sets boundary conditions
    void SetRandomGrainOrientations(PhaseField& Phase, const int seed = 0);     ///< Sets random grains orientations
    void WriteVTK(const int tStep, const Settings& locSettings,
                  const int precision=16);                                      ///< Writes rotations to the VTK file
    void WriteTotalVTK(const PhaseField& Phase, const int tStep,
                       const Settings& locSettings,
                       const int precision = 16) const;                         ///< Writes total rotations (including grains orientations) to the VTK file
    void Write(const int tStep) const;                                          ///< Writes orientations storage to binary file
    void Read(const int tStep);                                                 ///< Reads orientations from binary file
    void WriteGrainEBSDDataQuaternions(const PhaseField& Phase, const int tStep);///< Writes EBSD file (for MTex). In interface, majority phase is used.
    void WriteMisorientationsVTK(const int tStep,
                                 const Settings& locSettings,
                                 const int precision,
                                 const std::string measure = "deg") const;      ///< Writes misorientations in VTK format
    void WriteRotated100Vector(const int tStep);
    static double getMisorientation(const dMatrix3x3 RotMatA, const dMatrix3x3 RotMatB); ///< Calculates missorientation between two matrices without consideration of symmetries
    static double getMisorientationCubic(const dMatrix3x3 RotMatA, const dMatrix3x3 RotMatB, const Crystallography& CR); /// Calculates missorientation between two matrices considering cubic symmetry
    static double getDisorientation(const Quaternion OrientationA, const Quaternion OrientationB); /// Calculates Disorientation between two quaternions considering cubic symmetry
    static EulerAngles RotationToEuler(const dMatrix3x3& Rot, const EulerConvention locConvention);
    dVector3 MillerConversion(const std::vector<double>& hkil, dVector3 hkl, bool PlaneNormal = true); // Hexagonal Miller plane Normal indices to Miller-Bravais // flase for directions --

    void PrintPointStatistics(int x, int y, int z);
    void Advect(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi,
                const BoundaryConditions& BC,
                const double dt, const double tStep) override;                  ///< Advects orientations
    int Nx;
    int Ny;
    int Nz;
    int dNx;
    int dNy;
    int dNz;

    double dx;

    int Nphases;

    Storage3D<Quaternion, 0>    Quaternions;
    Storage3D<Quaternion, 0>    QuaternionsDot;

    std::string VTKDir;                                                         ///< Directory-path added in front of VTK files
    std::string RawDataDir;                                                     ///< Directory-path added in front of Restart files

    Orientations& operator= (const Orientations& rhs);

    dMatrix3x3 getEffectiveGrainRotation(const PhaseField& Phase,
                                    const int i, const int j, const int k) const///< Returns averaged grain rotations, averaging is done via quaternions
    {
        Quaternion locQuaternion;

        if(Phase.Fields(i,j,k).flag < 2)
        {
            int index = Phase.Fields(i,j,k).cbegin()->index;
            locQuaternion = Phase.FieldsStatistics[index].Orientation;
        }
        else
        {
            locQuaternion.set(0,0,0,0);

            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); alpha++)
            {
                locQuaternion += Phase.FieldsStatistics[alpha->index].Orientation*alpha->value;
            }
        }

        return locQuaternion.RotationMatrix;
    }

    dMatrix3x3 getTotalRotation(const PhaseField& Phase,
                                    const int i, const int j, const int k) const///< Returns the local rotation including grain orientation and deformation induced rotation
    {
        return Quaternions(i,j,k).RotationMatrix*getEffectiveGrainRotation(Phase, i,j,k);
    }

    dMatrix3x3 getTotalGrainRotation(PhaseField& Phase, int i, int j, int k,
                                                                      int alpha)///< Returns the total rotation seen by the specified grain including deformation induced rotation
    {
        return Quaternions(i,j,k).RotationMatrix*Phase.FieldsStatistics[alpha].Orientation.RotationMatrix;
    }

 protected:
 private:
};
}// namespace openphase
#endif
