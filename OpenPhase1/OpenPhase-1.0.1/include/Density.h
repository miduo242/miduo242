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
 *   Main contributors :   Oleg Shchyglo;
 *
 */

#ifndef DENSITY_H
#define DENSITY_H

#include "Base/Includes.h"

namespace openphase
{
class PhaseField;
class BoundaryConditions;
class Composition;
class Temperature;
class Settings;

class Density : public OPObject                                                 ///< Stores the mass densities of phases in each point.
{
 public:
    Density(){};                                                                ///< Empty constructor.
    Density(Settings& locSettings,
            const std::string InputFileName = DefaultInputFileName);
    void Initialize(Settings& locSettings) override;                            ///< Initializes storages, sets internal variables.
    void ReadInput(const std::string FileName) override;                        ///< Reads input values from file
    void ReadInput(std::stringstream& inp);                        ///< Reads input values from file
    void Set(PhaseField& Phase, Composition& Cx, Temperature& Tx);              ///< Sets the density depending on composition and temperature
    void SetInitial(PhaseField& Phase);                                         ///< Sets initial density

    void WriteVTK(const int tStep, const Settings& locSettings,
                  const int precision = 16);                                    ///< Writes density in VTK format
    void Write(const int tStep);                                                ///< Writes the raw composition into a file
    void Read(const int tStep);                                                 ///< Reads the raw composition from a file
    void SetBoundaryConditions(BoundaryConditions& BC);                         ///< Sets the boundary conditions

    void MoveFrame(int di, int dj, int dk, BoundaryConditions& BC);             ///< Shifts the data on the storage by di, dj and dk in x, y and z directions correspondingly.
    void Remesh(int newNx, int newNy, int newNz, BoundaryConditions& BC);       ///< Remesh the storage while keeping the data
    void PrintPointStatistics(int x, int y, int z);                             ///< Prints to screen density in a given point (x, y, z)

    Density& operator=(const Density& rhs);                                     ///< Copy operator for Density class

    size_t Nphases;                                                             ///< Number of thermodynamic phases

    int TotalNx;                                                                ///< X dimension of the system in MPI parallel mode
    int OffsetX;                                                                ///< X position of the current domain in MPI parallel mode
    int TotalNy;                                                                ///< Y dimension of the system in MPI parallel mode
    int OffsetY;                                                                ///< Y position of the current domain in MPI parallel mode
    int TotalNz;                                                                ///< Z dimension of the system in MPI parallel mode
    int OffsetZ;                                                                ///< Z position of the current domain in MPI parallel mode

    int Nx;                                                                     ///< Size of the inner calculation domain along X
    int Ny;                                                                     ///< Size of the inner calculation domain along Y
    int Nz;                                                                     ///< Size of the inner calculation domain along Z

    int dNx;                                                                    ///< Active X dimension
    int dNy;                                                                    ///< Active Y dimension
    int dNz;                                                                    ///< Active Z dimension

    double dx;                                                                  ///< Grid spacing

    Storage3D<double, 1> Phase;                                                 ///< Phase densities storage
    Storage3D<double, 0> Total;                                                 ///< Total density storage

    Tensor<double, 1> Initial;                                                  ///< Initial density of all phases

    std::string VTKDir;                                                         ///< Directory-path added in front of VTK files
    std::string RawDataDir;                                                     ///< Directory-path added in front of Restart files

 protected:
 private:
};

} // namespace openphase

#endif // DENSITY_H
