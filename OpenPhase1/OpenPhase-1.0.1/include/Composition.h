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
 *   File created :   2012
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Matthias Stratmann
 *
 */

#ifndef COMPOSITION_H
#define COMPOSITION_H

#include "Base/Includes.h"
#include "Chemistry/PeriodicTable.h"
#include "H5Interface.h"

namespace openphase
{
class PhaseField;
class BoundaryConditions;
class ElasticProperties;
class Velocities;
class PhaseField;
class Settings;

struct CompositionDelta
{
    double in;
    double out;

    CompositionDelta& operator*=(double m)
    {
        in  *= m;
        out *= m;
        return *this;
    }
    CompositionDelta& operator+=(CompositionDelta m)
    {
        in  += m.in;
        out += m.out;
        return *this;
    }
    CompositionDelta& operator-=(CompositionDelta m)
    {
        in  -= m.in;
        out -= m.out;
        return *this;
    }
    void pack(std::vector<double>& buffer)
    {
        buffer.push_back(in);
        buffer.push_back(out);
    }
    void unpack(std::vector<double>& buffer, size_t& it)
    {
        in = buffer[it]; ++it;
        out = buffer[it]; ++it;
    }
};

class OP_EXPORTS Composition : public OPObject                                             ///< Stores the composition as concentrations or molar densities or ... etc.
{
 public:
    Composition(){};
    Composition(Settings& locSettings,
                const std::string InputFileName = DefaultInputFileName)         ///< Initializes storages, sets internal variables.
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }
    void Initialize(Settings& locSettings) override;                                     ///< Initializes storages, sets internal variables.
    void ReadInput(const std::string InputFileName) override;                            ///<  Reads input parameters from a file
    void ReadInput(std::stringstream& inp) override;                                     ///<  Reads input parameters from a file

    void Remesh(int newNx, int newNy, int newNz,
                const BoundaryConditions& BC) override;                         ///< Remesh the storage while keeping the data

    void SetInitialMoleFractions(PhaseField& Phi, int mode = 0);
    void CalculateTotalMoleFractions(PhaseField& Phase);                        ///< Calculates total mole fractions from the phase mole fractions
    void WriteVTK(const int tStep, const Settings& locSettings,
            const int precision=16) const;                                      ///< Writes composition in VTK format (.vts file)

    void WriteDistortedVTK(int tStep, const Settings& locSettings, const ElasticProperties& EP);     ///< Writes composition in VTK format on a distorted grid (.vtk file)
    void Write(const int tStep);                   ///< Writes the raw composition into a file
    bool Read(BoundaryConditions& BC, const int tStep);    ///< Reads the raw composition from a file
    void WriteH5(int tStep, H5Interface& H5);
    bool ReadH5(int tStep, H5Interface& H5);

    void ReadSubset(BoundaryConditions& BC,
                    std::string FileName, int offsetX,
                    int offsetY,
                    int offsetZ,
                    bool legacy_format = false,
                    int inpNx = 1,
                    int inpNy = 1,
                    int inpNz = 1);                                             ///< Read subset of a raw (binary) composition from the file of possibly different (larger) size
    void SetBoundaryConditions(const BoundaryConditions& BC);                   ///< Sets the boundary conditions
    void SetLimitsBoundaryConditions(const BoundaryConditions& BC);             ///< Sets the boundary conditions for limits
    void MoveFrame(const int dx, const int dy, const int dz,
                   const BoundaryConditions& BC);                               ///< Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in x, y or z direction correspondingly.
    void ConsumePlane(const int dx, const int dy, const int dz,
                      const int x, const int y, const int z,
                      const BoundaryConditions& BC);                            ///< Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in the direction to/from (x, y, z) point.
    void PrintPointStatistics(int x, int y, int z);                             ///< Prints to screen composition at a given point (x, y, z)
    void WriteStatistics(int tStep, double dt);                                 ///< Writes composition statistics into file, input: time step

    void Advect(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi,
        const BoundaryConditions& BC, const double dt, const double tStep) override;

    size_t Nphases;                                                             ///< Number of thermodynamic phases
    size_t Ncomp;                                                               ///< Number of chemical components

    int TotalNx;                                                                ///< X dimension of the system in MPI parallel mode
    int OffsetX;                                                                ///< X position of the current domain in MPI parallel mode
    int TotalNy;                                                                ///< Y dimension of the system in MPI parallel mode
    int OffsetY;                                                                ///< Y position of the current domain in MPI parallel mode
    int TotalNz;                                                                ///< Z dimension of the system in MPI parallel mode
    int OffsetZ;                                                                ///< Z position of the current domain in MPI parallel mode

    int Nx;                                                                     ///< X dimension of the system
    int Ny;                                                                     ///< Y dimension of the system
    int Nz;                                                                     ///< Z dimension of the system

    int dNx;                                                                    ///< Active X dimension
    int dNy;                                                                    ///< Active Y dimension
    int dNz;                                                                    ///< Active Z dimension
    double dx;                                                                  ///< Grid spacing
    double Threshold;                                                           ///< Minimum phase fraction for which diffusion fluxes are still calculated
    double AtomicWeightMixture;                                                 ///< Atomic weight of the Mixture or Composition 
    double SpecificHeatCapacityMix;                                             ///< Atomic weight of the Mixture or Composition

    std::vector<std::string> ElementNames;
    std::vector<std::string> PhaseNames;

    std::string VTKDir;                                                         ///< Directory-path added in front of VTK files
    std::string RawDataDir;                                                     ///< Directory-path added in front of Restart files
    std::string TextDir;

    /* template specialization integer number stands for
     * the number of extra dimensions
     * order:
     *
     * Rank = 3: phase, phase, component
     * Rank = 2: phase, component
     * Rank = 1: phase or component
     */

    Storage3D<double, 2> MoleFractions;                                         ///< Mole fractions for each phase
    Storage3D<double, 1> MoleFractionsTotal;                                    ///< Total mole fractions

    Storage3D<CompositionDelta, 2> MoleFractionsDot;                            ///< Phase composition increments storage
    Storage3D<double, 1> MoleFractionsTotalDot;                                 ///< Total composition increments storage
    Storage3D<CompositionDelta, 2> Norm;                                        ///< Storage for the normalization of phase composition
    Tensor<double, 2> Initial;                                                  ///< Initial composition of components in all phases
    Tensor<double, 2> MoleFractionsAverage;                                     ///< Average composition of components in all phases
    Tensor<double, 1> MoleFractionsTotalAverage;                                ///< Average composition of components in the simulation domain
    Tensor<double, 3> MoleFractionsInterfaceAverage;                            ///< Average composition of components in the interfaces
    Tensor<double, 2> MoleFractionsMIN;                                         ///< Minimum mole fraction of a given component in a given phase
    Tensor<double, 2> MoleFractionsMAX;                                         ///< Minimum mole fraction of a given component in a given phase

    Tensor<double, 1> Reference;                                                ///< Reference composition
    Storage3D<int, 0> Limiting;                                                 ///< Storage for the limiting of the composition increments

    std::vector<double> TotInitial;                                             ///< Initial amount of each component
    double TotalMolarVolume;                                                    ///< Total molar volume of the system

    Tensor<double, 1> WeightFractionsTotal(int x, int y, int z) const;          ///< Returns local weight fractions of all elements
    void CalculateMoleFractionsTotalAverage(void);
    void CalculateMoleFractionsAverage(PhaseField& Phase);
    void CalculateTotalMolarVolume(void);
    void CollectStatistics(PhaseField& Phase)
    {
        CalculateMoleFractionsTotalAverage();
        CalculateMoleFractionsAverage(Phase);
    }

    std::vector<double> TotalAverage;

    PeriodicTable PT;                                                           ///< Periodic table of elements

 protected:
    int AtStart;
 private:
    double GetMolarVolume(std::string Element);                                 ///< Returns the molar volume of a given element
    double GetAtomicWeight(std::string Element);                                ///< Returns the atomic weight of a given element
};

} // namespace openphase

#endif // COMPOSITION_H
