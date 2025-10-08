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
 *   File created :   2011
 *   Main contributors :   Efim Borukhovich; Oleg Shchyglo
 *
 */

#ifndef ELECTRICALPOTENTIAL_H
#define ELECTRICALPOTENTIAL_H

#include "fftw3.h"
#include "Base/Includes.h"

namespace openphase
{

class Settings;
class PhaseField;
class Composition;
class BoundaryConditions;

struct electrode
{
    int                 index;
    double              charge;                                                 ///< charge of the electrode [Coulomb]
    int                 position[3];                                            ///< position of the electrode given in index numbers
    int                 counterElectrode;                                       ///< index of the counterelectrode
    double              interfaceVolume;
};

class ElectricalPotential : public OPObject
{
 public:

    ElectricalPotential(){};
    ElectricalPotential(Settings& locSettings,
                        const std::string InputFileName = DefaultInputFileName);

    void Initialize(Settings& locSettings) override;                            /// allocating memory and reading parameters from the Settings object
    void ReadInput(const std::string InputFileName) override;                   /// reading the input parameters
    void ReadInput(std::stringstream& inp);                   /// reading the input parameters
    void QXYZ(void);                                                            /// generates the wave vector, copied from the spectral solver
    void Solve(PhaseField& Phase, Composition& Cx, BoundaryConditions& BC);     /// solving the Poisson equation
//        void CalculateChargeDensity(MolarDensity& Nu, Diffusion& DF);         /// calculating the charge density
    void SetElectrodeCharges(PhaseField& Phase, Composition& Cx);
    inline double ChargeDensity(int i, int j, int k, const PhaseField& Phase,
                                                     const Composition& Cx);    /// assumes the electrode charge sitting in the middle of the electrode
    void CalculateElectrodeInterfaceVolumes(const PhaseField& Phase);           /// calculates the volume of the electrodes' interfaces
    void WriteChargeDensityVTK(const int tStep, const Settings& locSettings,
                               const PhaseField& Phase, const Composition& Cx); /// Writes the charge density into a VTK file
    void WritePotentialVTK(const int tStep, const Settings& locSettings) const; /// Writes the electrical potential into a VTK file

    Storage3D<double, 0>   Potential;

    std::vector<electrode>   electrodes;
    int                      Nelectrodes;

    double              Faraday;                                                /// Faraday constant: elementarcharge times Avogadro constant
    double              EpsInv;                                                 /// inverse of the (electrostatic) permittivity in vacuum

    std::string VTKDir;                                                         ///< Directory-path added in front of VTK files
    std::string RawDataDir;                                                     ///< Directory-path added in front of Restart files

 protected:
//        Storage3D<double>   Rho;                                              /// Storage for the charge density

    ///objects for fourier methode
    double*                 Q[3];                                               /// Wave vector
    double*                 RHS;                                                /// right hand side of the Poisson equation
    std::complex<double>*   ftRHS;                                              /// fourier transormed RHS
    double*                 rlPotential;                                        /// fourier transformed electrical potential
    std::complex<double>*   ftPotential;                                        /// fourier transformed electrical potential
    fftw_plan               ForwardPlan;                                        /// Forward FFT Plans
    fftw_plan               BackwardPlan;                                       /// Backward FFT Plans

    //read from the input file by readinput()
    std::vector<int>        ElementarCharges;
    std::vector<double>     MolarCharge;                                        /// Molar charges for the components

    //read from an OPSettings object by initiate()
    int                 Nx;
    int                 Ny;
    int                 Nz;

    int                 dNx;
    int                 dNy;
    int                 dNz;

    int                 Nz2;
    int                 Size, ftSize;
    double              DPi_nX, DPi_nY,DPi_nZ;

    double              dx;
    int                 Ncomp;
    int                 Nphases;

    double              CathodeCharge;
    double              AnodeCharge;

 private:
};

}
#endif // ELECTRICALPOTENTIAL_H
