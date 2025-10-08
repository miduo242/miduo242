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
 *   File created :   2010
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */

#ifndef PARABOLICDIFFUSION_H
#define PARABOLICDIFFUSION_H

#include <vector>
#include <string>
#include "Base/Tensor.h"
#include "Base/Storage3D.h"

namespace openphase
{

class Settings;
class PhaseField;
class DrivingForce;
class Composition;
class Temperature;
class BoundaryConditions;

class ParabolicDiffusion: public OPObject
{
 public:
    ParabolicDiffusion(){}
    ParabolicDiffusion(Settings& locSettings,
            std::string InputFileName = DefaultInputFileName);                  ///< Calls Initialize and ReadInput
    void Initialize(Settings& locSettings) override;                            ///< Initializes global settings
    void ReadInput(const std::string InputFileName) override;                   ///< Reads input parameters from a file
    void GetDrivingForce(PhaseField& Phase, Composition& Cx, Temperature& Tx,
            DrivingForce& dGab);                                                ///< Calculates the driving force for each point
    void CalculateMobility(PhaseField& Phasem, size_t comp);                    ///< Sets chemical mobility for each component in each point
    void Solve(PhaseField& Phase, Composition& Cx, Temperature& Tx,
            BoundaryConditions& BC, double dt);                                 ///< Calculates the change of total concentrations in one time step taking into account cross terms
    void ApplyIncrements(PhaseField& Phase, Composition& Elements, double dt,
            size_t comp);                                                       ///< Calculates and applies composition increments
    void ApplyIncrementsConstantIsotropicMobility(PhaseField& Phase,
            Composition& Elements, double dt, size_t comp);                     ///< Calculates and applies composition increments
    void CalculateChemicalPotential(PhaseField& Phase, Composition& Elements,
            Temperature& Tx);                                                   ///< Calculates the chemical potential of each component
    void CaluculateDiffusionFlux(size_t comp);                                  ///< Calculates the flux of each component
    void CalculateReferenceElementMoleFractions(Composition& Elements);         ///< Calculates mole fractions of the reference composition
    void SetPhaseMolarFractions(Composition& Elements);                         ///< Calculates phase molar factions
    void ResetChemicalPotential();                                              ///< Sets chemical potential to zero so that external contributions can be added
    void WriteVTK(int tStep, const Settings& locSettings,
            int precision = 16) const;                                          ///< Writes chemical potential and flux into VTK file

    int TotalNx;                                                                ///< Size of the inner calculation domain along X
    int Nx;                                                                     ///< Size of the inner calculation domain along X
    int Ny;                                                                     ///< Size of the inner calculation domain along Y
    int Nz;                                                                     ///< Size of the inner calculation domain along Z
    int dNx;                                                                    ///< Active X dimension
    int dNy;                                                                    ///< Active Y dimension
    int dNz;                                                                    ///< Active Z dimension
    int Bcells;                                                                 ///< Number of boundary cells
    size_t Ncomp;                                                               ///< Number of components
    size_t Nphases;                                                             ///< Number of thermodynamic phases
    double dx;                                                                  ///< Grid spacing
    size_t RefComp;                                                             ///< Index of the reference chemical component
    bool UseLaplacian;                                                          ///< Uses Laplacian stencil to solve diffusion equation if possible
    LaplacianStencil  LStencil;                                                 ///< Laplace stencil
    GradientStencil   GStencil;                                                 ///< Gradient stencil
    //DivergenceStencil DStencil;                                                 ///< Divergence stencil

    std::vector<std::string> ElementNames;                                      ///< Names of corresponding chemical components (e.g Au, Cu, Na, Cl etc.)
    std::vector<std::string> PhaseNames;                                        ///< Names of corresponding thermodynamic phases
    std::string VTKDir;                                                         ///< Directory-path added in front of VTK files
    std::string RawDataDir;                                                     ///< Directory-path added in front of Restart files
    std::string TextDir;                                                        ///< Directory-path added in front of text files

    Tensor <double, 2> C0;                                                      ///< Equilibrium/Solvus concentration
    Tensor <double, 2> EPS2;                                                    ///< Energy coefficient [J/m^3] for component in each phase (NOTE: EPS=RT/Vm)
    Tensor <double, 2> EPS1;                                                    ///< Energy coefficient [J/m^3] for component in each phase
    Tensor <double, 2> EPS0;                                                    ///< Energy coefficient [J/m^3] for component in each phase
    Tensor <double, 2> PhaseDiffusionCoefficient;                               ///< Solute diffusion coefficients for each thermodynamic phase
    Tensor <double, 2> PhaseMobility;                                           ///< Solute mobility coefficients for each thermodynamic phase
    Tensor <double, 1> ConstantIsotropicMobility;                               ///< Non-zero if mobility is constant and isotropic

    Storage3D <double,   1> ChemicalPotentialOld;                               ///< Chemical Potential of last time step (used for vtk file)
    Storage3D <double,   1> ChemicalPotential;                                  ///< Chemical Potential of each component in each point
    Storage3D <double,   1> Mobility;                                           ///< Chemical mobility of each component in each point
    Storage3D <dVector3, 1> Flux;                                               ///< Mobility coefficients of each component in each point
 protected:
 private:
};

} // namespace openphase
#endif
