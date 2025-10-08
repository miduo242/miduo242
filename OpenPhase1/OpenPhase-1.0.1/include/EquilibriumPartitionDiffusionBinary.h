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

#ifndef EQUILIBRIUMPARTITIONDIFFUSIONBINARY_H
#define EQUILIBRIUMPARTITIONDIFFUSIONBINARY_H

#include "Base/Includes.h"

namespace openphase
{

class Settings;
class PhaseField;
class DrivingForce;
class Composition;
class Temperature;
class InterfaceProperties;
class BoundaryConditions;

class EquilibriumPartitionDiffusionBinary : public OPObject                     ///<  Solver for diffusion using binary linearized phase diagrams
{
 public:
    EquilibriumPartitionDiffusionBinary(){}
    EquilibriumPartitionDiffusionBinary(Settings& locSettings,
            const std::string InputFileName = DefaultInputFileName)             ///<  Initializes global settings
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    };
    void Initialize(Settings& locSettings) override;                            ///<  Initializes global settings
    void ReadInput(const std::string InputFileName) override;                   ///<  Reads input parameters from a file
	void ReadInput(std::stringstream& inp) override;                            ///<  Reads input parameters from a file
    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC) override; ///<  Remesh the storage while keeping the data

    void GetDrivingForce(PhaseField& Phase, Composition& Cx, Temperature& Tx,
                                                            DrivingForce& dGab);///<  Calculates the driving force for each point
    double GetDrivingForceAlpha(PhaseField& Phase, Composition& Elements,
                                                   Temperature& Tx, size_t alpha,
                                                         int i, int j, int k);  ///<  Calculates the driving force for phase alpha in a given point

    void CalculateLocalPhaseConcentrations(PhaseField& Phase, Temperature& Tx,
                                       Composition& Cx, int i, int j, int k);   ///<  Distributes the total concentrations into concentrations in each phase for a given point
    void CalculatePhaseConcentrations(PhaseField& Phase,
                                              Temperature& Tx, Composition& Cx);///<  Distributes the total concentrations into concentrations in each phase
    void SetDiffusionCoefficients(PhaseField& Phase, Temperature& Tx);          ///<  Sets diffusion coefficients in each point
    void Solve(PhaseField& Phase, Composition& Cx,
                            Temperature& Tx, BoundaryConditions& BC, double dt,
                            bool InternalDiffusivities = true);                 ///<  Calculates the change of total concentrations in one time step taking into account cross terms
    void CalculateInterfaceMobility(PhaseField& Phase, Composition& Cx, Temperature& Tx,
                              BoundaryConditions& BC, InterfaceProperties& IProperties,
                                             bool InternalDiffusivities = true);///<  Calculates concentration-dependent mobility
    void RestoreStoichiometric(PhaseField& Phase, Temperature& Tx,
                                                            Composition& Cx);   ///<  Recovers the exact stoichiometric composition in the fully grown stoichiometric cell
    void RestoreStoichiometricThreadSafe(PhaseField& Phase, Temperature& Tx,
                                                            Composition& Cx);   ///<  Recovers the exact stoichiometric composition in the fully grown stoichiometric cell
    void CalculateDiffusionIncrements(PhaseField& Phase, Composition& Elements,
                                                               Temperature& Tx);///<  Calculates Fick's diffusion concentration increments
    void CalculateAntitrappingIncrements(PhaseField& Phase,
                                                      Composition& Elements);   ///<  Calculates antitrapping concentration increments
    void LimitAntitrappingIncrements(PhaseField& Phase,
                                     Composition& Elements, Temperature& Tx);   ///<  Limits antitrapping concentration increments
    void LimitAntitrappingIncrementsThreadSafe(PhaseField& Phase,
                                     Composition& Elements, Temperature& Tx);   ///<  Limits antitrapping concentration increments
    void CalculateLimits(PhaseField& Phase, Composition& Elements, double dt);  ///<  Calculated limits for increments in order to keep concentrations within physical range
    void LimitDiffusionIncrements(PhaseField& Phase,
                                     Composition& Elements, Temperature& Tx);   ///<  Limits Fick's diffusion concentration increments
    void ApplyIncrements(PhaseField& Phase, Composition& Elements, double dt);  ///<  Apply the limited increments

    void CalculateReferenceElementMoleFractions(Composition& Elements);         ///<  Calculates mole fractions of the reference chemical component

    double ReportMaximumTimeStep(void)                                          ///<  Returns maximum time step for diffusion
    {
        return 0.25*dx*dx/maxDC;
    }
    void PrintPointStatistics(int x, int y, int z);                             ///<  Prints diffusion coefficients in the point
    int Nx;                                                                     ///<  Size of the inner calculation domain along X
    int Ny;                                                                     ///<  Size of the inner calculation domain along Y
    int Nz;                                                                     ///<  Size of the inner calculation domain along Z
    int dNx;                                                                    ///<  Active X dimension
    int dNy;                                                                    ///<  Active Y dimension
    int dNz;                                                                    ///<  Active Z dimension
    bool LimitingNeeded;
    size_t RefComp;                                                             ///<  Index of the reference chemical component
    size_t Comp;                                                                ///<  Index of the chemical component for a given instance of the diffusion solver
    size_t Nphases;                                                             ///<  Number of thermodynamic phases
    size_t Ncomp;                                                               ///<  Number of chemical components

    std::vector<std::string> ElementNames;                                      ///<  Names of chemical elements

    double Precision;                                                           ///<  Precision of composition evaluation
    double Eta;                                                                 ///<  Interface width
    double dx;                                                                  ///<  Grid spacing
    double R;                                                                   ///<  Universal gas constant
    double TotalMass;                                                           ///<  Total amount of the component handled by a given instance of the solver
    double maxDC;                                                               ///<  Maximum diffusion coefficient in the simulation.

    std::vector<int> Stoichiometric;                                            ///<  Stoichiometry flags for each thermodynamic phase
    std::vector<double> DC0;                                                    ///<  Solute diffusion coefficients for each thermodynamic phase
    std::vector<double> AE;                                                     ///<  Diffusion activation energies for each thermodynamic phase
    Tensor<double, 2> Pt;                                                       ///<  Pairwise partition coefficients for a linear phase diagram
    Tensor<double, 2> mL;                                                       ///<  Pairwise liquidus slopes for a linear phase diagram
    Tensor<double, 2> mL_1;                                                     ///<  Inverse of the pairwise liquidus slopes for a linear phase diagram
    Tensor<double, 2> Ts;                                                       ///<  Pairwise temperatures of the liquidus-solidus or solvus lines intersections for a linear phase diagram
    Tensor<double, 2> Cs;                                                       ///<  Pairwise concentrations of the liquidus-solidus or solvus lines intersections for a linear phase diagram

    std::vector<double> S;                                                      ///<  Entropies of different phases. Their differences are used in the driving force

    int ThereAreStoichiometricPhases;                                           ///<  Indicates if there are stoichiometric phases
    Storage3D<double, 1> dMu;                                                   ///<  External contribution to the chemical potential.
    Storage3D<double, 1> DC;                                                    ///<  Diffusion coefficients in each point

    LaplacianStencil DStencil;                                                  ///<  Diffusion stencil. Uses Laplacian stencil as the basis
    bool EnableAntiTrapping;                                                    ///<  True if antitrapping current is enabled, false othewise
 protected:
 private:
};

} // namespace openphase
#endif
