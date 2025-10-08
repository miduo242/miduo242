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
 *   File created :   2021
 *   Main contributors :   Raphael Schiedung
 *
 */

#ifndef GRANDCANONICALDIFFUSION_H
#define GRANDCANONICALDIFFUSION_H

#include "Base/Storage3D.h"
#include "Base/dVector3.h"

namespace openphase
{
class BoundaryConditions;
class InterfaceProperties;
class PhaseField;
class Settings;

namespace GrandPotential
{
class Density;
class Solver
{
 public:
    Solver(){}
    Solver(const Settings& locSettings);                                        ///<  Initializes global settings
    Solver(const Settings& locSettings, const std::string InputFileName);       ///<  Initializes global settings

    void Initialize(const Settings& locSettings);                               ///<  Initializes global settings
    void ReadInput(const std::string InputFileName);                            ///<  Reads input parameters from a file

    // Restart:
    void Write(const char* filename) const;                                     ///< Writes the raw composition into a file
    void Write(const std::string& filename) const;                              ///< Writes the raw composition longo a file
    void Write(long tStep) const;                                               ///< Writes the raw composition longo a file named with class name and time step
    void Write() const;                                                         ///< Writes the raw composition longo a file named with class name
    void Read(const BoundaryConditions& BC, const PhaseField& Phase, const Density& omega, const char* filename);        ///< Reads the raw composition from a file
    void Read(const BoundaryConditions& BC, const PhaseField& Phase, const Density& omega, const std::string& filename); ///< Reads the raw composition from a file
    void Read(const BoundaryConditions& BC, const PhaseField& Phase, const Density& omega, long tStep);                  ///< Reads the raw composition from a file named with class name and time step
    void Read(const BoundaryConditions& BC, const PhaseField& Phase, const Density& omega);                              ///< Reads the raw composition from a file named with class name

    // Algorithm:
    void SetInitialConcentration (const PhaseField& Phase, const Density& omega); ///< Sets initial concentrations and chemical potentials
    void SetInitialPressure      (const PhaseField& Phase, const Density& omega); ///< Sets initial concentrations and chemical potentials
    void SetInitial (const PhaseField& Phase, const Density& omega, const BoundaryConditions& BC, double Temp = -DBL_EPSILON); ///< Sets initial concentrations and chemical potentials
    void Solve      (      PhaseField& Phase, const Density& omega, const BoundaryConditions& BC, const InterfaceProperties& IP, double dt, double Temp = -DBL_EPSILON); ///< Solves diffusion equation

    // Fields:
    Storage3D< double, 1> ChemicalPotential;                                    ///< Chemical potential of component in J/mol
    double PhasePressure  (long i, long j, long k, size_t PhaseIdx, const PhaseField& Phase, const Density& omega) const; ///< Calculates local pressure
    double MolarVolume    (long i, long j, long k, const PhaseField& Phase, const Density& omega) const; ///< Volume of component comp at (i,j,k) in m^3
    double Pressure       (long i, long j, long k, const PhaseField& Phase, const Density& omega) const; ///< Calculates local pressure
    double MassDensity    (long i, long j, long k, const PhaseField& Phase, const Density& omega) const; ///< Calculates local mass density
    double MoleFraction   (long i, long j, long k, size_t comp, const PhaseField& Phase, const Density& omega) const; ///< Fraction of component comp at (i,j,k) in %
    double Concentration  (long i, long j, long k, size_t comp, const PhaseField& Phase, const Density& omega) const; ///< Calculates local concentration
    double Susceptibility (long i, long j, long k, size_t comp, const PhaseField& Phase, const Density& omega) const; ///< Calculates local susceptibility

    // VTK Output:
    void WriteVTKChemicalPotential     (long tStep, const Settings& locSettings, long precision=16) const; ///< Writes chemical Potential in VTK format (.vts file)
    void WriteVTKDiffusionFlux         (long tStep, const Settings& locSettings, long precision=16) const; ///< Writes diffusion flux in VTK format (.vts file) 
    void WriteVTKMobility              (long tStep, const Settings& locSettings, long precision=16) const; ///< Writes mobility in VTK format (.vts file)
    void WriteVTKVelocity              (long tStep, const Settings& locSettings, long precision=16) const; ///< Writes velocity of element in VTK format (.vts file)
    void WriteVTK                      (long tStep, const Settings& locSettings, const PhaseField& Phase, const Density& omega, long precision=16) const; ///< Writes every possible output longo a VTK format file
    void WriteVTKConcentration         (long tStep, const Settings& locSettings, const PhaseField& Phase, const Density& omega, long precision=16) const; ///< Writes chemical composition in VTK format (.vts file)
    void WriteVTKDiffusivity           (long tStep, const Settings& locSettings, const PhaseField& Phase, const Density& omega, long precision=16) const; ///< Writes molar volume in VTK format (.vts file)
    void WriteVTKMoleFraction          (long tStep, const Settings& locSettings, const PhaseField& Phase, const Density& omega, long precision=16) const; ///< Writes mole fraction in VTK format (.vts file)
    void WriteVTKPressure              (long tStep, const Settings& locSettings, const PhaseField& Phase, const Density& omega, long precision=16) const; ///< Writes pressure in VTK format (.vts file)
    void WriteVTKSusceptibility        (long tStep, const Settings& locSettings, const PhaseField& Phase, const Density& omega, long precision=16) const; ///< Writes molar volume in VTK format (.vts file)
    void WriteVTKMolarVolume           (long tStep, const Settings& locSettings, const PhaseField& Phase, const Density& omega, long precision=16) const; ///< Writes molar volume in VTK format (.vts file)

    // Diagnostics and helper functions:
    double TotalAmountOfComponent(size_t comp) const;                           ///< Calculates total amount of comp
    double GrandPotential(const PhaseField& Phase, const Density& omega) const; ///< Calculates total Grand Potential
    double EquilibriumPhaseConcentration(size_t PhaseIdx, size_t comp) const;   ///< Returns the equilibrium concentration of component comp in phase PhaseIdx (Without Curvature Effect)

 protected:
 private:
    void SetTemperature(double Temp);                                           ///< Updates material parameters to new temp
    void CaclulateLocalDiffusionFlux       (long i, long j, long k);            ///< Calculates local diffusion flux
    void MergeLocalIncrements              (long i, long j, long k, double dt); ///< Merge local chemical potential increments
    void MergeLocalIncrements1             (long i, long j, long k, double dt); ///< Merge local chemical potential increments
    void MergeLocalIncrements2             (long i, long j, long k, double dt); ///< Merge local chemical potential increments
    void CalculateLocalConcentrations      (long i, long j, long k,              const PhaseField& Phase, const Density& omega); ///< Calculates local concentration and phase concentration
    void CalculateLocalConcentrations      (long i, long j, long k, size_t comp, const PhaseField& Phase, const Density& omega); ///< calculates local concentration and phase concentration
    void CalculateLocalMobilities          (long i, long j, long k, const PhaseField& Phase); ///< Calculates local mobility and phase phase mobilities
    void CalculateLocalIncrements          (long i, long j, long k, const PhaseField& Phase, const Density& omega, double dt); ///< Calculates local chemical potential increments
    void CalculateLocalIncrements1         (long i, long j, long k, const PhaseField& Phase, const Density& omega, double dt); ///< Calculates local chemical potential increments
    void CalculateLocalIncrements2         (long i, long j, long k, const PhaseField& Phase, const Density& omega, double dt); ///< Calculates local chemical potential increments
    void CalculateLocalPhaseFieldIncrements(long i, long j, long k, PhaseField& Phase, const Density& omega, const InterfaceProperties& IP); ///< Calculates local driving force
    void MergeLocalIncrements2Implicit     (long i, long j, long k, PhaseField& Phase, const Density& omega, const InterfaceProperties& IP, double dt); ///< Merge local chemical potential increments
    void SolveExplicit(PhaseField& Phase, const Density& omega, const BoundaryConditions& BC, const InterfaceProperties& IP, double dt, double Temp = -DBL_EPSILON); ///< Solves diffusion equation
    void SolveImplicit(PhaseField& Phase, const Density& omega, const BoundaryConditions& BC, const InterfaceProperties& IP, double dt, double Temp = -DBL_EPSILON); ///< Solves diffusion equation
    void EnforceConservationOfTOC(const PhaseField& Phase, const Density& omega);

    const std::string thisclassname = "GrandPotentialSolver";

    double dx;                                                                  ///< Grid spacing
    long TotalNx;                                                               ///< Size of the inner calculation domain along X
    long Nx;                                                                    ///< Size of the inner calculation domain along X
    long Ny;                                                                    ///< Size of the inner calculation domain along Y
    long Nz;                                                                    ///< Size of the inner calculation domain along Z
    int dNx;                                                                    ///< Active X dimension
    int dNy;                                                                    ///< Active Y dimension
    int dNz;                                                                    ///< Active Z dimension
    int boundary;                                                               ///< Number of boundary cells
    size_t Ncomp;                                                               ///< Number of chemical components
    size_t Nphases;                                                             ///< Number of thermodynamic phases

    std::vector<double>      ElementMasses;                                     ///< Names of corresponding chemical components (e.g Au, Cu, Na, Cl etc.)
    std::vector<std::string> ElementNames;                                      ///< Names of corresponding chemical components (e.g Au, Cu, Na, Cl etc.)
    std::vector<std::string> PhaseNames;                                        ///< Names of corresponding thermodynamic phases

    std::string VTKDir;                                                         ///< Directory-path added in front of VTK files
    std::string RawDataDir;                                                     ///< Directory-path added in front of Restart files
    std::string TextDir;

    Storage3D< double,   1> Concentrations;                                     ///< concentration of component in mol/m^3
    Storage3D< double,   1> ChemicalPotentialOld;                               ///< Old Chemical potential before time step
    Storage3D< double,   1> ChemicalPotentialDot;                               ///< Change of chemical with time
    Storage3D< double,   1> ChemicalPotentialDot2;                              ///< Change of chemical with time
    Storage3D< dVector3, 1> DiffusionFlux;                                      ///< Stores local diffusion flux

    double InitialPressure;
    Tensor<double, 2> CInitial;                                                 ///< Initial molar concentration
    bool ConserveTOC;                                                           ///< Enforce conservation of  total amount of components
    bool UseInitialPressure;                                                    ///< Set initial constants pressure
    bool UseImplicitSolver;                                                     ///< True if implicit Euler method is used
    double ChemicalPotentialAccuracy;                                           ///< The implicit solver will iterate until the residual is smaller than MaxResidual;
    double TOCAccuracy;                                                         ///< Accuracy of conservation of total amount of components
    double Temperature;                                                         ///< Global Temperature
    size_t MaxIterations;                                                       ///< Maximum number of implicit solver iterations
    size_t TOCMaxIterations;                                                    ///< Maximum number of iterations to conserve total amount of components
    std::vector<double> TOC0;                                                   ///< Total Amount of Component
    std::vector<double> InitialChemicalPotential;                               ///< Initial chemical potential

    // Mobility related storages:
    Storage3D< double, 1> Mobilities;                                           ///< Stores local mobility
    bool Use_InterfaceMobilities;                                               ///< True if interface mobilities is considered                                          
    bool Use_dPhaseMobility_dConcentration;                                     ///< True if phase mobilities depend on concentration
    //bool InterPhaseMobilityConcentrationCoupling;                             ///< True if interface mobilities depend on concentration
    Tensor<double, 2> PhaseMobilities;                                          ///< Bulk mobilities
    Tensor<double, 3> dPhaseMobilities_dConcentration;                          ///< Derivative of phase mobilities with respect to concentration
    Tensor<double, 3> InterfaceMobilities;                                      ///< Interface mobilities
    //Tensor<double, 4> dInterfaceMobilities_dConcentration;                    ///< Derivative of interface mobilities with respect to concentration

    template<typename field_t>
    double CalculateVolumeIntegral(field_t field, double dx) const
    {
        double local_value  = 0.0;
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,0,reduction(+:local_value))
        {
            local_value += field(i,j,k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        local_value*=dx*dx*dx;
        #ifdef MPI_PARALLEL
            double global_value = 0.0;
            MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            return global_value;
        #else
            return local_value;
        #endif
    }
};
} // namespace GrandPotential
} // namespace openphase
#endif
