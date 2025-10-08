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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Johannes Goerler
 *
 */

#ifndef ELASTICPROPERTIES_H
#define ELASTICPROPERTIES_H

#include "Base/Includes.h"
#include "Base/NodeV3.h"
#include "Mechanics/ElasticityModels/ElasticityKhachaturyan.h"
#include "Mechanics/ElasticityModels/ElasticitySteinbach.h"
#include "SymmetryVariants.h"
#include "Tools.h"

namespace openphase
{

class Advection;
class AdvectionHR;
class BoundaryConditions;
class Composition;
class DrivingForce;
class EquilibriumPartitionDiffusionBinary;
class InterfaceProperties;
class NodeV;
class Orientations;
class PhaseField;
class Settings;
class Temperature;
class Velocities;

enum class ElasticityModels                                                     ///< Elasticity homogenization modes and materials models
{
    Khachaturyan,                                                               ///< Khachaturyan's material model
    Steinbach,                                                                  ///< Steinbach's material model
    Voigt,                                                                      ///< Voigt/Taylor homogenization mode
    Reuss,                                                                      ///< Reuss/Sachs homogenization mode
    Rank1
};

struct NeuberParameters                                                         ///< Nueber plasticity correction parameters
{
    double YoungsModulus;
    double YieldStrength;
    double Hardening;
    bool   Active;                                                              ///< Neuber correction active/inactive indicator
};

class OP_EXPORTS ElasticProperties : public OPObject                            ///< Module which stores and handles elastic properties
{
 public:
    ElasticProperties(){};                                                      ///< Constructor (does nothing)
    ElasticProperties(Settings& locSettings,
                      const std::string InputFileName = DefaultInputFileName);  ///< Constructs and initializes the elastic properties

    void Initialize(Settings& locSettings) override;                            ///< Allocates storages, initializes internal parameters
    void ReadInput(const std::string InputFileName) override;                   ///< Reads elastic properties from the input file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads elastic properties from the input file
    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC) override;///< Remeshes the elasticity storages

    void SetBoundaryConditions(const BoundaryConditions& BC);                   ///< Sets boundary conditions for strains, stresses and displacements
    void SetGrainsProperties(PhaseField& Phase);                                ///< Sets elastic properties for each grain according to its orientation

    void SetEffectiveTransformationStretches(PhaseField& Phase);                  ///< Sets initial effective transformation stretches
    void SetEffectiveTransformationStretches(PhaseField& Phase, InterfaceProperties& IP); ///< Sets initial effective transformation stretches considering Vegard's contribution
    void SetEffectiveTransformationStretches(PhaseField& Phase, Composition& Cx); ///< Sets initial effective transformation stretches considering Vegard's contribution
    void SetEffectiveTransformationStretches(PhaseField& Phase, Temperature& Tx); ///< Sets initial effective transformation stretches considering thermal contribution
    void SetEffectiveTransformationStretches(PhaseField& Phase, Composition& Cx, Temperature& Tx);///< Sets initial effective transformation stretches considering Vegard's and thermal contribution

    void CalculateThermalExpansion(PhaseField& Phase, Temperature& Tx);         ///< Calculates thermal expansion contribution to the transformation stretches
    void CalculateVegardsExpansion(PhaseField& Phase, Composition& Cx);         ///< Calculates Vegard's expansion contribution to the transformation stretches
    void CalculateInterfaceStress(PhaseField& Phase, InterfaceProperties& IP);  ///< Calculates interface contribution to the transformation stretches

    void CalculateAverageElasticConstants(void);                                ///< Calculates AverageElasticConstants

    void SetEffectiveElasticConstants(PhaseField& Phase);                       ///< Sets effective elastic constants
    void SetEffectiveElasticConstants(PhaseField& Phase, Composition& Cx);      ///< Sets effective elastic constants considering chemo-mechanical coupling
    void SetEffectiveElasticConstants(PhaseField& Phase, Temperature& Tx);      ///< Sets effective elastic constants considering thermal coupling
    void SetEffectiveElasticConstants(PhaseField& Phase, Composition& Cx, Temperature& Tx);///< Sets effective elastic constants considering thermal and chemo-mechanical coupling

    void CalculateDrivingForce(PhaseField& Phase, DrivingForce& dGab) const;    ///< Calculates elastic driving force
    void CalculateDrivingForce(PhaseField& Phase, Composition& Cx, DrivingForce& dGab) const;    ///< Calculates elastic driving force
    void CalculateDrivingForce(PhaseField& Phase, Temperature& Tx, DrivingForce& dGab) const;    ///< Calculates elastic driving force
    vStrain CalculateNeuberCorrection(const vStrain ElasticStrains, const dMatrix6x6 locCompliance,
                                      const int i, const int j, const int k,
                                      const size_t pIndexA, const size_t pIndexB) const;/// Calculates Neuber corrected local elastic strains
    void CalculateDeformationJumps(const PhaseField& Phase);                    ///< Calculates deformation jumps from the Hadamard jump condition for all phase-field pairs

    void CalculateInterfaceEnergyContribution(PhaseField& Phase, InterfaceProperties& IP) const;///< Calculates interface energy contribution

    double Energy(void) const;                                                  ///< Elastic energy density in simulation domain [Joule]
    double EnergyDensity(int i, int j, int k) const;                            ///< Returns elastic energy density in a given point [Joule/m^3]
    double AverageEnergyDensity(void) const;                                    ///< Average elastic energy density in [Joule/<simulation cell>]

    void Advect(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi, const BoundaryConditions& BC, const double dt, const double tStep) override;

    // Binary input/
    void Read(const BoundaryConditions& BC, const int tSep);
    void ReadStresses(const int tStep);                                         ///< Reads restart input in binary format

    void ReadDeformationGradientsTotal(const int tStep);                        ///< Reads total deformation gradients in binary format
    void ReadDeformationGradientsEigen(const int tStep);                        ///< Reads eigen deformation gradients in binary format

    // Binary output
    void Write(const int tSep) const;                                           ///< Writes restart output in binary format
    void WriteStresses(const int tStep) const;                                  ///< Writes stresses in binary format

    void WriteDeformationGradientsTotal(const int tStep) const;                 ///< Writes total deformation gradients in binary format
    void WriteDeformationGradientsEigen(const int tStep) const;                 ///< Writes eigen deformation gradients in binary format

    // VTK output
    void WriteCauchyStressesVTK(const int tStep, const Settings& locSettings,
            const int precision=16) const;                                      ///< Writes Cauchy stresses in VTK format
    void WriteEffectiveElasticConstantsVTK(const int tStep, const Settings& locSettings,
            const int precision=16) const;                                      ///< Writes effective elastic constants in VTK format
    void WriteEigenStrainsVTK(const int tStep, const Settings& locSettings,
            const int precision=16) const;                                      ///< Writes effective eigenstrains using Voigt index notation, there is a factor 2 difference to the tensor notation
    void WriteDeformationGradientsTotalVTK(const int tStep, const Settings& locSettings,
            const int precision=20) const;                                      ///< Writes EffectiveEigenstrains using Voigt index notation, there is a factor 2 difference to the tensor notation
    void WriteElasticStrainsVTK(const int tStep, const Settings& locSettings,
            const int precision=16) const;                                      ///< Writes elastic strains in VTK format
    void WriteEnergyDensityVTK(const int tStep, const Settings& locSetings,
            const int precision=16) const;                                      ///< Writes local elastic energy in VTK format
    void WriteStressesVTK(const int tStep, const Settings& locSettings,
            const int precision=16) const;                                      ///< Writes stresses in VTK format
    void WriteTotalStrainsVTK(const int tStep, const Settings& locSettings,
            const int precision=16) const;                                      ///< Writes total strains in VTK format
    void WriteForceDensityVTK(const int tStep, const Settings& locSettings,
            const int precision=16) const;                                      ///< Writes force densities in VTK format
    void WriteTotalRotationsVTK(const int tStep, const Settings& locSettings,
            const int precision=16) const;
    void WriteDisplacementsVTK(const int tStep, const Settings& locSettings,
            const int precision=16) const;

    // ASCII output
    vStress WriteStressStrainData(std::string filename, std::string LDflag = "SD") const;

    void PrintPointStatistics(const int x, const int y, const int z) const;     ///< Prints various properties at a given point (x, y, z) to screen

    void GetAverageDeformationGradient();                                       ///< Calculates average deformation gradient

    ElasticProperties& operator= (const ElasticProperties& rhs);

    vStress AverageStress;                                                      ///< Average stress calculated by the spectral solver.
    vStrain AverageStrain;                                                      ///< Average strain calculated by the spectral solver.
    vStrain StrainToRemesh;                                                     ///< Accumulates the strain that has to be applied by remeshing
    vStrain RemeshedStrain;                                                     ///< Accumulates the strain that has been applied by remeshing
    vStress AppliedStress;                                                      ///< Applied stress tensor
    vStrain AppliedStrain;                                                      ///< Applied strain tensor
    vStrain AppliedStrainRate;                                                  ///< Applied strain rate tensor
    vStrain AppliedStrainOLD;                                                   ///< Applied strain tensor from the previous time step
    vStrain AppliedStressOLD;                                                   ///< Applied strain tensor from the previous time step
    dMatrix3x3 AverageDeformationGradient;                                      ///< Average deformation gradient

    //Elasticity Tensors:
    dMatrix6x6   AverageElasticConstants;                                       ///< Average elastic constants of the whole system
    dMatrix6x6   MAXElasticConstants;                                           ///< Maximum elastic constants of the system

    dVector6     AppliedStrainMask;                                             ///< Applied strain markers for individual strain components
    dVector6     AppliedStressMask;                                             ///< Applied stress markers for individual stress components

    Storage<dMatrix3x3>     TransformationStretches;                            ///< Reference transformation stretches for each phase field
    Storage<dMatrix6x6>     ElasticConstants;                                   ///< Reference elastic constants for each phase field
    Storage<dMatrix6x6>     Compliances;                                        ///< Inverse of the ElasticConstants

    Storage<dMatrix3x3>     PhaseTransformationStretches;                       ///< Reference transformation stretches gradients of each phase (relative to Cref), grain orientation not considered
    Storage<dMatrix6x6>     PhaseElasticConstants;                              ///< Reference elastic constants of each phase (relative to Cref), grain orientation not considered
    Storage<dMatrix6x6>     PhaseCompliences;                                   ///< Inverse of the PhaseElasticConstants, grain orientation not considered

    bool                    ThermoMechanicalCoupling;                           ///< Thermo-mechanical coupling flag
    Storage<dMatrix3x3>     PhaseAlpha;                                         ///< Linear thermal expansion coefficients for each phase
    Storage<dMatrix3x3>     Alpha;                                              ///< Linear thermal expansion coefficients for each phase field
    Storage<dMatrix6x6>     Gamma;                                              ///< Linear temperature dependence coefficients for the elasticity parameters for each phase field
    Storage<dMatrix6x6>     PhaseGamma;                                         ///< Linear temperature dependence coefficients for the elasticity parameters for each phase
    Storage<double>         Tref;                                               ///< Reference temperature for temperature dependence parameters for each phase
    Storage<double>         PoissonRatio;                                       ///< Poisson Ratio of each phase
    /* template specialization integer number stands for
     * the number of extra dimensions
     * order:
     * Rank = 2: phase, component
     * Rank = 1: phase or component
     */

    bool ChemoMechanicalCoupling;                                               ///< Chemo-mechanical coupling flag
    std::vector < std::string > Names;                                          ///< Names of corresponding chemical components (e.g Au, Cu, Na, Cl etc.)

    Tensor<dMatrix6x6, 2>   Kappa;                                              ///< Linear composition dependence coefficients for the elasticity parameters for each phase field
    Tensor<dMatrix6x6, 2>   PhaseKappa;                                         ///< Linear composition dependence coefficients for the elasticity parameters for each phase
    Tensor<dMatrix3x3, 2>   Lambda;                                             ///< Linear composition dependence coefficients for the transformation stretches for each phase field
    Tensor<dMatrix3x3, 2>   PhaseLambda;                                        ///< Linear composition dependence coefficients for the transformation stretches for each phase
    Tensor<double, 2>       Cref;                                               ///< Reference composition for composition dependence parameters for each phase

    Storage3D<dMatrix3x3, 0>    DeformationGradientsTotal;                      ///< Storage for total deformation gradients
    Storage3D<dMatrix3x3, 0>    DeformationGradientsEigen;                      ///< Storage for transformation induced deformation gradients
    Storage3D<NodeV3, 0>        DeformationJumps;                               ///< Storage for deformation jumps in the interfaces

    Storage3D<vStress, 0>       Stresses;                                       ///< Storage for stresses

    Storage3D<dMatrix6x6, 0>    EffectiveElasticConstants;                      ///< Storage for effective elastic constants

    Storage3D<dVector3,0>       ForceDensity;                                   ///< External force density
    Storage3D<dVector3,0>       Displacements;                                  ///< Displacements

    SymmetryVariants            Variants;                                       ///< Symmetry variants storage

    int TotalNx;                                                                ///< X dimension of the system in MPI parallel mode
    int OffsetX;                                                                ///< X position of the current domain in MPI parallel mode
    int TotalNy;                                                                ///< Y dimension of the system in MPI parallel mode
    int OffsetY;                                                                ///< Y position of the current domain in MPI parallel mode
    int TotalNz;                                                                ///< Z dimension of the system in MPI parallel mode
    int OffsetZ;                                                                ///< Z position of the current domain in MPI parallel mode
    int Nx;                                                                     ///< Size of domain along X dimension
    int Ny;                                                                     ///< Size of domain along Y dimension
    int Nz;                                                                     ///< Size of domain along Z dimension
    int dNx;                                                                    ///< Active X dimension
    int dNy;                                                                    ///< Active Y dimension
    int dNz;                                                                    ///< Active Z dimension

    double dx;                                                                  ///< Grid spacing
    size_t Nphases;                                                             ///< Number of thermodynamic phases
    size_t Ncomp;                                                               ///< Number of chemical components

    double StrainAccuracy;                                                      ///< Convergence parameter for the elasticity solver
    int MAXIterations;                                                          ///< Maximum iterations per time step for the elasticity solver

    ElasticityModels EModel;                                                    ///< Elasticity model selector

    bool KeepAspectRatio;                                                       ///< Restricts mechanical conditions to only volume relaxation
    bool PreventShear;                                                          ///< Prevents shear deformation of the simulation domain
    std::vector<NeuberParameters> NeuberCorrection;                             ///< Neuber correction parameters for each phase
    bool ConsiderExternalForces;                                                ///< Enables external force density consideration

    std::string VTKDir;
    std::string RawDataDir;

    /* These short methods are inlined and implemented
     * outside of the class declaration below
     */
    vStrain StrainSmall(const dMatrix3x3& locDefGrad) const;                    ///< Outputs strain considering small strain operation mode
    vStrain TotalStrains(int i, int j, int k) const;                            ///< Outputs strain considering strain operation mode (small, finite)
    vStrain EigenStrains(int i, int j, int k) const;                            ///< Outputs strain considering strain operation mode (small, finite)
    vStrain StressFreeStrains(int i, int j, int k) const;                       ///< Outputs strain considering strain operation mode (small, finite)
    vStrain ElasticStrains(int i, int j, int k) const;                          ///< Outputs strain considering strain operation mode (small, finite)
    vStress ShearStresses(int i, int j, int k) const;                           ///< Outputs shear stress
    vStress CauchyStress(int i, int j, int k) const;                            ///< Outputs Cauchy stress

    template<class T>
    void CalculateChemicalPotentialContribution(const PhaseField& Phase, T& DF) const ///< Calculates chemical potential contribution
    {
        /** This function calculates the first partial derivative of the mechanical
        free-energy density with respect to the phase-composition.*/

        if(ChemoMechanicalCoupling)
        {
            switch(EModel)
            {
                case ElasticityModels::Khachaturyan:
                {
                    ElasticityKhachaturyan::CalculateChemicalPotentialContribution(Phase, *this, DF);
                    break;
                }
                case ElasticityModels::Steinbach:
                {
                    ElasticitySteinbach::CalculateChemicalPotentialContribution(Phase, *this, DF);
                    break;
                }
                default:
                {
                    Info::WriteExit("Nonexistent elasticity model selected", thisclassname, "CalculateChemicalPotentialContribution()");
                    exit(1);
                }
            }
        }
    };
 protected:
 private:

};

// Inline methods implementation
inline vStrain ElasticProperties::StrainSmall(const dMatrix3x3& locDefGrad) const
{
    dMatrix3x3 locStrain = (locDefGrad.transposed() + locDefGrad)*0.5 - dMatrix3x3::UnitTensor();
    return VoigtStrain(locStrain);
};

inline vStrain ElasticProperties::TotalStrains(int i, int j, int k) const
{
    return StrainSmall(DeformationGradientsTotal(i,j,k));
};
inline vStrain ElasticProperties::EigenStrains(int i, int j, int k) const
{
    return StrainSmall(DeformationGradientsEigen(i,j,k));
};
inline vStrain ElasticProperties::StressFreeStrains(int i, int j, int k) const
{
    vStrain locStrain;

    locStrain = StrainSmall(DeformationGradientsEigen(i,j,k));

    return locStrain;
};
inline vStrain ElasticProperties::ElasticStrains(int i, int j, int k) const
{
    vStrain locStrain;
    locStrain = StrainSmall(DeformationGradientsTotal(i,j,k));
    locStrain -= StrainSmall(DeformationGradientsEigen(i,j,k));

    return locStrain;
};
inline vStress ElasticProperties::ShearStresses(int i, int j, int k) const
{
    vStress tmpStress = Stresses(i,j,k);

    double locTrace = tmpStress.trace()/(dNx+dNy+dNz);
    if(dNx) tmpStress[0] -= locTrace;
    if(dNy) tmpStress[1] -= locTrace;
    if(dNz) tmpStress[2] -= locTrace;

    return tmpStress;
}
inline vStress ElasticProperties::CauchyStress(int i, int j, int k) const
{
    dMatrix3x3 locDefGrad = DeformationGradientsTotal(i,j,k);
    double J_1 = 1.0/locDefGrad.determinant();
    dMatrix3x3 locCauchyStressTensor = (locDefGrad*Stresses(i,j,k).tensor()*
                                        locDefGrad.transposed())*J_1;
    return VoigtStress(locCauchyStressTensor);
}
}// namespace openphase
#endif
