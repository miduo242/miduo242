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
 *   File created :   2009
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev;
 *                         Reza Darvishi Kamachali; Raphael Schiedung;
 *                         Johannes Goerler; Marvin Tegeler
 *
 */

#ifndef PHASEFIELD_H
#define PHASEFIELD_H

#include "Base/Includes.h"
#include "Base/UserInterface.h"
#include "Base/NodeA.h"
#include "Base/NodePF.h"
#include "Base/NodeAB.h"
#include "Base/NodeV3.h"
#include "Base/Tensor.h"
#include "GrainInfo.h"
#include "H5Interface.h"

struct GrainPairStorage
{
    std::vector<std::vector<std::pair<size_t,double> > > content;               ///< Store Limits to apply for FieldsDot for each grain pair. Only upper triangle.
};

namespace openphase
{

class ElasticProperties;
class BoundaryConditions;
class ChemicalProperties;
class Velocities;
class FlowSolverLBM;
class Settings;
class H5Interface;

class OP_EXPORTS PhaseField : public OPObject                                   ///< Phase field class. It stores the phase fields, performs basic operations on them.
{
 public:

    PhaseField(){};
    PhaseField(Settings& locSettings, std::string InputFileName = DefaultInputFileName)
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    };
    void Initialize(Settings& locSettings) override;                            ///< Initializes the storage and initial variables of the driving force class.
    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC) override; ///< Changes the mesh size while keeping the data.

    void Coarsen(void);                                                         ///< Calculates single resolution phase-fields from the double resolution ones
    void CoarsenDot(void);                                                      ///< Calculates single resolution phase-field increments from the double resolution ones
    void Refine(void);                                                          ///< Calculates double resolution phase-fields from the single resolution ones

    dVector3 Normal(const int i, const int j, const int k,
                    const size_t alpha, const size_t beta) const;               ///< Returns interface normal for a given pair of phase fields at a given position
    dVector3 NormalDR(const int i, const int j, const int k,
                      const size_t alpha, const size_t beta) const;             ///< Returns interface normal for a given pair of phase fields at a given position in double resolution
    NodeV3 Normals(const int i, const int j, const int k) const;                ///< Returns interface normals for all pairs of phase fields at a given position
    NodeV3 NormalsDR(const int i, const int j, const int k) const;              ///< Returns interface normals for all pairs of phase fields at a given position in double resolution
    NodeV3 NormalsPhase(const int i, const int j, const int k) const;           ///< Returns interface normals for all pairs of phase fields at a given position
    NodeV3 NormalsPhaseDR(const int i, const int j, const int k) const;         ///< Returns interface normals for all pairs of phase fields at a given position
    bool Interface(const int x, const int y, const int z) const                 ///< Indicates the location of the interface.
    {
        return (Fields(x,y,z).flag > 1);
    };
    bool InterfaceDR(const int x, const int y, const int z) const               ///< Indicates the location of the interface.
    {
        return (FieldsDR(x,y,z).flag > 1);
    };

    double LocalNumberOfPhaseFieldsSR(int i, int j, int k)
    {
        if (!ConsiderNucleusVolume) return Fields(i,j,k).size();
        double norm = 0.0;
        for(auto alpha = Fields(i,j,k).cbegin();
                 alpha != Fields(i,j,k).cend(); ++alpha)
        {
            norm += FieldsStatistics[alpha->index].VolumeRatio;
        }
        return norm;
    };

    double LocalNumberOfPhaseFieldsDR(int i, int j, int k)
    {
        if (!ConsiderNucleusVolume) return FieldsDR(i,j,k).size();
        double norm = 0.0;
        for(auto alpha = FieldsDR(i,j,k).cbegin();
                 alpha != FieldsDR(i,j,k).cend(); ++alpha)
        {
            norm += FieldsStatistics[alpha->index].VolumeRatio;
        }
        return norm;
    };

    bool PhaseFieldPresent(const int i, const int j, const int k,
            const size_t Index) const;                                          ///< Returns true if phase field is present in a grid cell with coordinates (i,j,k)
    bool ThermodynamicPhasePresent(size_t alpha);
    bool ThermodynamicPhasePairPresent(size_t alpha, size_t beta);
    Storage3D<double,1> NewFractions(double dt);
    NodeA OrdinaryCurvatures(const int i, const int j, const int k) const;      ///< Returns ordinary curvatures of all phase fields at a given point
    NodeA OrdinaryCurvaturesDR(const int i, const int j, const int k) const;    ///< Returns ordinary curvatures of all phase fields at a given point
    double Curvature(const int i, const int j, const int k,
            const size_t phase) const;                                          ///< Curvature for plotting VTK
    double Interfaces(const int i, const int j, const int k) const;             ///< Indicates interfaces for plotting VTK
    double InterfacesDR(const int i, const int j, const int k) const;           ///< Indicates interfaces for plotting VTK
    double Junctions(const int i, const int j, const int k) const;              ///< Junctions for plotting VTK
    double JunctionsDR(const int i, const int j, const int k) const;            ///< Junctions for plotting VTK
    double LocalIndex(const int i, const int j, const int k) const;             ///< Phase-field index for plotting VTK
    double LocalIndexDR(const int i, const int j, const int k) const;           ///< Phase-field index for plotting VTK
    double Variants(const int i, const int j, const int k) const;               ///< Variants for plotting VTK
    double VariantsDR(const int i, const int j, const int k) const;             ///< Variants for plotting VTK
    size_t AddGrainInfo(size_t PhaseIndex);                                     ///< Adds new grain information for a phase "PhaseIndex", returns resulting phase field index
    size_t PlantGrainNucleus(size_t PhaseIndex, int x, int y, int z);           ///< Plants grain nucleus at position (x,y,z) of a phase "PhaseIndex", returns resulting phase field index
    std::array<double,2> PrincipalCurvatures(const int i, const int j,
            const int k, const size_t phase) const;                             ///< Principle Curvatures for plotting VTK
    std::pair<NodeA,NodeA> PrincipalCurvatures(const int i, const int j,
            const int k) const;                                                 ///< Returns principal curvatures of all phase fields at a given position
    std::vector<int> GetPresentPhaseFields() const;                             ///< Returns an ordered list with indices of present phase fields
    std::vector<int> ReturnVicinityPhaseFields(const int i, const int j,
            const int k) const;                                                 ///< Returns list of phase fields at point (i,j,k) and neighboring nodes
    void Advect(AdvectionHR& Adv, Velocities& Vel, BoundaryConditions& BC,
                FlowSolverLBM& LBM, double dt, double tStep);                   ///< Advects phase-fields for LBM cases
    void Advect(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi,
                const BoundaryConditions& BC,
                const double dt, const double tStep) override;                  ///< Advects phase-fields
    void CalculateDerivatives(void);                                            ///< Calculates Laplacians and gradients and stores them in the phase-fields nodes
    void CalculateDerivativesDR(void);                                          ///< Calculates Laplacians and gradients and stores them in the phase-fields nodes (Double resolution)
    void CalculateVolumes();                                                    ///< Collects volume for each phase field.
    void Clear();                                                               ///< Clears the phase field storage
    void ConsumePlane(const int dx, const int dy, const int dz, const int x,
                      const int y, const int z, const BoundaryConditions& BC);  ///< Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in the direction to/from (x, y, z) point.
    void FinalizeSR(const BoundaryConditions& BC, const bool finalize = true);  ///< Finalizing the phase fields calculations in this time step
    void FinalizeDR(const BoundaryConditions& BC, const bool finalize = true);  ///< Finalizing the phase fields calculations in this time step in double resolutions
    void Finalize  (const BoundaryConditions& BC, const bool finalize = true)
    {
        switch(Resolution)
        {
            case Resolutions::Single:
            {
                FinalizeSR(BC, finalize);
                break;
            }
            case Resolutions::Double:
            {
                FinalizeDR(BC, finalize);
                break;
            }
        }
    }

    double ColorScale(const int i, const int j, const int k) const;             ///< Variants for plotting VTK
    double ColorScaleDR(const int i, const int j, const int k) const;           ///< Variants for plotting VTK

    void FixSpreading(BoundaryConditions& BC, double cutoff);                   ///< Removes phase-fields with values below the cutoff (used for reducing parasitic diffusion in advection)
    void MergeIncrementsSR(const BoundaryConditions& BC, const double dt,
            const bool finalize = true);                                        ///< Merges the increments into the phase fields
    void MergeIncrementsDR(const BoundaryConditions& BC, const double dt,
               const bool finalize = true);                                     ///< Merges the increments into the phase fields
    void MergeIncrements(const BoundaryConditions& BC, const double dt,
                const bool finalize = true)                                     ///< Merges the increments into the phase fields
    {
        switch(Resolution)
        {
            case Resolutions::Single:
            {
                MergeIncrementsSR(BC, dt, finalize);
                break;
            }
            case Resolutions::Double:
            {
                MergeIncrementsDR(BC, dt, finalize);
                break;
            }
        }
    }

    void MoveFrame(const int dx, const int dy, const int dz,
            const BoundaryConditions& BC);                                      ///< Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in x, y and or z directions correspondingly.
    void MoveFrameDR(const int dx, const int dy, const int dz,
            const BoundaryConditions& BC);                                      ///< Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in x, y and or z directions correspondingly.
    void NormalizeIncrementsSR(const BoundaryConditions& BC, const double dt);  ///< Normalizes the interface fields such that resulting phase fields after merge do not escape interval [0,1] and sum up to 1.
    void NormalizeIncrementsDR(const BoundaryConditions& BC, const double dt);  ///< Normalizes the interface fields such that resulting phase fields after merge do not escape interval [0,1] and sum up to 1.
    void NormalizeIncrements(const BoundaryConditions& BC, const double dt)     ///< Normalizes the interface fields such that resulting phase fields after merge do not escape interval [0,1] and sum up to 1.
    {
        switch(Resolution)
        {
            case Resolutions::Single:
            {
                NormalizeIncrementsSR(BC, dt);
                break;
            }
            case Resolutions::Double:
            {
                NormalizeIncrementsDR(BC, dt);
                break;
            }
        }
    }

    void PrintPFVolumes() const;                                                ///< Prints volume of all phase fields
    void PrintPointStatistics(const int x, const int y, const int z) const;     ///< Prints the locally present phase fields and their values to the screen.
    void PrintVolumeFractions(void);                                            ///< Prints volume fraction of all thermodynamic phases
    bool Read(std::string FileName);                                            ///< Read raw (binary) phase fields from the file named FileName
    bool Read(const BoundaryConditions& BC, const int tStep,
            const bool finalize = true);                                        ///< Read raw (binary) phase fields from the file of specific time step
    bool Read(const BoundaryConditions& BC, const bool finalize = true);        ///< Read raw (binary) phase fields from the file thisclassname.dat

    bool ReadH5(const BoundaryConditions& BC, int tStep, H5Interface& H5);
    void SetBoundaryConditions(const BoundaryConditions& BC);                   ///< Set boundary conditions
    void SetBoundaryConditionsDR(const BoundaryConditions& BC);                 ///< Set boundary conditions in double resolution case
    void SetFlags();                                                            ///< Sets the flags which mark interfaces
    void SetFlagsDR();                                                          ///< Sets the flags which mark interfaces in double resolution case

    void KeepPhaseFieldsVolume(void);                                           ///< Keeps phase fields volume constant by allowing only grain shape change. Should be called before NormalizeIncrements()
    void KeepPhaseVolume(Tensor<bool,2> AllowedTransitions);                    ///< Keeps phases volume constant by allowing only grain shape change and grain transformations between the grains of the same phase (emulates coexistence of immiscible phases). Should be called before NormalizeIncrements()

    void SetIncrementsBoundaryConditions(const BoundaryConditions& BC);         ///< Set boundary conditions for phase field increments
    void SetIncrementsBoundaryConditionsDR(const BoundaryConditions& BC);       ///< Set boundary conditions for phase field increments in double resolution case
    void Write(const std::string& FileName) const;                              ///< Write raw (binary) phase fields to the file FileName

    void Write(const int tStep) const;                                          ///< Write raw (binary) phase fields to the file PhaseField_tStep.dat
    void Write() const;                                                         ///< Write raw (binary) phase fields to the file PhaseField.dat
    void WriteH5(const int tStep, H5Interface& H5);
    void WriteAverageVolume(const int tStep, const size_t PhaseIndex) const;    ///< Writes to the file the average volume of a given phase
    void WriteDistortedVTK(const int tStep,
            const Settings& locSettings,
            const ElasticProperties& EP,
            const bool CurvatureOutput = false,
            const int precision = 16) const;                                    ///< Write phase fields on the distorted grid to the file in VTK format in MPI environment

    void WriteGrainsStatistics(const int tStep);                                ///< Writes to the file full statistics on the grains, junctions etc. at the current time step.
    void WriteIndividualPhaseFieldValuesVTK(const int tStep,
                            const Settings& locSettings,
                            const std::initializer_list<size_t> FieldIndices,
                            const int precision = 16) const;                    ///< Writes values of individual phase fields to file in VTK format
    void WriteLaplacianVTK(const int tStep, const Settings& locSettings,
                            size_t PhiIndex, const int precision = 16);         ///< Write Laplacian of a given phase field to the file in VTK format
    void WriteVTK(const int tStep, const Settings& locSettings,
                  const bool CurvatureOutput = false,
                  const int precision = 16) const;                              ///< Write phase fields to the file in VTK format

    Storage3D< NodePF, 0 > Fields;                                              ///< Phase-field storage
    Storage3D< NodeAB, 0 > FieldsDot;                                           ///< Phase-field increments storage
    Storage3D< double, 1 > Fractions;                                           ///< Thermodynamic phase fractions storage

    Storage3D< NodePF, 0 > FieldsDR;                                            ///< Phase-field storage for double resolution mode
    Storage3D< NodeAB, 0 > FieldsDotDR;                                         ///< Phase-field increments storage for double resolution mode

    NodePF Dot (const int i, const int j, const int k, const double dt) const;  ///< Calculates the upcoming derivate of Fields with respect to time including finalization! (NOTE: Dot = Dot1 + Dot2)
    NodePF Dot1(const int i, const int j, const int k, const double dt) const;  ///< (used for semi-implicit solver Grand Potential solver!)
    NodePF Dot2(const int i, const int j, const int k, const double dt) const;  ///< (used for semi-implicit solver Grand Potential solver!)

    void CalculateFractions();                                                  ///< Calculates thermodynamic phase fractions
    Tensor<double,1> CalculateNewFractions(Tensor<double,1> oldFractions,
                                           NodeAB& FieldsDot, double locdt);    ///< Calculate new phase field values for each thermodynamic phase as if Phi.MergeIncrements would have been applied.
    Tensor<double,2> CalculatePsi(NodeAB& FieldsDot, double locdt);             ///< Calculate Phi.FieldsDot for thermodynamic phases.

    double dx;                                                                  ///< Grid spacing
    int    TotalNx;                                                             ///< X dimension of the system in MPI parallel mode
    int    OffsetX;                                                             ///< X position of the current domain in MPI parallel mode
    int    TotalNy;                                                             ///< Y dimension of the system in MPI parallel mode
    int    OffsetY;                                                             ///< Y position of the current domain in MPI parallel mode
    int    TotalNz;                                                             ///< Z dimension of the system in MPI parallel mode
    int    OffsetZ;                                                             ///< Z position of the current domain in MPI parallel mode

    int    Nx;                                                                  ///< X dimension of the system
    int    Ny;                                                                  ///< Y dimension of the system
    int    Nz;                                                                  ///< Z dimension of the system

    int    dNx;                                                                 ///< Active X dimension
    int    dNy;                                                                 ///< Active Y dimension
    int    dNz;                                                                 ///< Active Z dimension

    size_t Nphases;                                                             ///< Number of thermodynamic phases
    double Eta;                                                                 ///< Interface width in physical units
    double iWidth;                                                              ///< Interface width in grid points
    std::vector<std::string> PhaseNames;                                        ///< Names of phases for subset selection
    std::vector<AggregateStates> PhaseAggregateStates;                          ///< Aggregate states of all phases
    std::vector<double> FractionsTotal;                                         ///< Total phase fractions in the entire simulation domain
    bool NucleationPresent;                                                     ///< True if there are nuclei of any phase, false otherwise
    double RefVolume;                                                           ///< Reference volume for the nucleation
    bool ConsiderNucleusVolume;
    LaplacianStencil LStencil;                                                  ///< Laplacian stencil. Uses user specified stencil as the basis
    GradientStencil GStencil;                                                   ///< Gradient stencil. Uses user specified stencil as the basis

    GrainInfo FieldsStatistics;                                                 ///< Phase fields statistics. Contains information about location, velocity, orientation etc. of each phase field

    PhaseField& operator= (const PhaseField& rhs);                              ///< Copy operator for PhaseField class

    void CombinePhaseFields(const size_t PhaseIndex);                           ///< Merge phase fields of same phase to phase field with index PhaseIndex
    void SelectiveCombinePhaseFields(BoundaryConditions& BC,                    ///< Merge phase fields of sourceIndex to phase field with index targetPhaseIndex
            const size_t TargetPFIndex, const size_t SourcePFIndex)
    {
        switch(Resolution)
        {
            case Resolutions::Single:
            {
                SelectiveCombinePhaseFieldsSR(BC,TargetPFIndex,SourcePFIndex);
                break;
            }
            case Resolutions::Double:
            {
                SelectiveCombinePhaseFieldsDR(BC,TargetPFIndex,SourcePFIndex);
                break;
            }
        }
    }
    void SelectiveCombinePhaseFieldsSR(BoundaryConditions& BC,                  ///< Merge phase fields of sourceIndex to phase field with index targetPhaseIndex
            const size_t TargetPFIndex, const size_t SourcePFIndex);
    void SelectiveCombinePhaseFieldsDR(BoundaryConditions& BC,                  ///< Merge phase fields of sourceIndex to phase field with index targetPhaseIndex
            const size_t TargetPFIndex, const size_t SourcePFIndex);
    std::vector<int> GetMaxPhaseFieldOverlap(const size_t thPhase1,             ///< Returns a vector containing the phase indices thPhase1 and thPhase2 that have the biggest overlap, also the number of interface points between the phasefields. vector contains 1. phasefield of thPhase1, 2. phasefield of thePhase2, 3. number of overlap points. All values are -1 if no overlap is found.
            const size_t thPhase2);
    Matrix<int> GetPhaseFieldOverlap(const size_t thPhase1,                     ///< Returns a Matrix containing overlap of phase indices denoted in the columns and rows of the matrix. Only values above the diagonal are set.
            const size_t thPhase2);
    Resolutions Resolution;                                                     ///< Grid resolution for phase field. Can be single (1) or double (2)

    GrainPairStorage GrainPairLimits;

    std::string VTKDir;                                                         ///< Directory-path added in front of VTK files
    std::string RawDataDir;                                                     ///< Directory-path added in front of Restart files
 protected:
 private:
};

}// namespace openphase
#endif

