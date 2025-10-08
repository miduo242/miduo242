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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */

#ifndef DRIVINGFORCE_H
#define DRIVINGFORCE_H

#include "Base/Includes.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"

namespace openphase
{

class Settings;
class NodeAB;
class PhaseField;
class BoundaryConditions;
class InterfaceProperties;
/***************************************************************/
enum class NoiseModes{Zero, Relative, Absolute};

class OP_EXPORTS DrivingForce : public OPObject                                            ///< The module that handles the driving force. Provides the storage and manipulation methods.
{
 public:
    DrivingForce() : Averaging(false) {};
    DrivingForce(Settings& locSettings,
                 const std::string InputFileName = DefaultInputFileName);       ///< Initializes the storage and internal variables of the driving force class.
    void Initialize(Settings& locSettings) override;                            ///< Initializes the storage and internal variables of the driving force class.
    void ReadInput(const std::string InputFileName) override;                   ///< Reads driving force settings
    void ReadInput(std::stringstream& inp) override;                            ///< Reads driving force settings
    void Refine(PhaseField& Phase);
    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC) override;                      ///< Remeshs the storage while keeping the data
    void ClearSR();                                                             ///< Empties the storage
    void ClearDR();                                                             ///< Empties the double resolution storage
    void Average(const PhaseField& Phase, const BoundaryConditions& BC);        ///< Average driving forces across the interface
    void CollectAverage(const PhaseField& Phase);                               ///< First part of Average
    void DistributeAverage(const PhaseField& Phase);                            ///< Second part of Average
    void Clear(void);                                                           ///< Deletes driving forces in the storage. Needs to be called at the end/beginning of each time step!
    void SetBoundaryConditions(const BoundaryConditions& BC);                   ///< Sets the boundary conditions
    void MergePhaseFieldIncrementsSR(PhaseField& Phase, InterfaceProperties& IP);///< Merges the driving force into the phase field increments.
    void MergePhaseFieldIncrementsDR(PhaseField& Phase, InterfaceProperties& IP);///< Merges the driving force into the phase field increments.
    void MergePhaseFieldIncrements(PhaseField& Phase, InterfaceProperties& IP); ///< Merges the driving force into the phase field increments.

    NodeAB CalcUnified(const PhaseField& Phase, const BoundaryConditions& BC);  ///< Reduces raw driving force to a single number
    void   Unify(const PhaseField& Phase, const BoundaryConditions& BC);        ///< Reduces raw driving force to a single number

    double MaxTimeStep(PhaseField& Phase, InterfaceProperties& IP,
                       Settings& OPSettings, double TheorLimit,
                       double NumLimit);                                        ///< Calculates minimum number of phase field iterations for a given time step

    void PrintDiagnostics(void);                                                ///< Prints driving force statistics to screen. Indicates if driving force is overshooting.
    void PrintPointStatistics(const int i, const int j, const int k) const;
  
	ClearingModes dGClearingMode;
    int TotalNx;                                                                ///< X dimension of the system in MPI parallel mode
    int OffsetX;                                                                ///< X position of the current domain in MPI parallel mode
    int TotalNy;                                                                ///< Y dimension of the system in MPI parallel mode
    int OffsetY;                                                                ///< Y position of the current domain in MPI parallel mode
    int TotalNz;                                                                ///< Z dimension of the system in MPI parallel mode
    int OffsetZ;                                                                ///< Z position of the current domain in MPI parallel mode

    int Nx;                                                                     ///< X dimension of the system.
    int Ny;                                                                     ///< Y dimension of the system.
    int Nz;                                                                     ///< Z dimension of the system.

    int dNx;                                                                    ///< Active X dimension
    int dNy;                                                                    ///< Active Y dimension
    int dNz;                                                                    ///< Active Z dimension

    size_t Nphases;                                                             ///< Number of thermodynamic phases
    int Range;                                                                  ///< Radius of a sphere around a grid point over which the driving force is averaged
    double PhiThreshold;                                                        ///< Outlines the inner part of the interface for driving force averaging.
    void WriteVTK(const int tStep, const Settings& locSettings,
                  const size_t indexA, const size_t indexB,
                  const int precision = 16) const;                              ///< Writes the driving force acting from phase field with indexB onto phase field with indexA
    void WriteVTKforPhases(const int tStep, const Settings& locSettings,
                           const PhaseField& Phi,
                           const int precision = 16) const;                     ///< Writes the average driving force acting between all thermodynamic phases, not between individual grains
    Storage3D<NodeAB, 0> Raw;                                                   ///< Raw driving force storage
    Storage3D<NodeAB, 0> RawDR;                                                 ///< Raw driving force storage in double resolution

    DrivingForce& operator= (const DrivingForce& rhs);                          ///< Copy operator for DrivingForce class

    long int OvershootCounter;                                                  ///< Number of driving force overshooting events
    Matrix<double> MAXOvershootPOS;                                             ///< Maximum positive driving force overshoot for each phase pair
    Matrix<double> MAXOvershootNEG;                                             ///< Maximum negative driving force overshoot for each phase pair
    Matrix<double> MAXDrivingForcePOS;                                          ///< Maximum positive driving force value for each phase pair
    Matrix<double> MAXDrivingForceNEG;                                          ///< Maximum negative driving force value for each phase pair

    bool Averaging;                                                             ///< Control parameter for averaging the driving force over the interface.
    bool Unifying;                                                              ///< Control parameter for unification of the driving force over the interface.
    bool dGcut;                                                                 ///< If set to true enables the driving force cutoff to stabilize interface profile
    const double CutOffDefault = 0.95;                                          ///< Default driving force cutoff value
    Matrix<double> CutOff;                                                      ///< Driving force cutoff for pairs of phases

    NoiseModes NoiseMode;                                                       ///< Noise mode
    Matrix<double> Noise;                                                       ///< Noise amplitude for pairs of phases
    Resolutions Resolution;                                                     ///< Single or double resolution switch

    std::string VTKDir;                                                         ///< Directory-path added in front of VTK files
    std::string RawDataDir;                                                     ///< Directory-path added in front of Restart files

    NodeAB Raw_at(const double x, const double y, const double z) const;        ///< Arbitrary point access operator for driving force. Uses tri-linear interpolation

    inline double CalculatePrefactorSR(PhaseField& Phase, int i, int j, int k)  ///< Returns driving force prefactor in single resolution
    {
        return 2.0*Pi/(Phase.Eta*Phase.LocalNumberOfPhaseFieldsSR(i,j,k));
    };

    inline double CalculatePrefactorDR(PhaseField& Phase, int i, int j, int k)  ///< Returns driving force prefactor in double resolution
    {
        return 2.0*Pi/(Phase.Eta*Phase.LocalNumberOfPhaseFieldsDR(i,j,k));
    };
    double GetDrivingForce(PhaseField& Phi, const int i, const int j, const int k, const size_t alpha, const size_t beta) const
    {
    double tempdG = 0.0;
        for(auto it1 = Phi.Fields(i,j,k).cbegin();
                 it1 < Phi.Fields(i,j,k).cend(); ++it1)
        for(auto it2 = Phi.Fields(i,j,k).cbegin();
                 it2 < Phi.Fields(i,j,k).cend(); ++it2)
        if((Phi.FieldsStatistics[it1->index].Phase == alpha)
        and(Phi.FieldsStatistics[it2->index].Phase == beta))
        {
            tempdG += Raw(i,j,k).get_asym1(it1->index, it2->index);
        }
        return tempdG;
    }
    void AverageGlobal(PhaseField& Phase, double time)
    {
        AverageDG2.clear();
        AverageDG2 = AverageDG1;
        AverageDG1.clear();
        AverageDG1 = AverageDG;
        AverageDG.clear();
        STORAGE_LOOP_BEGIN(i,j,k,Raw,0)
        {
            if (Phase.Interface(i,j,k))
            for(auto it = Raw(i,j,k).begin();
                     it < Raw(i,j,k).end(); ++it)
            {
                AverageDG.add_asym1(it->indexA,it->indexB,Raw(i, j, k).get_asym1(it->indexA, it->indexB));
                AverageDG.add_sym2(it->indexA,it->indexB,1.);
            }
        }
        STORAGE_LOOP_END
        std::ofstream outfile("AvDg.dat", std::ofstream::app);
        outfile << time << " ";
        for(auto it = AverageDG.begin();
                     it < AverageDG.end(); ++it)
        {
            it->value1 /= it->value2;

        }
        outfile << AverageDG.get_asym1(0,1) << " ";
        outfile << AverageDG.get_asym1(0,2) << " ";
        outfile << AverageDG.get_asym1(0,3) << " ";
        outfile << AverageDG.get_asym1(1,3) << " ";
        outfile << AverageDG.get_asym1(2,3) << " ";
        outfile << std::endl;
        outfile.close();
   }
   NodeAB AverageDG;
   NodeAB AverageDG1;
   NodeAB AverageDG2;
 protected:
 private:
    double maxPsi;
};
} // namespace openphase
#endif
