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
 *   Main contributors :   Oleg Shchyglo; Philipp Engels; Raphael Schiedung;
 *                         Marvin Tegeler; Helge Schaar
 *
 */

#ifndef HEATDIFFUSION_H
#define HEATDIFFUSION_H

#include "Base/Includes.h"
#include "BoundaryConditions.h"
#include "Temperature.h"

namespace openphase
{
class Settings;
class GrainInfo;
class PhaseField;
class FlowSolverLBM;
class Temperature;
class BoundaryConditions;
class Composition;

class OP_EXPORTS HeatDiffusion : public OPObject                                ///< Heat equation solver class
{
 public:
    HeatDiffusion(){};                                                          ///< Constructor
    HeatDiffusion(Settings& locSettings,
                  const std::string InputFileName = DefaultInputFileName)       ///< Initializes storages, sets internal variables.
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }
    void Initialize(Settings& locSettings) override;                            ///< Allocates memory, initializes the settings);
    void ReadInput(const std::string FileName) override;                        ///< Reads input parameters
    void ReadInput(std::stringstream& inp);                                     ///< Reads input parameters
    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC) override; ///< Changes system size while keeping the data

    void SetEffectiveProperties(const PhaseField& Phase,
                                const Temperature& Tx);                         ///< Set effective thermal properties

    void SetEffectiveProperties(const PhaseField& Phase,
                                const Temperature& Tx,
                                const FlowSolverLBM& FL);                       ///< Set effective thermal properties

    void SetLatentHeat(const PhaseField& Phase);                                ///< Sets local latent heat sources

    int  SolveImplicit(const PhaseField& Phase, const BoundaryConditions& BC,   ///< Implicit Solver
                       Temperature& Tx, const double dt);

    double MaxThermalConductivity;                                              ///< Maximum thermal conductivity in the system

    Storage<double> PhaseThermalConductivity;                                   ///< Thermal conductivity for each phase
    Storage<double> PhaseVolumetricHeatCapacity;                                ///< Volumetric heat capacity for each phase
    Storage<double> PhaseSpecificHeatCapacity;                                  ///< Specific heat capacity for each phase

    Storage3D<double,0> EffectiveThermalConductivity;                           ///< Thermal conductivity
    Storage3D<double,0> EffectiveHeatCapacity;                                  ///< Volumetric heat capacity
    Storage3D<double,0> Qdot;                                                   ///< Heat source in the bulk of the simulation domain

    Storage3D<double,0> TxOld;                                                  ///< Temporary temperature storage for the iterative solver
    Storage3D<double,0> dTx;                                                    ///< Temperature increments

    Matrix<double> LatentHeat;                                                  ///< Latent heat values for each phase pair

    double Tolerance;                                                           ///< Solver convergence tolerance
    int MaxIterations;                                                          ///< Maximum number of implicit solver iterations

    int TotalNx;                                                                ///< X dimension of the system in MPI parallel mode
    int OffsetX;                                                                ///< X position of the current domain in MPI parallel mode
    int TotalNy;                                                                ///< Y dimension of the system in MPI parallel mode
    int OffsetY;                                                                ///< Y position of the current domain in MPI parallel mode
    int TotalNz;                                                                ///< Z dimension of the system in MPI parallel mode
    int OffsetZ;                                                                ///< Z position of the current domain in MPI parallel mode

    int Nx;                                                                     ///< Size of the inner calculation domain along X
    int Ny;                                                                     ///< Size of the inner calculation domain along Y
    int Nz;                                                                     ///< Size of the inner calculation domain along Y

    int dNx;                                                                    ///< Active X dimension
    int dNy;                                                                    ///< Active Y dimension
    int dNz;                                                                    ///< Active Z dimension

    size_t Nphases;                                                             ///< Number of phases
    double dx;                                                                  ///< Grid spacing

 protected:
 private:
};

} // namespace openphase
#endif
