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
 *   File created :   2016
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Johannes Goerler
 *
 */

#ifndef ELASTICITYSOLVERSPECTRAL_H
#define ELASTICITYSOLVERSPECTRAL_H

/******************************************************************************
*  This module is based on the iterative algorithm of [S.Y.Hu, L.Q.Chen,      *
*  Acta. Mater. 49 (2001) 1879 - 1890] with modifications allowing to use     *
*  inhomogeneous elasticity parameters not limited to composition dependence. *
*  The boundary condition allowing free volume expansion is also introduced   *
*  following the approach outlined in [I.Steinbach, M.Apel, Physica D 217     *
*  (2006) 153 - 160].                                                         *
*  Finite strain extension considers Green-Lagrangian strain within the       *
*  St. Venant-Kirchhoff hyper-elastic model. The resulting solver algorithm   *
*  is similar to [P. Eisenlohr et al., International Journal of Plasticity 46 *
*  (2013) 37–53]                                                              *
******************************************************************************/

#ifdef MPI_PARALLEL
#include "fftw3-mpi.h"
#else
#include "fftw3.h"
#endif
#include "Base/Includes.h"

namespace openphase
{
class Settings;
class ElasticProperties;

class OP_EXPORTS ElasticitySolverSpectral : public OPObject                                ///< Elastic problem solver based on spectral method.
{
 public:
    ElasticitySolverSpectral(){};
    ElasticitySolverSpectral(Settings& locSettings);                            ///< Constructor calls Initialize and ReadInput
    ElasticitySolverSpectral(Settings& locSettings,
                             const BoundaryConditions& BC);                     ///< Constructor calls Initialize and ReadInput
    ~ElasticitySolverSpectral(void);
    void Initialize(Settings& locSettings) override;                            ///< Named constructor, allocates internal storages and initialized variables
    void Initialize(Settings& locSettings, const BoundaryConditions& BC);       ///< Named constructor (considering non periodic boundary conditions)
    void Reinitialize(ElasticProperties& EP);                                   ///< Needed if the system dimensions have been altered
    void Reinitialize(ElasticProperties& EP, const BoundaryConditions& BC);     ///< Needed if the system dimensions have been altered (considering non periodic boundary conditions)

    int  Solve(ElasticProperties& EP, BoundaryConditions& BC, double dt,
               std::function<bool()> EndSolve = [](){return false;});           ///< Solves mechanical equilibrium problem

    int TotalNx;                                                                ///< X dimension of the system in MPI parallel mode
    int OffsetX;                                                                ///< X position of the current domain in MPI parallel mode
    int TotalNy;                                                                ///< Y dimension of the system in MPI parallel mode
    int OffsetY;                                                                ///< Y position of the current domain in MPI parallel mode
    int TotalNz;                                                                ///< Z dimension of the system in MPI parallel mode
    int OffsetZ;                                                                ///< Z position of the current domain in MPI parallel mode

    int Nx;                                                                     ///< System size along X direction
    int Ny;                                                                     ///< System size along Y direction
    int Nz;                                                                     ///< System size along Z direction
    int dNx;                                                                    ///< Active X dimension
    int dNy;                                                                    ///< Active Y dimension
    int dNz;                                                                    ///< Active Z dimension
    double dx;                                                                  ///< Grid spacing

    bool Smooth;                                                                ///< Dampen high-frequency modes in Fourier space
    double relaxation_parameter;                                                ///< Convergence parameter for the homogeneous strain relaxation
 private:
    long int  Nz2;                                                              ///< Half of the system size along Z direction

    double *  UandForce[3];                                                     ///< Force and displacements in real and reciprocal space
    double *  RHSandDefGrad[9];                                                 ///< RHS and deformation gradient in real and reciprocal space

    fftw_plan ForwardPlanRHS[9];                                                ///< Forward FFT plans for the RHSide
    fftw_plan ForwardPlanForce[3];                                              ///< Forward FFT plans for the force
    fftw_plan BackwardPlanDefGrad[9];                                           ///< Backward FFT plans for the deformation gradients
    fftw_plan BackwardPlanU[3];                                                 ///< Backward FFT plans for the displacements

    void CopyForceDensity(ElasticProperties& EP);                               ///< Copies force density into the internal storage of the solver
    void CalculateRHS(ElasticProperties& EP);                                   ///< Calculates right hand side of the mechanical equilibrium equation
    void ExecuteForwardFFT(ElasticProperties& EP);                              ///< Executes forward FFT of the RHS
    void CalculateFourierSolution(ElasticProperties& EP);                       ///< Calculates solution in Fourier space
    void ExecuteBackwardFFT(void);                                              ///< Executes backward FFT of the deformation gradients
    void ExecuteBackwardFFTdisplacements(void);                                 ///< Execute forward FFT of the displacements
    void ExecuteForwardFFTforces(void);                                         ///< Execute forward FFT of the external forces
    void CalculateTargetStress(ElasticProperties& EP, double& MAXStrainDeviation,
                                                      vStress& TargetStress);   ///< Calculates average stress
    void CopyDisplacements(ElasticProperties& EP);                              ///< Copies displacements to the ElasticProperties
    void SetElasticProperties(ElasticProperties& EP);                           ///< Sets stresses and deformation gradients

    void ApplyMechanicalBC(ElasticProperties& EP, vStress TargetStress);        ///< Applies mechanical boundary conditions
};
} // namespace openphase
#endif //SpectralElasticitySolver
