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

#ifdef MPI_PARALLEL
#include "fftw3-mpi.h"
#else
#include "fftw3.h"
#endif
#include "Info.h"
#include "Settings.h"
#include "Base/UserInterface.h"
#include "BoundaryConditions.h"
#include "Mechanics/ElasticProperties.h"
#include "Mechanics/ElasticitySolverSpectral.h"
#include "VTK.h"

namespace openphase
{
using namespace std;

ElasticitySolverSpectral::ElasticitySolverSpectral(Settings& locSettings)
{
    Initialize(locSettings);
}

ElasticitySolverSpectral::ElasticitySolverSpectral(Settings& locSettings,
        const BoundaryConditions& BC)
{
    Initialize(locSettings,BC);
}

void ElasticitySolverSpectral::Initialize(Settings& locSettings)
{
    thisclassname = "ElasticitySolverSpectral";

    Smooth = false;
    relaxation_parameter  = 0.25;

    //Initialize FFTW OpenMP threads
#ifdef _OPENMP
    fftw_init_threads();
#endif

    // Set system dimensions
    TotalNx = locSettings.TotalNx;
    OffsetX = locSettings.OffsetX;
    TotalNy = locSettings.TotalNy;
    OffsetY = locSettings.OffsetY;
    TotalNz = locSettings.TotalNz;
    OffsetZ = locSettings.OffsetZ;

    Nx  = locSettings.Nx;
    Ny  = locSettings.Ny;
    Nz  = locSettings.Nz;

    dNx = locSettings.dNx;
    dNy = locSettings.dNy;
    dNz = locSettings.dNz;

    Nz2 = Nz/2 + 1;

    dx = locSettings.dx;

    // Initialize FFTW MPI parallelism
#ifdef MPI_PARALLEL
    fftw_mpi_init();

    ptrdiff_t local_n0      = 0;
    ptrdiff_t local_0_start = 0;

    size_t SIZE = 2*fftw_mpi_local_size_3d(TotalNx, TotalNy, TotalNz, MPI_COMM_WORLD,
                                  &local_n0, &local_0_start);
#else
    size_t SIZE = Nx*Ny*Nz2*2;
#endif

    // Allocate arrays
    for(int n = 0; n < 9; n++)
    {
        RHSandDefGrad[n] = (double *)fftw_malloc(sizeof(double)*SIZE);

    }
    for(int n = 0; n < 3; n++)
    {
        UandForce[n] = (double *)fftw_malloc(sizeof(double)*SIZE);
    }

    // Set number of FFTW OpenMP threads
#ifdef _OPENMP
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    // Create FFT plans
    for(int n = 0; n < 9; n++)
    {
        #ifdef MPI_PARALLEL
        ForwardPlanRHS[n] = fftw_mpi_plan_dft_r2c_3d (TotalNx, TotalNy, TotalNz, RHSandDefGrad[n],
                        reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                        MPI_COMM_WORLD,
                        FFTW_ESTIMATE);
        #else
        ForwardPlanRHS[n] = fftw_plan_dft_r2c_3d
                        (Nx, Ny, Nz, RHSandDefGrad[n],
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         FFTW_ESTIMATE);
        #endif
    }

    for(int n = 0; n < 9; n++)
    {
        #ifdef MPI_PARALLEL
        BackwardPlanDefGrad[n] = fftw_mpi_plan_dft_c2r_3d
                        (TotalNx, TotalNy, TotalNz,
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         RHSandDefGrad[n],
                         MPI_COMM_WORLD,
                         FFTW_ESTIMATE);
        #else
        BackwardPlanDefGrad[n] = fftw_plan_dft_c2r_3d
                        (Nx, Ny, Nz,
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         RHSandDefGrad[n],
                         FFTW_ESTIMATE);
        #endif
    }

    for(int n = 0; n < 3; n++)
    {
        #ifdef MPI_PARALLEL
        BackwardPlanU[n] = fftw_mpi_plan_dft_c2r_3d
                            (TotalNx, TotalNy, TotalNz,
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             UandForce[n],
                             MPI_COMM_WORLD,
                             FFTW_ESTIMATE);

        ForwardPlanForce[n] = fftw_mpi_plan_dft_r2c_3d
                           (TotalNx, TotalNy, TotalNz, UandForce[n],
                            reinterpret_cast<fftw_complex*> (UandForce[n]),
                            MPI_COMM_WORLD,
                            FFTW_ESTIMATE);
        #else
        BackwardPlanU[n] = fftw_plan_dft_c2r_3d
                            (Nx, Ny, Nz,
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             UandForce[n],
                             FFTW_ESTIMATE);

        ForwardPlanForce[n] = fftw_plan_dft_r2c_3d
                            (Nx, Ny, Nz, UandForce[n],
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             FFTW_ESTIMATE);
        #endif
    }

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void ElasticitySolverSpectral::Initialize(Settings& locSettings, const BoundaryConditions& BC)
{
    thisclassname = "ElasticitySolverSpectral";

    Smooth = false;
    relaxation_parameter  = 0.25;

    //Initialize FFTW OpenMP threads
#ifdef _OPENMP
    fftw_init_threads();
#endif

    // Set system dimensions
    TotalNx = locSettings.TotalNx;
    OffsetX = locSettings.OffsetX;
    TotalNy = locSettings.TotalNy;
    OffsetY = locSettings.OffsetY;
    TotalNz = locSettings.TotalNz;
    OffsetZ = locSettings.OffsetZ;

    Nx      = locSettings.Nx;
    Ny      = locSettings.Ny;
    Nz      = locSettings.Nz;

    if(BC.BC0X != BoundaryConditionTypes::Periodic or
       BC.BCNX != BoundaryConditionTypes::Periodic)
    {
#ifdef MPI_PARALLEL
        std::cerr << "ElasticitySolverSpectral::Initialize()\n"
                  << "\tNonperiodic boundary conditions in X-direction are not permitted in MPI parallel mode!" << std::endl;
        exit(1);
#else
        TotalNx *= 2;
        Nx *= 2;
#endif
    }
    if(BC.BC0Y != BoundaryConditionTypes::Periodic or
       BC.BCNY != BoundaryConditionTypes::Periodic)
    {
        TotalNy *= 2;
        Ny *= 2;
    }
    if(BC.BC0Z != BoundaryConditionTypes::Periodic or
       BC.BCNZ != BoundaryConditionTypes::Periodic)
    {
        TotalNz *= 2;
        Nz *= 2;
    }

    dNx = locSettings.dNx;
    dNy = locSettings.dNy;
    dNz = locSettings.dNz;

    Nz2 = Nz/2+1;

    dx = locSettings.dx;

    // Initialize FFTW MPI parallelism
#ifdef MPI_PARALLEL
    fftw_mpi_init();

    ptrdiff_t local_n0      = 0;
    ptrdiff_t local_0_start = 0;

    size_t SIZE = 2*fftw_mpi_local_size_3d(TotalNx, TotalNy, TotalNz, MPI_COMM_WORLD,
                                  &local_n0, &local_0_start);
#else
    size_t SIZE = Nx*Ny*Nz2*2;
#endif

    // Allocate arrays
    for(int n = 0; n < 9; n++)
    {
        RHSandDefGrad[n] = (double *)fftw_malloc(sizeof(double)*SIZE);
    }
    for(int n = 0; n < 3; n++)
    {
        UandForce[n] = (double *)fftw_malloc(sizeof(double)*SIZE);
    }

    // Set number of FFTW OpenMP threads
#ifdef _OPENMP
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    // Create FFT plans
    for(int n = 0; n < 9; n++)
    {
        #ifdef MPI_PARALLEL
        ForwardPlanRHS[n] = fftw_mpi_plan_dft_r2c_3d
                       (TotalNx, TotalNy, TotalNz, RHSandDefGrad[n],
                        reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                        MPI_COMM_WORLD,
                        FFTW_ESTIMATE);
        #else
        ForwardPlanRHS[n] = fftw_plan_dft_r2c_3d
                        (Nx, Ny, Nz, RHSandDefGrad[n],
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         FFTW_ESTIMATE);
        #endif
    }

    for(int n = 0; n < 9; n++)
    {
        #ifdef MPI_PARALLEL
        BackwardPlanDefGrad[n] = fftw_mpi_plan_dft_c2r_3d
                        (TotalNx, TotalNy, TotalNz,
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         RHSandDefGrad[n],
                         MPI_COMM_WORLD,
                         FFTW_ESTIMATE);
        #else
        BackwardPlanDefGrad[n] = fftw_plan_dft_c2r_3d
                        (Nx, Ny, Nz,
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         RHSandDefGrad[n],
                         FFTW_ESTIMATE);
        #endif
    }

    for(int n = 0; n < 3; n++)
    {
        #ifdef MPI_PARALLEL
        BackwardPlanU[n] = fftw_mpi_plan_dft_c2r_3d
                            (TotalNx, TotalNy, TotalNz,
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             UandForce[n],
                             MPI_COMM_WORLD,
                             FFTW_ESTIMATE);

        ForwardPlanForce[n] = fftw_mpi_plan_dft_r2c_3d
                           (TotalNx, TotalNy, TotalNz, UandForce[n],
                            reinterpret_cast<fftw_complex*> (UandForce[n]),
                            MPI_COMM_WORLD,
                            FFTW_ESTIMATE);
        #else
        BackwardPlanU[n] = fftw_plan_dft_c2r_3d
                            (Nx, Ny, Nz,
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             UandForce[n],
                             FFTW_ESTIMATE);

        ForwardPlanForce[n] = fftw_plan_dft_r2c_3d
                            (Nx, Ny, Nz, UandForce[n],
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             FFTW_ESTIMATE);
        #endif
    }
    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void ElasticitySolverSpectral::Reinitialize(ElasticProperties& EP)
{
    // Set new system dimensions
    TotalNx = EP.TotalNx;
    OffsetX = EP.OffsetX;
    TotalNy = EP.TotalNy;
    OffsetY = EP.OffsetY;
    TotalNz = EP.TotalNz;
    OffsetZ = EP.OffsetZ;

    Nx  = EP.Nx;
    Ny  = EP.Ny;
    Nz  = EP.Nz;

    dNx = EP.dNx;
    dNy = EP.dNy;
    dNz = EP.dNz;

    Nz2 = Nz/2 + 1;

    dx = EP.dx;

#ifdef MPI_PARALLEL

    ptrdiff_t local_n0      = 0;
    ptrdiff_t local_0_start = 0;

    size_t SIZE = 2*fftw_mpi_local_size_3d(TotalNx, TotalNy, TotalNz, MPI_COMM_WORLD,
                                  &local_n0, &local_0_start);
#else
    size_t SIZE = Nx*Ny*Nz2*2;
#endif

    // Destroy old FFT plans and free allocated memory
    for(int n = 0; n < 9; n++)
    {
        fftw_destroy_plan(ForwardPlanRHS[n]);
        fftw_destroy_plan(BackwardPlanDefGrad[n]);

        fftw_free(RHSandDefGrad[n]);
    }
    for(int n = 0; n < 3; n++)
    {
        fftw_destroy_plan(BackwardPlanU[n]);
        fftw_destroy_plan(ForwardPlanForce[n]);

        fftw_free(UandForce[n]);
    }

    // Reallocate arrays
    for(int n = 0; n < 9; n++)
    {
        RHSandDefGrad[n] = (double *)fftw_malloc(sizeof(double)*SIZE);
    }
    for(int n = 0; n < 3; n++)
    {
        UandForce[n] = (double *)fftw_malloc(sizeof(double)*SIZE);
    }

    // Create new FFT plans
    for(int n = 0; n < 9; n++)
    {
        #ifdef MPI_PARALLEL
        ForwardPlanRHS[n] = fftw_mpi_plan_dft_r2c_3d
                       (TotalNx, TotalNy, TotalNz, RHSandDefGrad[n],
                        reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                        MPI_COMM_WORLD,
                        FFTW_ESTIMATE);
        #else
        ForwardPlanRHS[n] = fftw_plan_dft_r2c_3d
                        (Nx, Ny, Nz, RHSandDefGrad[n],
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         FFTW_ESTIMATE);
        #endif
    }

    for(int n = 0; n < 9; n++)
    {
        #ifdef MPI_PARALLEL
        BackwardPlanDefGrad[n] = fftw_mpi_plan_dft_c2r_3d
                        (TotalNx, TotalNy, TotalNz,
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         RHSandDefGrad[n],
                         MPI_COMM_WORLD,
                         FFTW_ESTIMATE);
        #else
        BackwardPlanDefGrad[n] = fftw_plan_dft_c2r_3d
                        (Nx, Ny, Nz,
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         RHSandDefGrad[n],
                         FFTW_ESTIMATE);
        #endif
    }

    for(int n = 0; n < 3; n++)
    {
        #ifdef MPI_PARALLEL
        BackwardPlanU[n] = fftw_mpi_plan_dft_c2r_3d
                            (TotalNx, TotalNy, TotalNz,
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             UandForce[n],
                             MPI_COMM_WORLD,
                             FFTW_ESTIMATE);

        ForwardPlanForce[n] = fftw_mpi_plan_dft_r2c_3d
                           (TotalNx, TotalNy, TotalNz, UandForce[n],
                            reinterpret_cast<fftw_complex*> (UandForce[n]),
                            MPI_COMM_WORLD,
                            FFTW_ESTIMATE);
        #else
        BackwardPlanU[n] = fftw_plan_dft_c2r_3d
                            (Nx, Ny, Nz,
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             UandForce[n],
                             FFTW_ESTIMATE);

        ForwardPlanForce[n] = fftw_plan_dft_r2c_3d
                            (Nx, Ny, Nz, UandForce[n],
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             FFTW_ESTIMATE);
        #endif
    }
    Info::WriteStandard(thisclassname, "Reinitialized");
}

void ElasticitySolverSpectral::Reinitialize(ElasticProperties& EP, const BoundaryConditions& BC)
{
    // Set new system dimensions
    TotalNx = EP.TotalNx;
    OffsetX = EP.OffsetX;
    TotalNy = EP.TotalNy;
    OffsetY = EP.OffsetY;
    TotalNz = EP.TotalNz;
    OffsetZ = EP.OffsetZ;

    Nx  = EP.Nx;
    Ny  = EP.Ny;
    Nz  = EP.Nz;

    if(BC.BC0X != BoundaryConditionTypes::Periodic or
       BC.BCNX != BoundaryConditionTypes::Periodic)
    {
#ifdef MPI_PARALLEL
        std::cerr << "ElasticitySolverSpectral::Initialize()\n"
                  << "\tNonperiodic boundary conditions in X-direction are not permitted in MPI parallel mode!" << std::endl;
        exit(1);
#else
        TotalNx *= 2;
        Nx *= 2;
#endif
    }
    if(BC.BC0Y != BoundaryConditionTypes::Periodic or
       BC.BCNY != BoundaryConditionTypes::Periodic)
    {
        TotalNy *= 2;
        Ny *= 2;
    }
    if(BC.BC0Z != BoundaryConditionTypes::Periodic or
       BC.BCNZ != BoundaryConditionTypes::Periodic)
    {
        TotalNz *= 2;
        Nz *= 2;
    }

    dNx = EP.dNx;
    dNy = EP.dNy;
    dNz = EP.dNz;

    Nz2 = Nz/2 + 1;

    dx = EP.dx;

#ifdef MPI_PARALLEL

    ptrdiff_t local_n0      = 0;
    ptrdiff_t local_0_start = 0;

    size_t SIZE = 2*fftw_mpi_local_size_3d(TotalNx, TotalNy, TotalNz, MPI_COMM_WORLD,
                                  &local_n0, &local_0_start);
#else
    size_t SIZE = Nx*Ny*Nz2*2;
#endif

    // Destroy old FFT plans and free allocated memory
    for(int n = 0; n < 9; n++)
    {
        fftw_destroy_plan(ForwardPlanRHS[n]);
        fftw_destroy_plan(BackwardPlanDefGrad[n]);

        fftw_free(RHSandDefGrad[n]);
    }
    for(int n = 0; n < 3; n++)
    {
        fftw_destroy_plan(BackwardPlanU[n]);
        fftw_destroy_plan(ForwardPlanForce[n]);

        fftw_free(UandForce[n]);
    }

    // Reallocate arrays
    for(int n = 0; n < 9; n++)
    {
        RHSandDefGrad[n] = (double *)fftw_malloc(sizeof(double)*SIZE);
    }
    for(int n = 0; n < 3; n++)
    {
        UandForce[n] = (double *)fftw_malloc(sizeof(double)*SIZE);
    }

    // Create new FFT plans
    for(int n = 0; n < 9; n++)
    {
        #ifdef MPI_PARALLEL
        ForwardPlanRHS[n] = fftw_mpi_plan_dft_r2c_3d
                       (TotalNx, TotalNy, TotalNz, RHSandDefGrad[n],
                        reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                        MPI_COMM_WORLD,
                        FFTW_ESTIMATE);
        #else
        ForwardPlanRHS[n] = fftw_plan_dft_r2c_3d
                        (Nx, Ny, Nz, RHSandDefGrad[n],
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         FFTW_ESTIMATE);
        #endif
    }

    for(int n = 0; n < 9; n++)
    {
        #ifdef MPI_PARALLEL
        BackwardPlanDefGrad[n] = fftw_mpi_plan_dft_c2r_3d
                        (TotalNx, TotalNy, TotalNz,
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         RHSandDefGrad[n],
                         MPI_COMM_WORLD,
                         FFTW_ESTIMATE);
        #else
        BackwardPlanDefGrad[n] = fftw_plan_dft_c2r_3d
                        (Nx, Ny, Nz,
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         RHSandDefGrad[n],
                         FFTW_ESTIMATE);
        #endif
    }

    for(int n = 0; n < 3; n++)
    {
        #ifdef MPI_PARALLEL
        BackwardPlanU[n] = fftw_mpi_plan_dft_c2r_3d
                            (TotalNx, TotalNy, TotalNz,
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             UandForce[n],
                             MPI_COMM_WORLD,
                             FFTW_ESTIMATE);

        ForwardPlanForce[n] = fftw_mpi_plan_dft_r2c_3d
                           (TotalNx, TotalNy, TotalNz, UandForce[n],
                            reinterpret_cast<fftw_complex*> (UandForce[n]),
                            MPI_COMM_WORLD,
                            FFTW_ESTIMATE);
        #else
        BackwardPlanU[n] = fftw_plan_dft_c2r_3d
                            (Nx, Ny, Nz,
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             UandForce[n],
                             FFTW_ESTIMATE);

        ForwardPlanForce[n] = fftw_plan_dft_r2c_3d
                            (Nx, Ny, Nz, UandForce[n],
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             FFTW_ESTIMATE);
        #endif
    }
    Info::WriteStandard(thisclassname, "Reinitialized");
}

//++++++++++++ Destructor +++++++++++++++++++++++++

ElasticitySolverSpectral::~ElasticitySolverSpectral(void)
{
    for(int n = 0; n < 9; n++)
    {
        fftw_destroy_plan(ForwardPlanRHS[n]);
        fftw_destroy_plan(BackwardPlanDefGrad[n]);

        fftw_free(RHSandDefGrad[n]);
    }

    for(int n = 0; n < 3; n++)
    {
        fftw_destroy_plan(BackwardPlanU[n]);
        fftw_destroy_plan(ForwardPlanForce[n]);

        fftw_free(UandForce[n]);
    }

#ifdef MPI_PARALLEL
    fftw_mpi_cleanup();
#endif

#ifdef _OPENMP
    fftw_cleanup_threads();
#endif
}

int ElasticitySolverSpectral::Solve(ElasticProperties& EP,
        BoundaryConditions& BC, double dt,
        std::function<bool()> EndSolve)
{
    vStress TargetStress;
    vStrain oldAverageStrain;

    int    IterationsCount = 0;
    double MAXStrainDifference = 0.0;
    double MAXTargetStrainDifference = 0.0;

    EP.AppliedStrain += EP.AppliedStrainRate*dt;

    do // Iteration loop begin
    {
        IterationsCount++;

        MAXStrainDifference = 0.0;
        MAXTargetStrainDifference = 0.0;

        if(EP.ConsiderExternalForces)
        {
            CopyForceDensity(EP);
            ExecuteForwardFFTforces();
        }

        CalculateRHS(EP);
        ExecuteForwardFFT(EP);
        CalculateFourierSolution(EP);
        ExecuteBackwardFFT();
        CalculateTargetStress(EP, MAXStrainDifference, TargetStress);
        oldAverageStrain = EP.AverageStrain;

        ApplyMechanicalBC(EP, TargetStress);
        SetElasticProperties(EP);

        MAXTargetStrainDifference = max(MAXTargetStrainDifference,(EP.AverageStrain - oldAverageStrain).norm());

        #ifdef MPI_PARALLEL
        double sMAXTargetStrainDifference = MAXTargetStrainDifference;
        MPI_Allreduce(&sMAXTargetStrainDifference,&MAXTargetStrainDifference,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        double sMAXStrainDifference = MAXStrainDifference;
        MPI_Allreduce(&sMAXStrainDifference,&MAXStrainDifference,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        #endif
        if(IterationsCount > EP.MAXIterations)
        {
            std::string message = "Maximum number of iterations (" + std::to_string(EP.MAXIterations) + ") reached\n";
            message += Info::GetStandard("Strains converged to",
                       Info::to_string_with_precision(max(MAXTargetStrainDifference,MAXStrainDifference)));
            Info::WriteWarning(message, thisclassname, "Solve()");
            break;
        }
#ifdef DEBUG
        //std::cout << "Iteration: " << IterationsCount << "\tMAXStrainDifference: " << max(MAXTargetStrainDifference,MAXStrainDifference) << "\n";
#endif
    } // Iteration loop end
    while(not EndSolve() and (MAXStrainDifference > EP.StrainAccuracy or MAXTargetStrainDifference > 10.0*EP.StrainAccuracy));

    ExecuteBackwardFFTdisplacements();
    CopyDisplacements(EP);

    EP.SetBoundaryConditions(BC);

    return IterationsCount;
}

void ElasticitySolverSpectral::CalculateRHS(ElasticProperties& EP)
{
    dMatrix6x6 Cij = EP.MAXElasticConstants;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.DeformationGradientsTotal, 0,)
    {

        dMatrix3x3 locRHStensor;
        vStress locStress = EP.EffectiveElasticConstants(i,j,k)*EP.ElasticStrains(i,j,k);

        locRHStensor = (Cij*EP.StrainSmall(EP.DeformationGradientsTotal(i,j,k))).tensor()
                     -  locStress.tensor();

        const long int xyz = k + 2*Nz2*(j + Ny*i);

        RHSandDefGrad[0][xyz] = locRHStensor(0,0);
        RHSandDefGrad[1][xyz] = locRHStensor(0,1);
        RHSandDefGrad[2][xyz] = locRHStensor(0,2);

        RHSandDefGrad[3][xyz] = locRHStensor(1,0);
        RHSandDefGrad[4][xyz] = locRHStensor(1,1);
        RHSandDefGrad[5][xyz] = locRHStensor(1,2);

        RHSandDefGrad[6][xyz] = locRHStensor(2,0);
        RHSandDefGrad[7][xyz] = locRHStensor(2,1);
        RHSandDefGrad[8][xyz] = locRHStensor(2,2);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    if(Ny > EP.Ny)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.DeformationGradientsTotal, 0,)
        for(int n = 0; n < 9; n++)
        {
            int y = Ny - j - 1;
            RHSandDefGrad[n][k + 2*Nz2*(y + Ny*i)] = RHSandDefGrad[n][k + 2*Nz2*(j + Ny*i)];
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    if(Nz > EP.Nz)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.DeformationGradientsTotal, 0,)
        for(int n = 0; n < 9; n++)
        {
            int z = Nz - k - 1;
            RHSandDefGrad[n][z + 2*Nz2*(j + Ny*i)] = RHSandDefGrad[n][k + 2*Nz2*(j + Ny*i)];
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
#else
    if(Nx > EP.Nx)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.DeformationGradientsTotal, 0,)
        for(int n = 0; n < 9; n++)
        {
            int x = Nx - i - 1;
            RHSandDefGrad[n][k + 2*Nz2*(j + Ny*x)] = RHSandDefGrad[n][k + 2*Nz2*(j + Ny*i)];
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    if(Ny > EP.Ny)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.DeformationGradientsTotal, 0,)
        for(int n = 0; n < 9; n++)
        {
            int y = Ny - j - 1;
            RHSandDefGrad[n][k + 2*Nz2*(y + Ny*i)] = RHSandDefGrad[n][k + 2*Nz2*(j + Ny*i)];
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    if(Nz > EP.Nz)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.DeformationGradientsTotal, 0,)
        for(int n = 0; n < 9; n++)
        {
            int z = Nz - k - 1;
            RHSandDefGrad[n][z + 2*Nz2*(j + Ny*i)] = RHSandDefGrad[n][k + 2*Nz2*(j + Ny*i)];
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
#endif
}

void ElasticitySolverSpectral::CalculateFourierSolution(ElasticProperties& EP)
{
    dMatrix6x6 Cij = EP.MAXElasticConstants;

    double Norm = 1.0/double(TotalNx*TotalNy*TotalNz);

    double DPi_Nx = 2.0*Pi/double(TotalNx*dx);
    double DPi_Ny = 2.0*Pi/double(TotalNy*dx);
    double DPi_Nz = 2.0*Pi/double(TotalNz*dx);

    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Nx ; i++)
    for(int j = 0; j < Ny ; j++)
    for(int k = 0; k < Nz2; k++)
    {
        long int XYZ = k + Nz2*(j + Ny*i);
#ifdef MPI_PARALLEL
        long int ii = i + OffsetX;
        double Qx = DPi_Nx*(ii*(ii <= TotalNx/2) - (TotalNx-ii)*(ii > TotalNx/2));
        long int jj = j + OffsetY;
        double Qy = DPi_Ny*(jj*(jj <= TotalNy/2) - (TotalNy-jj)*(jj > TotalNy/2));
        long int kk = k + OffsetZ;
        double Qz = DPi_Nz*(kk*(kk <= TotalNz/2) - (TotalNz-kk)*(kk > TotalNz/2));
#else
        double Qx = DPi_Nx*(i*(i <= Nx/2) - (Nx-i)*(i > Nx/2));
        double Qy = DPi_Ny*(j*(j <= Ny/2) - (Ny-j)*(j > Ny/2));
        double Qz = DPi_Nz*(k*(k <= Nz/2) - (Nz-k)*(k > Nz/2));
#endif

        if(Smooth)
        {
            Qx *= (1.0 + cos(fabs(Qx*dx)))/2.0;
            Qy *= (1.0 + cos(fabs(Qy*dx)))/2.0;
            Qz *= (1.0 + cos(fabs(Qz*dx)))/2.0;
        }

        complex<double> rhsX = -I*(Qx*reinterpret_cast<complex<double>*>(RHSandDefGrad[0])[XYZ] +
                                   Qy*reinterpret_cast<complex<double>*>(RHSandDefGrad[1])[XYZ] +
                                   Qz*reinterpret_cast<complex<double>*>(RHSandDefGrad[2])[XYZ]);
        complex<double> rhsY = -I*(Qx*reinterpret_cast<complex<double>*>(RHSandDefGrad[3])[XYZ] +
                                   Qy*reinterpret_cast<complex<double>*>(RHSandDefGrad[4])[XYZ] +
                                   Qz*reinterpret_cast<complex<double>*>(RHSandDefGrad[5])[XYZ]);
        complex<double> rhsZ = -I*(Qx*reinterpret_cast<complex<double>*>(RHSandDefGrad[6])[XYZ] +
                                   Qy*reinterpret_cast<complex<double>*>(RHSandDefGrad[7])[XYZ] +
                                   Qz*reinterpret_cast<complex<double>*>(RHSandDefGrad[8])[XYZ]);

        if(EP.ConsiderExternalForces)
        {
            rhsX += reinterpret_cast<complex<double>*>(UandForce[0])[XYZ];
            rhsY += reinterpret_cast<complex<double>*>(UandForce[1])[XYZ];
            rhsZ += reinterpret_cast<complex<double>*>(UandForce[2])[XYZ];
        }

        double a11 = (Cij(0,0)*Qx*Qx + 2.0*Cij(0,5)*Qx*Qy + Cij(5,5)*Qy*Qy +
                  2.0*Cij(0,4)*Qx*Qz + 2.0*Cij(4,5)*Qy*Qz + Cij(4,4)*Qz*Qz);

        double a21 = (Cij(0,5)*Qx*Qx + Cij(0,1)*Qx*Qy + Cij(5,5)*Qx*Qy +
                      Cij(1,5)*Qy*Qy + Cij(0,3)*Qx*Qz + Cij(4,5)*Qx*Qz +
                      Cij(1,4)*Qy*Qz + Cij(3,5)*Qy*Qz + Cij(3,4)*Qz*Qz);

        double a31 = (Cij(0,4)*Qx*Qx + Cij(0,3)*Qx*Qy + Cij(4,5)*Qx*Qy +
                      Cij(3,5)*Qy*Qy + Cij(0,2)*Qx*Qz + Cij(4,4)*Qx*Qz +
                      Cij(2,5)*Qy*Qz + Cij(3,4)*Qy*Qz + Cij(2,4)*Qz*Qz);

        double a12 = a21;

        double a22 = (Cij(5,5)*Qx*Qx + 2.0*Cij(1,5)*Qx*Qy + Cij(1,1)*Qy*Qy +
                  2.0*Cij(3,5)*Qx*Qz + 2.0*Cij(1,3)*Qy*Qz + Cij(3,3)*Qz*Qz);

        double a32 = (Cij(4,5)*Qx*Qx + Cij(1,4)*Qx*Qy + Cij(3,5)*Qx*Qy +
                      Cij(1,3)*Qy*Qy + Cij(2,5)*Qx*Qz + Cij(3,4)*Qx*Qz +
                      Cij(1,2)*Qy*Qz + Cij(3,3)*Qy*Qz + Cij(2,3)*Qz*Qz);

        double a13 = a31;

        double a23 = a32;

        double a33 = (Cij(4,4)*Qx*Qx + 2.0*Cij(3,4)*Qx*Qy + Cij(3,3)*Qy*Qy +
                  2.0*Cij(2,4)*Qx*Qz + 2.0*Cij(2,3)*Qy*Qz + Cij(2,2)*Qz*Qz);

        double denominator = (-a13*a22*a31 + a12*a23*a31 + a13*a21*a32 -
                               a11*a23*a32 - a12*a21*a33 + a11*a22*a33);

        if(std::abs(denominator) > DBL_EPSILON and std::abs(denominator) < DBL_MAX)
        {
            denominator = 1.0/denominator;
        }
        else
        {
            denominator = 0.0;
        }

        complex<double> locUrcX = (-a23*a32*rhsX + a22*a33*rhsX + a13*a32*rhsY -
                                    a12*a33*rhsY - a13*a22*rhsZ + a12*a23*rhsZ)*denominator*Norm;

        complex<double> locUrcY = ( a23*a31*rhsX - a21*a33*rhsX - a13*a31*rhsY +
                                    a11*a33*rhsY + a13*a21*rhsZ - a11*a23*rhsZ)*denominator*Norm;

        complex<double> locUrcZ = (-a22*a31*rhsX + a21*a32*rhsX + a12*a31*rhsY -
                                    a11*a32*rhsY - a12*a21*rhsZ + a11*a22*rhsZ)*denominator*Norm;

        reinterpret_cast<complex<double>*>(UandForce[0])[XYZ] = locUrcX;
        reinterpret_cast<complex<double>*>(UandForce[1])[XYZ] = locUrcY;
        reinterpret_cast<complex<double>*>(UandForce[2])[XYZ] = locUrcZ;

        reinterpret_cast<complex<double>*>(RHSandDefGrad[0])[XYZ] = I*(Qx*locUrcX);
        reinterpret_cast<complex<double>*>(RHSandDefGrad[1])[XYZ] = I*(Qy*locUrcX);
        reinterpret_cast<complex<double>*>(RHSandDefGrad[2])[XYZ] = I*(Qz*locUrcX);
        reinterpret_cast<complex<double>*>(RHSandDefGrad[3])[XYZ] = I*(Qx*locUrcY);
        reinterpret_cast<complex<double>*>(RHSandDefGrad[4])[XYZ] = I*(Qy*locUrcY);
        reinterpret_cast<complex<double>*>(RHSandDefGrad[5])[XYZ] = I*(Qz*locUrcY);
        reinterpret_cast<complex<double>*>(RHSandDefGrad[6])[XYZ] = I*(Qx*locUrcZ);
        reinterpret_cast<complex<double>*>(RHSandDefGrad[7])[XYZ] = I*(Qy*locUrcZ);
        reinterpret_cast<complex<double>*>(RHSandDefGrad[8])[XYZ] = I*(Qz*locUrcZ);
    }
}

void ElasticitySolverSpectral::CalculateTargetStress(ElasticProperties& EP,
                            double& MAXStrainDeviation, vStress& TargetStress)
{
    vStress locAverageStress;

    EP.GetAverageDeformationGradient();

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.DeformationGradientsTotal,0,reduction(max:MAXStrainDeviation) reduction(vStressSUM:locAverageStress))
    {
        const long int xyz = k + 2*Nz2*(j + Ny*i);

        dMatrix3x3 locDefGrad = EP.AverageDeformationGradient;
        locDefGrad(0,0) += RHSandDefGrad[0][xyz];
        locDefGrad(0,1) += RHSandDefGrad[1][xyz];
        locDefGrad(0,2) += RHSandDefGrad[2][xyz];

        locDefGrad(1,0) += RHSandDefGrad[3][xyz];
        locDefGrad(1,1) += RHSandDefGrad[4][xyz];
        locDefGrad(1,2) += RHSandDefGrad[5][xyz];

        locDefGrad(2,0) += RHSandDefGrad[6][xyz];
        locDefGrad(2,1) += RHSandDefGrad[7][xyz];
        locDefGrad(2,2) += RHSandDefGrad[8][xyz];

        vStrain locStrain;

        locStrain = EP.StrainSmall(locDefGrad);
        locAverageStress += EP.EffectiveElasticConstants(i, j, k)*
                            (locStrain - EP.StressFreeStrains(i, j, k));

        MAXStrainDeviation = max(MAXStrainDeviation, (EP.TotalStrains(i,j,k) - locStrain).norm());
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    EP.AverageStress = locAverageStress;
#ifdef MPI_PARALLEL
    MPI_Allreduce(locAverageStress.data(),EP.AverageStress.data(),6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double sMAXStrainDeviation = MAXStrainDeviation;
    MPI_Allreduce(&sMAXStrainDeviation,&MAXStrainDeviation,1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    EP.AverageStress /= double(TotalNx*TotalNy*TotalNz);

    TargetStress[0] = -EP.AverageStress[0]*fabs(1.0 - EP.AppliedStressMask[0]) + (EP.AppliedStress[0]-EP.AverageStress[0])*fabs(EP.AppliedStressMask[0]);
    TargetStress[1] = -EP.AverageStress[1]*fabs(1.0 - EP.AppliedStressMask[1]) + (EP.AppliedStress[1]-EP.AverageStress[1])*fabs(EP.AppliedStressMask[1]);
    TargetStress[2] = -EP.AverageStress[2]*fabs(1.0 - EP.AppliedStressMask[2]) + (EP.AppliedStress[2]-EP.AverageStress[2])*fabs(EP.AppliedStressMask[2]);
    TargetStress[3] = -EP.AverageStress[3]*fabs(1.0 - EP.AppliedStressMask[3]) + (EP.AppliedStress[3]-EP.AverageStress[3])*fabs(EP.AppliedStressMask[3]);
    TargetStress[4] = -EP.AverageStress[4]*fabs(1.0 - EP.AppliedStressMask[4]) + (EP.AppliedStress[4]-EP.AverageStress[4])*fabs(EP.AppliedStressMask[4]);
    TargetStress[5] = -EP.AverageStress[5]*fabs(1.0 - EP.AppliedStressMask[5]) + (EP.AppliedStress[5]-EP.AverageStress[5])*fabs(EP.AppliedStressMask[5]);
}

void ElasticitySolverSpectral::SetElasticProperties(ElasticProperties& EP)
{
    EP.GetAverageDeformationGradient();

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.DeformationGradientsTotal,0,)
    {
        const long int xyz = k + 2*Nz2*(j + Ny*i);

        EP.DeformationGradientsTotal(i,j,k) = EP.AverageDeformationGradient;

        EP.DeformationGradientsTotal(i,j,k)(0,0) += RHSandDefGrad[0][xyz];
        EP.DeformationGradientsTotal(i,j,k)(0,1) += RHSandDefGrad[1][xyz];
        EP.DeformationGradientsTotal(i,j,k)(0,2) += RHSandDefGrad[2][xyz];

        EP.DeformationGradientsTotal(i,j,k)(1,0) += RHSandDefGrad[3][xyz];
        EP.DeformationGradientsTotal(i,j,k)(1,1) += RHSandDefGrad[4][xyz];
        EP.DeformationGradientsTotal(i,j,k)(1,2) += RHSandDefGrad[5][xyz];

        EP.DeformationGradientsTotal(i,j,k)(2,0) += RHSandDefGrad[6][xyz];
        EP.DeformationGradientsTotal(i,j,k)(2,1) += RHSandDefGrad[7][xyz];
        EP.DeformationGradientsTotal(i,j,k)(2,2) += RHSandDefGrad[8][xyz];

        vStress locStress = EP.EffectiveElasticConstants(i,j,k)*EP.ElasticStrains(i,j,k);

        EP.Stresses(i,j,k) = locStress;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySolverSpectral::ExecuteForwardFFT(ElasticProperties& EP)
{
    fftw_execute(ForwardPlanRHS[0]);
    fftw_execute(ForwardPlanRHS[1]);
    fftw_execute(ForwardPlanRHS[2]);
    fftw_execute(ForwardPlanRHS[3]);
    fftw_execute(ForwardPlanRHS[4]);
    fftw_execute(ForwardPlanRHS[5]);
    fftw_execute(ForwardPlanRHS[6]);
    fftw_execute(ForwardPlanRHS[7]);
    fftw_execute(ForwardPlanRHS[8]);
}

void ElasticitySolverSpectral::ExecuteBackwardFFT(void)
{
    fftw_execute(BackwardPlanDefGrad[0]);
    fftw_execute(BackwardPlanDefGrad[1]);
    fftw_execute(BackwardPlanDefGrad[2]);
    fftw_execute(BackwardPlanDefGrad[3]);
    fftw_execute(BackwardPlanDefGrad[4]);
    fftw_execute(BackwardPlanDefGrad[5]);
    fftw_execute(BackwardPlanDefGrad[6]);
    fftw_execute(BackwardPlanDefGrad[7]);
    fftw_execute(BackwardPlanDefGrad[8]);
}

void ElasticitySolverSpectral::ExecuteBackwardFFTdisplacements(void)
{
    fftw_execute(BackwardPlanU[0]);
    fftw_execute(BackwardPlanU[1]);
    fftw_execute(BackwardPlanU[2]);
}

void ElasticitySolverSpectral::ExecuteForwardFFTforces(void)
{
    fftw_execute(ForwardPlanForce[0]);
    fftw_execute(ForwardPlanForce[1]);
    fftw_execute(ForwardPlanForce[2]);
}

void ElasticitySolverSpectral::CopyForceDensity(ElasticProperties& EP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.ForceDensity, 0,)
    {
        const long int xyz = k + 2*Nz2*(j + Ny*i);

        UandForce[0][xyz] = EP.ForceDensity(i,j,k)[0];
        UandForce[1][xyz] = EP.ForceDensity(i,j,k)[1];
        UandForce[2][xyz] = EP.ForceDensity(i,j,k)[2];
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySolverSpectral::CopyDisplacements(ElasticProperties& EP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.Displacements,0,)
    {
        const long int xyz = k + 2*Nz2*(j + Ny*i);

        EP.Displacements(i,j,k)[0] = UandForce[0][xyz];
        EP.Displacements(i,j,k)[1] = UandForce[1][xyz];
        EP.Displacements(i,j,k)[2] = UandForce[2][xyz];
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySolverSpectral::ApplyMechanicalBC(ElasticProperties& EP, vStress TargetStress)
{
    size_t size = 6;
    for(int n = 0; n < 6; n++)
    if(EP.AppliedStrainMask[n])
    {
        size--;
    }

    dVectorN locRHS(size);
    dMatrixNxN locLHS(size);

    int x = 0;
    for(int n = 0; n < 6; n++)
    if(!EP.AppliedStrainMask[n])
    {
        locRHS[x] = TargetStress[n];
        int y = 0;
        for(int m = 0; m < 6; m++)
        if(!EP.AppliedStrainMask[m])
        {
            locLHS(x,y) = EP.AverageElasticConstants(n,m);
            y++;
        }
        x++;
    }

    dVectorN Solution = locLHS.inverted()*locRHS;

    x = 0;
    for(int n = 0; n < 6; n++)
    if(!EP.AppliedStrainMask[n])
    {
        EP.AverageStrain[n] += Solution[x];
        x++;
    }
    else
    {
        EP.AverageStrain[n] = EP.AppliedStrain[n];
    }

    if(EP.KeepAspectRatio)
    {
        double trace = (1.0/(dNx+dNy+dNz))*(EP.AverageStrain[0]*dNx + EP.AverageStrain[1]*dNy + EP.AverageStrain[2]*dNz);
        EP.AverageStrain[0] = trace * dNx;
        EP.AverageStrain[1] = trace * dNy;
        EP.AverageStrain[2] = trace * dNz;
    }

    if(EP.KeepAspectRatio or EP.PreventShear)
    {
        EP.AverageStrain[3] = 0.0;
        EP.AverageStrain[4] = 0.0;
        EP.AverageStrain[5] = 0.0;
    }
}

} // namespace openphase
