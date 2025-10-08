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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Philipp Engels;
 *                         Muhammad Adil Ali; Hesham Salama
 *
 */

#ifndef TOOLS_H
#define TOOLS_H

#include "Base/Includes.h"
#include "Crystallography.h"

namespace openphase
{

class PhaseField;
class Orientations;

class OP_EXPORTS Tools
{
 public:
    /*static void WriteAbaqusInput(const PhaseField& Phi,                       /// Writes Simulia Abaqus input file (Maintainer: PE)
                                    const Orientations& OR);
    */
    static void Decompose(dMatrix3x3& Def, dMatrix3x3& Rot, dMatrix3x3& Stretch);/// Decomposes deformation gradient tensor into pure Rotation and pure Stretch
    static void Eigensystem(dMatrix3x3& M, dMatrix3x3& Eigenvectors,
                                                       dMatrix3x3& Eigenvalues);/// M MUST BE _SYMMETRIC_, calculates eigenvalues and (column) eigenvectors of M
    static void Eigensystem2(dMatrix3x3& M, dMatrix3x3& Eigenvectors,
                                                       dMatrix3x3& Eigenvalues,
                                             bool calculateEigenVectors = true);/// Non-iterative! M MUST BE _SYMMETRIC_, calculates eigenvalues and (column) eigenvectors of M
    static void Eigensystem3(dMatrix3x3& M, dMatrix3x3& Eigenvectors,
                                                       dMatrix3x3& Eigenvalues); /// M MUST BE _SYMMETRIC_. Jacobian method. Calculates eigenvalues and (column) eigenvectors of M.
    static dMatrix3x3 sqrtM3x3(dMatrix3x3 M);                                   /// M MUST BE _SYMMETRIC_, calculates sqrt(M)
    static dMatrix3x3 qurtM3x3(dMatrix3x3 M);                                   /// M MUST BE _SYMMETRIC_, calculates fourth order (quartic) root of M
    static dMatrix3x3 logM3x3(dMatrix3x3 M);                                    /// M MUST BE _SYMMETRIC_, calculates log(M)
    static dMatrix3x3 expM3x3(dMatrix3x3 M);                                    /// M MUST BE _SYMMETRIC_, calculates exp(M)
    static dMatrix3x3 AlignBaseAxes(dMatrix3x3 M);                              /// Aligns X, Y and Z axes of the deformation gradient tensor with the base vectors.

    /////////// Extract Rotation Matrix from Total Deformation ////////
    static void jacobiRotate(dMatrix3x3 &A, dMatrix3x3 &R, int p, int q);
    static void eigenDecomposition(const dMatrix3x3 &A, dMatrix3x3 &eigenVecs, dVector3 &eigenVals);
    static void rotationMatrixIrving(const dMatrix3x3 &A, dMatrix3x3 &R);
    static void GetAxisAngleFromRotationMatrix(const dMatrix3x3 RotMatrix, double& Angle, dVector3& Axis);
    static void getAxisAngle(const dMatrix3x3& TransformationMatrix, dVector3& Axis, double &Angle);
    static dVector3 IPFColor(size_t sd, EulerAngles& tempEuler, Crystallography& CR);

    static void ExtractRotation(dMatrix3x3 &Mat,  Quaternion& Quat, const size_t maxIter);
    template <typename T>
    static void Smooth(Storage3D<T, 0> &Field, size_t SmoothIterations = 1)
    {
        int Nx = Field.sizeX();
        int Ny = Field.sizeY();
        int Nz = Field.sizeZ();
        int dNx = Field.dNx();
        int dNy = Field.dNy();
        int dNz = Field.dNz();
        double Stencil[3][3][3] = {{{1.0/64.0,   1.0/32.0, 1.0/64.0},
                                   {1.0/32.0,   1.0/16.0, 1.0/32.0},
                                   {1.0/64.0,   1.0/32.0, 1.0/64.0}},

                                  {{1.0/32.0,   1.0/16.0, 1.0/32.0},
                                   {1.0/16.0,   1.0/8.0, 1.0/16.0},
                                   {1.0/32.0,   1.0/16.0, 1.0/32.0}},

                                  {{1.0/64.0,   1.0/32.0, 1.0/64.0},
                                   {1.0/32.0,   1.0/16.0, 1.0/32.0},
                                   {1.0/64.0,   1.0/32.0, 1.0/64.0}}};

        std::function<bool(int,int,int)> condition = [Nx,Ny,Nz](int u,int v,int w)
        {
            return u < Nx and v < Ny and w < Nz and u >= 0 and v >= 0 and w >= 0;
        };
        Storage3D<T, 0> FieldTmp;
        FieldTmp.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, 1);
        for(size_t start = 0; start < SmoothIterations; start++)
        {
            if(start)
            {
                FieldTmp.Reallocate(Nx, Ny, Nz);
            }
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, Field,Field.Bcells(),)
            {
                        double tmpStencil[3][3][3] = {{{0.0,   0.0,  0.0},
                                                     {0.0,   0.0,  0.0},
                                                     {0.0,   0.0,  0.0}},

                                                     {{0.0,   0.0,  0.0},
                                                      {0.0,   0.0,  0.0},
                                                      {0.0,   0.0,  0.0}},

                                                     {{0.0,   0.0,  0.0},
                                                      {0.0,   0.0,  0.0},
                                                      {0.0,   0.0,  0.0}}};
                        double WS = 0;
                        for(int a = -dNx; a <= +dNx; a++)
                        for(int b = -dNy; b <= +dNy; b++)
                        for(int c = -dNz; c <= +dNz; c++)
                        {
                            if (condition(i+a,j+b,k+c))
                            {
                                WS += Stencil[a+1][b+1][c+1];
                            }
                        }
                        for(int a = -dNx; a <= +dNx; a++)
                        for(int b = -dNy; b <= +dNy; b++)
                        for(int c = -dNz; c <= +dNz; c++)
                        {
                            if (condition(i+a,j+b,k+c))
                            {
                                tmpStencil[a+1][b+1][c+1] = Stencil[a+1][b+1][c+1]/WS;
                            }
                            else
                            {
                                tmpStencil[a+1][b+1][c+1] = Stencil[a+1][b+1][c+1];
                            }
                        }

                        for(int a = -dNx; a <= +dNx; a++)
                        for(int b = -dNy; b <= +dNy; b++)
                        for(int c = -dNz; c <= +dNz; c++)
                        {
                            if(condition(i+a,j+b,k+c))
                            {
                                FieldTmp(i,j,k) += Field(i+a,j+b,k+c) * tmpStencil[a+1][b+1][c+1];
                            }
                        }
            }
            OMP_PARALLEL_STORAGE_LOOP_END
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,0,)
        {
                Field(i,j,k) = FieldTmp(i,j,k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        }
    }

    template <typename T>
    static T Average(Storage3D<T, 0> &Field)
    {
        int Nx = Field.sizeX();
        int Ny = Field.sizeY();
        int Nz = Field.sizeZ();

        size_t Nthreads = 1;

        #ifdef _OPENMP
        Nthreads = omp_get_max_threads();
        #endif

        std::vector<T> avgThreads;
        avgThreads.resize(Nthreads);
        T AverageField;

        double Norm = 1.0/double(Nx*Ny*Nz);

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,0,)
        {
            size_t thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif
            avgThreads[thread] += Field(i,j,k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        for(size_t thread = 0; thread < Nthreads; thread++)
            AverageField += avgThreads[thread];

        #ifdef MPI_PARALLEL
        for(int m = 0; m < 6; m++)
        {
            double tmp = AverageField[m];
            MPI_Reduce(&tmp,&AverageField[m],1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
        }
        #endif
        return AverageField*Norm;
    }
 protected:
 private:
};
}
#endif
