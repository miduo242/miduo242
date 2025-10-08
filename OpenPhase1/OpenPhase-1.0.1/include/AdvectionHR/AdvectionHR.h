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
 *   File created :   2015
 *   Main contributors :   Philipp Engels; Marvin Tegeler; Raphael Schiedung
 *
 */

#ifndef ADVECTIONHR_H
#define ADVECTIONHR_H

#include "Settings.h"
#include "Base/CommonFunctions.h"
#include "Base/Includes.h"
#include "Base/UserInterface.h"
#include "BoundaryConditions.h"
#include "Composition.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "Info.h"
#include "Orientations.h"
#include "PhaseField.h"
#include "Temperature.h"
#include "Temperature.h"
#include "Velocities.h"

namespace openphase
{

/*enum class AdvectionModes
{
    Minmod,
    MonotonizedCentral,
    Superbee,
    //LaxWendroff,
    Upwind
};*/

template <typename T>
struct has_size
{
    static const int value = 0;
};

template <>
struct has_size<double>
{
    static const int value = 1;
};

template <>
struct has_size<dVector3>
{
    static const int value = 3;
};

template <>
struct has_size<Quaternion>
{
    static const int value = 4;
};

template <>
struct has_size<vStress>
{
    static const int value = 6;
};

template <>
struct has_size<dMatrix3x3>
{
    static const int value = 9;
};

inline double Slope_Limiter_Minmod(const double a, const double b)
{
    return 0.5*((a > 0) - (a < 0)+(b > 0) - (b < 0))*std::min(std::fabs(a),std::fabs(b));
}

inline double Slope_Limiter_MC(const double a, const double b)
{
    return 0.5*((a > 0) - (a < 0)+(b > 0) - (b < 0))*std::min(0.5*std::fabs(a+b),std::min(2.0*std::fabs(a),2.0*std::fabs(b)));
}

inline double Slope_Limiter_Superbee(const double a, const double b)
{
    return 0.5*((a > 0) - (a < 0)+(b > 0) - (b < 0))*std::max(std::min(2.0*std::fabs(a),std::fabs(b)),std::min(std::fabs(a),2.0*std::fabs(b)));
}

inline double Slope_Limiter_Upwind(const double a, const double b)
{
    return 0.0;
}
/*inline double Slope_Limiter_LaxWendroff(const double a, const double b)
{
    return 1.0;
}*/

template<class T, int Rank>
class AdvectionMethod
{
 public:
    static void Advection(Storage3D<T, Rank>& Field,
                   const BoundaryConditions& BC, const Velocities& Vel, const int direction,
                   const double dt, const double dx, double (*Limiter)(const double,const double), bool compressive = false)
    {
        Info::WriteExit("Advection of template Storage3D<T, rank> not yet implemented.",
                        "AdvectionMethod", "Advection()");
        raise(SIGABRT);
    }
};

template<class T>
class AdvectionMethod<T,0>
{
 public:
    static void Advection(Storage3D<T, 0>& Field,
                          const BoundaryConditions& BC,
                          const Velocities& Vel,
                          const int direction,
                          const double dt,
                          const double dx,
                          double (*Limiter)(const double,const double),
                          bool compressible = false)
    {
        Storage3D<T, 0> FieldBackup(Field);
        Storage3D<T, 0> FieldUpdate;
        FieldUpdate.Allocate(Field);

        int CFL = 0;
        double locdt = dt;
        std::vector<int> dir(3);
        dir[0] = (direction == 0)?1:0;
        dir[1] = (direction == 1)?1:0;
        dir[2] = (direction == 2)?1:0;
        const int offset = Field.Bcells()-2;
        int its = 1;
        do
        {
            CFL = 0;
            for (int it = 0; it < its; ++it)
            {
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,offset,reduction(+:CFL))
                for (int n = 0; n < has_size<T>::value; ++n)
                {
                    double q   = Field(i,j,k)[n];
                    double qp  = Field(i+dir[0],j+dir[1],k+dir[2])[n];
                    double qm  = Field(i-dir[0],j-dir[1],k-dir[2])[n];
                    double qpp = Field(i+2*dir[0],j+2*dir[1],k+2*dir[2])[n];
                    double qmm = Field(i-2*dir[0],j-2*dir[1],k-2*dir[2])[n];
                    double v   = Vel.Average(i,j,k)[direction];
                    double vp  = v;
                    double vm  = v;
                    if(compressible)
                    {
                        vp = Vel.Average(i+dir[0],j+dir[1],k+dir[2])[direction];
                        vm = Vel.Average(i-dir[0],j-dir[1],k-dir[2])[direction];
                    }
                    double L0  = 0.5*Limiter(q-qm,qp-q);
                    double Lp1 = 0.5*Limiter(qp-q,qpp-qp);
                    double Lm1 = 0.5*Limiter(qm-qmm,q-qm);
                    double ap  = 0.5*(vp+v);
                    double am  = 0.5*(vm+v);
                    double Fp  = std::max(ap,0.0)*(q+L0)
                               + std::min(ap,0.0)*(qp-Lp1);
                    double Fm  = std::max(am,0.0)*(qm+Lm1)
                               + std::min(am,0.0)*(q-L0);

                    if (std::fabs(ap)*locdt > 0.49*dx or std::fabs(am)*locdt > 0.49*dx)
                    {
                        CFL = 1;
                    }
                    FieldUpdate(i,j,k)[n] = Field(i,j,k)[n] - (Fp-Fm)*locdt/dx;
                }
                OMP_PARALLEL_STORAGE_LOOP_END
                if (CFL == 0)
                {
                    Field = FieldUpdate;
                    BC.SetX(Field);
                    BC.SetY(Field);
                    BC.SetZ(Field);
                }
                else
                {
                    locdt *= 0.5;
                    Field = FieldBackup;
                }
            }
            its *= 2;
        }
        while (CFL > 0);
    }
};

template<int Rank>
class AdvectionMethod<double,Rank>
{
 public:
    static void Advection(Storage3D<double, Rank>& Field,
                          const BoundaryConditions& BC,
                          const Velocities& Vel,
                          const int direction,
                          const double dt,
                          const double dx,
                          double (*Limiter)(const double,const double),
                          bool compressible = false)
    {
        const size_t totsize = Field(0,0,0).size();

        Storage3D<double, Rank> FieldBackup(Field);
        Storage3D<double, Rank> FieldUpdate;
        FieldUpdate.Allocate(Field);

        int CFL = 0;
        double locdt = dt;
        std::vector<int> dir(3);
        dir[0] = (direction == 0)?1:0;
        dir[1] = (direction == 1)?1:0;
        dir[2] = (direction == 2)?1:0;
        const int offset = Field.Bcells()-2;
        int its = 1;
        do
        {
            CFL = 0;
            for (int it = 0; it < its; ++it)
            {
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,offset,reduction(+:CFL))
                for (size_t n = 0; n < totsize; ++n)
                {
                    double q   = Field(i,j,k)[n];
                    double qp  = Field(i+dir[0],j+dir[1],k+dir[2])[n];
                    double qm  = Field(i-dir[0],j-dir[1],k-dir[2])[n];
                    double qpp = Field(i+2*dir[0],j+2*dir[1],k+2*dir[2])[n];
                    double qmm = Field(i-2*dir[0],j-2*dir[1],k-2*dir[2])[n];
                    double v   = Vel.Average(i,j,k)[direction];
                    double vp  = v;
                    double vm  = v;
                    if(compressible)
                    {
                        vp = Vel.Average(i+dir[0],j+dir[1],k+dir[2])[direction];
                        vm = Vel.Average(i-dir[0],j-dir[1],k-dir[2])[direction];
                    }
                    double L0  = 0.5*Limiter(q-qm,qp-q);
                    double Lp1 = 0.5*Limiter(qp-q,qpp-qp);
                    double Lm1 = 0.5*Limiter(qm-qmm,q-qm);
                    double ap  = 0.5*(vp+v);
                    double am  = 0.5*(vm+v);
                    double Fp  = std::max(ap,0.0)*(q+L0)
                               + std::min(ap,0.0)*(qp-Lp1);
                    double Fm  = std::max(am,0.0)*(qm+Lm1)
                               + std::min(am,0.0)*(q-L0);

                    if (std::fabs(ap)*locdt > 0.49*dx or std::fabs(am)*locdt > 0.49*dx)
                    {
                        CFL = 1;
                    }
                    FieldUpdate(i,j,k)[n] = Field(i,j,k)[n] - (Fp-Fm)*locdt/dx;
                }
                OMP_PARALLEL_STORAGE_LOOP_END
                if (CFL == 0)
                {
                    Field = FieldUpdate;
                    BC.SetX(Field);
                    BC.SetY(Field);
                    BC.SetZ(Field);
                }
                else
                {
                    locdt *= 0.5;
                    Field = FieldBackup;
                }
            }
            its *= 2;
        }
        while (CFL>0);
    }
};

template<>
class AdvectionMethod<double,0>
{
 public:
    static void Advection(Storage3D<double, 0>& Field,
                          const BoundaryConditions& BC,
                          const Velocities& Vel,
                          const int direction,
                          const double dt,
                          const double dx,
                          double (*Limiter)(const double,const double),
                          bool compressible = false)
    {
        Storage3D<double, 0> FieldBackup(Field);
        Storage3D<double, 0> FieldUpdate;
        FieldUpdate.Allocate(Field);

        int CFL = 0;
        double locdt = dt;
        std::vector<int> dir(3);
        dir[0] = (direction == 0)?1:0;
        dir[1] = (direction == 1)?1:0;
        dir[2] = (direction == 2)?1:0;
        const int offset = Field.Bcells()-2;
        int its = 1;
        do
        {
            CFL = 0;
            for (int it = 0; it < its; ++it)
            {
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,offset,reduction(+:CFL))
                {
                    double q   = Field(i,j,k);
                    double qp  = Field(i+dir[0],j+dir[1],k+dir[2]);
                    double qm  = Field(i-dir[0],j-dir[1],k-dir[2]);
                    double qpp = Field(i+2*dir[0],j+2*dir[1],k+2*dir[2]);
                    double qmm = Field(i-2*dir[0],j-2*dir[1],k-2*dir[2]);
                    double v   = Vel.Average(i,j,k)[direction];
                    double vp  = v;
                    double vm  = v;
                    if(compressible)
                    {
                        vp = Vel.Average(i+dir[0],j+dir[1],k+dir[2])[direction];
                        vm = Vel.Average(i-dir[0],j-dir[1],k-dir[2])[direction];
                    }
                    double L0  = 0.5*Limiter(q-qm,qp-q);
                    double Lp1 = 0.5*Limiter(qp-q,qpp-qp);
                    double Lm1 = 0.5*Limiter(qm-qmm,q-qm);
                    double ap  = 0.5*(vp+v);
                    double am  = 0.5*(vm+v);
                    double Fp  = std::max(ap,0.0)*(q+L0)
                               + std::min(ap,0.0)*(qp-Lp1);
                    double Fm  = std::max(am,0.0)*(qm+Lm1)
                               + std::min(am,0.0)*(q-L0);

                    if (std::fabs(ap)*locdt > 0.49*dx or std::fabs(am)*locdt > 0.49*dx)
                    {
                        CFL = 1;
                    }
                    FieldUpdate(i,j,k) = Field(i,j,k) - (Fp-Fm)*locdt/dx;
                }
                OMP_PARALLEL_STORAGE_LOOP_END
                if (CFL == 0)
                {
                    Field = FieldUpdate;
                    BC.SetX(Field);
                    BC.SetY(Field);
                    BC.SetZ(Field);
                }
                else
                {
                    locdt *= 0.5;
                    Field = FieldBackup;
                }
            }
            its *= 2;
        }
        while (CFL>0);
    }
};

class OP_EXPORTS AdvectionHR : OPObject
{
 public:

    AdvectionHR(){};
    AdvectionHR(Settings& locSettings, std::string InputFileName = DefaultInputFileName);

    void Initialize(Settings& Settings) override;
    void ReadInput(const std::string InputFileName) override;
    void ReadInput(std::stringstream& inp) override;
    void AdvectPhaseField(
            PhaseField& Phase,
            const Velocities& Vel,
            const BoundaryConditions& BC,
            const double dx,
            const double dt,
            const int tStep,
            const bool finalize = true);

    template<class T, int Rank>
    void AdvectField(
            Storage3D<T, Rank>& Field,
            const Velocities& Vel,
            const BoundaryConditions& BC,
            const double dx,
            const double dt,
            const int tStep,
            bool compressible = false)
    {
        if (Field.Bcells() < 2)
        {
            Info::WriteExit("Number of Bcells for storage needs to be 2 or higher.",
                            "AdvectionHR", "Advection(Field)");
            raise(SIGABRT);
        }

        BC.SetX(Field);
        BC.SetY(Field);
        BC.SetZ(Field);
        double (*Limiter)(const double,const double);
        switch(Scheme)
        {
            case AdvectionSchemes::Minmod:
            {
                Limiter = &Slope_Limiter_Minmod;
                break;
            }
            case AdvectionSchemes::MonotonizedCentral:
            {
                Limiter = &Slope_Limiter_MC;
                break;
            }
            case AdvectionSchemes::Superbee:
            {
                Limiter = &Slope_Limiter_Superbee;
                break;
            }
            /*case AdvectionSchemes::LaxWendroff:
            {
                Limiter = &Slope_Limiter_LaxWendroff;
                break;
            }*/
            case AdvectionSchemes::Upwind:
            default:
            {
                Limiter = &Slope_Limiter_Upwind;
                break;
            }
        }
        if (tStep % 2 == 0)
        {
            if(dNx) AdvectionMethod<T,Rank>::Advection(Field, BC, Vel, 0, dt, dx, Limiter, compressible);
            if(dNy) AdvectionMethod<T,Rank>::Advection(Field, BC, Vel, 1, dt, dx, Limiter, compressible);
            if(dNz) AdvectionMethod<T,Rank>::Advection(Field, BC, Vel, 2, dt, dx, Limiter, compressible);
        }
        else
        {
            if(dNz) AdvectionMethod<T,Rank>::Advection(Field, BC, Vel, 2, dt, dx, Limiter, compressible);
            if(dNy) AdvectionMethod<T,Rank>::Advection(Field, BC, Vel, 1, dt, dx, Limiter, compressible);
            if(dNx) AdvectionMethod<T,Rank>::Advection(Field, BC, Vel, 0, dt, dx, Limiter, compressible);
        }
    }

 private:
    AdvectionSchemes Scheme;
    int dNx;
    int dNy;
    int dNz;

    void Advection(
            PhaseField& Phase,
            const BoundaryConditions& BC,
            const Velocities& Vel,
            const int direction,
            const double dt,
            const double dx,
            double (*Limiter)(const double,const double));
};

}// namespace openphase
#endif
