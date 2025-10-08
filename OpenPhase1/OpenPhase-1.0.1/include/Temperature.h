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
 *   Main contributors :   Oleg Shchyglo; Raphael Schiedung; Marvin Tegeler;
 *                         Matthias Stratmann
 *
 */

#ifndef TEMPERATURE_H
#define TEMPERATURE_H

#include "Base/Includes.h"
#include "H5Interface.h"
#include "Info.h"
#include "Temperature1Dextension.h"
#include "HeatDiffusion.h"

namespace openphase
{
class AdvectionHR;
class BoundaryConditions;
class Composition;
class FlowSolverLBM;
class GrainInfo;
class H5Interface;
class PhaseField;
class Settings;
class Velocities;

class Extension1D                                                               ///< 1D temperature field extension storage
{
 public:
    void Initialize(size_t size)
    {
        Data.resize(size+2,0.0);
        DataOld.resize(size+2,0.0);
    }
    Extension1D& operator=(Extension1D& RHS)
    {
        Data = RHS.Data;
        DataOld = RHS.DataOld;
        return *this;
    }
    void setBC()
    {
        Data[size()-1] = Data[size()-2];
    }
    void setBC(double value)
    {
        Data[size()-1] = value;
    }
    bool isActive()
    {
        return Data.size() > 0;
    };
    size_t size() const
    {
        return Data.size();
    }
    void store_temporary(void)
    {
        DataOld = Data;
    }
    void read(std::ifstream& out)
    {
        out.read(reinterpret_cast<char*>(Data.data()),Data.size()*sizeof(double));
    }
    void write(std::ofstream& out)
    {
        out.write(reinterpret_cast<char*>(Data.data()),Data.size()*sizeof(double));
    }
    std::vector<double> Data;                                                   ///< Data storage array
    std::vector<double> DataOld;                                                ///< Temporary data storage for iterative heat diffusion solver
    };

class OP_EXPORTS Temperature : public OPObject                                             ///< Storage for the temperature
{
 public:
    Temperature(){};
    Temperature(Settings& locSettings, const std::string InputFileName = DefaultInputFileName);
    void Initialize(Settings& locSettings) override;                            ///< Allocates internal storages, initializes internal parameters
    void ReadInput(const std::string InputFileName) override;                   ///< Reads initial values for internal parameters.
    void ReadInput(std::stringstream& inp) override;                            ///< Reads initial values for internal parameters.
    void Remesh(const int newNx, const int newNy, const int newNz,
                                        const BoundaryConditions& BC) override; ///< Changes system size while keeping the data
    void MoveFrame(const int dx, const int dy, const int dz,
                   const BoundaryConditions& BC);                               ///< Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in x, y or z direction correspondigly.
    void ConsumePlane(const int dx, const int dy, const int dz,
                   const int x, const int y, const int z,
                   const BoundaryConditions& BC);                               ///< Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in the direction to/from (x, y, z) point.

    void SetBoundaryConditions(const BoundaryConditions& BC);                   ///< Sets boundary conditions for the temperature by setting the appropriate values in ghost nodes
    void Write(const int tStep);                                                ///< Writes the raw temperature into a file
    void WriteH5(const int tStep, H5Interface& H5);
    bool Read(BoundaryConditions& BC, const int tStep);                         ///< Reads the raw temperature from a file
    bool ReadH5(const int tStep, H5Interface& H5);
    void WriteVTK(const int tStep, Settings& locSettings) const;                ///< Writes temperature in the VTK format into a file
    void WriteGradientVTK(const int tStep, Settings& locSettings) const;        ///< Writes temperature gradient in the VTK format into a file

    void PrintPointStatistics(const int x, const int y, const int z) const;     ///< Prints temperature at a given point (x, y, z) to screen
    void PrintStatistics() const;                                               ///< Prints min, max and average temperature to screen
    void Set(const BoundaryConditions& BC, const PhaseField& Phase, const double dt);///< Sets temperature according to the cooling rate
    void SetMinMaxAvg();                                                        ///< Evaluates min, max and average temperatures in the simulation domain

    void SetInitial(const BoundaryConditions& BC);                              ///< Sets initial temperature according to the starting temperature and temperature gradient

    double& operator()(const int x, const int y, const int z)                   ///< Bidirectional access operator
    {
        return Tx(x, y, z);
    };
    double const& operator()(const int x, const int y, const int z) const       ///< Constant reference access operator
    {
        return Tx(x, y, z);
    };
    double at(const double x, const double y, const double z) const             ///< Interpolating access operator
    {
        return Tx.at(x, y, z);
    };
    double Average(void);                                                       ///< Averages temperature in the entire domain
    double Average(PhaseField& Phi, size_t n, size_t m);                        ///< Averages temperature in all interfaces of a phase pair (n,m)
    Temperature& operator=(const Temperature& rhs);                             ///< Copy operator for Temperature class
    void Advect(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi, const BoundaryConditions& BC,
                         const double dt, const double tStep) override;         ///< Advects temperature

    double T0;                                                                  ///< Starting temperature
    double TBC0X;                                                               ///< x0 Boundary starting temperature
    double TBCNX;                                                               ///< xN Boundary starting temperature
    double TBC0Y;                                                               ///< y0 Boundary starting temperature
    double TBCNY;                                                               ///< yN Boundary starting temperature
    double TBC0Z;                                                               ///< z0 Boundary starting temperature
    double TBCNZ;                                                               ///< zN Boundary starting temperature

    double TSphere;                                                             ///< Starting temperature
    double Pr;                                                                  ///< Prandth Number

    dVector3 dT_dr;                                                             ///< Initial Temperature gradient
    dVector3 r0;                                                                ///< Initial position with T0 temperature gradient (dT_dr will be applied with respect to this position)
    double dT_dt;                                                               ///< Cooling rate

    int TotalNx;                                                                ///< X dimension of the system in MPI parallel mode
    int OffsetX;                                                                ///< X position of the current domain in MPI parallel mode
    int TotalNy;                                                                ///< Y dimension of the system in MPI parallel mode
    int OffsetY;                                                                ///< Y position of the current domain in MPI parallel mode
    int TotalNz;                                                                ///< Z dimension of the system in MPI parallel mode
    int OffsetZ;                                                                ///< Z position of the current domain in MPI parallel mode

    int Nx;                                                                     ///< Size of the inner calculation domain along X
    int Ny;                                                                     ///< Size of the inner calculation domain along Y
    int Nz;                                                                     ///< Size of the inner calculation domain along Z

    int dNx;                                                                    ///< Active X dimension
    int dNy;                                                                    ///< Active Y dimension
    int dNz;                                                                    ///< Active Z dimension

    size_t Nphases;                                                             ///< Number of thermodynamic phases
    double dx;                                                                  ///< Grid spacing

    Storage3D<double, 0> Tx;                                                    ///< Array of temperature
    Storage3D<double, 0> TxDot;                                                 ///< Array of temperature increments needed for advection
    Storage3D<double, 0> TxOld;                                                 ///< Temporary temperature storage

    double Tmin;                                                                ///< Minimum temperature in the simulation domain
    double Tmax;                                                                ///< Maximum temperature in the simulation domain
    double Tavg;                                                                ///< Average temperature in the simulation domain
    bool   Tneg;                                                                ///< Negative temperature control: "true" if negative temperature is allowed

    std::string VTKDir;                                                         ///< Directory-path added in front of VTK files
    std::string RawDataDir;                                                     ///< Directory-path added in front of Restart files

    Temperature1Dextension ExtensionX0;                                         ///< 1D temperature field extension at the lower X boundary
    Temperature1Dextension ExtensionXN;                                         ///< 1D temperature field extension at the upper X boundary
    Temperature1Dextension ExtensionY0;                                         ///< 1D temperature field extension at the lower Y boundary
    Temperature1Dextension ExtensionYN;                                         ///< 1D temperature field extension at the upper Y boundary
    Temperature1Dextension ExtensionZ0;                                         ///< 1D temperature field extension at the lower Z boundary
    Temperature1Dextension ExtensionZN;                                         ///< 1D temperature field extension at the upper Z boundary
    bool ExtensionsActive;                                                      ///< True if at least one 1D extension is active, false otherwise

 protected:
 private:
    void SetToValue(double value);                                              ///< Sets the temperature in the system to a given value
    void SetToValueWithGrad(double value, dVector3 dT_dr, dVector3 r0);         ///< Sets the temperature in the system to a given value considering temperature gradient
    void IncrementWithValue(double value);                                      ///< Increments the temperature in the system with a given value
};

} // namespace openphase
#endif
