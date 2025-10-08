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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitri Medvedev
 *
 */

#ifndef VELOCITIES_H
#define VELOCITIES_H

#include "Base/Includes.h"
namespace openphase
{

class Settings;
class PhaseField;
class BoundaryConditions;
class ElasticProperties;
class RunTimeControl;

class OP_EXPORTS Velocities : public OPObject                                              ///< Storage for the advection velocities
{
 public:
    Velocities(){};
    Velocities(Settings& locSettings);
    void Initialize(Settings& locSettings) override;                            ///<  Allocates memory, initializes the settings
    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC) override; ///<  Remeshes the system while keeping the data

    double GetMaxVelocity();                                                    ///<  returns maximal velocity component
    void Advect(BoundaryConditions& BC, double dt, AdvectionSchemes scheme = AdvectionSchemes::Upwind);        ///<  Advects average velocities
    void CalculateAverage(const PhaseField& Phase);                             ///<  Calculates the average velocity in each grid point using phase velocities
    void Clear(void);                                                           ///<  Clears the velocity storage
    void MoveFrame(int dx, int dy, int dz, BoundaryConditions& BC);             ///<  Moves frame by (dx,dy,dz).
    void PrescribePhaseVelocities(PhaseField& Phase);                           ///<  Sets phase velocities equal to average velocities
    void Read(BoundaryConditions& BC, const int tStep);
    void SetAllPhases(dVector3& value);                                         ///<  Set phase velocities to a fixed value
    void SetAverage(dVector3& value);                                           ///<  Set average velocity to a fixed value
    void SetAverage(ElasticProperties& EP, double dt);                          ///<  Set average velocity from EP.Displacements
    void SetBoundaryConditions(const BoundaryConditions& BC);                   ///<  Sets boundary conditions for the velocities
    void Write(int tStep) const;
    void WriteVTK(int tStep, Settings& locSettings) const;                      ///<  Writes velocities in VTK format into a file

    void WriteStatistics(RunTimeControl& RTC);                                  ///<  Writes velocity statistics into a file, input: time step

    Storage3D < dVector3, 0 > Average;                                          ///<  Average velocities storage
    Storage3D < dVector3, 0 > AverageDot;                                       ///<  Average velocities increments storage
    Storage3D < dVector3, 1 > Phase;                                            ///<  Phase velocities storage
    double dx;                                                                  ///<  Grid spacing
    size_t Nphases;                                                             ///<  Number of thermodynamic phases
    int Nx;                                                                     ///<  Size of the inner calculation domain along X
    int Ny;                                                                     ///<  Size of the inner calculation domain along Y
    int Nz;                                                                     ///<  Size of the inner calculation domain along Z
    int dNx;                                                                    ///<  Active X dimension
    int dNy;                                                                    ///<  Active Y dimension
    int dNz;                                                                    ///<  Active Z dimension

    void PrintPointStatistics(int x, int y, int z);                             ///<  Prints local velocity values to screen.

    std::string VTKDir;                                                         ///< Directory-path added in front of VTK files
    std::string RawDataDir;                                                     ///< Directory-path added in front of Restart files

    Velocities& operator= (const Velocities& rhs);                              ///< Copy operator for Velocities class

 protected:
 private:
};

}
#endif
