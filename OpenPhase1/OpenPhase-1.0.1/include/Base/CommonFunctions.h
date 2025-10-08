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
 *   File created :   2017
 *   Main contributors :   Marvin Tegeler, Raphael Schiedung
 *
 */

#ifndef COMMONFUNCTIONS_H
#define COMMONFUNCTIONS_H

#include "Base/Storages.h"
#include "BoundaryConditions.h"
#include "Settings.h"
namespace openphase
{
class CommonFunctions
{
 public:
    static void CalculateDistancePeriodic(dVector3 A, dVector3 B,
            dVector3 &dist, double Nx, double Ny, double Nz)
    {
        dist[0] = A[0]-B[0];
        if (std::abs(dist[0]) > std::abs(dist[0]+Nx)) dist[0] += Nx;
        else if (std::abs(dist[0]) > std::abs(dist[0]-Nx)) dist[0] -= Nx;
        dist[1] = A[1]-B[1];
        if (std::abs(dist[1]) > std::abs(dist[1]+Ny)) dist[1] += Ny;
        else if (std::abs(dist[1]) > std::abs(dist[1]-Ny)) dist[1] -= Ny;
        dist[2] = A[2]-B[2];
        if (std::abs(dist[2]) > std::abs(dist[2]+Nz)) dist[2] += Nz;
        else if (std::abs(dist[2]) > std::abs(dist[2]-Nz)) dist[2] -= Nz;
    };

    template<typename vector>
    static vector Distance( const vector A, const vector B,
            const BoundaryConditions& BC, const Settings& locSettings)
    {
        vector dist = A-B;
#ifdef MPI_PARALLEL
        if (BC.MPIperiodicX or BC.BC0X == BoundaryConditionTypes::Periodic)
#else
        if (BC.BC0X == BoundaryConditionTypes::Periodic)
#endif
        {
            double val = dist[0] - locSettings.TotalNx;
            if (std::abs(val) < std::abs(dist[0])) dist[0] = val;
        }
#ifdef MPI_PARALLEL
        if (BC.MPIperiodicX or BC.BCNX == BoundaryConditionTypes::Periodic)
#else
        if (BC.BCNX == BoundaryConditionTypes::Periodic)
#endif
        {
            double val = dist[0] + locSettings.TotalNx;
            if (std::abs(val) < std::abs(dist[0])) dist[0] = val;
        }
#ifdef MPI_PARALLEL
        if (BC.MPIperiodicY or BC.BC0Y == BoundaryConditionTypes::Periodic)
#else
        if (BC.BC0Y == BoundaryConditionTypes::Periodic)
#endif
        {
            double val = dist[1] - locSettings.TotalNy;
            if (std::abs(val) < std::abs(dist[1])) dist[1] = val;
        }
#ifdef MPI_PARALLEL
        if (BC.MPIperiodicY or BC.BCNY == BoundaryConditionTypes::Periodic)
#else
        if (BC.BCNY == BoundaryConditionTypes::Periodic)
#endif
        {
            double val = dist[1] + locSettings.TotalNy;
            if (std::abs(val) < std::abs(dist[1])) dist[1] = val;
        }
#ifdef MPI_PARALLEL
        if (BC.MPIperiodicZ or BC.BC0Z == BoundaryConditionTypes::Periodic)
#else
        if (BC.BC0Z == BoundaryConditionTypes::Periodic)
#endif
        {
            double val = dist[2] - locSettings.TotalNz;
            if (std::abs(val) < std::abs(dist[2])) dist[2] = val;
        }
#ifdef MPI_PARALLEL
        if (BC.MPIperiodicZ or BC.BCNZ == BoundaryConditionTypes::Periodic)
#else
        if (BC.BCNZ == BoundaryConditionTypes::Periodic)
#endif
        {
            double val = dist[2] + locSettings.TotalNz;
            if (std::abs(val) < std::abs(dist[2])) dist[2] = val;
        }

        return dist;
    };
};

}// namespace openphase
#endif
