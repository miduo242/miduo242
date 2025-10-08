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
 *   File created :   2021
 *   Main contributors :   Oleg Shchyglo; Raphael Schiedung; Marvin Tegeler;
 *                         Matthias Stratmann
 *
 */

#include "Temperature1Dextension.h"
#include "Temperature.h"
#include "HeatDiffusion.h"

namespace openphase
{
using namespace std;

pair<int, int> three_state_bounds_selector(int selector, int lower_bound, int upper_bound)
{
    pair<int, int> result;
    switch (selector)
    {
        case -1: // lower bound only
        {
            result = make_pair(lower_bound, lower_bound + 1);
            break;
        }
        case 0: // both bounds
        {
            result = make_pair(lower_bound, upper_bound);
            break;
        }
        case 1: // upper bound only
        {
            result = make_pair(upper_bound - 1, upper_bound);
            break;
        }
    }
    return result;
}

void Temperature1Dextension::SetInitial(const Temperature& Tx)
{
    //Get the average temperature at the boundary of interest
    Xbounds = three_state_bounds_selector(Direction[0], 0, Tx.Nx);
    Ybounds = three_state_bounds_selector(Direction[1], 0, Tx.Ny);
    Zbounds = three_state_bounds_selector(Direction[2], 0, Tx.Nz);

    double boundaryValue = 0.0;

    double area          = 0.0;

    for (int i = Xbounds.first; i < Xbounds.second; ++i)
    for (int j = Ybounds.first; j < Ybounds.second; ++j)
    for (int k = Zbounds.first; k < Zbounds.second; ++k)
    {
        boundaryValue += Tx(i,j,k);
        area++;
    }

    Data[0] = boundaryValue/area;
    //Done getting the average temperature at the boundary of interest

    int offset_x = 0.5*Tx.Nx*Direction[1]*Direction[2];
    int offset_y = 0.5*Tx.Ny*Direction[0]*Direction[2];
    int offset_z = 0.5*Tx.Nz*Direction[0]*Direction[1];

    double dx = Tx.dx;
    bool limit = false;

    for (size_t i = 0; i < Data.size(); ++i)
    {
        Data[i] = Data[0] + (Tx.dT_dr[0]*(offset_x + i*Direction[0] - Tx.r0[0]) +
                           Tx.dT_dr[1]*(offset_y + i*Direction[1] - Tx.r0[1]) +
                           Tx.dT_dr[2]*(offset_z + i*Direction[2] - Tx.r0[2]))*dx;
        if(Data[i] < 0.0)
        {
            limit = true;
        }
    }
    if(limit and !Tx.Tneg)
    {
        stringstream message;
        message << "Negative temperature detected!\n"
                << "If it is an intended behavior use <$Tneg : Yes>\n"
                << "in the temperature input to allow negative temperatures\n"
                << "or adjust initial temperature input parameters.\n";
        Info::WriteExit(message.str(), Tx.thisclassname, "1Dextension::SetInitial()");
        exit(1);
    }
}

void Temperature1Dextension::PerformImplicitIteration(Temperature& Tx, HeatDiffusion& HD, double& residual, double dt)
{
    /** Calculation one iteration of the heat diffusion using a Jacobi implicit
        algorithm using the simplified equation for different 1D-Extensions:

        double dT = (RhoCp*TxOld[x]+Lambda*dt/(dx*dx)*(Tx[x-1]+Tx[x+1]))
                   /(RhoCp + 2.0*Lambda*dt/(dx*dx)) - Tx[x];
    */
    Xbounds = three_state_bounds_selector(Direction[0], 0, Tx.Nx);
    Ybounds = three_state_bounds_selector(Direction[1], 0, Tx.Ny);
    Zbounds = three_state_bounds_selector(Direction[2], 0, Tx.Nz);

    double dt_dx2 = dt/(HD.dx*HD.dx);

    double boundaryValue = 0.0;
    double averageRhoCP  = 0.0;
    double averageLambda = 0.0;
    double area          = 0.0;

    for (int i = Xbounds.first; i < Xbounds.second; ++i)
    for (int j = Ybounds.first; j < Ybounds.second; ++j)
    for (int k = Zbounds.first; k < Zbounds.second; ++k)
    {
        averageRhoCP  += HD.EffectiveHeatCapacity(i,j,k);                       // Average volumetric heat capacity [J/(m^3 K)]
        averageLambda += HD.EffectiveThermalConductivity(i,j,k);                // Average thermal conductivity [J/(m s K)]
        boundaryValue += Tx(i,j,k);
        area++;
    }

#ifdef MPI_PARALLEL
    if(Direction[0] == 0)
    {
        double loc_averageRhoCP = averageRhoCP;
        MPI_Allreduce(&loc_averageRhoCP, &averageRhoCP, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double loc_averageLambda = averageLambda;
        MPI_Allreduce(&loc_averageLambda, &averageLambda, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double loc_boundaryValue = boundaryValue;
        MPI_Allreduce(&loc_boundaryValue, &boundaryValue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double loc_area = area;
        MPI_Allreduce(&loc_area, &area, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
#endif

    averageRhoCP  /= area;
    averageLambda /= area;

    Data[0] = boundaryValue/area;

    for (size_t i = 1; i < size()-1; ++i)
    {
        double locQ = 0.0;
        if(i == size()-2) locQ = Qdot *dt;

        double loc_deltaT = (averageRhoCP*DataOld[i] +
                             averageLambda*dt_dx2*(Data[i+1] + Data[i-1]) + locQ)
                           /(averageRhoCP + 2.0*averageLambda*dt_dx2) - Data[i];

        residual = max(residual,loc_deltaT*loc_deltaT);
        Data[i] += loc_deltaT;
    }

    for (int i = Xbounds.first; i < Xbounds.second; ++i)
    for (int j = Ybounds.first; j < Ybounds.second; ++j)
    for (int k = Zbounds.first; k < Zbounds.second; ++k)
    {
        Tx(i+Direction[0],j+Direction[1],k+Direction[2]) = Data[1];
    }
}

}// namespace openphase
