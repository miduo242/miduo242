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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev;
 *                         Reza Darvishi Kamachali; Philipp Engels;
 *                         Raphael Schiedung
 *
 */

#include "Base/CommonFunctions.h"
#include "Base/Includes.h"
#include "BoundaryConditions.h"
#include "DoubleObstacle.h"
#include "DrivingForce.h"
#include "GrainInfo.h"
#include "Info.h"
#include "Initializations.h"
#include "InterfaceProperties.h"
#include "Mechanics/ElasticProperties.h"
#include "Orientations.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Tools.h"
#include "UserDrivingForce.h"
#include "CoordinationShells.h"

/***************************************************************/

namespace openphase
{

double rnd()
{
    return double(rand())/RAND_MAX;
}

using namespace std;

vector<iVector3> Initializations::QuasiRandomNuclei(PhaseField& Phase, Settings& locSettings, size_t phaseIndex, int dist, int seed, int offset)
{
    vector<iVector3> result;

    if (seed == -1)
    {
        srand(time(NULL));
    }
    else
    {
        srand(seed);
    }

#ifdef MPI_PARALLEL
    int Nx = locSettings.TotalNx;
    int Ny = locSettings.TotalNy;
    int Nz = locSettings.TotalNz;
#else
    int Nx = locSettings.Nx;
    int Ny = locSettings.Ny;
    int Nz = locSettings.Nz;
#endif

    double threshold = 0.05;
    int distx = (dist < Nx) ? dist : 0;
    int disty = (dist < Ny) ? dist : 0;
    int distz = (dist < Nz) ? dist : 0;

    for (int i = distx; i < Nx; i += 2 * dist)
    for (int j = disty; j < Ny; j += 2 * dist)
    for (int k = distz; k < Nz; k += 2 * dist)
    {
        double chance = double(rand()) / double(RAND_MAX);
#ifdef MPI_PARALLEL
MPI_Bcast(&(chance), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
        if (chance > threshold)
        {
            int di = (dist >= Nx or offset == 0) ? 0 : (rand() % (2 * offset) - offset);
            int dj = (dist >= Ny or offset == 0) ? 0 : (rand() % (2 * offset) - offset);
            int dk = (dist >= Nz or offset == 0) ? 0 : (rand() % (2 * offset) - offset);
#ifdef MPI_PARALLEL
MPI_Bcast(&(di), 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&(dj), 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&(dk), 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
            if (i + di >= 0 and i + di < Nx and
                j + dj >= 0 and j + dj < Ny and
                k + dk >= 0 and k + dk < Nz)
            {
                Phase.PlantGrainNucleus(phaseIndex, i + di, j + dj, k + dk);

                iVector3 temp;
                temp[0] = i + di;
                temp[1] = j + dj;
                temp[2] = k + dk;
                result.push_back(temp);
            }
        }
    }

    return result;
}

vector<iVector3> Initializations::QuasiRandomSpheres(PhaseField& Phase,
                        BoundaryConditions& BC, Settings& locSettings,
                        int phaseIndex1, int phaseIndex2, int dist,
                        double radius1, double radius2,
                        double probabilityPhase1, int seed, int offset)
{
    vector<iVector3> result;

    if (seed == -1)
    {
        srand(time(NULL));
    }
    else
    {
        srand(seed);
    }

    int distx = (dist < locSettings.Nx) ? dist : 0;
    int disty = (dist < locSettings.Ny) ? dist : 0;
    int distz = (dist < locSettings.Nz) ? dist : 0;
    double threshold = 0.05;

    for (int i = distx; i < locSettings.Nx; i += 2 * dist)
    for (int j = disty; j < locSettings.Ny; j += 2 * dist)
    for (int k = distz; k < locSettings.Nz; k += 2 * dist)
    {
        double chance = double(rand()) / double(RAND_MAX);

        if (chance > threshold)
        {
            int di = (distx == 0 or offset == 0) ? 0 : (rand() % (2 * offset) - offset);
            int dj = (disty == 0 or offset == 0) ? 0 : (rand() % (2 * offset) - offset);
            int dk = (distz == 0 or offset == 0) ? 0 : (rand() % (2 * offset) - offset);

            if (i + di >= 0 and i + di < locSettings.Nx and
                j + dj >= 0 and j + dj < locSettings.Ny and
                k + dk >= 0 and k + dk < locSettings.Nz)
            {
                double randPhase = double(rand()) / double(RAND_MAX);
                if (randPhase < probabilityPhase1)
                {
                    Initializations::Sphere(Phase, phaseIndex1, radius1, i+di, j+dj, k+dk, BC, locSettings);
                }
                else
                {
                    Initializations::Sphere(Phase, phaseIndex2, radius2, i+di, j+dj, k+dk, BC, locSettings);
                }

                iVector3 temp;
                temp[0] = i + di;
                temp[1] = j + dj;
                temp[2] = k + dk;

                result.push_back(temp);
            }
        }
    }

    return result;
}

void Initializations::RandomNucleiOnPlane(PhaseField& Phase,
                      Settings& locSettings, size_t phaseIndex, size_t Nparticles,
                      int seed, string axis, string position)
{
    if (seed == -1){ srand(time(NULL)); }
    else{ srand(seed); }

    transform(axis.begin(), axis.end(),axis.begin(), ::toupper);
    transform(position.begin(), position.end(),position.begin(), ::toupper);

    enum AlongAxis { X, Y, Z };
    enum OnPosition { Top, Bottom };

    AlongAxis along;
    OnPosition at;

    if(axis == "X"){ along = X;}
    else if(axis == "y"){ along = Y;}
    else{ along = Z;}

    if(position == "TOP"){ at = Top;}
    else{ at = Bottom;}

    for (size_t n = 0; n < Nparticles; n++)
    {
        bool planted = false;
        int iterations = 0;
        while(!planted and iterations < 1000)
        {
            iterations++;
            int di = rand() % locSettings.Nx;
            int dj = rand() % locSettings.Ny;
            int dk = rand() % locSettings.Nz*0;

            switch(along)
            {
                case X :
                {
                    switch(at){ case Top    : di = 0;   break;
                                case Bottom : di = locSettings.Nx; break;}
                } break;

                case Y :
                {
                    switch(at){ case Top    : dj = 0;   break;
                                case Bottom : dj = locSettings.Ny; break; }
                } break;

                case Z :
                {
                    switch(at){ case Top    : dj = 0;   break;
                                case Bottom : dj = locSettings.Nz; break;
                    }
                } break;
            }

            bool freecell = true;

            for (int ii = - locSettings.iWidth; ii <= locSettings.iWidth; ii++)
            for (int jj = - locSettings.iWidth; jj <= locSettings.iWidth; jj++)
            for (int kk = - locSettings.iWidth; kk <= locSettings.iWidth; kk++)
            if ((di+ii > -Phase.Fields.Bcells()
             and di+ii < locSettings.Nx + Phase.Fields.Bcells() and
                 dj+jj > -Phase.Fields.Bcells()
             and dj+jj < locSettings.Ny + Phase.Fields.Bcells() and
                 dk+kk > -Phase.Fields.Bcells()
             and dk+kk < locSettings.Nz + Phase.Fields.Bcells()) and
                Phase.Fields(di+ii,dj+jj,dk+kk).flag)
            {
                freecell = false;
                break;
            }
            if(freecell)
            {
                //int locIndex = Phase.PlantGrainNucleus(phaseIndex, di, dj, dk);
                Phase.PlantGrainNucleus(phaseIndex, di, dj, dk);
                planted = true;
            }
        }
    }
}

void Initializations::RandomNuclei(PhaseField& Phase, Settings& locSettings, size_t phaseIndex, size_t Nparticles, int seed)
{
    if (seed == -1)
    {
        srand(time(NULL));
    }
    else
    {
        srand(seed);
    }
    for (size_t n = 0; n < Nparticles; n++)
    {
        bool planted = false;
        int iterations = 0;
        while(!planted and iterations < 1000)
        {
            iterations++;
            int di = rand() % locSettings.Nx;
            int dj = rand() % locSettings.Ny;
            int dk = rand() % locSettings.Nz*0;

            bool freecell = true;

            for (int ii = - locSettings.iWidth; ii <= locSettings.iWidth; ii++)
            for (int jj = - locSettings.iWidth; jj <= locSettings.iWidth; jj++)
            for (int kk = - locSettings.iWidth; kk <= locSettings.iWidth; kk++)
            if ((di+ii > -Phase.Fields.Bcells() and di+ii < locSettings.Nx + Phase.Fields.Bcells() and
                 dj+jj > -Phase.Fields.Bcells() and dj+jj < locSettings.Ny + Phase.Fields.Bcells() and
                 dk+kk > -Phase.Fields.Bcells() and dk+kk < locSettings.Nz + Phase.Fields.Bcells()) and
                Phase.Fields(di+ii,dj+jj,dk+kk).flag)
            {
                freecell = false;
                break;
            }
            if(freecell)
            {
                int locIndex = Phase.PlantGrainNucleus(phaseIndex, di, dj, dk);
                int Q1 = rand() % 180;
                int Q2 = rand() % 180;
                int Q3 = rand() % 180;
                EulerAngles locAngles({Q1*Pi/180, Q2*Pi/180, Q3*Pi/180}, XYZ);

                Phase.FieldsStatistics[locIndex].Orientation = locAngles.getQuaternion();
                planted = true;
            }
        }
    }
}
/**
 * Initializes a new phase on the negative side of a plane
 * @param Phase       Phase field
 * @param Point       Point on plane
 * @param Orientation Orientation of the plane
 * @param PhaseIndex  Phase index
 * @param BC          Boundary conditions
 * @param locSettings Project settings
 */
size_t Initializations::SectionalPlane(PhaseField& Phase, const size_t PhaseIndex,
        const dVector3 Point, const dVector3 Orientation,
        const BoundaryConditions& BC, const Settings& locSettings, const bool NewGrain, const size_t FieldIndex)
{

    const dVector3 n        = Orientation.normalized();
    const double   iWidth   = locSettings.iWidth;

    int locIndex;
    if (NewGrain)
    {
         locIndex = Phase.AddGrainInfo(PhaseIndex);
    }
    else
    {
         locIndex = FieldIndex;
    }

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        const dVector3 x = {double(i),double(j),double(k)}; // position in space
        const double distance = (x-Point)*n;   // distance from plane

        if (distance < - 0.5*iWidth)
        {
            Phase.Fields(i, j, k).clear();
            Phase.Fields(i, j, k).set_value(locIndex, 1.0);
        }
        else if (distance <= 0.5*iWidth)
        {
            const double IntProf = 0.5 - 0.5*sin(Pi*(distance)/iWidth);
            for(auto alpha = Phase.Fields(i,j,k).begin();
                     alpha < Phase.Fields(i,j,k).end(); alpha++)
            {
                alpha->value *= 1.0 - IntProf;
            }
            Phase.Fields(i, j, k).set_value(locIndex, IntProf);
        }
    }
    STORAGE_LOOP_END

    Phase.Finalize(BC);
    return locIndex;
}

size_t Initializations::Sphere(PhaseField& Phase, const size_t PhaseIndex,
        const double Radius, const double x0, const double y0, const double z0,
        const BoundaryConditions& BC, const Settings& locSettings,
        const bool Finalize)
{
    //TODO assertion should by Info::WriteExit because the should be warned about it
    assert(x0 >= 0);
    assert(y0 >= 0);
    assert(z0 >= 0);
    assert(x0 < locSettings.TotalNx);
    assert(y0 < locSettings.TotalNy);
    assert(z0 < locSettings.TotalNz);
    assert(Radius < std::max(locSettings.TotalNx,std::max(locSettings.TotalNy,locSettings.TotalNz)));

    const double iWidth = (locSettings.Resolution == Resolutions::Double) ? 0.5*locSettings.iWidth : locSettings.iWidth;
    const size_t locIndex = Phase.AddGrainInfo(PhaseIndex);
    auto set_sphere = [&iWidth,&Radius,&locIndex,&Phase](long int i,long int j,long int k, double rad)
    {
        if (rad < Radius - iWidth*0.5)
        {
            Phase.Fields(i,j,k).clear();
            Phase.Fields(i,j,k).set_value(locIndex, 1.0);
        }
        else if (rad <= Radius + iWidth*0.5)
        {
            const double IntProf = 0.5 - 0.5*std::sin(Pi*(rad - Radius)/iWidth);
            for(auto alpha = Phase.Fields(i,j,k).begin();
                     alpha < Phase.Fields(i,j,k).end(); alpha++)
            {
                alpha->value *= 1.0 - IntProf;
            }
            Phase.Fields(i,j,k).set_value(locIndex, IntProf);
        }
        return false;
    };
    const double loop_radius = Radius + 1.1*iWidth;
    loop_sphere(Phase.Fields, set_sphere, x0, y0, z0, loop_radius, locSettings, BC);

    if (Finalize)
    {
        Phase.FinalizeSR(BC);
        if(locSettings.Resolution == Resolutions::Double)
        {
            Phase.FinalizeDR(BC);
        }
    }
    return locIndex;
}

void Initializations::Young4Periodic(PhaseField& Phase, size_t PhaseIndex,
        BoundaryConditions& BC, Settings& locSettings)
{
    int Nx = locSettings.Nx;
    int Ny = locSettings.Ny;
    int Nz = locSettings.Nz;

    size_t locIndex0 = Phase.AddGrainInfo(PhaseIndex);
    size_t locIndex1 = Phase.AddGrainInfo(PhaseIndex+1);
    size_t locIndex2 = Phase.AddGrainInfo(PhaseIndex+2);
    size_t locIndex3 = Phase.AddGrainInfo(PhaseIndex+3);
    size_t locIndex4 = Phase.AddGrainInfo(PhaseIndex+4);
    size_t locIndex5 = Phase.AddGrainInfo(PhaseIndex+5);
    size_t locIndex6 = Phase.AddGrainInfo(PhaseIndex+6);
    size_t locIndex7 = Phase.AddGrainInfo(PhaseIndex+7);

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i,j,k).flag = 2;

        if (j <= Ny/2)
        {
            if      (i <= Nx/2 and k >  Nz/4 and k <= Nz*3./4)
                Phase.Fields(i, j, k).set_value(locIndex1,1);
            else if (i >  Nx/2 and k >  Nz/4 and k <= Nz*3./4)
                Phase.Fields(i, j, k).set_value(locIndex2,1);
            else if (i >  Nx/4 and i <= Nx*3./4 and k <= Nz/4)
                Phase.Fields(i, j, k).set_value(locIndex3,1);
            else if (i >  Nx/4 and i <= Nx*3./4 and k >  Nz*3./4)
                Phase.Fields(i, j, k).set_value(locIndex3,1);
            else Phase.Fields(i,j,k).set_value(locIndex0,1);
        }
        else
        {
            if      (i <= Nx/8 and k >  Nz/8 and k <= Nz*5./8)
                Phase.Fields(i, j, k).set_value(locIndex5,1);
            else if (i >  Nx/8 and i <= Nx*5./8 and k >  Nz/8 and k <= Nz*5./8)
                Phase.Fields(i, j, k).set_value(locIndex4,1);
            else if (i >  Nx*5./8 and k >  Nz/8 and k <= Nz*5./8)
                Phase.Fields(i, j, k).set_value(locIndex5,1);
            else if (i >  Nx*3./8 and i <= Nx*7./8 and k <= Nz/8)
                Phase.Fields(i, j, k).set_value(locIndex6,1);
            else if (i >  Nx*3./8 and i <= Nx*7./8 and k >  Nz*5./8)
                Phase.Fields(i, j, k).set_value(locIndex6,1);
            else Phase.Fields(i,j,k).set_value(locIndex7,1);
        }
    }
    STORAGE_LOOP_END

    Phase.SetBoundaryConditions(BC);
    Phase.CalculateDerivatives();
    Phase.CalculateFractions();
}

void Initializations::TripleJunction(PhaseField& Phase, size_t PhaseIndex, BoundaryConditions& BC, Settings& locSettings)
{
    int Nx = locSettings.Nx;
    int Nz = locSettings.Nz;

    Rectangular(Phase, PhaseIndex  , 2.0*Nx/3.0, 0, Nz/2.0, 1.0*Nx/3.0, 0,     Nz/4.0,  BC, locSettings, false);
    Rectangular(Phase, PhaseIndex+1, 2.0*Nx/3.0, 0, Nz/2.0, 1.0*Nx/3.0, 0, 3.0*Nz/4.0,  BC, locSettings, false);
    Rectangular(Phase, PhaseIndex+2,     Nx/3.0, 0, Nz    , 5.0*Nx/6.0, 0,     Nz/2.0,  BC, locSettings, false);

    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }
}

vector<size_t> Initializations::Young4(PhaseField& Phase, size_t PhaseIndex,
        BoundaryConditions& BC, Settings& locSettings)
{
    vector<size_t> result;

    int Nx = locSettings.Nx;
    int Ny = locSettings.Ny;
    int Nz = locSettings.Nz;

    result.push_back(Single    (Phase,PhaseIndex, BC,locSettings));
    result.push_back(Rectangular(Phase,PhaseIndex, 3.0*Nx/4.0+Phase.iWidth,    Ny/3.0+Phase.iWidth,Nz+2*Phase.iWidth,5.0*Nx/8.0,    Ny/6.0-Phase.iWidth,Nz/2, BC,locSettings, false));
    result.push_back(Rectangular(Phase,PhaseIndex, 3.0*Nx/4.0+Phase.iWidth,2.0*Ny/3.0+Phase.iWidth,Nz/2+Phase.iWidth,5.0*Nx/8.0,4.0*Ny/6.0,             Nz/4, BC,locSettings, false));
    result.push_back(Rectangular(Phase,PhaseIndex, 3.0*Nx/4.0+Phase.iWidth,2.0*Ny/3.0+Phase.iWidth,Nz/2+Phase.iWidth,5.0*Nx/8.0,4.0*Ny/6.0,         3.0*Nz/4, BC,locSettings, false));

    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }
    return result;
}

vector<size_t> Initializations::Young3(PhaseField& Phase, size_t alpha, size_t beta,
                             size_t gamma, size_t delta, BoundaryConditions& BC,
                             Settings& locSettings)
{
    vector<size_t> result;

    int Nx = locSettings.Nx;
    int Nz = locSettings.Nz;

    double x = Nx/2.0;
    double y = 0;
    double z = Nz/2.0;

    result.push_back(Single    (Phase,alpha,                   BC,locSettings));
    result.push_back(Rectangular(Phase,beta,  x,y,z, 0.5*x,0,z, BC,locSettings,false));
    result.push_back(Rectangular(Phase,gamma, x,y,z,     x,0,0, BC,locSettings,false));
    result.push_back(Rectangular(Phase,delta, x,y,z, 1.5*x,0,z, BC,locSettings,false));


    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }

    return result;
}

size_t Initializations::Single(PhaseField& Phase, size_t PhaseIndex,
        BoundaryConditions& BC, Settings& locSettings)
{
    size_t locIndex = Phase.AddGrainInfo(PhaseIndex);
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();
        Phase.Fields(i, j, k).set_value(locIndex, 1.0);
    }
    STORAGE_LOOP_END

    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }

    return locIndex;
}

size_t Initializations::Zlayer(PhaseField& Phase, size_t PhaseIndex, int Position,
                            int LayerThickness, BoundaryConditions& BC, Settings& locSettings)
{
    double iWidth = locSettings.iWidth;

    size_t index = Phase.AddGrainInfo(PhaseIndex);

    int Zmin = Position - LayerThickness/2;
    int Zmax = Position + LayerThickness/2;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        if (k > Zmin + iWidth/2 and k < Zmax - iWidth/2)
        {
            Phase.Fields(i, j, k).clear();
            Phase.Fields(i, j, k).set_value(index, 1);
        }
        if (k >= Zmin - iWidth/2 and k <= Zmin + iWidth/2)
        {
            double IntProfile = 0.5 - 0.5*sin(Pi*(k - Zmin)/iWidth);
            for(auto beta = Phase.Fields(i,j,k).begin();
                     beta != Phase.Fields(i,j,k).end(); beta++)
            {
                beta->value *= (IntProfile);
            }
            Phase.Fields(i, j, k).set_value(index, 1.0 - IntProfile);
            //Phase.Fields(i,j,k).flag = 2;
        }
        if (k >= Zmax - iWidth/2 and k <= Zmax + iWidth/2)
        {
            double IntProfile = 0.5 - 0.5*sin(Pi*(k - Zmax)/iWidth);
            for(auto beta = Phase.Fields(i,j,k).begin();
                     beta != Phase.Fields(i,j,k).end(); beta++)
            {
                beta->value *= 1.0 - IntProfile;
            }
            Phase.Fields(i, j, k).set_value(index, IntProfile);
            //Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }

    return index;
}

vector<size_t> Initializations::Fractional(PhaseField& Phase,
        size_t MajorityPhaseIndex, size_t MinorityPhaseIndex,
                                 double MinorityPhaseLayerThickness,
                                 BoundaryConditions& BC, Settings& locSettings)
{

    if( BC.BC0Z != BoundaryConditionTypes::NoFlux || BC.BCNZ != BoundaryConditionTypes::NoFlux )
    {
        Info::WriteWarning(
                "Consider NoFlux Boundary condition at least along \n"
                "Z-direction for Fractional initialization", thisclassname, "Fractional");
    }

    double iWidth = locSettings.iWidth;
    if(locSettings.Resolution == Resolutions::Double)
    {
        iWidth *= 0.5;
    }

    size_t index1 = Phase.AddGrainInfo(MajorityPhaseIndex);
    size_t index2 = Phase.AddGrainInfo(MinorityPhaseIndex);

    int offset = MinorityPhaseLayerThickness;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();

        if (k > offset + iWidth/2)
        {
            Phase.Fields(i, j, k).set_value(index1, 1);
        }
        else if (k < offset - iWidth/2)
        {
            Phase.Fields(i, j, k).set_value(index2, 1);
        }
        else
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(k - offset)/iWidth);
            Phase.Fields(i, j, k).set_value(index1, 1.0 - IntProf);
            Phase.Fields(i, j, k).set_value(index2, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }
    return vector<size_t>({index1,index2});
}

vector<size_t> Initializations::ThreeFractionals(PhaseField& Phase,
                                 size_t MajorityPhaseIndex,
                                 double MajorityPhaseLayerThickness,
                                 size_t MinorityPhaseIndex1,
                                 double MinorityPhaseLayerThickness1,
                                 size_t MinorityPhaseIndex2,
                                 BoundaryConditions& BC, Settings& locSettings)
{

    if( BC.BC0Z != BoundaryConditionTypes::NoFlux || BC.BCNZ != BoundaryConditionTypes::NoFlux )
    {
        Info::WriteWarning(
                "Consider NoFlux Boundary condition atleast along \n"
                "Z-direction for Fractional initialization", thisclassname, "Fractional");
    }

    double iWidth = locSettings.iWidth;

    size_t index1 = Phase.AddGrainInfo(MajorityPhaseIndex);
    size_t index2 = Phase.AddGrainInfo(MinorityPhaseIndex1);
    size_t index3 = Phase.AddGrainInfo(MinorityPhaseIndex2);

    int offset1 = MajorityPhaseLayerThickness;
    int offset2 = MinorityPhaseLayerThickness1;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();
        Phase.Fields(i,j,k).flag = 0;

        if (k < offset1 - iWidth/2)
        {
            Phase.Fields(i, j, k).set_value(index1, 1);
        }
        else if ((k > offset1 + iWidth/2) and (k < offset1+offset2 - iWidth/2))
        {
            Phase.Fields(i, j, k).set_value(index2, 1);
        }
        else if (k > offset1+offset2 + iWidth/2)
        {
            Phase.Fields(i, j, k).set_value(index3, 1);
        }
        else if ((k <= offset1 + iWidth/2) and (k >= offset1 - iWidth/2))
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(k - offset1)/iWidth);
            Phase.Fields(i, j, k).set_value(index2, 1.0 - IntProf);
            Phase.Fields(i, j, k).set_value(index1, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
        else if ((k <= offset1+offset2 + iWidth/2) and
                 (k >= offset1+offset2 - iWidth/2))
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(k - (offset1+offset2))/iWidth);
            Phase.Fields(i, j, k).set_value(index3, 1.0 - IntProf);
            Phase.Fields(i, j, k).set_value(index2, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }
    return vector<size_t>({index1,index2,index3});
}

vector<size_t> Initializations::TwoWalls(PhaseField& Phase, size_t ChannelPhaseIndex, size_t WallsPhaseIndex,
                                 double WallsThickness, BoundaryConditions& BC, Settings& locSettings)
{
    int Nz = locSettings.TotalNz;
    double iWidth = locSettings.iWidth;

    size_t index1 = Phase.AddGrainInfo(ChannelPhaseIndex);
    size_t index2 = Phase.AddGrainInfo(WallsPhaseIndex);

    int offset = WallsThickness;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();
        Phase.Fields(i,j,k).flag = 0;

        if (k + locSettings.OffsetZ > offset + iWidth/2 && k + locSettings.OffsetZ < Nz - offset - 1 - iWidth/2)
        {
            Phase.Fields(i, j, k).set_value(index1, 1);
        }
        else if (k+ locSettings.OffsetZ < offset - iWidth/2 || k+ locSettings.OffsetZ > Nz - offset - 1 + iWidth/2)
        {
            Phase.Fields(i, j, k).set_value(index2, 1);
        }
        else if (k+ locSettings.OffsetZ >= offset - iWidth/2 && k+ locSettings.OffsetZ <= offset + iWidth/2)
        {
            double IntProf = 0.5 - 0.5*sin(Pi*((k+ locSettings.OffsetZ) - offset)/iWidth);
            Phase.Fields(i, j, k).set_value(index1, 1.0 - IntProf);
            Phase.Fields(i, j, k).set_value(index2, IntProf);

            Phase.Fields(i,j,k).flag = 2;
        }
        else if (k+ locSettings.OffsetZ >= Nz - offset - 1 - iWidth/2 && k+ locSettings.OffsetZ <= Nz - offset - 1 + iWidth/2)
        {
            double IntProf = 0.5 + 0.5*sin(Pi*(Nz - (k + locSettings.OffsetZ)- offset - 1)/iWidth);
            Phase.Fields(i, j, k).set_value(index1, IntProf);
            Phase.Fields(i, j, k).set_value(index2, 1.0 - IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }
    return vector<size_t>({index1,index2});
}

vector<size_t> Initializations::TwoDifferentWalls(PhaseField& Phase, size_t ChannelPhaseIndex, size_t WallsPhaseIndex,
                                 double WallsThickness, BoundaryConditions& BC, Settings& locSettings)
{
    int Nz = locSettings.Nz;
    double iWidth = locSettings.iWidth;

    size_t index1 = Phase.AddGrainInfo(ChannelPhaseIndex);
    size_t index2 = Phase.AddGrainInfo(WallsPhaseIndex);
    size_t index3 = Phase.AddGrainInfo(WallsPhaseIndex);

    int offset = WallsThickness;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();

        if (k > offset + iWidth/2 && k < (Nz) - offset - iWidth/2 - 1)
        {
            Phase.Fields(i, j, k).set_value(index1, 1);
        }
        else if (k < offset - iWidth/2)
        {
            Phase.Fields(i, j, k).set_value(index2, 1);
        }
        else if (k > (Nz) - offset + iWidth/2 - 1)
        {
            Phase.Fields(i, j, k).set_value(index3, 1);
        }
        else if (k >= offset - iWidth/2 && k <= offset + iWidth/2)
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(k - offset)/iWidth);
            Phase.Fields(i, j, k).set_value(index1, 1.0 - IntProf);
            Phase.Fields(i, j, k).set_value(index2, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
        else if (k >= (Nz) - offset - iWidth/2 - 1 && k <= (Nz) - offset + iWidth/2 - 1)
        {
            double IntProf = 0.5 + 0.5*sin(Pi*((Nz) - k - offset - 1)/iWidth);
            Phase.Fields(i, j, k).set_value(index1, IntProf);
            Phase.Fields(i, j, k).set_value(index3, 1.0 - IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }
    return vector<size_t>({index1,index2,index3});
}

double PinE(const double x, const double y, const double z, const double a, const double b, const double c) /// Auxiliary function, check whether a point (x,y,z) lies within ellipsoid with axes a, b, c
{
  return (x*x/a/a + y*y/b/b + z*z/c/c - 1.0);
}

size_t Initializations::Ellipsoid(PhaseField& Phase,
                                size_t PhaseIndex,
                                double RadiusX, double RadiusY, double RadiusZ,
                                double x0, double y0, double z0,
                                BoundaryConditions& BC, Settings& locSettings)
{
    double iWidth = locSettings.iWidth;
    size_t index = Phase.AddGrainInfo(PhaseIndex);

    ///< Inner ellipsoid
    double Rxi = RadiusX - iWidth*0.5;
    double Ryi = RadiusY - iWidth*0.5;
    double Rzi = RadiusZ - iWidth*0.5;
    ///< Outer ellipsoid
    double Rxo = RadiusX + iWidth*0.5;
    double Ryo = RadiusY + iWidth*0.5;
    double Rzo = RadiusZ + iWidth*0.5;

    double Rxi2 = Rxi*Rxi;
    double Ryi2 = Ryi*Ryi;
    double Rzi2 = Rzi*Rzi;

    double Rxo2 = Rxo*Rxo;
    double Ryo2 = Ryo*Ryo;
    double Rzo2 = Rzo*Rzo;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        double radX = (i+Phase.OffsetX-x0);
        double radY = (j+Phase.OffsetY-y0);
        double radZ = (k+Phase.OffsetZ-z0);

        if ((radX*radX/Rxi2 + radY*radY/Ryi2 + radZ*radZ/Rzi2) <= 1.0)  ///< Point in the pure phase
        {
            Phase.Fields(i, j, k).clear();     ///< remove all phase fields that may already be present
            Phase.Fields(i, j, k).set_value(index, 1.0);
        }
        else if ((radX*radX/Rxo2 + radY*radY/Ryo2 + radZ*radZ/Rzo2) <= 1.0)  ///< Point in the interface
        {
            // Find coordinates inside the interface
            double r1 = -iWidth*0.5;
            double r2 =  iWidth*0.5;
            double rr = (r1 + r2)*0.5;
            double tolerance = 1e-6;
            while (fabs(r2 - r1) > tolerance)
            {
                if(PinE(radX, radY, radZ, RadiusX+r1, RadiusY+r1, RadiusZ+r1) > 0.0)
                   r1 += 0.25*(r2 - r1);
                else r1 -= 0.25*(r2 - r1);
                if(PinE(radX, radY, radZ, RadiusX+r2, RadiusY+r2, RadiusZ+r2) < 0.0)
                    r2 -= 0.25*(r2 - r1);
                else r2 += 0.25*(r2 - r1);

                rr = (r1 + r2)*0.5;
            }

            double IntProf = 0.5 - 0.5*sin(Pi*rr/iWidth);

            for(auto alpha = Phase.Fields(i,j,k).begin();
                     alpha < Phase.Fields(i,j,k).end(); ++alpha)
            {
                alpha->value *= 1.0 - IntProf;
            }
            Phase.Fields(i, j, k).set_value(index, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }
    return index;
}

void Initializations::Read(PhaseField& Phase, string FileName, BoundaryConditions& BC, Settings& locSettings)
{
    Phase.Read(FileName);

    if(Phase.Resolution == Resolutions::Single)
    {
        STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
        {
            for(auto alpha = Phase.Fields(i,j,k).begin();
                     alpha != Phase.Fields(i,j,k).end(); alpha++)
            {
               if(alpha->index >= Phase.FieldsStatistics.size())
               {
                   Phase.FieldsStatistics.Reallocate(alpha->index+1);
               }
               Phase.FieldsStatistics[alpha->index].Exist = true;
           }
        }
        STORAGE_LOOP_END
    }

    if(Phase.Resolution == Resolutions::Double)
    {
        STORAGE_LOOP_BEGIN(i,j,k,Phase.FieldsDR,Phase.FieldsDR.Bcells())
        {
            for(auto alpha = Phase.FieldsDR(i,j,k).begin();
                     alpha != Phase.FieldsDR(i,j,k).end(); alpha++)
            {
                if(alpha->index >= Phase.FieldsStatistics.size())
                {
                    Phase.FieldsStatistics.Reallocate(alpha->index+1);
                }
                Phase.FieldsStatistics[alpha->index].Exist = true;
            }
        }
        STORAGE_LOOP_END
    }
    Phase.Finalize(BC);
}

size_t Initializations::Rectangular(PhaseField& Phase, size_t PhaseIndex,
        double Lx, double Ly, double Lz,
        double x0, double y0, double z0,
        const BoundaryConditions& BC, const Settings& locSettings, const bool Finalize)
{
    const double iWidth = locSettings.iWidth;
    // Determine the area where the phase-field is one
    const double SizeX = (Lx < iWidth) ? 1: (Lx - iWidth)/2.0;
    const double SizeY = (Ly < iWidth) ? 1: (Ly - iWidth)/2.0;
    const double SizeZ = (Lz < iWidth) ? 1: (Lz - iWidth)/2.0;

    const size_t locIndex = Phase.AddGrainInfo(PhaseIndex);
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        // Determine relative coordinates
        const double x = double(i+locSettings.OffsetX);
        const double y = double(j+locSettings.OffsetY);
        const double z = double(k+locSettings.OffsetZ);
        const dVector3 pos  ({x , y, z});
        const dVector3 pos0 ({x0,y0,z0});
        dVector3 dist = CommonFunctions::Distance(pos, pos0, BC, locSettings);
        dist[0] = std::abs(dist[0]) - SizeX;
        dist[1] = std::abs(dist[1]) - SizeY;
        dist[2] = std::abs(dist[2]) - SizeZ;

        // Set cubic phase field
        if (dist[0] <= 0 and dist[1] <= 0 and dist[2] <= 0)
        {
            Phase.Fields(i, j, k).clear();
            Phase.Fields(i, j, k).set_value(locIndex, 1);
        }
        else if (dist[0] <= iWidth and dist[1] <= iWidth and dist[2] <= iWidth)
        {
            const double xx = std::max(dist[0], std::max(dist[1], dist[2]));
            const double IntProf = 0.5 + 0.5 * std::cos(Pi * xx / iWidth);
            Phase.Fields(i,j,k).flag = 2;

            double tmp = IntProf;
            for(auto alpha = Phase.Fields(i,j,k).begin();
                     alpha != Phase.Fields(i,j,k).end(); alpha++)
            if (alpha->value > tmp)
            {
                alpha->value -= tmp;
                tmp = 0;
            }
            else
            {
                tmp -= alpha->value;
                alpha->value = 0;
            }
            Phase.Fields(i,j,k).set_value(locIndex, IntProf);
        }
    }
    STORAGE_LOOP_END

    if (Finalize)
    {
        Phase.Finalize(BC);

        if(locSettings.Resolution == Resolutions::Double)
        {
            Phase.Refine();
        }
    }
    return locIndex;
}

size_t Initializations::CylinderSimple(PhaseField& Phase,
        const size_t PhaseIndex, const double Radius, const double length, const int Axis,
        const double x0, const double y0, const double z0,
        const BoundaryConditions& BC, const Settings& locSettings)
{
    double iWidth = locSettings.iWidth;
    if(locSettings.Resolution == Resolutions::Double)
    {
        iWidth *= 0.5;
    }
    const size_t locIndex = Phase.AddGrainInfo(PhaseIndex);

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        dVector3 pos__plane ({double(i+locSettings.OffsetX),double(j),double(k)});
        dVector3 pos0_plane ({double(x0),double(y0),double(z0)});
        dVector3 pos__axis  ({0,0,0});
        dVector3 pos0_axis  ({0,0,0});
        pos__axis [Axis] = pos__plane[Axis];
        pos0_axis [Axis] = pos0_plane[Axis];
        pos__plane[Axis] = 0;
        pos0_plane[Axis] = 0;

        dVector3 dist_plane = CommonFunctions::Distance(pos__plane, pos0_plane, BC, locSettings);
        dVector3 dist_axis  = CommonFunctions::Distance(pos__axis,  pos0_axis,  BC, locSettings);

        const double rad = dist_plane.abs();
        const double z   = dist_axis.abs();

        if ((rad < Radius - iWidth*0.5) and (z < length/2 - iWidth*0.5))
        {
            Phase.Fields(i,j,k).clear();
            Phase.Fields(i,j,k).set_value(locIndex, 1.0);
        }
        else if ((rad <= Radius + iWidth*0.5) and (z <= length/2 + iWidth*0.5))
        {
            double IntProf = 0;
            if ((rad > Radius - iWidth*0.5) and (z > length/2 - iWidth*0.5))
            {
                const double dist = std::sqrt(rad*rad + z*z);
                const double Dist = std::sqrt(Radius*Radius + length*length/4);
                IntProf = 0.5 - 0.5*std::sin(Pi*(dist - Dist)/iWidth);
            }
            else if (z > length/2 - iWidth*0.5)
            {
                IntProf = 0.5 - 0.5*std::sin(Pi*(z- length/2)/iWidth);
            }
            else
            {
                IntProf = 0.5 - 0.5*std::sin(Pi*(rad - Radius)/iWidth);
            }

            for(auto alpha = Phase.Fields(i,j,k).begin();
                     alpha < Phase.Fields(i,j,k).end(); alpha++)
            {
                alpha->value *= 1.0 - IntProf;
            }
            Phase.Fields(i,j,k).set_value(locIndex, IntProf);
        }
    }
    STORAGE_LOOP_END
    if(locSettings.Resolution == Resolutions::Double)
    {
        Phase.Refine();
    }
    Phase.Finalize(BC);

    return locIndex;
}


size_t Initializations::Cylinder(PhaseField& Phase, size_t PhaseIndex,
                               double Radius, double length, int Axis,
                               double x0, double y0, double z0,
                               const BoundaryConditions& BC,
                               Settings& locSettings)
{
    double iWidth = locSettings.iWidth;
    int halfIWidth = iWidth/2;

    size_t PhaseFieldIndex = Phase.AddGrainInfo(PhaseIndex);                    //TODO: TO CHECK

    switch(Axis)
    {
        case 0:
        {
    //        SphereFixedIdx(Phase, PhaseIndex, iWidth*0.5, x0-length*0.5+halfIWidth, y0, z0,
    //                locSettings);
    //        SphereFixedIdx(Phase, PhaseIndex, iWidth*0.5, x0+length*0.5-halfIWidth, y0, z0,
    //                locSettings);

            for(int i = x0-length*0.5 + halfIWidth; i < x0+length*0.5 - halfIWidth; i++)
            {
                Disc(Phase, PhaseFieldIndex, Radius, Axis, i,y0,z0, locSettings);
            }
            break;
        }
        case 1:
        {
    //        SphereFixedIdx(Phase, PhaseIndex, iWidth*0.5,
    //                          x0, y0-length*0.5, z0, locSettings);
    //        SphereFixedIdx(Phase, PhaseIndex, iWidth*0.5,
    //                          x0, y0+length*0.5, z0, locSettings);

            for(int j = y0-length*0.5 + halfIWidth; j < y0+length*0.5 - halfIWidth; j++)
            {
                Disc(Phase, PhaseFieldIndex, Radius, Axis, x0,j,z0, locSettings);
            }
            break;
        }
        case 2:
        {
    //        SphereFixedIdx(Phase, PhaseIndex, iWidth*0.5,
    //                          x0, y0, z0-length*0.5, locSettings);
    //        SphereFixedIdx(Phase, PhaseIndex, iWidth*0.5,
    //                          x0, y0, z0+length*0.5, locSettings);

            for(int k = z0-length*0.5 + halfIWidth; k < z0+length*0.5 - halfIWidth; k++)
            {
                Disc(Phase, PhaseFieldIndex, Radius, Axis, x0, y0, k, locSettings);
            }

            break;
        }
        default:
        {
            stringstream message;
            message<<"Axis = "<<Axis<<"! Choose a value between 0 and 2!\n";
            Info::WriteExit(message.str(), "Initialisations", "Cylinder()");
            exit(13);
        }
    }
    if(locSettings.Resolution == Resolutions::Double)
    {
        Phase.Refine();
    }
    Phase.Finalize(BC);
    return PhaseFieldIndex;
}

void Initializations::Disc(PhaseField& Phase, size_t PhaseFieldIndex, double Radius,
         int NormalAxis, double x0, double y0, double z0, Settings& locSettings)
{
    double iWidth = locSettings.iWidth;

//    int locIndex = Phase.AddGrainInfo(PhaseIndex);

    switch(NormalAxis)
    {
        case 0:
        {
            for(int j = y0-Radius-iWidth/2-1; j < y0+Radius+iWidth/2+1; ++j)
            for(int k = z0-Radius-iWidth/2-1; k < z0+Radius+iWidth/2+1; ++k)
            {
                double rad = sqrt((j-y0)*(j-y0)+(k-z0)*(k-z0));
                if (rad < Radius - iWidth*0.5)
                {
                    Phase.Fields(x0, j, k).clear();
                    Phase.Fields(x0, j, k).set_value(PhaseFieldIndex, 1.0);
                }
                else if (rad < Radius + iWidth*0.5)
                {
                    double IntProf = 0.5 - 0.5*sin(Pi*(rad - Radius)/iWidth);
                    if(IntProf > Phase.Fields(x0, j, k)[PhaseFieldIndex])
                    {
                        for(auto alpha = Phase.Fields(x0,j,k).begin();
                                 alpha < Phase.Fields(x0,j,k).end(); alpha++)
                        {
                            alpha->value *= 1.0 - IntProf;
                        }
                        Phase.Fields(x0, j, k).set_value(PhaseFieldIndex, IntProf);
                        Phase.Fields(x0, j, k).flag = 2;
                    }
                }
            }
            break;
        }
        case 1:
        {
            for(int i = x0-Radius-iWidth/2-1; i < x0+Radius+iWidth/2+1; ++i)
            for(int k = z0-Radius-iWidth/2-1; k < z0+Radius+iWidth/2+1; ++k)
            {
                double rad = sqrt((i-x0)*(i-x0)+(k-z0)*(k-z0));
                if (rad < Radius - iWidth*0.5)
                {
                    Phase.Fields(i, y0, k).clear();
                    Phase.Fields(i, y0, k).set_value(PhaseFieldIndex, 1.0);
                }
                else if (rad < Radius + iWidth*0.5)
                {
                    double IntProf = 0.5 - 0.5*sin(Pi*(rad - Radius)/iWidth);
                    if(IntProf > Phase.Fields(i, y0, k)[PhaseFieldIndex])
                    {
                        for(auto alpha = Phase.Fields(i,y0,k).begin();
                                 alpha < Phase.Fields(i,y0,k).end(); alpha++)
                        {
                            alpha->value *= 1.0 - IntProf;
                        }
                        Phase.Fields(i, y0, k).set_value(PhaseFieldIndex, IntProf);
                        Phase.Fields(i, y0, k).flag = 2;
                    }
                }
            }
            break;
        }
        case 2:
        {
            for(int i = x0-Radius-iWidth/2-1; i < x0+Radius+iWidth/2+1; ++i)
            for(int j = y0-Radius-iWidth/2-1; j < y0+Radius+iWidth/2+1; ++j)
            {
                double rad = sqrt((i-x0)*(i-x0)+(j-y0)*(j-y0));
                if (rad < Radius - iWidth*0.5)
                {
                    Phase.Fields(i, j, z0).clear();
                    Phase.Fields(i, j, z0).set_value(PhaseFieldIndex, 1.0);
                }
                else if (rad < Radius + iWidth*0.5)
                {
                    double IntProf = 0.5 - 0.5*sin(Pi*(rad - Radius)/iWidth);
                    if(IntProf > Phase.Fields(i, j, z0)[PhaseFieldIndex])
                    {
                        for(auto alpha = Phase.Fields(i,j,z0).begin();
                                 alpha < Phase.Fields(i,j,z0).end(); alpha++)
                        {
                            alpha->value *= 1.0 - IntProf;
                        }
                        Phase.Fields(i, j, z0).set_value(PhaseFieldIndex, IntProf);
                        Phase.Fields(i, j, z0).flag = 2;
                    }
                }
            }
            break;
        }
        default:
        {
            stringstream message;
            message << "NormalAxis = " << NormalAxis
                    << "! Choose a value between 0 and 2!\n";
            Info::WriteExit(message.str(), thisclassname, "Disc()");
            exit(13);
        }
    }
}

void Initializations::SphereFixedIdx(PhaseField& Phase, size_t PhaseIndex,
                             double Radius, double x0, double y0, double z0,
                             BoundaryConditions& BC, Settings& locSettings)
{
    const double iWidth = (locSettings.Resolution == Resolutions::Double) ? 0.5*locSettings.iWidth : locSettings.iWidth;
    Phase.FieldsStatistics[PhaseIndex].Exist = true;
    Phase.FieldsStatistics[PhaseIndex].Stage = 0;

    auto set_sphere = [&iWidth,&Radius,&PhaseIndex,&Phase,&z0](long int i,long int j,long int k, double rad)
    {
       if (rad < Radius - iWidth*0.5)
        {
            Phase.Fields(i, j, k).clear();
            Phase.Fields(i, j, k).set_value(PhaseIndex, 1.0);
        }
        else if (rad < Radius + iWidth*0.5)
        {
            const double IntProf = 0.5 - 0.5*sin(Pi*(rad - Radius)/iWidth);
            if(IntProf > Phase.Fields(i, j, z0)[PhaseIndex])
            {
                for(auto alpha = Phase.Fields(i,j,k).begin();
                         alpha < Phase.Fields(i,j,k).end(); alpha++)
                {
                    alpha->value *= 1.0 - IntProf;
                }
                Phase.Fields(i, j, k).set_value(PhaseIndex, IntProf);
                Phase.Fields(i,j,k).flag = 2;
            }
        }
        return false;
    };
    const double loop_radius = Radius + 1.1*iWidth;
    loop_sphere(Phase.Fields, set_sphere, x0, y0, z0, loop_radius, locSettings, BC);

    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }
}

int Initializations::TwoDimEBSD(std::string filename, std::vector<int> columns,
        PhaseField& Phase, BoundaryConditions& BC, Settings& locSettings)
{
    // Method to read in 2D EBSD data from input file with the following format:
    //    x-coord     y-coordinate     phase-field index (this line is not required in the file)
    // "    0            0                    1       "
    // "    0.1          0                    2       "
    // "                     etc                      "

    // The columns are seperated by whitespaces, the exact column numbers can be given
    // via the columns vector using the following correlations:
    // x coordinate -> columns[0]
    // y coordinate -> columns[1]
    // phase field index -> columns[2]
    //
    // Example call:
    //    Initializations::TwoDimEBSDWithOrientations("EBSDmap.dat",
    //                      {4, 5, 6, 1, 2, 3}, Phi, EP, BC, OPSettings);
    //
    // Note: For hexagonal grids run /scripts/EBSDHexRegParser.m first (Matlab)

    if(columns.size() != 3)
    {
        columns.clear(); columns = {1, 2, 3};
        Info::WriteWarning("Column list incomplete. Set to default {1, 2, 3}.", thisclassname, "TwoDimEBSD");
    }

    if(locSettings.Nz != 1)
    {
        Info::WriteExit("Only 2d supported", thisclassname, "TwoDimEBSD");
        exit(1);
    }

    ifstream inp(filename.c_str(), ios::out);

    if (!inp)
    {
        Info::WriteExit("File " + filename + " could not be opened",
                thisclassname, "TwoDimEBSD()");
        exit(1);
    };

    Info::WriteBlankLine();
    Info::WriteLineInsert("EBSD reader");
    Info::WriteStandard("Source", filename);
    Info::WriteStandard("Read x coordinates from column", std::to_string(columns[0]));
    Info::WriteStandard("Read y coordinates from column", std::to_string(columns[1]));
    Info::WriteStandard("Read phase indices from column", std::to_string(columns[2]));

    // Count number of lines
    int nol = 0;
    std::string line;
    while(getline (inp,line)) nol++;
    inp.clear(); // back to begin
    inp.seekg(0, ios::beg);
    Info::WriteStandard("Number of lines", std::to_string(nol));

    // Read input file

    vector<size_t> phaseindex;
    vector<size_t> individualphaseindices;
    vector<double> xcoord;
    vector<double> individualxcoord;
    vector<double> ycoord;
    vector<double> individualycoord;

    while(true)
    {
        std::string line;
        std::getline( inp, line );
        if( !inp ) break;
        std::istringstream iline( line );
//        if (line.size() > 0) cout << line << endl;

        int spos = 1;

        while(true)
        {
            std::string str;
            std::getline(iline, str, ' ' );

            if(spos == columns[0] and str.size() > 0) // Read x coordinates
            {
                double xx = std::stod(str);

                // Identify individual number of xcoordinates
                if(std::find(individualxcoord.begin(), individualxcoord.end(), xx) == individualxcoord.end())
                {
                    individualxcoord.push_back(xx);
    //                cout << xx << endl;
                }
                xcoord.push_back(xx);
            }

            if(spos == columns[1] and str.size() > 0)  // Read y coordinates
            {
                double yy = std::stod(str);
                // Identify individual number of ycoordinates
                if(std::find(individualycoord.begin(), individualycoord.end(), yy) == individualycoord.end())
                {
                    individualycoord.push_back(yy);
                }
                ycoord.push_back(yy);
            }

            if(spos == columns[2] and str.size() > 0) // Read phase field index
            {
                int lindex = std::stoi(str);
                // Identify individual number of phases
                if (std::find(individualphaseindices.begin(), individualphaseindices.end(), lindex) == individualphaseindices.end())
                {
                    individualphaseindices.push_back(lindex);
                }
                phaseindex.push_back(lindex);
            }
            if (str.size() > 0) spos++;
            if( !iline ) break;
        }
    }

    size_t individualgrainnum = individualphaseindices.size();
    Info::WriteStandard("Number of individual phases", std::to_string(individualgrainnum));
    size_t maxphaseindex = *max_element(individualphaseindices.begin(), individualphaseindices.end());
    Info::WriteStandard("Largest phase-field index", std::to_string(maxphaseindex));
    size_t pfsize = phaseindex.size();
    Info::WriteStandard("Phasefield storage size", std::to_string(pfsize));
    size_t NxEBSD = individualxcoord.size();
    Info::WriteStandard("Number of (individual) x coordinates", std::to_string(NxEBSD));
    double maxX = xcoord.back();
    Info::WriteStandard("Biggest X Coord", std::to_string(maxX));
    double distX = maxX/double(NxEBSD);
    Info::WriteStandard("Rastering X", std::to_string(distX));
    size_t NyEBSD = individualycoord.size();
    Info::WriteStandard("Number of (individual) y coordinates", std::to_string(NyEBSD));
    double maxY = ycoord.back();
    Info::WriteStandard("Biggest Y Coord", std::to_string(maxY));
    double distY = maxY/double(NyEBSD);
    Info::WriteStandard("Rastering Y", std::to_string(distY));

//    for(size_t i = 0 ; i < phaseindex.size(); ++i) cout << phaseindex[i] << " ";
//    cout << endl;
//
//    for(size_t i = 0 ; i < xcoord.size(); ++i) cout << xcoord[i] << " ";
//    cout << endl;
//
    stringstream outPhaseIndeces;
    for(size_t i = 0 ; i < individualphaseindices.size(); ++i) outPhaseIndeces << individualphaseindices[i] << " ";
    Info::WriteSimple(outPhaseIndeces.str());
    Info::WriteBlankLine();

    const int Nx = locSettings.Nx;
    const int Ny = locSettings.Ny;

    for (size_t idx = 0; idx < maxphaseindex+1; idx++)
    {
        if(std::find(individualphaseindices.begin(), individualphaseindices.end(), idx) != individualphaseindices.end())
        {
            size_t locIndex = Phase.AddGrainInfo(0) + idx%locSettings.Nphases; // eigene zv

            Info::WriteStandard("Add grain/phase", std::to_string(locIndex));

            int incx = 1; // Skip points in x-direction
            int incy = 1; // Skip points in y-direction

            int itery = 1;
            int colstep = 0;
            for(int j = 0; j < Ny; ++j)
            {
                int iterx = 1;

                for(int i = 0; i < Nx; ++i)
                {
                    if (idx + 1 == phaseindex[-1 + iterx + colstep])
    //                if ((locIndex)%locSettings.Nphases + 1 + idx == phaseindex[-1 + iterx + colstep])
                    {
    //                    cout << i << " " << j << " " << -1 + iterx + (itery-1)*Nx << endl;
                        Phase.Fields(i, j, 0).clear();
                        Phase.Fields(i, j, 0).flag = 2;
                        Phase.Fields(i, j, 0).set_value(locIndex, 1.0);
                    }
                    iterx += incx;
                }
                if(incx == 1 and incy == 1)
                {
                    colstep += Nx;
                }
                else
                {
                // Odd columns
                    if(j%2 == 1) colstep += Nx;
                    else colstep += (Nx-1);
                }
                itery += incy;
            }
        }
    }

    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }
    Info::WriteStandard("EBSD", "done");
    Info::WriteLine();

    return 0;
}

int Initializations::TwoDimEBSDWithOrientations(std::string filename, std::vector<int> columns, std::string anglerepresentation,
        PhaseField& Phase, ElasticProperties& EP, Orientations& OR, BoundaryConditions& BC, Settings& locSettings)
{
    // Method to read in 2D EBSD data from input file with the following format:
    //    x-coord     y-coordinate     phase-field index (this line is not required in the file)
    // "    0            0                    1       "
    // "    0.1          0                    2       "
    // "                     etc                      "

    // The columns are seperated by whitespaces, the exact column numbers can be given
    // via the columns vector using the following correlations:
    // x coordinate -> columns[0]
    // y coordinate -> columns[1]
    // phase field index -> columns[2]
    // Euler angle 1 -> columns[3]
    // Euler angle 2 -> columns[4]
    // Euler angle 3 -> columns[5]
    //
    // Example call:
    //    Initializations::TwoDimEBSDWithOrientations("EBSDmap.dat",
    //                      {4, 5, 6, 1, 2, 3}, Phi, EP, BC, OPSettings);
    //
    // Note: For hexagonal grids run /scripts/EBSDHexRegParser.m first (Matlab)

    const int Nx = locSettings.Nx;
    const int Ny = locSettings.Ny;
    const int Nz = locSettings.Nz;

    if(columns.size() != 6)
    {
        columns.clear(); columns = {1, 2, 3, 4, 5, 6};
        Info::WriteWarning("Column list incomplete. Set to default {1, 2, 3, 4, 5, 6}.",
                thisclassname, "TwoDimEBSD");
    }

    if(locSettings.Nz != 1)
    {
        Info::WriteExit("Only 2d supported", thisclassname, "TwoDimEBSD");
        exit(1);
    }

    ifstream inp(filename.c_str(), ios::out);

    if (!inp)
    {
        Info::WriteExit("File " + filename + " could not be opened",
                thisclassname, "TwoDimEBSD()");
        exit(1);
    };

    Info::WriteBlankLine();
    Info::WriteLineInsert("EBSD reader");
    Info::WriteStandard("Source", filename);
    Info::WriteStandard("Read x coordinates from column", std::to_string(columns[0]));
    Info::WriteStandard("Read y coordinates from column", std::to_string(columns[1]));
    Info::WriteStandard("Read phase indices from column", std::to_string(columns[2]));
    Info::WriteStandard("Read Euler angle 1 from column", std::to_string(columns[3]));
    Info::WriteStandard("Read Euler angle 2 from column", std::to_string(columns[4]));
    Info::WriteStandard("Read Euler angle 3 from column", std::to_string(columns[5]));

    // Count number of lines
    int nol = 0;
    std::string line;
    while(getline (inp,line)) nol++;
    inp.clear(); // back to begin
    inp.seekg(0, ios::beg);
    Info::WriteStandard("Number of lines", std::to_string(nol));

    // Read input file

    vector<size_t> phaseindex;
    vector<size_t> individualphaseindices;
    vector<double> xcoord;
    vector<double> individualxcoord;
    vector<double> ycoord;
    vector<double> individualycoord;
    vector<EulerAngles> Eang;

    while(true)
    {
        std::string line;
        std::getline( inp, line );
        if( !inp ) break;
        std::istringstream iline( line );
//        if (line.size() > 0) cout << line << endl;

        EulerAngles tempEang;
        tempEang.set_to_zero();
        int spos = 1;

        while(true)
        {
            std::string str;
            std::getline(iline, str, ' ' );

            if(spos == columns[0] and str.size() > 0) // Read x coordinates
            {
                double xx = std::stod(str);

                // Identify individual number of xcoordinates
                if(std::find(individualxcoord.begin(), individualxcoord.end(), xx) == individualxcoord.end())
                {
                    individualxcoord.push_back(xx);
    //                cout << xx << endl;
                }
                xcoord.push_back(xx);
            }

            if(spos == columns[1] and str.size() > 0)  // Read y coordinates
            {
                double yy = std::stod(str);
                // Identify individual number of ycoordinates
                if(std::find(individualycoord.begin(), individualycoord.end(), yy) == individualycoord.end())
                {
                    individualycoord.push_back(yy);
                }
                ycoord.push_back(yy);
            }

            if(spos == columns[2] and str.size() > 0) // Read phase field index
            {
                int lindex = std::stoi(str);
                // Identify individual number of phases
                if (std::find(individualphaseindices.begin(), individualphaseindices.end(), lindex) == individualphaseindices.end())
                {
                    individualphaseindices.push_back(lindex);
                }
                phaseindex.push_back(lindex);
            }

            double anglefactor = 180.0/Pi;
            if(!anglerepresentation.compare("degree") or !anglerepresentation.compare("deg")
                    or !anglerepresentation.compare("Degree") or !anglerepresentation.compare("Deg"))
            {
                anglefactor = 1.0;
            }

            if(spos == columns[3] and str.size() > 0) // Read Euler angle 1
            {
                tempEang.Q[0] = std::stod(str)*anglefactor;
            }

            if(spos == columns[4] and str.size() > 0) // Read Euler angle 2
            {
                tempEang.Q[1] = std::stod(str)*anglefactor;
            }

            if(spos == columns[5] and str.size() > 0) // Read Euler angle 3
            {
                tempEang.Q[2] = std::stod(str)*anglefactor;
            }

            if (str.size() > 0) spos++;
            if( !iline ) break;
        }

        tempEang.set_convention(ZXZ);

        Info::WriteWarning("Using ZXZ angle convention", thisclassname, "TwoDimEBSDWithOrientations()");

        tempEang.setTrigonometricFunctions();
        Eang.push_back(tempEang);
    }

    size_t individualgrainnum = individualphaseindices.size();
    Info::WriteStandard("Number of individual phases", std::to_string(individualgrainnum));
    size_t pfsize = phaseindex.size();
    Info::WriteStandard("Phasefield storage size", std::to_string(pfsize));
    size_t NxEBSD = individualxcoord.size();
    Info::WriteStandard("Number of (individual) x coordinates", std::to_string(NxEBSD));
    double maxX = xcoord.back();
    Info::WriteStandard("Biggest X Coord", std::to_string(maxX));
    double distX = maxX/double(NxEBSD);
    Info::WriteStandard("Rastering X", std::to_string(distX));
    size_t NyEBSD = individualycoord.size();
    Info::WriteStandard("Number of (individual) y coordinates", std::to_string(NyEBSD));
    double maxY = ycoord.back();
    Info::WriteStandard("Biggest Y Coord", std::to_string(maxY));
    double distY = maxY/double(NyEBSD);
    Info::WriteStandard("Rastering Y", std::to_string(distY));

    Storage3D <EulerAngles, 0> StorageEulerAngles;
    StorageEulerAngles.Allocate(Nx, Ny, 1, 1,1,0, 0);

    for (size_t idx = 0; idx < individualgrainnum; idx++)
    {
        size_t locIndex = Phase.AddGrainInfo(0) + idx%locSettings.Nphases;

        Info::WriteStandard("Add grain/phase", std::to_string(locIndex));

        int incx = 1; // Skip points in x-direction
        int incy = 1; // Skip points in y-direction

        int itery = 1;
        int colstep = 0;
        for(int j = 0; j < Ny; ++j)
        {
            int iterx = 1;

            for(int i = 0; i < Nx; ++i)
            {
                // SetInitialOrientations
                EulerAngles EangLocal = Eang[-1 + iterx + colstep];
                OR.Quaternions(i,j,0) = EangLocal.getQuaternion();
                StorageEulerAngles(i,j,0) = EangLocal;

                // Set phase fields
                if (idx + 1 == phaseindex[-1 + iterx + colstep])
//                if ((locIndex)%locSettings.Nphases + 1 + idx == phaseindex[-1 + iterx + colstep])
                {
//                    cout << i << " " << j << " " << -1 + iterx + (itery-1)*Nx << endl;
                    Phase.Fields(i, j, 0).clear();
                    Phase.Fields(i, j, 0).flag = 2;
                    Phase.Fields(i, j, 0).set_value(locIndex, 1.0);
                }
                iterx += incx;
            }
            if(incx == 1 and incy == 1)
            {
                colstep += Nx;
            }
            else
            {
            // Odd columns
                if(j%2 == 1) colstep += Nx;
                else colstep += (Nx-1);
            }
            itery += incy;
        }
    }

    // Create output to VTK

    stringstream outbuffer;

    outbuffer << "# vtk DataFile Version 3.0\n";
    outbuffer << "InitialEulerAngles\n";
    outbuffer << "ASCII\n";
    outbuffer << "DATASET STRUCTURED_GRID\n";
    outbuffer << "DIMENSIONS " << Nx << " " << Ny << " " << 1 << "\n";
    outbuffer << "POINTS " <<  Nx*Ny << " int\n";

    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        outbuffer << i << " " << j << " " << k << "\n";
    }
    outbuffer << " \n";
    outbuffer << "POINT_DATA " << Nx*Ny << " \n";

    outbuffer << "SCALARS Eang_" << 1 << " double\n";
    outbuffer << "LOOKUP_TABLE default\n";
    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        outbuffer << StorageEulerAngles(i,j,k).Q[0] << " ";
    }
    outbuffer << " \n";
    outbuffer << "SCALARS Eang_" << 2 << " double\n";
    outbuffer << "LOOKUP_TABLE default\n";
    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        outbuffer << StorageEulerAngles(i,j,k).Q[1] << " ";
    }
    outbuffer << " \n";
    outbuffer << "SCALARS Eang_" << 3 << " double\n";
    outbuffer << "LOOKUP_TABLE default\n";
    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        outbuffer << StorageEulerAngles(i,j,k).Q[2] << " ";
    }

    string FileName = DefaultVTKDir + "InitializedEulerAngles.vtk";

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbuffer.rdbuf();
    vtk_file.close();

    // Write table of orientations

    outbuffer.str("");

    for(size_t i = 0 ; i < Eang.size(); ++i)
    {
        outbuffer << i << "  " << Eang[i].Q[0] << "  "
                               << Eang[i].Q[1] << "  "
                               << Eang[i].Q[2] << endl;
    }
    FileName = "EBSDorientations.dat";

    ofstream orientations_file(FileName.c_str());
    orientations_file << outbuffer.rdbuf();
    orientations_file.close();

    // End write Euler angles

    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }
    Info::WriteStandard("EBSD", "done");
    Info::WriteLine();

    return 0;
}

size_t Initializations::SphereInGrain(PhaseField& Phase,
        const size_t ParrentPhaseFieldIndex, const size_t PhaseIndex,
        const double Radius, const double x0, const double y0, const double z0,
        const BoundaryConditions& BC, const Settings& locSettings,
        const bool Finalize)
{
    const double iWidth = (locSettings.Resolution == Resolutions::Double) ? 0.5*locSettings.iWidth : locSettings.iWidth;
    const size_t locIndex = Phase.AddGrainInfo(PhaseIndex);
    auto set_sphere = [&iWidth,&Radius,&locIndex,&Phase,&ParrentPhaseFieldIndex](long int i,long int j,long int k, double rad)
    {
        //const double ParrentPhaseAmount = Phase.Fields(i, j, k)[ParrentPhaseFieldIndex];

        if (rad < (Radius - iWidth*0.5))
        {
            Phase.Fields(i,j,k).clear();
            Phase.Fields(i,j,k).set_value(locIndex, 1.0);
            Phase.Fields(i,j,k).flag = 2;
        }
        else if (rad < (Radius + iWidth*0.5))
        {
            const double Profile = (0.5 - 0.5*sin(Pi*(rad - Radius)/iWidth));
            Phase.Fields(i, j, k).set_value(ParrentPhaseFieldIndex, /*ParrentPhaseAmount**/(1.0-Profile));
            Phase.Fields(i, j, k).set_value(locIndex, /*ParrentPhaseAmount**/Profile);
            Phase.Fields(i,j,k).flag = 2;
        }
        return false;
    };
    const double loop_radius = Radius + 1.1*iWidth;
    loop_sphere(Phase.Fields, set_sphere, x0, y0, z0, loop_radius, locSettings, BC);
    if (Finalize)
    {
        if(locSettings.Resolution == Resolutions::Double)
        {
            Phase.Refine();
        }
        Phase.Finalize(BC);
    }
    Phase.FinalizeSR(BC);
    return locIndex;
}

size_t  Initializations::FillGrainWithSpheres(PhaseField& Phase,
        size_t ParrentPhaseFieldIndex, size_t SpheresPhaseIndex,
        double MinR, double MaxR, BoundaryConditions& BC, Settings& locSettings,
        size_t Nspheres, double MinDistance)
{
    if (MinDistance < 0.0) MinDistance = MinR;

    vector<iVector3> spheres;
    vector<double> Radii;

    int counter = 0;
    int globalCounter = 0;
    const int MaxNumberOfTries = 10.0*(locSettings.TotalNx)*
                                      (locSettings.TotalNy)*
                                      (locSettings.TotalNz);

    //default_random_engine generator;
    std::random_device generator;
    //std::mt19937_64 generator(rd());
    uniform_int_distribution<int> distributionX(0, (locSettings.TotalNx-1));
    uniform_int_distribution<int> distributionY(0, (locSettings.TotalNy-1));
    uniform_int_distribution<int> distributionZ(0, (locSettings.TotalNz-1));
    uniform_real_distribution<double> distributionR(MinR, MaxR);

    while(counter < MaxNumberOfTries)
    {
        int i = distributionX(generator);
        int j = distributionY(generator);
        int k = distributionZ(generator);
        double Radius  = distributionR(generator);
        #ifdef MPI_PARALLEL
        MPI_Bcast(&(i), 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(j), 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(k), 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(Radius), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        #endif

        bool overlapping = false;
        for(size_t n = 0; n < spheres.size(); n++)
        {
            double distance = CommonFunctions::Distance({i,j,k}, spheres[n], BC, locSettings).abs();
            if(distance < Radius + Radii[n] + MinDistance) overlapping = true;
        }

        #ifdef MPI_PARALLEL
        double locParentPhaseFraction = 0;
        double ParentPhaseFraction = 0;
        if (i > locSettings.OffsetX and i < locSettings.OffsetX+locSettings.Nx)
        {
            locParentPhaseFraction = Phase.Fields(i - locSettings.OffsetX, j, k)[ParrentPhaseFieldIndex];
        }
        MPI_Allreduce(&(locParentPhaseFraction), &(ParentPhaseFraction), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        #else
            double ParentPhaseFraction = Phase.Fields(i,j,k)[ParrentPhaseFieldIndex];
        #endif

        if(overlapping or ParentPhaseFraction < 1.0)
        {
            counter++;
        }
        else
        {
            size_t sphereIndex =
            SphereInGrain(Phase, ParrentPhaseFieldIndex, SpheresPhaseIndex,
                          Radius, i, j, k, BC, locSettings, false);
            spheres.push_back({i,j,k});
            Radii.push_back(Radius);
            stringstream message;
            message << "*********************************************************\n";
            message << "after " << counter << " tries sphere nr: "
                    << spheres.size() << " at (" << i << "," << j << "," << k <<")"
                    << " with R = " << Radius<< " and phase field index = "
                    << sphereIndex << " set." << "\n";
            //Phase.PrintPointStatistics(i,j,k); //NOTE: not mpi-parallel
            message << "*********************************************************\n";
            Info::WriteSimple(message.str());
            globalCounter += counter;
            counter = 0;

            if(Nspheres != 0 and spheres.size() == Nspheres)
            {
                break;
            }
        }
    }
    cout << "after " << globalCounter << " tries " << spheres.size()
         << " spheres could be initialized!\n";

    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }
    return spheres.size();
}

std::vector<size_t> Initializations::FillRectangularWithSpheres(
        PhaseField& Phase,
        const BoundaryConditions& BC, const Settings& locSettings,
        std::function<size_t(long int, long int, long int)> PhaseIndex,
        double MeanRadius, double StdRadius,
        long int x_min, long int x_max,
        long int y_min, long int y_max,
        long int z_min, long int z_max,
        double MinDistance)
{
    const int Nx = locSettings.TotalNx;
    const int Ny = locSettings.TotalNy;
    const int Nz = locSettings.TotalNz;
    const double MinRadius = locSettings.iWidth/2;
    const double MaxRadius = std::max(Nx/2, std::max(Ny/2, Nz/2));

    if (StdRadius < 0.0)
    {
        Info::WriteExit("Standard deviation needs to be positive!", thisclassname, "FillWithSpheres");
        std::exit(EXIT_FAILURE);
    }

    if (MeanRadius < MinRadius)
    {
        Info::WriteExit("The mean radius needs to be at least half the interface width!", thisclassname, "FillWithSpheres");
        std::exit(EXIT_FAILURE);
    }
    if (MeanRadius > MaxRadius)
    {
        Info::WriteExit("The mean radius needs to be smaller than the simulation box!", thisclassname, "FillWithSpheres");
        std::exit(EXIT_FAILURE);
    }

    std::vector<size_t> SpheresIdx;

    std::random_device rd;
    std::mt19937_64 generator(rd());
    std::normal_distribution<double> radius_gen(MeanRadius,StdRadius);
    std::uniform_real_distribution<double> dist_01(0.0,1.0);
    std::uniform_real_distribution<double> phi_gen(0,2.*M_PI);
    std::uniform_real_distribution<double> theta_gen(0,M_PI);

    auto rand01 = std::bind(dist_01,generator);

    Storage3D<bool,0> Shielded;
    Shielded.Allocate(locSettings.Nx, locSettings.Ny, locSettings.Nz, locSettings.dNx, locSettings.dNy, locSettings.dNz, 0);
    struct Planted_t {int x; int y; int z; double radius;};
    std::vector<Planted_t> Planted;
    {
        // Die solid radius
        double radius = radius_gen(generator);
        #ifdef MPI_PARALLEL
        MPI_Bcast(&(radius), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        #endif

        while(radius < MinRadius or radius > MaxRadius) radius = radius_gen(generator);

        Info::Write("Plant initial Grain");
        assert(radius < MaxRadius);

        // Select random position inside boundaries
        Planted_t PInit;
        PInit.x = (Nx > 1 ) ? x_min + rand01()*(x_max-x_min) : 0;
        PInit.y = (Ny > 1 ) ? y_min + rand01()*(y_max-y_min) : 0;
        PInit.z = (Nz > 1 ) ? z_min + rand01()*(z_max-z_min) : 0;
        PInit.radius = radius;

        // Check if simulation box is 2D
        while (PInit.x >= Nx) PInit.x -= Nx;
        while (PInit.y >= Ny) PInit.y -= Ny;
        while (PInit.z >= Nz) PInit.z -= Nz;
        while (PInit.x <   0) PInit.x += Nx;
        while (PInit.y <   0) PInit.y += Ny;
        while (PInit.z <   0) PInit.z += Nz;

        // Initialize phase-field
        #ifdef MPI_PARALLEL
        MPI_Bcast(&(PInit.x), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(PInit.y), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(PInit.z), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        #endif
        Planted.push_back(PInit);

        const size_t phaseIdx =  PhaseIndex(PInit.x, PInit.y, PInit.z);
        const size_t idx = Initializations::Sphere(Phase, phaseIdx, radius, PInit.x, PInit.y, PInit.z, BC, locSettings, false);
        SpheresIdx.push_back(idx);

        // Shield surrounding against places of other solid particles
        auto set_shield = [&Shielded](long int i,long int j,long int k, double rad){Shielded(i,j,k) = true; return false;};
        Initializations::loop_sphere(Shielded, set_shield, PInit.x, PInit.y, PInit.z, radius + MinDistance/2.0, locSettings, BC);

        Info::Write("Planted Solid number ", Planted.size());

        size_t i = 0;
        while(i < Planted.size())
        {
            size_t N = 1000;//4.0/3.0*M_PI*(std::pow(2*radius,3) - std::pow(radius,3));
            for (size_t j = 0; j < N; j++)
            {
                radius = radius_gen(generator);
                #ifdef MPI_PARALLEL
                //MPI_Bcast(&(i),      1, MPI_INT, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&(j),      1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Bcast(&(radius), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                #endif
                if (radius < MinRadius) continue;

                // Place new solid next to the old one
                double phi   = phi_gen(generator);
                double theta = theta_gen(generator);
                #ifdef MPI_PARALLEL
                MPI_Bcast(&(phi),   1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Bcast(&(theta), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                #endif

                Planted_t PNew;
                PNew.x = (Nx > 1) ? Planted[i].x + (Planted[i].radius+MinDistance+radius+3)*cos(phi)   : 0;//*sin(theta);
                PNew.y = (Ny > 1) ? Planted[i].y + (Planted[i].radius+MinDistance+radius+3)*sin(phi)   : 0;//*sin(theta);
                PNew.z = (Nz > 1) ? Planted[i].z + (Planted[i].radius+MinDistance+radius+3)*cos(theta) : 0;//TODO z is too large
                PNew.radius = radius;
                // Check if new solid is in simulation box
                if (PNew.x >= Nx) continue;
                if (PNew.y >= Ny) continue;
                if (PNew.z >= Nz) continue;
                if (PNew.x <   0) continue;
                if (PNew.y <   0) continue;
                if (PNew.z <   0) continue;

                // Check if floor and ceiling are not touched
                if (x_min != 0  and PNew.x-PNew.radius < x_min) continue;
                if (x_max != Nx and PNew.x+PNew.radius > x_max) continue;
                if (y_min != 0  and PNew.y-PNew.radius < y_min) continue;
                if (y_max != Ny and PNew.y+PNew.radius > y_max) continue;
                if (z_min != 0  and PNew.z-PNew.radius < z_min) continue;
                if (z_max != Nz and PNew.z+PNew.radius > z_max) continue;

                // Check if surrounding shielded against planting
                auto check_shield = [&Shielded](long int ii,long int jj,long int kk, double rad){return (Shielded(ii,jj,kk) == true) ? true : false;};
                int  IsShielded = Initializations::loop_sphere_with_exit(Shielded, check_shield, PNew.x, PNew.y, PNew.z, radius + MinDistance/2.0, locSettings, BC);
                #ifdef MPI_PARALLEL
                int tmpIsShielded = IsShielded;
                MPI_Allreduce(&tmpIsShielded, &(IsShielded), 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                #endif
                if (IsShielded) continue;

                // Shield surrounding against places of other solid particles
                Initializations::loop_sphere(Shielded, set_shield, PNew.x, PNew.y, PNew.z, radius + MinDistance/2.0, locSettings, BC);

                // Initialize spherical phase field
                Planted.push_back(PNew);
                int phaseIdx2 = PhaseIndex(PNew.x, PNew.y, PNew.z);
                #ifdef MPI_PARALLEL
                MPI_Bcast(&(phaseIdx2), 1, MPI_INT, 0, MPI_COMM_WORLD);
                #endif
                const size_t idx = Initializations::Sphere(Phase, phaseIdx2, radius, PNew.x, PNew.y, PNew.z, BC, locSettings, false);
                SpheresIdx.push_back(idx);

                Info::Write("Planted Solid number ", Planted.size());
            }
            ++i;
        }
    }

    // Finalise initialization
    if(locSettings.Resolution == Resolutions::Double)
    {
        Phase.Refine();
    }
    Phase.Finalize(BC);
    return SpheresIdx;
}

void RemoveParrentGrain(PhaseField& Phase, int ParrentPhaseFieldIndex,
                        BoundaryConditions& BC, Settings& locSettings)
{
    /*
     * Removes phase number ParrentPhaseFieldIndex and increases phase number 0
     * by the removed amount.
     */
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        double ParrentPhaseAmount = Phase.Fields(i, j, k)[ParrentPhaseFieldIndex];

        if(ParrentPhaseAmount != 0.0)
        {
            Phase.Fields(i, j, k).add_value(0, ParrentPhaseAmount);
            Phase.Fields(i, j, k).set_value(ParrentPhaseFieldIndex, 0);
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }
}

std::vector<size_t> Initializations::ThermalGrooving(PhaseField& Phase,
                                                     size_t alpha, size_t beta,
                                                     BoundaryConditions& BC,
                                                     Settings& locSettings)
{
    size_t Nx = locSettings.Nx;
    size_t Nz = locSettings.Nz;

    int x = Nx/2.0;
    int y = 0;
    int z = Nz+1;
    int dx = 3.0*Nx/4.0+1;
    int dy = 0;
    int dz = Nz/2.0;

    vector<size_t> locIndex = Fractional(Phase,beta,beta,dz,BC,locSettings);
    locIndex.push_back(Rectangular(Phase,alpha,x,y,z,dx,dy,dz,BC,locSettings));

    Phase.SetBoundaryConditions(BC);
    Phase.CalculateDerivatives();
    Phase.CalculateFractions();

    return locIndex;
}

void Initializations::VoronoiTesselation(PhaseField& Phase,
                      BoundaryConditions& BC, const Settings& locSettings,
                      const size_t Ngrains, const size_t GrainsPhase)
{
    Info::WriteLineInsert("Voronoi creator");
    Info::WriteStandard("Number of grains", std::to_string(Ngrains));
    Info::WriteStandard("Thermodynamic phase", std::to_string(GrainsPhase));
    Info::WriteLineInsert("Generating grain seeds");

    Phase.Clear();

    for(size_t n = 0; n < Ngrains; n++)
    {
        Phase.AddGrainInfo(GrainsPhase);
    }

    int TotalNx = locSettings.TotalNx;
    int TotalNy = locSettings.TotalNy;
    int TotalNz = locSettings.TotalNz;

    int Nx = locSettings.Nx;
    int Ny = locSettings.Ny;
    int Nz = locSettings.Nz;

    size_t SeedX = 253;
    size_t SeedY = 4958;
    size_t SeedZ = 54861;

    mt19937_64 xPosGenerator(SeedX);
    mt19937_64 yPosGenerator(SeedY);
    mt19937_64 zPosGenerator(SeedZ);

    uniform_int_distribution <int> xPosDistribution(0, TotalNx - 1);
    uniform_int_distribution <int> yPosDistribution(0, TotalNy - 1);
    uniform_int_distribution <int> zPosDistribution(0, TotalNz - 1);

    for (size_t n = 0; n < Phase.FieldsStatistics.size(); ++n)
    {
        Phase.FieldsStatistics[n].Rcm[0] = xPosDistribution(xPosGenerator);
        Phase.FieldsStatistics[n].Rcm[1] = yPosDistribution(yPosGenerator);
        Phase.FieldsStatistics[n].Rcm[2] = zPosDistribution(zPosGenerator);
    }

    size_t NpointsTotal = TotalNx*TotalNy*TotalNz;
    size_t Npoints = 0;
    size_t NpointsRun = 0;
    int OldFractionComplete = 0;
    Info::WriteSimple("Creating coordination shells...");

    CoordinationShells CS(locSettings);
    Info::WriteSimple("Creating Voronoi grain structure...");
    Info::WriteSimple("0-------50-------100%");

    for(size_t s = 0; s < CS.Nshells(); s++)
    {
        vector <iVector3> locShell = CS[s];
        for (size_t n = 0; n < Phase.FieldsStatistics.size(); ++n)
        {
            int Rx = Phase.FieldsStatistics[n].Rcm[0];
            int Ry = Phase.FieldsStatistics[n].Rcm[1];
            int Rz = Phase.FieldsStatistics[n].Rcm[2];

            for(size_t m = 0; m < CS[s].size(); m++)
            {
                int x = Rx + locShell[m][0];
                int y = Ry + locShell[m][1];
                int z = Rz + locShell[m][2];
#ifdef MPI_PARALLEL
                if(BC.BC0X == BoundaryConditionTypes::Periodic or BC.MPIperiodicX) x = (x + TotalNx) % TotalNx;
                if(BC.BC0Y == BoundaryConditionTypes::Periodic or BC.MPIperiodicY) y = (y + TotalNy) % TotalNy;
                if(BC.BC0Z == BoundaryConditionTypes::Periodic or BC.MPIperiodicZ) z = (z + TotalNz) % TotalNz;

                x -= locSettings.OffsetX;
                y -= locSettings.OffsetY;
                z -= locSettings.OffsetZ;
#else
                if(BC.BC0X == BoundaryConditionTypes::Periodic) x = (x + Nx) % Nx;
                if(BC.BC0Y == BoundaryConditionTypes::Periodic) y = (y + Ny) % Ny;
                if(BC.BC0Z == BoundaryConditionTypes::Periodic) z = (z + Nz) % Nz;
#endif
                if(x >= 0 and x < Nx and
                   y >= 0 and y < Ny and
                   z >= 0 and z < Nz and
                   Phase.Fields(x,y,z).size() == 0)
                {
                    Phase.Fields(x,y,z).set_value(n, 1.0);
                    Npoints ++;
                }
            }
            NpointsRun = Npoints;
#ifdef MPI_PARALLEL
            MPI_Allreduce(&Npoints, &NpointsRun, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
            if(MPI_RANK == 0)
            {
#endif
            int FractionComplete = ((NpointsRun + 5)*100)/(5*NpointsTotal);

            if(FractionComplete > OldFractionComplete)
            {
                OldFractionComplete = FractionComplete;
                std::cout << "+" << flush;
            }
#ifdef MPI_PARALLEL
            }
#endif
            if(NpointsRun == NpointsTotal) break;
        }
        if(NpointsRun == NpointsTotal) break;
    }
    Phase.SetBoundaryConditions(BC);

#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
#endif
    std::cout << endl;

    // Spread interfaces by overlapping the grains near the interface
    Info::WriteSimple("Creating grain boundaries...");
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    if(Phase.Fields(i,j,k).size() == 1)
    {
        size_t alpha_index = Phase.Fields(i,j,k).begin()->index;

        for(int di = -1; di <= 1; di++)
        for(int dj = -1; dj <= 1; dj++)
        for(int dk = -1; dk <= 1; dk++)
        if(i+di >= 0 and i+di < Nx and
           j+dj >= 0 and j+dj < Ny and
           k+dk >= 0 and k+dk < Nz and
           (di != 0 or dj != 0 or dk != 0))
        {
            bool grain_not_present = true;
            for(auto beta = Phase.Fields(i+di,j+dj,k+dk).begin();
                     beta < Phase.Fields(i+di,j+dj,k+dk).end(); ++beta)
            if(beta->index == alpha_index)
            {
                grain_not_present = false;
                break;
            }
            else
            {
                break;
            }
            if(grain_not_present)
            {
                for(auto beta = Phase.Fields(i+di,j+dj,k+dk).begin();
                         beta < Phase.Fields(i+di,j+dj,k+dk).end(); ++beta)
                {
                    Phase.Fields(i,j,k).add_value(beta->index, 1.0);
                }
                Phase.Fields(i+di,j+dj,k+dk).add_value(alpha_index, 1.0);
            }
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }

    Info::WriteSimple("Generating random orientations...");

    mt19937_64 OrientGenerator1(45);
    mt19937_64 OrientGenerator2(697);
    mt19937_64 OrientGenerator3(255);

    uniform_int_distribution <int> Q1Distribution(0, 360);
    uniform_int_distribution <int> Q2Distribution(0, 180);
    uniform_int_distribution <int> Q3Distribution(0, 360);

    for(size_t n = 0; n < Phase.FieldsStatistics.size(); n++)
    {
        double Q1 = Q1Distribution(OrientGenerator1) * Pi/180.0;
        double Q2 = Q2Distribution(OrientGenerator2) * Pi/180.0;
        double Q3 = Q3Distribution(OrientGenerator3) * Pi/180.0;
        EulerAngles locAngles({Q1, Q2, Q3}, XYZ);
        Phase.FieldsStatistics[n].Orientation = locAngles.getQuaternion();
    }

    Info::WriteSimple("Done!");
}

void Initializations::ReadSubset(PhaseField& Phase, string PFFileName, string GSFileName,
                                             int offsetX,
                                             int offsetY,
                                             int offsetZ,
                                             const BoundaryConditions& BC,
                                             int inpNx,
                                             int inpNy,
                                             int inpNz)
{
    ReadSubsetGeneral(Phase, PFFileName, GSFileName, offsetX, offsetY, offsetZ,
            inpNx-offsetX, inpNy-offsetY, inpNz-offsetZ, inpNx, inpNy, inpNz,
            0, 0, 0, {false}, BC);
}

void Initializations::ReadSubsetGeneral(PhaseField& Phase, string PFFileName, string GSFileName,
                                             int offsetInpX,
                                             int offsetInpY,
                                             int offsetInpZ,
                                             int sizeInpX,
                                             int sizeInpY,
                                             int sizeInpZ,
                                             int totalSizeInpNx,
                                             int totalSizeInpNy,
                                             int totalSizeInpNz,
                                             int offsetLocalX,
                                             int offsetLocalY,
                                             int offsetLocalZ,
                                             std::initializer_list<bool> newGrain,
                                             const BoundaryConditions& BC)
{
    fstream inp(PFFileName.c_str(), ios::in | ios::binary);

    GrainInfo FieldsStatisticsInp;
    FieldsStatisticsInp.Read(GSFileName);
    int FSInpSize = FieldsStatisticsInp.size();
    int allocateNewGrains = Phase.FieldsStatistics.size();
    vector<int> targetPF(FSInpSize);                                            //Target phase field for every phase field in the source file

    for (int i = 0; i < FSInpSize; i++)
    {
        if ((unsigned)i > newGrain.size() -1)
        {
            Phase.FieldsStatistics.add_grain(FieldsStatisticsInp[i].Phase);
            targetPF[i] = Phase.FieldsStatistics.size()-1;                            //add the PF index of the current Grain to target PF
            allocateNewGrains++;
        }
        else if (newGrain.begin()[i])
        {
            Phase.FieldsStatistics.add_grain(FieldsStatisticsInp[i].Phase);
            targetPF[i] = Phase.FieldsStatistics.size()-1;                            //add the PF index of the current Grain to target PF
            allocateNewGrains++;
        }
        else
        {
            if (!Phase.FieldsStatistics[i].Exist or
                 Phase.FieldsStatistics[i].Phase != FieldsStatisticsInp[i].Phase)
            {
                stringstream message;
                message << "Target phase-field does not exist!\n"
                        << "Grain No.: " << i << "\n";
                Info::WriteExit(message.str(),
                        thisclassname, "ReadSubsetGeneral()");
                exit(1);
            }
            targetPF[i] = i;                                                    //do not change the PF index
        }
    }

    if (!inp)
    {
        Info::WriteExit(PFFileName + " could not be opened",
                thisclassname, "ReadSubsetGeneral()");
        exit(1);
    };

    int locNx = Phase.Nx;
    int locNy = Phase.Ny;
    int locNz = Phase.Nz;

    inp.read(reinterpret_cast<char*>(&locNx), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNy), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNz), sizeof(int));

    if(Phase.Nx - offsetLocalX < sizeInpX or
       Phase.Ny - offsetLocalY < sizeInpY or
       Phase.Nz - offsetLocalZ < sizeInpZ)
    {
        stringstream message;
        message << "Inconsistent system dimensions!\n"
                << "Input data dimensions: (" << locNx << ", "
                << locNy << ", " << locNz << ") grid points.\n"
                << "Local data dimensions: ("
                << Phase.Nx << ", " << Phase.Ny << ", " << Phase.Nz
                << ") grid points.\n"
                << "Local data offset relative to the input data: ("
                << offsetInpX << ", " << offsetInpY << ", " << offsetInpZ
                << ") grid points.\n";
        Info::WriteExit(message.str(), thisclassname, "ReadSubsetGeneral()");
        exit(1);
    }

    for(int i = 0; i < locNx; i++)
    for(int j = 0; j < locNy; j++)
    for(int k = 0; k < locNz; k++)
    {
        NodePF locPF;
        int   num = 0;
        int   idx = 0;
        double val = 0.0;

        inp.read(reinterpret_cast<char*>(&num), sizeof(int));

        for(int n = 0; n < num; n++)
        {
            inp.read(reinterpret_cast<char*>(&idx), sizeof(int));
            inp.read(reinterpret_cast<char*>(&val), sizeof(double));
            locPF.set_value(targetPF[idx], val);
        }

        if((i >= offsetInpX and j >= offsetInpY and k >= offsetInpZ) and
           (i <= offsetInpX + sizeInpX and
            j <= offsetInpY + sizeInpY and
            k <= offsetInpZ + sizeInpZ))
        {
            int x = i - offsetInpX + offsetLocalX;
            int y = j - offsetInpY + offsetLocalY;
            int z = k - offsetInpZ + offsetLocalZ;

            Phase.Fields(x,y,z) = locPF;
            Phase.Fields(x,y,z).flag = Phase.Fields(x,y,z).finalize();
        }
    }

    inp.close();

    Phase.FinalizeSR(BC);

    if(Phase.Resolution == Resolutions::Double)
    {
        Phase.Refine();
        Phase.FinalizeDR(BC);
    }
    Info::WriteStandard(thisclassname, "Binary input loaded");
}

}// namespace openphase

