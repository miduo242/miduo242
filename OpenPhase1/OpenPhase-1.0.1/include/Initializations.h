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

#ifndef INITIALIZATIONS_H
#define INITIALIZATIONS_H

#include "Base/CommonFunctions.h"
#include "Base/Includes.h"
#include "BoundaryConditions.h"
#include "Settings.h"
#include <cassert>
#include <vector>

namespace openphase
{

class BoundaryConditions;
class ElasticProperties;
class Orientations;
class PhaseField;
class Settings;
class ChemicalProperties;

class OP_EXPORTS Initializations : public OPObject
{
 public:
    static std::vector<iVector3> QuasiRandomNuclei(PhaseField& Phase,
            Settings& locSettings, size_t phaseIndex, int dist, int seed = -1,
            int offset = 2);
    static void RandomNuclei(PhaseField& Phase, Settings& locSettings,
            size_t phaseIndex, size_t Nparticles, int seed = -1);
    static size_t Ellipsoid(PhaseField& Phase, size_t PhaseIndex,
            double RadiusX, double RadiusY, double RadiusZ,
            double x0, double y0, double z0,
            BoundaryConditions& BC, Settings& locSettings);                     ///< Initializes a spherical grain
    static size_t SectionalPlane(PhaseField& Phase, const size_t PhaseIndex,
            const dVector3 Point, const dVector3 Orientation,
            const BoundaryConditions& BC, const Settings& locSettings,
            const bool NewGrain = true, const size_t FieldIndex = 0);           ///< Initializes a new phase on the negative side of a sectional plane
    static size_t Single(PhaseField& Phase, size_t PhaseIndex,
            BoundaryConditions& BC, Settings& locSettings);                     ///< Initializes a single space filling phase-field
    static size_t Zlayer(PhaseField& Phase, size_t PhaseIndex, int Position,
            int Thickness, BoundaryConditions& BC, Settings& locSettings);      ///< ?? TODO
    static size_t Sphere(PhaseField& Phase, const size_t PhaseIndex,
            const double Radius, const double x0, const double y0,
            const double z0, const BoundaryConditions& BC,
            const Settings& locSettings, const bool Finalize = true);           ///< Initializes a spherical grain
    static std::vector<size_t> Fractional(PhaseField& Phase,
            size_t MajorityPhaseIndex, size_t MinorityPhaseIndex,
            double MinorityPhaseLayerThickness, BoundaryConditions& BC,
            Settings& locSettings);                                             ///< ?? TODO
    static std::vector<size_t> ThreeFractionals(PhaseField& Phase,
            size_t MajorityPhaseIndex, double MajorityPhaseLayerThickness,
            size_t MinorityPhaseIndex1, double MinorityPhaseLayerThickness1,
            size_t MinorityPhaseIndex2, BoundaryConditions& BC,
            Settings& locSettings);                                             ///< ?? TODO
    static std::vector<size_t> TwoDifferentWalls(PhaseField& Phase,
            size_t ChannelPhaseIndex, size_t WallsPhaseIndex,
            double WallsThickness, BoundaryConditions& BC,
            Settings& locSettings);
    static std::vector<size_t> TwoWalls(PhaseField& Phase,
            size_t ChannelPhaseIndex, size_t WallsPhaseIndex,
            double WallsThickness, BoundaryConditions& BC,
            Settings& locSettings);
    static size_t Rectangular(PhaseField& Phase, const size_t PhaseIndex,
            const double Lx, const double Ly, const double Lz,
            const double x0, const double y0, const double z0,
            const BoundaryConditions& BC, const Settings& locSettings,
            const bool Finalize = true);                                        ///< Initializes a rectangular grain
    static size_t CylinderSimple(PhaseField& Phase, const size_t PhaseIndex,
            const double Radius, const double length, const int Axis,
            const double x0, const double y0, const double z0,
            const BoundaryConditions& BC, const Settings& locSettings);         ///< Initializes a cylindrical grain
    static size_t Cylinder(PhaseField& Phase, size_t PhaseFieldIndex,
            double Radius, double length, int Axis,
            double x0, double y0, double z0,
            const BoundaryConditions& BC, Settings& locSettings);               ///< Initializes a cylindrical grain
    static void RandomNucleiOnPlane(PhaseField& Phase, Settings& locSettings,
            size_t phaseIndex, size_t Nparticles, int seed, std::string axis,
            std::string position);
    // Special initialization functions:
    static size_t  FillGrainWithSpheres(PhaseField& Phase,
            size_t ParrentPhaseFieldIndex, size_t SpheresPhaseIndex,
            double MinR, double MaxR, BoundaryConditions& BC,
            Settings& locSettings, size_t Nspheres = 0,
            double MinDistance = -1);                                           ///< Fills grain "loosely" with spheres spaced MinDistance apart
    static std::vector<size_t> FillRectangularWithSpheres(
            PhaseField& Phase,
            const BoundaryConditions& BC,
            const Settings& locSettings,
            std::function<size_t(long int, long int, long int)> PhaseIndex,     ///< function that determines the phases of the grain which may depend on the position (i,j,k).
            double MeanRadius,                                                  ///< Mean radius of grains
            double StdRadius,                                                   ///< Standard deviation grain radius
            long int x_min, long int x_max,
            long int y_min, long int y_max,
            long int z_min, long int z_max,
            double MinDistance = 0.0                                            ///< Minimal distance of grains which can be negative for even denser microstructure!
            );                                                                  ///< Fills rectangular domain densely with spherical grains of normal distributed size.
    static size_t SphereInGrain(PhaseField& Phase,
            const size_t ParrentPhaseFieldIndex,
            const size_t PhaseIndex, const double Radius,
            const double x0, const double y0, const double z0,
            const BoundaryConditions& BC, const Settings& locSettings,
            const bool Finalize = true);
    static int TwoDimEBSD(std::string filename, std::vector<int> columns,
            PhaseField& Phase, BoundaryConditions& BC, Settings& locSettings);  ///< Reads microstructure from EBSD result files
    static int TwoDimEBSDWithOrientations(std::string filename,
            std::vector<int> columns, std::string anglerepresentation,
            PhaseField& Phase, ElasticProperties& EP, Orientations& OR,
            BoundaryConditions& BC, Settings& locSettings);                     ///< Reads microstructure from EBSD result files
    static void VoronoiTesselation(PhaseField& Phase, BoundaryConditions& BC,
            const Settings& OPSettings, const size_t seedpoints,
            const size_t matrixphase);                                          ///< Creates Voronoi grain structures using the algorithm based on coordination shells.
    static void TripleJunction(PhaseField& Phase, size_t PhaseIndex,
                BoundaryConditions& BC, Settings& locSettings);
    static std::vector<size_t> Young3(PhaseField& Phase, size_t alpha,
            size_t beta, size_t gamma, size_t delta, BoundaryConditions& BC,
            Settings& locSettings);
    static std::vector<size_t> Young4(PhaseField& Phase, size_t PhaseIndex,
            BoundaryConditions& BC, Settings& locSettings);
    static void Young4Periodic(PhaseField& Phase, size_t PhaseIndex,
            BoundaryConditions& BC, Settings& locSettings);
    static std::vector<size_t> ThermalGrooving(PhaseField& Phase,
            size_t PhaseIndex1, size_t PhaseIndex2, BoundaryConditions& BC,
            Settings& locSettings);
    static std::vector<iVector3> QuasiRandomSpheres(PhaseField& Phase,
            BoundaryConditions& BC, Settings& locSettings, int phaseIndex1,
            int phaseIndex2, int dist, double radius1, double radius2,
            double probabilityPhase1, int seed, int offset);

    static void Read(PhaseField& Phase, std::string FileInputName,
            BoundaryConditions& BC, Settings& locSettings);
    static void ReadSubset(PhaseField& Phase,
            std::string PFFileName,
            std::string GSFileName,
            int offsetX,
            int offsetY,
            int offsetZ,
            const BoundaryConditions& BC,
            int inpNx = 1,
            int inpNy = 1,
            int inpNz = 1);                                                     ///< Read subset of a raw (binary) phase fields from the file of possibly different (larger) size
    static void ReadSubsetGeneral(PhaseField& Phase,
            std::string PFFileName,                                             ///< Phase field input file name
            std::string GSFileName,                                             ///< Grains statistics input file name
            int offsetInpX,                                                     ///< Starting point for input reading
            int offsetInpY,
            int offsetInpZ,
            int sizeInpX,                                                       ///< size of data read from input starting from offsetInpX/Y/Z
            int sizeInpY,
            int sizeInpZ,
            int totalSizeInpNx,                                                 ///< total size of input data, for legacy format
            int totalSizeInpNy,
            int totalSizeInpNz,
            int offsetLocalX,                                                   ///< starting point for pasting to local data
            int offsetLocalY,
            int offsetLocalZ,
            std::initializer_list<bool> newGrain,                               ///< true for create new grain, false for copy to same grain in the target, if no value is specified, true is assumed
            const BoundaryConditions& BC);                                      ///< Read subset of a raw (binary) phase fields from the file of possibly different size, paste to custom location

    /// Calculates global coordinates from local mpi coordinates
    template<typename T>
    static inline void ApplyPeriodicBoundaryConditionsOnCoordinates
        (T& i, T& j, T& k, const BoundaryConditions& BC, const Settings& locSettings)
    {
#ifdef MPI_PARALLEL
        if (BC.MPIperiodicX or BC.BCNX == BoundaryConditionTypes::Periodic) while (i >= locSettings.TotalNx) i -= locSettings.TotalNx;
        if (BC.MPIperiodicY or BC.BCNY == BoundaryConditionTypes::Periodic) while (j >= locSettings.TotalNy) j -= locSettings.TotalNy;
        if (BC.MPIperiodicZ or BC.BCNZ == BoundaryConditionTypes::Periodic) while (k >= locSettings.TotalNz) k -= locSettings.TotalNz;
        if (BC.MPIperiodicX or BC.BC0X == BoundaryConditionTypes::Periodic) while (i <  0) i += locSettings.TotalNx;
        if (BC.MPIperiodicY or BC.BC0Y == BoundaryConditionTypes::Periodic) while (j <  0) j += locSettings.TotalNy;
        if (BC.MPIperiodicZ or BC.BC0Z == BoundaryConditionTypes::Periodic) while (k <  0) k += locSettings.TotalNz;
#else
        if (BC.BCNX == BoundaryConditionTypes::Periodic) while (i >= locSettings.Nx) i -= locSettings.Nx;
        if (BC.BCNY == BoundaryConditionTypes::Periodic) while (j >= locSettings.Ny) j -= locSettings.Ny;
        if (BC.BCNZ == BoundaryConditionTypes::Periodic) while (k >= locSettings.Nz) k -= locSettings.Nz;
        if (BC.BC0X == BoundaryConditionTypes::Periodic) while (i <  0) i += locSettings.Nx;
        if (BC.BC0Y == BoundaryConditionTypes::Periodic) while (j <  0) j += locSettings.Ny;
        if (BC.BC0Z == BoundaryConditionTypes::Periodic) while (k <  0) k += locSettings.Nz;
#endif
    };

    /// Checks if (i,j,k) is the
    template<typename T, typename field_t>
    static inline bool CoordinatesInBoundaries(T& i, T& j, T& k, const field_t& field, const Settings& locSettings)
    {
#ifdef MPI_PARALLEL
        if ((i < locSettings.OffsetX) or (i >= locSettings.OffsetX + field.sizeX())) return false;
        if ((j < locSettings.OffsetY) or (j >= locSettings.OffsetY + field.sizeY())) return false;
        if ((k < locSettings.OffsetZ) or (k >= locSettings.OffsetZ + field.sizeZ())) return false;
#else
        if ((i < 0) or (i >= field.sizeX())) return false;
        if ((j < 0) or (j >= field.sizeY())) return false;
        if ((k < 0) or (k >= field.sizeZ())) return false;
#endif
        return true;
    }

    /// Loops over all points (i,j,k) of a sphere with radius at (x0,y0,z0) and executes func(i,j,k,radius)
    template <class T, int Rank>
    static void loop_sphere(Storage3D<T,Rank>& field,
            const std::function<bool(long int, long int, long int, double)>& func,
            const double x0, const double y0, const double z0,
            const double radius, const Settings& locSettings, const BoundaryConditions& BC)
    {
        assert(x0 >= 0);
        assert(y0 >= 0);
        assert(z0 >= 0);
        assert(x0 < locSettings.TotalNx);
        assert(y0 < locSettings.TotalNy);
        assert(z0 < locSettings.TotalNz);
        assert(radius < std::max(locSettings.TotalNx,std::max(locSettings.TotalNy,locSettings.TotalNz)));

        dVector3 pos0 ({double(x0),double(y0),double(z0)});
        const long int iRadius = std::ceil(radius);
        const long int ii_min = (field.dNx() != 0) ? -iRadius : 0;
        const long int ii_max = (field.dNx() != 0) ?  iRadius : 0;
        const long int jj_min = (field.dNy() != 0) ? -iRadius : 0;
        const long int jj_max = (field.dNy() != 0) ?  iRadius : 0;
        const long int kk_min = (field.dNz() != 0) ? -iRadius : 0;
        const long int kk_max = (field.dNz() != 0) ?  iRadius : 0;

        #pragma omp parallel for collapse(3) //NOTE has crashes MPI calculation in older version of loop!
        for (long int ii = ii_min; ii <= ii_max; ++ii)
        for (long int jj = jj_min; jj <= jj_max; ++jj)
        for (long int kk = kk_min; kk <= kk_max; ++kk)
        {
            // Calculate global coordinates (i,j,k)
            long int i = std::round(x0 + ii);
            long int j = std::round(y0 + jj);
            long int k = std::round(z0 + kk);

            ApplyPeriodicBoundaryConditionsOnCoordinates(i,j,k,BC,locSettings);
            if(not CoordinatesInBoundaries(i,j,k,field,locSettings)) continue;

            // Calculate distance between (i,j,k) and (x0,y0,z0)
            dVector3 pos ({double(i),double(j),double(k)});
            double rad = CommonFunctions::Distance(pos, pos0, BC, locSettings).abs();
            if (rad > radius) continue;

            // Convert global coordinates to local block coordinates
#ifdef MPI_PARALLEL
            i -= locSettings.OffsetX;
            j -= locSettings.OffsetY;
            k -= locSettings.OffsetZ;
#endif
            func(i,j,k,rad);
        }
    }

    /// Loops over all points (i,j,k) of a sphere with radius at (x0,y0,z0) and executes func(i,j,k,radius)
    /// Exits loop when func(i,j,k,radius) return true
    template <class T, int Rank>
    static bool loop_sphere_with_exit(Storage3D<T,Rank>& field,
            const std::function<bool(long int, long int, long int, double)>& func,
            const double x0, const double y0, const double z0,
            const double radius, const Settings& locSettings, const BoundaryConditions& BC)
    {
        assert(x0 >= 0);
        assert(y0 >= 0);
        assert(z0 >= 0);
        assert(x0 < locSettings.TotalNx);
        assert(y0 < locSettings.TotalNy);
        assert(z0 < locSettings.TotalNz);
        assert(radius < std::max(locSettings.TotalNx,std::max(locSettings.TotalNy,locSettings.TotalNz)));

        dVector3 pos0 ({double(x0),double(y0),double(z0)});
        const long int iRadius = std::ceil(radius);
        const long int ii_min = (field.dNx() != 0) ? -iRadius : 0;
        const long int ii_max = (field.dNx() != 0) ?  iRadius : 0;
        const long int jj_min = (field.dNy() != 0) ? -iRadius : 0;
        const long int jj_max = (field.dNy() != 0) ?  iRadius : 0;
        const long int kk_min = (field.dNz() != 0) ? -iRadius : 0;
        const long int kk_max = (field.dNz() != 0) ?  iRadius : 0;

        int loc_exit_loop = 0; // is used to exit loop if a certain condition is met e.g. the presence of an existing phase
        #pragma omp parallel for shared(loc_exit_loop) collapse(3) //NOTE has crashes MPI calculation in older version of loop!
        for (long int ii = ii_min; ii <= ii_max; ++ii)
        for (long int jj = jj_min; jj <= jj_max; ++jj)
        for (long int kk = kk_min; kk <= kk_max; ++kk)
        {
            if (loc_exit_loop) continue;
            // Calculate global coordinates (i,j,k)
            long int i = std::round(x0 + ii);
            long int j = std::round(y0 + jj);
            long int k = std::round(z0 + kk);

            ApplyPeriodicBoundaryConditionsOnCoordinates(i,j,k,BC,locSettings);
            if(not CoordinatesInBoundaries(i,j,k,field,locSettings)) continue;

            // Calculate distance between (i,j,k) and (x0,y0,z0)
            dVector3 pos ({double(i),double(j),double(k)});
            double rad = CommonFunctions::Distance(pos, pos0, BC, locSettings).abs();
            if (rad > radius) continue;

            // Convert global coordinates to local block coordinates
#ifdef MPI_PARALLEL
            i -= locSettings.OffsetX;
            j -= locSettings.OffsetY;
            k -= locSettings.OffsetZ;
#endif
            if (func(i,j,k,rad)) loc_exit_loop++;
        }
        int exit_loop = 0;
#ifdef MPI_PARALLEL
        MPI_Allreduce(&loc_exit_loop, &exit_loop, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
        exit_loop += loc_exit_loop;
#endif

        if (exit_loop) return true;
        else return false;
    }

    static constexpr auto thisclassname = "Initializations";

 private:
    static void Disc(PhaseField& Phase, size_t PhaseIndex, double Radius,
            int NormalAxis, double x0, double y0, double z0,
            Settings& locSettings);                                             ///< initializes a 2D disc, only to be used within other initialization functions, does NOT call Phi.Finalize() at the end
    static void SphereFixedIdx(PhaseField& Phase, size_t PhaseIndex,
            double Radius, double x0, double y0, double z0,
            BoundaryConditions& BC, Settings& locSettings);
   };
}// namespace openphase
#endif
