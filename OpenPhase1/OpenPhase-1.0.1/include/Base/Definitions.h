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
 *   File created :   2022
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *                         Raphael Schiedung
 *
 */

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <string>
#include <complex>

#ifdef SERIAL
    // overloading several OpenMP methods used in OpenPhase for serial operation
    inline int omp_get_num_threads()
    {
        return 1;
    }
    inline int omp_get_max_threads()
    {
        return 1;
    }
    inline int omp_get_thread_num()
    {
        return 0;
    }
    inline void omp_set_num_threads(int)
    {
        // Do nothing
    };
#endif

#ifdef _WIN32
#define NOMINMAX
#define and &&
#define or ||
#define M_PI 3.14159265359
#endif

#ifdef _WIN32
#define OP_EXPORTS __declspec(dllexport)
#define OP_IMPORTS __declspec(dllimport)
#ifdef OPENPHASE_EXPORTS
#define OPENPHASE_API __declspec(dllexport)
#else
#define OPENPHASE_API __declspec(dllimport)
#endif
#else
#define OP_EXPORTS
#define OP_IMPORTS
#ifdef OPENPHASE_EXPORTS
#define OPENPHASE_API
#else
#define OPENPHASE_API
#endif
#endif

#ifdef MPI_PARALLEL
    extern int MPI_RANK;                                                        ///< Local MPI RANK in MPI parallel mode
    extern int MPI_SIZE;                                                        ///< Total number of MPI RANKs in MPI parallel mode
    extern int MPI_CART_RANK[3];                                                ///< Cartesian coordinates of the current process if using MPI 3D domain decomposition
    extern int MPI_CART_SIZE[3];                                                ///< Number of processes used in each direction if using MPI 3D domain decomposition
    extern bool MPI_3D_DECOMPOSITION;                                           ///< "true" if MPI should decompose in 3 dimensions
#endif

namespace openphase
{

inline void simulation_end()
{
#ifdef _WIN32
    std::getchar();                                                             ///< Prevents closing of the terminal window at the end of the simulation
#endif
}

#ifdef _WIN32
    const std::string dirSeparator = "\\";                                      ///< Windows style directory separator
#else
    const std::string dirSeparator = "/";                                       ///< Unix/Linux style directory separator
#endif

static constexpr double Pi = 3.14159265358979323846;                            ///< Pi constant value

const std::complex< double > I(0.0, 1.0);                                       ///< sqrt(-1) declaration

// Definitions of special keywords
enum class Resolutions                                                          ///< Possible simulation resolutions
{
    Single,                                                                     ///< Single resolution
    Double                                                                      ///< Double resolution
};

enum class AdvectionSchemes                                                     ///< Available advection schemes
{
    Upwind,
    Minmod,
//    VanLeer,
    Superbee,
    //LaxWendroff,
    MonotonizedCentral
};

enum class AggregateStates                                                      ///< Aggregate state markers
{
    Solid,
    Liquid,
    Gas
};

enum class EventTriggers                                                        ///< Events triggers/conditions
{
    User,                                                                       ///< The event is manually controlled by the user
    Tmax,                                                                       ///< Trigger at a given maximum temperature approaching from below [K]
    Tmin,                                                                       ///< Trigger at a given minimum temperature approaching from above [K]
    Time,                                                                       ///< Trigger at a particular moment in time [s]
    TimeStep,                                                                   ///< Trigger at a particular time step
    PhaseFractionMax,                                                           ///< Trigger on particular phase fraction value approaching from below
    PhaseFractionMin                                                            ///< Trigger on particular phase fraction value approaching from above
};

enum class ClearingModes                                                        ///< Convenience feature to control major Clear() methods calls in various classes
{
    Automatic,                                                                  ///< Relevant storages are cleared automatically at the end of the time step
    Manual                                                                      ///< To clear the relevant storages Clear() method should be called explicitly
};

class InterfaceTypes                                                            ///< Types of heterogeneous interfaces (between different aggregate states)
{
 public:
    static bool SolidSolid(const AggregateStates InA, const AggregateStates InB)
    {
        return (InA == AggregateStates::Solid && InB == AggregateStates::Solid);
    }
    static bool SolidLiquid(const AggregateStates InA, const AggregateStates InB)
    {
        return (InA == AggregateStates::Solid  && InB == AggregateStates::Liquid) or
               (InA == AggregateStates::Liquid && InB == AggregateStates::Solid);
    }
    static bool LiquidLiquid(const AggregateStates InA, const AggregateStates InB)
    {
        return (InA == AggregateStates::Liquid && InB == AggregateStates::Liquid);
    }
    static bool SolidGas(const AggregateStates InA, const AggregateStates InB)
    {
        return (InA == AggregateStates::Solid && InB == AggregateStates::Gas) or
               (InA == AggregateStates::Gas   && InB == AggregateStates::Solid);
    }
    static bool LiquidGas(const AggregateStates InA, const AggregateStates InB)
    {
        return (InA == AggregateStates::Liquid && InB == AggregateStates::Gas) or
               (InA == AggregateStates::Gas    && InB == AggregateStates::Liquid);
    }
    static bool GasGas(const AggregateStates InA, const AggregateStates InB)
    {
        return (InA == AggregateStates::Gas && InB == AggregateStates::Gas);
    }
};

template<typename T>
inline void ignore_result(T /* unused result */) {}                             ///< Allows to suppress warnings on unused function return values

// Standard OpenPhase files and directories
const std::string DefaultInputFileName = "ProjectInput.opi";
const std::string DefaultVTKDir        = "VTK" + dirSeparator;
const std::string DefaultRawDataDir    = "RawData" + dirSeparator;
const std::string DefaultTextDir       = "TextData" + dirSeparator;

}// namespace openphase
#endif
