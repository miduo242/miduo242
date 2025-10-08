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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */

#ifndef SETTINGS_H
#define SETTINGS_H

#include "Base/Definitions.h"
#include "Base/OPObject.h"

namespace openphase
{
class BoundaryConditions;

class OP_EXPORTS Settings                                                       ///< System settings module. Reads and stores system settings
{
 public:
    std::string thisclassname;                                                  ///< Object's implementation class name
    std::string thisobjectname;                                                 ///< Object's name

    int Nx;                                                                     ///< Actual number of grid points in X direction
    int Ny;                                                                     ///< Actual number of grid points in Y direction
    int Nz;                                                                     ///< Actual number of grid points in Z direction

    int dNx;                                                                    ///< Stencil step in X direction for reduced dimensions simulations
    int dNy;                                                                    ///< Stencil step in Y direction for reduced dimensions simulations
    int dNz;                                                                    ///< Stencil step in Z direction for reduced dimensions simulations

    int newNx;                                                                  ///< New number of grid points in X direction
    int newNy;                                                                  ///< New number of grid points in Y direction
    int newNz;                                                                  ///< New number of grid points in Z direction

    int maxNx;                                                                  ///< Max number of grid points in X direction
    int maxNy;                                                                  ///< Max number of grid points in Y direction
    int maxNz;                                                                  ///< Max number of grid points in Z direction

    int TotalNx;                                                                ///< Total system size in X direction (in case of MPI parallelism)
    int OffsetX;                                                                ///< X coordinate of the local coordinate system origin (in case of MPI parallelism)
    int TotalNy;                                                                ///< Total system size in Y direction (in case of MPI parallelism)
    int OffsetY;                                                                ///< Y coordinate of the local coordinate system origin (in case of MPI parallelism)
    int TotalNz;                                                                ///< Total system size in Z direction (in case of MPI parallelism)
    int OffsetZ;                                                                ///< Z coordinate of the local coordinate system origin (in case of MPI parallelism)

    double dx;                                                                  ///< Grid spacing in true units
    double dx_2;                                                                ///< Grid spacing in true units in double resolution mode
    double iWidth;                                                              ///< Interface width in grid cells
    double Eta;                                                                 ///< Interface width in true units

    bool ConsiderNucleusVolume;                                                 ///< Enables weighted number of local phase-fields to support nucleation

    size_t Nphases;                                                             ///< Number of thermodynamic phases
    size_t Ncomp;                                                               ///< Number of chemical elements
    size_t Bcells;                                                              ///< Number of boundary cells for setting the boundary conditions
    std::vector<size_t> Nvariants;                                              ///< Number of crystallographic (symmetry/translation/...) variants

    Resolutions Resolution;                                                     ///< Phase field resolution (Single or Double)

    LaplacianStencils DiffusionStencil;                                         ///< Diffusion stencil selector
    LaplacianStencils PhaseFieldLaplacianStencil;                               ///< Phase-field Laplacian stencil selector
    GradientStencils  PhaseFieldGradientStencil;                                ///< Phase-field gradient stencil selector

    int  ActiveDimensions() const;                                              ///< Indicates dimensionality of the simulation, e.g. 1D, 2D or 3D

    bool RemeshingAllowed;                                                      ///< Remeshing allowed (Yes or No)
    bool initialized = false;                                                   ///< Indicates if the settings have been initialized

    Settings()                                                                  ///< Default constructor
    {
        Initialize();
    }

    Settings(const std::string InputFileName)                                   ///< Constructor
    {
        Initialize();
        ReadInput(InputFileName);
    }
    void Initialize(void);                                                      ///< Initializes the class
    void ReadInput(const std::string InputFileName = DefaultInputFileName);     ///< Reads input from the specified input file
    void ReadInput(std::stringstream& data);                                    ///< Reads input from the specified input file

    void Resize(int size_x, int size_y, int size_z, int tStep,
                const BoundaryConditions& BC);                                  ///< Resizes the system to new dimensions
    void ReadDimensions(int tStep);                                             ///< Reads system dimensions history from the file
    void WriteDimensions(int tStep);                                            ///< Writes system dimensions history to the file
    void AddForRemesh(OPObject& obj);                                           ///< Adds object to be remeshed to the ObjectsToRemesh
    void RemeshAll(const BoundaryConditions& BC);                               ///< Calls Remesh() on all objects in ObjectsToRemesh
    void Remesh(std::string ObjectName, int size_x, int size_y, int size_z,
                const BoundaryConditions& BC);                                  ///< Calls Remesh() on object(s) with the given name base
    void AddForAdvection(OPObject& obj);                                        ///< Adds object to be advected to the ObjectsToAdvect
    void AdvectAll(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi,
                   const BoundaryConditions& BC, RunTimeControl& RTC);          ///< Calls Advect() on all objects in ObjectsToAdvect
    void Advect(std::string ObjectName, AdvectionHR& Adv, const Velocities& Vel,
                PhaseField& Phi, const BoundaryConditions& BC,
                RunTimeControl& RTC);                                           ///< Calls Advect() on object(s) with the given name base

    std::vector<std::array<int, 4>> DimensionsHistory;                          ///< Stores dimensions evolution history over the simulation time. Syntax: {{tStep1, Nx1, Ny1, Nz1}, ... ,{tStepN, NxN, NyN, NzN}}
    std::vector<OPObject*> ObjectsToRemesh;                                     ///< Stores pointers to objects sensitive to remeshing
    std::vector<std::string> NamesOfObjectsToRemesh;                            ///< Stores the names of objects sensitive to remeshing
    std::vector<OPObject*> ObjectsToAdvect;                                     ///< Stores pointers to objects sensitive to advection
    std::vector<std::string> NamesOfObjectsToAdvect;                            ///< Stores the names of objects sensitive to advection

    std::vector<std::string> PhaseNames;                                        ///< Names of thermodynamic phases
    std::vector<std::string> ElementNames;                                      ///< Names of chemical elements
    std::vector<AggregateStates> PhaseAggregateStates;                          ///< Aggregate states of all phases

    std::string VTKDir;                                                         ///< Directory name for the VTK files
    std::string RawDataDir;                                                     ///< Directory name for the raw data files
    std::string TextDir;                                                        ///< Directory name for the text files

    std::string move_frame_direction;
    size_t move_frame_phase;
    size_t move_frame_pos;

    Settings& operator= (const Settings& rhs);                                  ///< Copy operator

 protected:
 private:
#ifdef MPI_PARALLEL
    void Setup_MPI();                                                           ///< Sets up MPI parallel environment in 1D decomposition mode
    void Setup_MPI3D();                                                         ///< Sets up MPI parallel environment in 3D decomposition mode
#endif
};

}// namespace openphase

#endif
