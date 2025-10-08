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
 *   Main contributors :   Oleg Shchyglo; Alexander Monas; Marvin Tegeler;
 *                         Matthias Stratmann
 *
 */

#ifndef NUCLEATION_H
#define NUCLEATION_H

#include "Base/Includes.h"
#include "H5Interface.h"

namespace openphase
{
class PhaseField;
class Settings;
class InterfaceProperties;
class DrivingForce;
class BoundaryConditions;
class Orientations;
class Temperature;
class SymmetryVariants;
class H5Interface;
class ThermodynamicProperties;


enum class NucleiGenerationModes
{
    Static,                                                                     ///< Nucleation seeds are generated once and should be regenerated if new seed positions are required
    Dynamic                                                                     ///< Nucleation seeds are generated on the fly for each nucleation event
};

enum class NucleiOrientationModes                                               ///< Possible orientation modes for new grains
{
    Reference,                                                                  ///< Use frame of reference orientation
    Parent,                                                                     ///< Use parent grain orientation
    Random,                                                                     ///< Generate random orientation
};

enum class NucleiLocationModes                                                  ///< Possible location modes for new grains
{
    Bulk,                                                                       ///< In the bulk of the parent phase
    BulkAndGB,                                                                  ///< In the bulk and grain boundaries of the parent phase
    GB,                                                                         ///< In the entire region of the grain boundaries of the parent phase
    Junctions,                                                                  ///< In the triple and higher order junctions containing 100% matrix phase
    Interfaces,                                                                 ///< In the interfaces of the matrix phase with any other phase
    XBottom,                                                                    ///< At the bottom of x-axis
    XTop,                                                                       ///< At the top of x-axis
    YBottom,                                                                    ///< At the bottom of y-axis
    YTop,                                                                       ///< At the top of y-axis
    ZBottom,                                                                    ///< At the bottom of z-axis
    ZTop                                                                        ///< At the top of z-axis
};

enum class NucleiSizeDistributions                                              ///< Available seed size distributions
{
    Normal,
    Cauchy,
    Uniform,
    None                                                                        ///< Fixed size seeds
};

class OP_EXPORTS Nucleation : public OPObject                                   ///< Handles the nucleation of new phases and grains
{
  public:

    Nucleation(){};
    Nucleation(Settings& locSettings, const std::string InputFileName = DefaultInputFileName);
    void Initialize(Settings& locSettings) override;                            ///< Initializes storages, sets internal variables.
    void ReadInput(const std::string InputFileName) override;                   ///< Reads input values from file
	void ReadInput(std::stringstream& inp) override;                            ///< Reads input values from file
    void GenerateNucleationSites(PhaseField& Phase, Temperature& Tx);           ///< Randomly generates nucleation sites and gives them weights according to the chosen distribution parameters
    void ReGenerateNucleationSites(PhaseField& Phase, Temperature& Tx);         ///< Clears existing seeds and generates new nucleation sites
    void Clear();                                                               ///< Clears the particles storage
    void GenerateRandomSeeds(void);                                             ///< Generates new seeds for random number generators

    void PlantNuclei(PhaseField& Phi, int tstep);                               ///< Plants generated nuclei according to their nucleation parameters
    void CheckNuclei(PhaseField& Phi, InterfaceProperties& IP, DrivingForce& dG, int tstep);///< Checks planted nuclei, removes unstable ones

    void WriteStatistics(int tstep, long int PFindex, size_t NucleatingPhase,
                         size_t MatrixPhase, size_t Variant,
                         int x, int y, int z, Quaternion Q, double dGnuc,
                         double dGmin, std::string status) const;               ///< Writes nulceation statistics to a file
    void Read(int tstep);                                                       ///< Read stored nucleation information from a file
    void Write(int tstep);                                                      ///< Write nucleation information to a file

    void ReadH5(int tstep, H5Interface& H5);                                    ///< Read stored nucleation information from a file
    void WriteH5(int tstep, H5Interface& H5);                                   ///< Write nucleation information to a file
    size_t Nphases;                                                             ///< Number of thermodynamic phases
    std::vector<size_t> Nvariants;                                              ///< Number of crystallographic (symmetry/translation/...) variants

    double iWidth;                                                              ///< Interface width in grid points
    double dx;                                                                  ///< Grid spacing

    int TotalNx;                                                                ///< X dimension of the system in MPI parallel mode
    int OffsetX;                                                                ///< X dimension offset of the current domain in MPI parallel mode
    int TotalNy;                                                                ///< Y dimension of the system in MPI parallel mode
    int OffsetY;                                                                ///< Y dimension offset of the current domain in MPI parallel mode
    int TotalNz;                                                                ///< Z dimension of the system in MPI parallel mode
    int OffsetZ;                                                                ///< Z dimension offset of the current domain in MPI parallel mode

    int Nx;                                                                     ///< X dimension of the system
    int Ny;                                                                     ///< Y dimension of the system
    int Nz;                                                                     ///< Z dimension of the system

    int dNx;                                                                    ///< Active X dimension of the system
    int dNy;                                                                    ///< Active Y dimension of the system
    int dNz;                                                                    ///< Active Z dimension of the system

    size_t SeedX;                                                               ///< Seed for X dimension random number generator
    size_t SeedY;                                                               ///< Seed for Y dimension random number generator
    size_t SeedZ;                                                               ///< Seed for Z dimension random number generator

    size_t SeedA;                                                               ///< Seed for orientation around X axis random number generator
    size_t SeedB;                                                               ///< Seed for orientation around Y axis random number generator
    size_t SeedC;                                                               ///< Seed for orientation around Z axis random number generator

    size_t SeedR;                                                               ///< Seed for nuclei radius random number generator


    Matrix<double> MobilityReduction;                                           ///< Reduction of mobilities for growing nuclei

    size_t SeedV;                                                               ///< Seed for variants random number generator

    std::mt19937_64 VariantsGenerator;                                          ///< Random number generator for variants selection

    int    NucleateEvery;                                                       ///< How frequently nucleation attempts should be performed
    int    NumberOfAttempts;                                                    ///< Maximum number of attempts to generate each nucleation site

    std::string VTKDir;                                                         ///< Directory-path added in front of VTK files
    std::string RawDataDir;                                                     ///< Directory-path added in front of Restart files
    std::string TextDir;                                                        ///< Directory-path for statistic output in ASCII format

    // Parameters used to set particles size distribution for heterogeneous nucleation.
    // At the moment only normal distribution is implemented:
    // p = 1.0/(DistSigma*sqrt(2*Pi))*exp{-(x - DistMu)^2/(2.0*DistSigma^2)}

    struct NucleationParameters                                                 ///< Pairwise nucleation settings ("phase alpha" in "phase beta")
    {
        bool   Allowed;                                                         ///< True if nucleation is allowed
        double Tmin;                                                            ///< Upper temperature limit for nucleation
        double Tmax;                                                            ///< Lower temperature limit for nucleation
        bool   Generated;                                                       ///< True if particles are already generated
        size_t Nucleated;                                                       ///< Number of already nucleated seeds
        double DistSigma;                                                       ///< Normal nuclei radius distribution standard deviation
        double DistMu;                                                          ///< Normal nuclei radius distribution mean value
        double Density;                                                         ///< Average nucleation density
        size_t Nsites;                                                          ///< Number of generated nucleation sites
        size_t part;
        NucleiLocationModes LocationMode;                                       ///< Nuclei location mode
        NucleiOrientationModes OrientationMode;                                 ///< Nuclei orientation mode
        NucleiSizeDistributions Distribution;                                   ///< Using normal distribution if "true" or user specified number of seeds otherwise
        size_t Nvariants;                                                       ///< Number of permitted crystallographic variants per nucleation site
        double Shielding;                                                       ///< Shielding radius for the nuclei to avoid numerical artifacts
        bool   RelativeDensity;                                                 ///< Relative nucleation density if "true" or absolute otherwise
        double RadiusMIN;                                                       ///< Minimum nuclei radius
        double RadiusMAX;                                                       ///< Maximum nuclei radius
    };

    Matrix< NucleationParameters > Parameters;                                  ///< Stores nucleation parameters for each pair of phases (A in B)

    struct NucSite                                                              ///< Particle properties
    {
        std::vector<size_t> PFindices;                                          ///< Phase-field indices after nucleation (considering multiple symmetry variants)
        /// position indexes
        int x;
        int y;
        int z;
        /// orientation
        Quaternion Q;
        /// particle radius
        double radius;
        bool is_active()
        {
            return (not planted);
        }
        bool planted;                                                           ///< True if nucleus is planted successfully (and growing)
        bool growing;                                                           ///< True if nucleus is growing
        size_t time_stamp;                                                      ///< Time step of nucleation event
    };

    std::vector<iVector3> GetNucleationSites(PhaseField& Phi);                  ///< Vector of coordinates of all active nucleation sites
    std::vector<std::vector<bool> > GetNucleationEvents(void);
    bool IsShielded(std::vector<NucSite> LocGrainStorage, const int i, const int j, const int k,
            const double shielding) const;                                      ///< Returns true if member of the "NucleatedGrains" is inside the "shielding"-radius, false otherwise.

    Matrix <std::vector <NucSite> > GeneratedParticles;                         ///< Generated particles storage
    Matrix <std::vector <NucSite> > NucleatedParticles;                         ///< Nucleated particles storage

  private:
    void SetNucleationSites(PhaseField& Phase, Temperature& Tx, size_t n, size_t m);
};

}// namespace openphase

#endif // NUCLEATION_H
