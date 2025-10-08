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

#include "Nucleation.h"
#include "Settings.h"
#include "PhaseField.h"
#include "Temperature.h"
#include "DrivingForce.h"
#include "InterfaceProperties.h"
#include "Base/UserInterface.h"
#include "Orientations.h"
#include "Info.h"
#include "Mechanics/SymmetryVariants.h"

namespace openphase
{
using namespace std;

Nucleation::Nucleation(Settings& locSettings, const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void Nucleation::Initialize(Settings& locSettings)
{
    thisclassname = "Nucleation";

    TotalNx = locSettings.TotalNx;
    OffsetX = locSettings.OffsetX;
    TotalNy = locSettings.TotalNy;
    OffsetY = locSettings.OffsetY;
    TotalNz = locSettings.TotalNz;
    OffsetZ = locSettings.OffsetZ;

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;

    dNx = locSettings.dNx;
    dNy = locSettings.dNy;
    dNz = locSettings.dNz;

    Nphases = locSettings.Nphases;
    Nvariants = locSettings.Nvariants;

    for(size_t n = 0; n < Nvariants.size(); n++)
    {
        Nvariants[n] = max(Nvariants[n], (size_t)1);
    }

    Parameters.Allocate(Nphases, Nphases);
    GeneratedParticles.Allocate(Nphases, Nphases);
    NucleatedParticles.Allocate(Nphases, Nphases);

    SeedX = 253;
    SeedY = 4958;
    SeedZ = 54861;

    SeedA = 45;
    SeedB = 697;
    SeedC = 2597;

    SeedR = 34784;

    MobilityReduction.Allocate(Nphases, Nphases);
    for (size_t i = 0; i < Nphases; ++i)
    for (size_t j = 0; j < Nphases; ++j)
    {
        MobilityReduction(i,j) = 1.0;
    }
    NucleateEvery = 1;
    NumberOfAttempts = 100;

    iWidth = locSettings.iWidth;
    dx = locSettings.dx;

    VTKDir = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;
    TextDir = locSettings.TextDir;

    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    {
        Parameters(n, m).Allowed = false;
    }
    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void Nucleation::ReadInput(const std::string InputFileName)
{
    Info::WriteLineInsert("Nucleation input");
    Info::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };
    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);

    inp.close();
}

void Nucleation::ReadInput(std::stringstream& inp)
{
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
    // Reading Nucleation Allowance Booleans and corresponding nucleation densities:
    if(moduleLocation == -1)
    {
        Info::WriteWarning("Input for module \"" + thisclassname + "\" is not found in the input file! Nucleation is disabled!",thisclassname, "ReadInput()");
    }
    else
    {
        NucleateEvery = UserInterface::ReadParameterI(inp, moduleLocation,
                                               "NucleateEvery", false, 1);

        NumberOfAttempts = UserInterface::ReadParameterI(inp, moduleLocation,
                           "NumberOfAttempts", false, 100);

        for(size_t n = 0; n < Nphases; n++)
        for(size_t m = 0; m < Nphases; m++)
        {
            stringstream converter;
            converter << n << "_" << m;
            string counter = converter.str();

            string userinput = UserInterface::ReadParameterK(inp, moduleLocation,
                                      string("Allowed_") + counter, false, "NO");
            if(userinput == "NO")
            {
                Parameters(n, m).OrientationMode    = NucleiOrientationModes::Reference;
                Parameters(n, m).Density            = 0.0;
                Parameters(n, m).Generated          = false;
                Parameters(n, m).Distribution       = NucleiSizeDistributions::None;
            }
            else
            {
                NucleationParameters locNP;
                bool validuserinput = false;

                // Reading nucleation sites location
                if(userinput == "YES")
                {
                    locNP.LocationMode = NucleiLocationModes::BulkAndGB;
                    validuserinput = true;
                }
                if(userinput == "BULK")
                {
                    locNP.LocationMode = NucleiLocationModes::Bulk;
                    validuserinput = true;
                }
                if(userinput == "GB")
                {
                    locNP.LocationMode = NucleiLocationModes::GB;
                    validuserinput = true;
                }
                if(userinput == "INTERFACE")
                {
                    locNP.LocationMode = NucleiLocationModes::Interfaces;
                    validuserinput = true;
                }
                if(userinput == "JUNCTIONS")
                {
                    locNP.LocationMode = NucleiLocationModes::Junctions;
                    validuserinput = true;
                }
                if(userinput == "XBOTTOM")
                {
                    locNP.LocationMode = NucleiLocationModes::XBottom;
                    validuserinput = true;
                }
                if(userinput == "XTOP")
                {
                    locNP.LocationMode = NucleiLocationModes::ZBottom;
                    validuserinput = true;
                }
                if(userinput == "YBOTTOM")
                {
                    locNP.LocationMode = NucleiLocationModes::YBottom;
                    validuserinput = true;
                }
                if(userinput == "YTOP")
                {
                    locNP.LocationMode = NucleiLocationModes::YTop;
                    validuserinput = true;
                }
                if(userinput == "ZBOTTOM")
                {
                    locNP.LocationMode = NucleiLocationModes::ZBottom;
                    validuserinput = true;
                }
                if(userinput == "ZTOP")
                {
                    locNP.LocationMode = NucleiLocationModes::ZTop;
                    validuserinput = true;
                }
                if(validuserinput)
                {
                    locNP.Allowed = true;
                    // Reading seeds number or density
                    locNP.Nsites  = UserInterface::ReadParameterI(inp, moduleLocation,
                                               string("Nsites_")+counter, false, 0);
                    if(locNP.Nsites == 0)
                    {
                        locNP.Density = UserInterface::ReadParameterD(inp, moduleLocation,
                                              string("Density_")+counter, true, 0.0);

                        locNP.RelativeDensity = UserInterface::ReadParameterB(inp,
                           moduleLocation,string("RelDensity_") + counter, false, true);
                    }
    //////////////////////////// NEEDS CARE ///////////////////////////////////////////////////////////
                    if(locNP.Density <= 0.0)
                    {
                        Info::WriteExit("The number of nucleation sites Nsites_" + counter +
                                        " or the nucleation density Density_" + counter +
                                        " for phase pair " + counter + " is not valid",
                                        thisclassname, "ReadInput()");
                        exit(16);
                    }
    //////////////////////////////////////////////////////////////////////////////////////////////////
                    // Reading nucleation temperature
                    locNP.Tmin = UserInterface::ReadParameterD(inp, moduleLocation,
                                                        string("Tmin_") + counter);
                    locNP.Tmax = UserInterface::ReadParameterD(inp, moduleLocation,
                                                        string("Tmax_") + counter);
                    locNP.Shielding = UserInterface::ReadParameterD(inp, moduleLocation,
                                     string("Shielding_") + counter, false, iWidth);

                    locNP.Generated  = false;

                    // Reading seeds size distribution
                    string dist = UserInterface::ReadParameterK(inp, moduleLocation,
                                  string("Distribution_") + counter, false, "NONE");
                    bool distribution_set = false;
                    if(dist == "NONE")
                    {
                        locNP.Distribution = NucleiSizeDistributions::None;
                        locNP.DistMu    = UserInterface::ReadParameterD(inp, moduleLocation,
                                       string("Radius_") + counter, false, 0.5*iWidth*dx);
                        locNP.RadiusMIN = 0.0;
                        locNP.RadiusMAX = 2.0*locNP.DistMu;
                        distribution_set = true;
                    }
                    if(dist == "NORMAL")
                    {
                        locNP.Distribution = NucleiSizeDistributions::Normal;
                        locNP.DistMu    = UserInterface::ReadParameterD(inp, moduleLocation,
                                                   string("Center_") + counter, true, 0.0);
                        locNP.DistSigma = UserInterface::ReadParameterD(inp, moduleLocation,
                                                  string("Deviation_") + counter, true, 0.0);
                        locNP.RadiusMIN = 0.0;
                        locNP.RadiusMAX = 3.0*locNP.DistSigma;
                        distribution_set = true;
                    }
                    if(dist == "CAUCHY")
                    {
                        locNP.Distribution = NucleiSizeDistributions::Cauchy;
                        locNP.DistMu    = UserInterface::ReadParameterD(inp, moduleLocation,
                                                   string("Center_") + counter, true, 0.0);
                        locNP.DistSigma = UserInterface::ReadParameterD(inp, moduleLocation,
                                                  string("HalfWidth_") + counter, true, 0.0);
                        locNP.RadiusMIN = 0.0;
                        locNP.RadiusMAX = 10.0*locNP.DistSigma;
                        distribution_set = true;
                    }
                    if(dist == "UNIFORM")
                    {
                        locNP.Distribution = NucleiSizeDistributions::Uniform;
                        locNP.RadiusMIN = UserInterface::ReadParameterD(inp,
                                moduleLocation,string("RadiusMIN_") + counter, true, 0.0);
                        locNP.RadiusMAX = UserInterface::ReadParameterD(inp,
                                moduleLocation,string("RadiusMAX_") + counter, true, 0.0);
                        distribution_set = true;
                    }

                    if(!distribution_set)
                    {
                        //Wrong seed size distribution input
                        string message  = "Wrong or no input is given for the seeds size distribution!";
                        Info::WriteExit(message, thisclassname, "ReadInput()");
                        exit(1);
                    }

                    // Reading Nuclei Orientation Mode
                    string orient = UserInterface::ReadParameterK(inp,moduleLocation,
                                                             string("Orientation_")+counter);
                    bool orientation_set = false;
                    if(orient == "RANDOM")
                    {
                        locNP.OrientationMode = NucleiOrientationModes::Random;
                        orientation_set = true;
                    }
                    if(orient == "PARENT")
                    {
                        locNP.OrientationMode = NucleiOrientationModes::Parent;
                        orientation_set = true;
                    }
                    if(orient == "REFERENCE")
                    {
                        locNP.OrientationMode = NucleiOrientationModes::Reference;
                        orientation_set = true;
                    }
                    if(!orientation_set)
                    {
                        //Wrong seed orientation input
                        string message  = "Wrong or no input is given for the seeds orientation!";
                        Info::WriteExit(message, thisclassname, "ReadInput()");
                        exit(1);
                    }

                    if(Nvariants[n] != 0)
                    {
                        locNP.Nvariants = UserInterface::ReadParameterI(inp, moduleLocation,
                                                        string("Variants_") + counter, false, 1);
                    }

                    MobilityReduction(n,m) = UserInterface::ReadParameterD(inp, moduleLocation,
                                             string("MobilityReduction") + counter, false, 1.0);

                    locNP.Nucleated = 0;
                    Parameters(n, m) = locNP;
                }
                else
                {
                    Info::WriteExit("Nucleation mode for phase pair " + counter +
                                    " could not be read", thisclassname, "ReadInput()");
                    exit(16);
                }
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_RANK == 0)
        {
#endif
        fstream NucFile;
        NucFile.open (TextDir + "NucleationStatistics.dat", ios::out);
        NucFile.setf(ios::left);
        NucFile << setw(10) << "TimeStep";
        NucFile.setf(ios::right);
        NucFile << setw(10) << "PFindex"
                << setw(10) << "Nphase"
                << setw(10) << "Mphase"
                << setw(10) << "Variant"
                << setw(10) << "X"
                << setw(10) << "Y"
                << setw(10) << "Z"
                << setw(10) << "Q1"
                << setw(10) << "Q2"
                << setw(10) << "Q3"
                << setw(10) << "Q4"
                << setw(14) << "dGnuc"
                << setw(14) << "dGmin"
                << setw(10) << "Status" << endl;
        NucFile.close();
#ifdef MPI_PARALLEL
        }
#endif
    }

    Info::WriteLine();
    Info::WriteBlankLine();
}

void Nucleation::GenerateRandomSeeds(void)
{
    int min = 1;
    int max = 100000;

    SeedX = min + (rand() % (int)(max - min + 1));
    SeedY = min + (rand() % (int)(max - min + 1));
    SeedZ = min + (rand() % (int)(max - min + 1));

    SeedA = min + (rand() % (int)(max - min + 1));
    SeedB = min + (rand() % (int)(max - min + 1));
    SeedC = min + (rand() % (int)(max - min + 1));
}

void Nucleation::Clear()
{
    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    {
        GeneratedParticles(n,m).clear();
    }
}

vector<vector<bool> > Nucleation::GetNucleationEvents(void)
{
    /** This function returns pairs of phases indicating if nucleation events
     *  have happened or a given pair. Condition is "planted = true"*/

    vector<bool> tmp(Nphases, false);
    vector<vector<bool> > result(Nphases,tmp);
    if (initialized)
    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    for(auto ind = GeneratedParticles(n, m).begin();
             ind != GeneratedParticles(n, m).end(); ind++)
    if(ind->planted)
    {
        result[n][m] = true;
        result[m][n] = true;
        result[n][n] = true;
        result[m][m] = true;
    }

    return result;
}

vector<iVector3> Nucleation::GetNucleationSites(PhaseField& Phi)
{
    /** This function returns a vector of coordinates for all points
    where nucleation events have happened. Condition is "planted = true" and
    MAXVolume < RefVolume*/

    vector<iVector3> result;

    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    for(auto ind = GeneratedParticles(n, m).begin();
             ind != GeneratedParticles(n, m).end(); ind++)
    {
        if(ind->planted)
        {
            for(size_t idx = 0; idx < ind->PFindices.size(); idx++)
            {
                if(Phi.FieldsStatistics[ind->PFindices[idx]].VolumeRatio < 1.0)
                {
                    iVector3 temp({ind->x,ind->y,ind->z});
                    result.push_back(temp);
                }
            }
        }
    }
    return result;
}

bool Nucleation::IsShielded(std::vector<NucSite> LocGrainStorage, const int x, const int y, const int z, const double shielding) const
{
    for (auto it = LocGrainStorage.begin(); it != LocGrainStorage.end(); ++it)
    {
        const double xdis = std::min(std::fabs(x - it->x), std::min( std::fabs(x - it->x + TotalNx), std::fabs(x - it->x - TotalNx)));
        const double ydis = std::min(std::fabs(y - it->y), std::min( std::fabs(y - it->y + TotalNy), std::fabs(y - it->y - TotalNy)));
        const double zdis = std::min(std::fabs(z - it->z), std::min( std::fabs(z - it->z + TotalNz), std::fabs(z - it->z - TotalNz)));
        const double distanceSquare = xdis*xdis + ydis*ydis + zdis*zdis;
        if(distanceSquare < shielding*shielding) return true;
    }
    return false;
}

void Nucleation::SetNucleationSites(PhaseField& Phase, Temperature& Tx, size_t n, size_t m)
{
    // Using mersenne-twister 64 bit pseudo-random generator engine:
    mt19937_64 SizeGenerator(SeedR*Parameters(n, m).Nsites);

    mt19937_64 xPosGenerator(SeedX*n);
    mt19937_64 yPosGenerator(SeedY*n);
    mt19937_64 zPosGenerator(SeedZ*n);

    mt19937_64 OrientationGenerator1(SeedA*m);
    mt19937_64 OrientationGenerator2(SeedB*m);
    mt19937_64 OrientationGenerator3(SeedC*m);

    double DistMu = Parameters(n, m).DistMu;
    double DistSigma = Parameters(n, m).DistSigma;

    normal_distribution <double> SizeDistributionNormal(DistMu, DistSigma);
    cauchy_distribution <double> SizeDistributionCauchy(DistMu, DistSigma);
    uniform_real_distribution <double> SizeDistributionUniform(Parameters(n, m).RadiusMIN,
                                                               Parameters(n, m).RadiusMAX);

    uniform_int_distribution <int> xPosDistribution(0, TotalNx - 1);
    uniform_int_distribution <int> yPosDistribution(0, TotalNy - 1);
    uniform_int_distribution <int> zPosDistribution(0, TotalNz - 1);

    uniform_real_distribution <double> Q1Distribution(0, 1);
    uniform_real_distribution <double> Q2Distribution(0, 1);
    uniform_real_distribution <double> Q3Distribution(0, 1);

    uniform_real_distribution <double> A1Distribution(0.0, 2.0*Pi);
    uniform_real_distribution <double> A2Distribution(0.0, 2.0*Pi);
    uniform_real_distribution <double> A3Distribution(0.0, 2.0*Pi);

    size_t attempts = 0;
    double shielding = Parameters(n, m).Shielding;

    bool PrintMsg2 = false;
    while (Parameters(n, m).part < Parameters(n, m).Nsites and
           attempts < Parameters(n, m).Nsites*NumberOfAttempts)
    {
        bool enable_seed = false;

        int increment_part = 0;
        int increment_attempts = 0;

        double radius = 0.0;

        switch(Parameters(n, m).Distribution)
        {
            case NucleiSizeDistributions::Normal:
            {
                radius = 0.5*SizeDistributionNormal(SizeGenerator);
                break;
            }
            case NucleiSizeDistributions::Cauchy:
            {
                radius = 0.5*SizeDistributionCauchy(SizeGenerator);
                break;
            }
            case NucleiSizeDistributions::Uniform:
            {
                radius = 0.5*SizeDistributionUniform(SizeGenerator);
                break;
            }
            case NucleiSizeDistributions::None:
            {
                radius = 0.5*Phase.Eta;
                break;
            }
        }

        if(radius >= Parameters(n,m).RadiusMIN and
           radius <= Parameters(n,m).RadiusMAX)
        {
            int x = 0;
            int y = 0;
            int z = 0;

            if(Parameters(n, m).LocationMode == NucleiLocationModes::XBottom)
            {
                x = 0;
            }
            else if(Parameters(n, m).LocationMode == NucleiLocationModes::XTop)
            {
                x = TotalNx - 1;
            }
            else
            {
                x = xPosDistribution(xPosGenerator);
            }

            if(Parameters(n, m).LocationMode == NucleiLocationModes::YBottom)
            {
                y = 0;
            }
            else if(Parameters(n, m).LocationMode == NucleiLocationModes::YTop)
            {
                y = TotalNy - 1;
            }
            else
            {
                y = yPosDistribution(yPosGenerator);
            }

            if(Parameters(n, m).LocationMode == NucleiLocationModes::ZBottom)
            {
                z = 0;
            }
            else if(Parameters(n, m).LocationMode == NucleiLocationModes::ZTop)
            {
                z = TotalNz - 1;
            }
            else
            {
                z = zPosDistribution(zPosGenerator);
            }

            if (x - OffsetX >= 0.0 && x - OffsetX < Nx and
                y - OffsetY >= 0.0 && y - OffsetY < Ny and
                z - OffsetZ >= 0.0 && z - OffsetZ < Nz and
                !IsShielded(GeneratedParticles(n, m),x,y,z,shielding))
            {
                int i = x - OffsetX;
                int j = y - OffsetY;
                int k = z - OffsetZ;

                switch (Parameters(n, m).LocationMode)
                {
                    case NucleiLocationModes::XBottom:
                    case NucleiLocationModes::XTop:
                    case NucleiLocationModes::YBottom:
                    case NucleiLocationModes::YTop:
                    case NucleiLocationModes::ZBottom:
                    case NucleiLocationModes::ZTop:
                    {
                        if(Phase.Fractions(i,j,k)({m}) > 0.25)
                        {
                            enable_seed = true;
                        }
                        break;
                    }
                    case NucleiLocationModes::Bulk:
                    {
                        if( !Phase.Interface(i,j,k) and
                            (Phase.Fractions(i,j,k)({m}) == 1.0))
                        {
                            enable_seed = true;
                        }
                        break;
                    }
                    case NucleiLocationModes::BulkAndGB:
                    {
                        if(Phase.Fractions(i,j,k)({m}) == 1.0)
                        {
                            enable_seed = true;
                        }
                        break;
                    }
                    case NucleiLocationModes::GB:
                    {
                        if( Phase.Interface(i,j,k) and
                           (Phase.Fractions(i,j,k)({m}) == 1.0))
                        {
                            enable_seed = true;
                        }
                        break;
                    }
                    case NucleiLocationModes::Junctions:
                    {
                        if((Phase.Fractions(i,j,k)({m}) == 1.0) and
                           (Phase.Fields(i,j,k).size() > 2))
                        {
                            enable_seed = true;
                        }
                        break;
                    }
                    case NucleiLocationModes::Interfaces:
                    {
                        if( Phase.Interface(i,j,k) and
                           (Phase.Fractions(i,j,k)({m}) < 1.0) and
                           (Phase.Fractions(i,j,k)({m}) > 0.25))
                        {
                            enable_seed = true;
                        }
                        break;
                    }
                }
            }

#ifdef MPI_PARALLEL
            int enable_seed_tmp = enable_seed;
            MPI_Allreduce(&enable_seed_tmp, &enable_seed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
            if(enable_seed)
            {
                NucSite locSite;
                locSite.x = x;
                locSite.y = y;
                locSite.z = z;
                locSite.radius = radius;
                locSite.planted = false;
                locSite.growing = false;

                if(Parameters(n, m).OrientationMode == NucleiOrientationModes::Random)
                {
                    Quaternion tempQuat;

                    switch(100*dNx + 10*dNy + dNz)
                    {
                        case 100:
                        case 010:
                        case 001: // No rotation in 1D
                        {
                            locSite.Q.set(1.0, 0.0, 0.0, 0.0);
                            break;
                        }
                        case 110:
                        case 101:
                        case 011: // Only in-plain rotation in 2D
                        {
                            double a1 = 0.0;
                            double a2 = 0.0;
                            double a3 = 0.0;

                            if(Phase.dNx == 0) a1 = A1Distribution(OrientationGenerator1);
                            if(Phase.dNy == 0) a2 = A2Distribution(OrientationGenerator2);
                            if(Phase.dNz == 0) a3 = A3Distribution(OrientationGenerator3);

                            EulerAngles ph1({a1,a2,a3},XYZ);

                            locSite.Q = ph1.getQuaternion().normalized();
                            break;
                        }
                        case 111: // Full rotation freedom in 3D
                        {
                            double u1 = Q1Distribution(OrientationGenerator1);
                            double u2 = Q2Distribution(OrientationGenerator2);
                            double u3 = Q3Distribution(OrientationGenerator3);

                            tempQuat.set(sqrt(1.0-u1)*sin(2.0*Pi*u2),
                                         sqrt(1.0-u1)*cos(2.0*Pi*u2),
                                         sqrt(u1)*sin(2.0*Pi*u3),
                                         sqrt(u1)*cos(2.0*Pi*u3));

                            locSite.Q = tempQuat.normalized();
                            break;
                        }
                        default: // No rotation in zero dimensions.
                        {
                            locSite.Q.set(1.0, 0.0, 0.0, 0.0);
                            break;
                        }
                    }
                }

                std::stringstream message;
                message << "Nucleation: Generated seed particle " << Parameters(n, m).part << " at ["
                        << x << ", " << y << ", " << z << "] and effective radius of "
                        << radius;
                Info::WriteSimple(message.str());

                GeneratedParticles(n, m).push_back(locSite);
                increment_part = 1;
                PrintMsg2 = true;
            }
            else
            {
                increment_attempts = 1;
            }
        }
        Parameters(n, m).part += increment_part;
        attempts += increment_attempts;
    }
    if (PrintMsg2)
    {
        std::stringstream message2;
        message2 << "Nucleation: Generated " << Parameters(n, m).part
                 << " nucleation sites (of " << Parameters(n, m).Nsites
                 << ") for phase " << n
                 << " in phase " << m << ".";
        Info::WriteSimple(message2.str());
    }
}

void Nucleation::GenerateNucleationSites(PhaseField& Phase, Temperature& Tx)
{
    size_t ngrains = Phase.FieldsStatistics.size();
    double TotalVolume = TotalNx*TotalNy*TotalNz;
    double RealUnitsVolume = TotalVolume * pow(Phase.dx,3);

    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    if(Parameters(n, m).Allowed and
       Parameters(n, m).Tmin <= Tx.Tmax and
       Parameters(n, m).Tmax >= Tx.Tmin)
    {
        double PhaseFractions_m = 0.0;

        if(Parameters(n, m).RelativeDensity)
        {
            for(size_t idx = 0; idx < ngrains; idx++)
            if(Phase.FieldsStatistics[idx].Phase == m)
            PhaseFractions_m += Phase.FieldsStatistics[idx].Volume/TotalVolume;
        }
        else
        {
            PhaseFractions_m = 1.0;
        }

        if(Parameters(n, m).Nsites == 0 and
           Parameters(n, m).Density != 0.0 and !Parameters(n, m).Generated)
        {
            Parameters(n, m).Nsites = Parameters(n, m).Density*RealUnitsVolume
                                                              *PhaseFractions_m;
        }

        if(Parameters(n, m).Nsites)
        {
            SetNucleationSites(Phase,Tx,n,m);
            Parameters(n, m).Generated = true;
        }
        else
        {
            std::stringstream message3;
            message3 << "Nucleation: Too low particles density! No nucleation "
                     << "sites were generated for phase " << n << " in phase "
                     << m << ".";
            Info::WriteSimple(message3.str());
        }
    }
}
///TODO: implement MPI parallelism
void Nucleation::ReGenerateNucleationSites(PhaseField& Phase, Temperature& Tx)
{
    GenerateRandomSeeds();

    size_t ngrains = Phase.FieldsStatistics.size();
    double TotalVolume = TotalNx*TotalNy*TotalNz;
    double RealUnitsVolume = TotalVolume * pow(Phase.dx,3);

    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    if(Parameters(n, m).Allowed and
       Parameters(n, m).Tmin <= Tx.Tmax and
       Parameters(n, m).Tmax >= Tx.Tmin)
    {
        double PhaseFractions_n = 0.0;
        double PhaseFractions_m = 0.0;

        for(size_t idx = 0; idx < ngrains; idx++)
        {
            if(Phase.FieldsStatistics[idx].Phase == m)
            {
                PhaseFractions_m += Phase.FieldsStatistics[idx].Volume;
            }
            if(Phase.FieldsStatistics[idx].Phase == n)
            {
                PhaseFractions_n += Phase.FieldsStatistics[idx].Volume;
            }
        }

        if(PhaseFractions_n < DBL_EPSILON)
        {
            Parameters(n, m).Nucleated = 0;
        }

        if(Parameters(n, m).RelativeDensity)
        {
            PhaseFractions_m /= TotalVolume;
        }
        else
        {
            PhaseFractions_m = 1.0;
        }

        Parameters(n,m).Nsites = Parameters(n, m).Density*RealUnitsVolume
                                                         *PhaseFractions_m;

        vector<NucSite> tempNucSites;
        for (size_t i = 0; i < GeneratedParticles(n, m).size(); i++)
        {
            if (GeneratedParticles(n, m)[i].growing)
            {
                int x = GeneratedParticles(n, m)[i].x - OffsetX;
                int y = GeneratedParticles(n, m)[i].y - OffsetY;
                int z = GeneratedParticles(n, m)[i].z - OffsetZ;

                if ((Phase.Fractions(x, y, z) ({ n }) > 0.0)
                and (Phase.Fractions(x, y, z) ({ m }) > 0.0))
                {
                    tempNucSites.push_back(GeneratedParticles(n, m)[i]);
                }
                else
                {
                    Parameters(n, m).Nucleated--;
                }
            }
        }
        GeneratedParticles(n, m).clear();
        GeneratedParticles(n, m) = tempNucSites;

        if(Parameters(n,m).Nsites > Parameters(n,m).Nucleated)
        {
            Parameters(n,m).Nsites = Parameters(n,m).Nsites -
                                               Parameters(n,m).Nucleated;
        }
        else
        {
            Parameters(n,m).Nsites = 0;
        }

        if(Parameters(n, m).Nsites)
        {
            SetNucleationSites(Phase,Tx,n,m);
            Parameters(n, m).Generated = true;
        }
        else
        {
            std::stringstream message3;
            message3 << "Nucleation: Too low particles density! No nucleation "
                     << "sites were generated for phase " << n << " in phase "
                     << m << ".";
            Info::WriteSimple(message3.str());
        }
    }
}

void Nucleation::WriteStatistics(int tStep, long int PFindex, size_t NucleatingPhase,
                                 size_t MatrixPhase, size_t Variant, int Xpos, int Ypos, int Zpos,
                                 Quaternion Q, double dGnuc,
                                 double dGmin, string status) const
{
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
    {
#endif

    fstream NucFile;
    NucFile.open (TextDir + "NucleationStatistics.dat", ios::app);
    NucFile.setf(ios::left);
    NucFile << setw(10)  << setprecision(6) << tStep;
    NucFile.setf(ios::right);
    NucFile << setw(10)  << setprecision(6) << PFindex
            << setw(10)  << setprecision(5) << NucleatingPhase
            << setw(10)  << setprecision(5) << MatrixPhase
            << setw(10)  << setprecision(5) << Variant
            << setw(10)  << setprecision(5) << Xpos
            << setw(10)  << setprecision(5) << Ypos
            << setw(10)  << setprecision(5) << Zpos
            << setw(10)  << setprecision(4) << Q[0]
            << setw(10)  << setprecision(4) << Q[1]
            << setw(10)  << setprecision(4) << Q[2]
            << setw(10)  << setprecision(4) << Q[3]
            << setw(14)  << setprecision(6) << dGnuc
            << setw(14)  << setprecision(6) << dGmin
            << setw(10)  << status          << endl;
    NucFile.close();
#ifdef MPI_PARALLEL
    }
#endif
}

void Nucleation::PlantNuclei(PhaseField& Phase, int tStep)
{
    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    for(auto ind = GeneratedParticles(n, m).begin();
             ind != GeneratedParticles(n, m).end(); ind++)
    {
        if(!ind->growing)
        {
            size_t ParentGrainIndex = 0;
            vector<int> locVariantsOfN(Nvariants[n], 0);
            if(ind->x - OffsetX >= 0 and ind->x - OffsetX < Nx and
               ind->y - OffsetY >= 0 and ind->y - OffsetY < Ny and
               ind->z - OffsetZ >= 0 and ind->z - OffsetZ < Nz)
            {
                double locMax = 0.0;
                for(auto alpha  = Phase.Fields(ind->x - OffsetX, ind->y - OffsetY, ind->z - OffsetZ).begin();
                         alpha != Phase.Fields(ind->x - OffsetX, ind->y - OffsetY, ind->z - OffsetZ).end(); ++alpha)
                {
                    if(Phase.FieldsStatistics[alpha->index].Phase == m and alpha->value > locMax)
                    {
                        locMax = alpha->value;
                        ParentGrainIndex = alpha->index;
                    }
                    if(Phase.FieldsStatistics[alpha->index].Phase == n)
                    {
                        size_t locVariant = Phase.FieldsStatistics[alpha->index].Variant;
                        locVariantsOfN[locVariant] = true;
                    }
                }
            }
#ifdef MPI_PARALLEL
            unsigned long loc_ParentGrainIndex = ParentGrainIndex;
            MPI_Allreduce(&loc_ParentGrainIndex, &ParentGrainIndex, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
            vector<int> tmpVariansOfN = locVariantsOfN;
            MPI_Allreduce(tmpVariansOfN.data(), locVariantsOfN.data(), Nvariants[n], MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
            for(size_t v = 0; v < Nvariants[n]; v++)
            if(!locVariantsOfN[v])
            {
                size_t locIndex = Phase.PlantGrainNucleus(n, ind->x, ind->y, ind->z);

//                double RefVolume = Pi;
//                double locRadius = ind->radius/Phase.dx;
//
//                if(Nx < locRadius) RefVolume *= (TotalNx+1)/2.0;
//                else RefVolume *= 1.1*locRadius;
//                if(Ny < locRadius) RefVolume *= (TotalNy+1)/2.0;
//                else RefVolume *= 1.1*locRadius;
//                if(Nz < locRadius) RefVolume *= (TotalNz+1)/2.0;
//                else RefVolume *= 1.1*locRadius;
//
//                Phase.FieldsStatistics[locIndex].RefVolume = max(RefVolume,
//                        Phase.FieldsStatistics[locIndex].RefVolume)/
//                        Parameters(n,m).Nvariants;
//                //Phase.FieldsStatistics[locIndex].RefVolume = RefVolume;

                if(Parameters(n, m).OrientationMode == NucleiOrientationModes::Parent)
                {
                    Phase.FieldsStatistics[locIndex].Orientation =
                           Phase.FieldsStatistics[ParentGrainIndex].Orientation;
                }
                else
                {
                    Phase.FieldsStatistics[locIndex].Orientation = ind->Q;
                }

                Phase.FieldsStatistics[locIndex].Variant = v;
                Phase.FieldsStatistics[locIndex].Parent  = ParentGrainIndex;

                ind->PFindices.push_back(locIndex);
            }
            ind->planted = true;
            ind->time_stamp = tStep;
            ind->growing = false;
        }
    }
}

void Nucleation::CheckNuclei(PhaseField& Phase, InterfaceProperties& IP, DrivingForce& dG, int tStep)
{
    // Using mersenne-twister 64 bit pseudo-random generator engine:
    double minProbability = 0.0;
    double maxProbability = 1.0;
    uniform_real_distribution <double> VariantSelector(minProbability, maxProbability);

    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    for(auto ind = GeneratedParticles(n, m).begin();
             ind != GeneratedParticles(n, m).end(); ind++)
    {
        if(ind->planted and !ind->growing)
        {
            /** Minimum driving force barrier the nucleus has to overcome.*/
            double dGmin = 2.0*IP.InterfaceEnergy(n, m).MaxEnergy/ind->radius;

            vector<double> dGloc(Nvariants[n], 0.0);
            if(ind->x - OffsetX >= 0 and ind->x - OffsetX < Nx and
               ind->y - OffsetY >= 0 and ind->y - OffsetY < Ny and
               ind->z - OffsetZ >= 0 and ind->z - OffsetZ < Nz)
            {
                for(auto it  = dG.Raw(ind->x - OffsetX, ind->y - OffsetY, ind->z - OffsetZ).begin();
                         it != dG.Raw(ind->x - OffsetX, ind->y - OffsetY, ind->z - OffsetZ).end(); ++it)
                {
                    for(size_t v = 0; v < ind->PFindices.size(); v++)
                    {
                        if(it->indexA == ind->PFindices[v])
                        {
                            dGloc[v] += it->value1;
                        }
                        if(it->indexB == ind->PFindices[v])
                        {
                            dGloc[v] -= it->value1;
                        }
                    }
                }
            }
#ifdef MPI_PARALLEL
            vector<double> loc_dGloc = dGloc;
            MPI_Allreduce(loc_dGloc.data(), dGloc.data(), Nvariants[n], MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
            vector<double> Probability(ind->PFindices.size(), 0);
            for(size_t v = 0; v < ind->PFindices.size(); v++)
            {
                if(dGloc[v] > dGmin)
                {
                    Probability[v] = VariantSelector(VariantsGenerator);
                }
            }

            vector<bool> ProbableVariants(Nvariants[n], false);

            for(size_t p = 0; p < Parameters(n,m).Nvariants; p++)
            {
                double loc_minProbability = minProbability;
                int loc_variant = -1;
                for(size_t v = 0; v < ind->PFindices.size(); v++)
                {
                    if(Probability[v] > loc_minProbability and !ProbableVariants[v])
                    {
                        loc_minProbability = Probability[v];
                        loc_variant = v;
                    }
                }
                if(loc_variant != -1)
                {
                    ProbableVariants[loc_variant] = true;
                }
            }

            size_t NucleiCounter = ind->PFindices.size();
            for(size_t v = 0; v < ind->PFindices.size(); v++)
            {
                if(!ProbableVariants[v])
                {
                    NucleiCounter--;

                    Phase.FieldsStatistics[ind->PFindices[v]].Exist = false;
                    Phase.FieldsStatistics[ind->PFindices[v]].Stage = 0;

                    if(ind->x - OffsetX >= 0 and ind->x - OffsetX < Nx and
                       ind->y - OffsetY >= 0 and ind->y - OffsetY < Ny and
                       ind->z - OffsetZ >= 0 and ind->z - OffsetZ < Nz)
                    {
                        for(auto it  = dG.Raw(ind->x - OffsetX, ind->y - OffsetY, ind->z - OffsetZ).begin();
                                 it != dG.Raw(ind->x - OffsetX, ind->y - OffsetY, ind->z - OffsetZ).end(); )
                        {
                            if(it->indexA == ind->PFindices[v] or
                               it->indexB == ind->PFindices[v])
                            {
                                it = dG.Raw(ind->x - OffsetX, ind->y - OffsetY, ind->z - OffsetZ).erase(it);
                            }
                            else
                            {
                                ++it;
                            }
                        }
                    }
                }
                else if (Probability[v] != 0.0)
                {
                    WriteStatistics(tStep, ind->PFindices[v], n, m, v,
                        ind->x, ind->y, ind->z, ind->Q, dGloc[v], dGmin, "Planted");
                    ind->growing = true;
                    Parameters(n, m).Nucleated++;
                }
            }

            if(NucleiCounter == 0)
            {
                ind->planted = false;
                ind->PFindices.resize(0);
            }
        }
    }
}

void Nucleation::Write(int tStep)
{
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
    {
#endif
    string FileName = UserInterface::MakeFileName(RawDataDir, "Nucleation_",
                                                             tStep, ".dat");
    fstream out(FileName.c_str(), ios::out);

    if (!out)
    {
        Info::WriteExit("File/" + FileName + "could not be created",
                                                     thisclassname);
        exit(1);
    }

    for(size_t n = 0; n < Nphases; ++n)
    {
        for(size_t m = 0; m < Nphases; ++m)
        {
            out << GeneratedParticles(n,m).size() << " " ;
        }
        out << endl;
    }
    for(size_t n = 0; n < Nphases; ++n)
    for(size_t m = 0; m < Nphases; ++m)
    for(size_t i = 0; i < GeneratedParticles(n,m).size(); ++i)
    {
        out << GeneratedParticles(n,m)[i].x  << " "
            << GeneratedParticles(n,m)[i].y  << " "
            << GeneratedParticles(n,m)[i].z  << " " << endl;
        out << GeneratedParticles(n,m)[i].Q[0] << " "
            << GeneratedParticles(n,m)[i].Q[1] << " "
            << GeneratedParticles(n,m)[i].Q[2] << " "
            << GeneratedParticles(n,m)[i].Q[3] << " " << endl;
        out << GeneratedParticles(n,m)[i].radius << endl;
    }
    out.close();
#ifdef MPI_PARALLEL
    }
#endif
}

void Nucleation::Read(int tStep)
{
    string FileName = UserInterface::MakeFileName(RawDataDir, "Nucleation_",
                                                             tStep, ".dat");
    fstream inp(FileName.c_str(), ios::in);

    if (!inp)
    {
        Info::WriteExit("File/" + FileName + "could not be opened",
                                                    thisclassname);
        exit(1);
    };

    for(size_t n = 0; n < Nphases; ++n)
    for(size_t m = 0; m < Nphases; ++m)
    {
        int tmp;
        inp >> tmp;

        if(tmp)
        {
            GeneratedParticles(n,m).resize(tmp);
        }
    }

    for(size_t n = 0; n < Nphases; ++n)
    for(size_t m = 0; m < Nphases; ++m)
    for(size_t i = 0; i < GeneratedParticles(n,m).size(); ++i)
    {
        inp >> GeneratedParticles(n,m)[i].x
            >> GeneratedParticles(n,m)[i].y
            >> GeneratedParticles(n,m)[i].z;
        inp >> GeneratedParticles(n,m)[i].Q[0]
            >> GeneratedParticles(n,m)[i].Q[1]
            >> GeneratedParticles(n,m)[i].Q[2]
            >> GeneratedParticles(n,m)[i].Q[3];
        inp >> GeneratedParticles(n,m)[i].radius;
    }
    inp.close();
    Info::WriteStandard(thisclassname, "Binary input loaded");

}

void Nucleation::WriteH5(int tStep, H5Interface& H5)
{
    #ifdef H5OP
    std::vector<double> dbuffer;

    for(size_t n = 0; n < Nphases; ++n)
    {
        for(size_t m = 0; m < Nphases; ++m)
        {
            dbuffer.push_back(GeneratedParticles(n,m).size());
        }
    }
    for(size_t n = 0; n < Nphases; ++n)
    for(size_t m = 0; m < Nphases; ++m)
    for(size_t i = 0; i < GeneratedParticles(n,m).size(); ++i)
    {
        dbuffer.push_back(GeneratedParticles(n,m)[i].x);
        dbuffer.push_back(GeneratedParticles(n,m)[i].y);
        dbuffer.push_back(GeneratedParticles(n,m)[i].z);
        dbuffer.push_back(GeneratedParticles(n,m)[i].Q[0]);
        dbuffer.push_back(GeneratedParticles(n,m)[i].Q[1]);
        dbuffer.push_back(GeneratedParticles(n,m)[i].Q[2]);
        dbuffer.push_back(GeneratedParticles(n,m)[i].Q[3]);
        dbuffer.push_back(GeneratedParticles(n,m)[i].radius);
    }
    H5.WriteCheckPoint(tStep, "Nucleation", dbuffer);
    #else
    std::cerr << "OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"" << std::endl;
    exit(5);
    #endif
}

void Nucleation::ReadH5(int tStep, H5Interface& H5)
{
    #ifdef H5OP
    std::vector<double> dbuffer;
    H5.ReadCheckPoint(tStep, "Nucleation", dbuffer);
    int i = 0;
    for(size_t n = 0; n < Nphases; ++n)
    for(size_t m = 0; m < Nphases; ++m)
    {
        int tmp;
        tmp = dbuffer[i]; ++i;
        if(tmp)
        {
            GeneratedParticles(n,m).resize(tmp);
        }
    }

    for(size_t n = 0; n < Nphases; ++n)
    for(size_t m = 0; m < Nphases; ++m)
    for(size_t i = 0; i < GeneratedParticles(n,m).size(); ++i)
    {
        GeneratedParticles(n,m)[i].x = dbuffer[i]; ++i;
        GeneratedParticles(n,m)[i].y = dbuffer[i]; ++i;
        GeneratedParticles(n,m)[i].z = dbuffer[i]; ++i;
        GeneratedParticles(n,m)[i].Q[0] = dbuffer[i]; ++i;
        GeneratedParticles(n,m)[i].Q[1] = dbuffer[i]; ++i;
        GeneratedParticles(n,m)[i].Q[2] = dbuffer[i]; ++i;
        GeneratedParticles(n,m)[i].Q[3] = dbuffer[i]; ++i;
        GeneratedParticles(n,m)[i].radius= dbuffer[i]; ++i;
    }
    #else
    std::cerr << "OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"" << std::endl;
    exit(5);
    #endif
}
} //namespace openphase
