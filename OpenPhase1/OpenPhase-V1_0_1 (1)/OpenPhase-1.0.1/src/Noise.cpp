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
 *   Main contributors :   Raphael Schiedung
 *
 */

#include "Base/Includes.h"
#include "DrivingForce.h"
#include "Info.h"
#include "Noise.h"
#include "Settings.h"
#include "Temperature.h"
#include "VTK.h"
#include "fftw3.h"
#include <chrono>

namespace openphase
{
using namespace std;

void Noise::Initialize(Settings& locSettings)
{
    thisclassname        = "Noise";

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;

    dNx = locSettings.dNx;
    dNy = locSettings.dNy;
    dNz = locSettings.dNz;

    dx = locSettings.dx;
    dy = locSettings.dx;
    dz = locSettings.dx;

    Nphases   = locSettings.Nphases;

    // Allocate three-dimensional Fourier space noise
    RandomFourier = new complex<double>[Nx*Ny*Nz]();
    RandomReal    = new complex<double>[Nx*Ny*Nz]();

    // Allocate raw Noise Storage
    Raw.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, 0);

    // Assign storages to fftw pointers
    fftw_In  = reinterpret_cast<fftw_complex*>(RandomFourier);
    fftw_Out = reinterpret_cast<fftw_complex*>(RandomReal);

    // Create FFT-Plan, which determines how the FFT will be executed
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
    FFTBackward = fftw_plan_dft_3d(Nx, Ny, Nz, fftw_In,fftw_Out,FFTW_BACKWARD,FFTW_ESTIMATE);

    // Initialize random number distribution
    distribution = normal_distribution<double>(0.0,0.5);

    VTKDir = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;

    size_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator = std::default_random_engine(seed);

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void Noise::ReadInput(const std::string InputFileName)
{
    Info::WriteLineInsert("Noise input");
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

void Noise::ReadInput(std::stringstream& inp)
{
    // Read Parameters
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    gamma            = UserInterface::ReadParameterD(inp, moduleLocation, string("dgamma"),false,0);
    TimeSteps        = UserInterface::ReadParameterI(inp, moduleLocation, string("iTime"));
    CutOffWaveLength = UserInterface::ReadParameterI(inp, moduleLocation, string("dWaveLength"));
    FixAmpl          = UserInterface::ReadParameterI(inp, moduleLocation, string("dFixAmpl"));

    // Check if the input TimeSteps makes sense
    if (TimeSteps < 1)
    {
        std::string message = "A input value TimeSteps < 1 does not make sense! "
            "Set TimeSteps to 1!";
        Info::WriteWarning( message, thisclassname, "ReadInput");
        TimeSteps = 0;
    }

    // Check if the input CutOffWaveLength makes sense
    if (CutOffWaveLength < 0)
    {
        std::string message = "A input value CutOffWaveLength < 0 does not make sense! "
            "Set CutOffWaveLength to 0!";
        Info::WriteWarning( message, thisclassname, "ReadInput");
        CutOffWaveLength = 0;
    }
    Info::WriteLine();
    Info::WriteBlankLine();
}

void Noise::Generate(void)
{
    double Lx = 2.0 * Pi / (double(Nx) * dx);
    double Ly = 2.0 * Pi / (double(Ny) * dy);
    double Lz = 2.0 * Pi / (double(Nz) * dz);


    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        dVector3 WaveVector;
        WaveVector[0] = Lx*(i*(i <= Nx/2) - (Nx-i)*(i > Nx/2));
        WaveVector[1] = Ly*(j*(j <= Ny/2) - (Ny-j)*(j > Ny/2));
        WaveVector[2] = Lz*(k*(k <= Nz/2) - (Nz-k)*(k > Nz/2));

        if(WaveVector.abs() < 2.0 * Pi/CutOffWaveLength)
        {
            // Generate Gaussian white noise coefficients
            complex<double> cijk{distribution(generator),distribution(generator)};
            //cijk /= abs(cijk);

            // Store computed coefficient
            RandomFourier[k + Nz*(j + Ny*i)] = cijk;
        }
        else
        {
            RandomFourier[k + Nz*(j + Ny*i)] = complex<double>(0.0,0.0);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    // Do Fourier transform
    fftw_execute(FFTBackward);

    // Calculate global maximum of RandomReal
    double Max = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,reduction(max:Max))
    {
        Max = max(Max, abs(real(RandomReal[k + Nz*(j + Ny*i)])));
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        RandomReal[k + Nz*(j + Ny*i)] /= Max;
        RandomReal[k + Nz*(j + Ny*i)] -= Raw(i,j,k);
        RandomReal[k + Nz*(j + Ny*i)] /= TimeSteps;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Noise::UpdateRaw(int tStep)
{
    if (tStep == 0)
    {
        // Generate initial noise distribution
        Generate();

        // Save initial noise distribution and add it to the driving force

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
        {
            // Set initial noise
            Raw(i,j,k) = TimeSteps * real(RandomReal[k + Nz*(j + Ny*i)]);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        // Generate new noise distribution in order to able to interpolate
        // between the old distribution and the new one.
        Generate();
    }
    else
    {
        // First check if a new noise Distribution has to be generated
        if (!(tStep%TimeSteps)) Generate();

        // Interpolate between old and new noise distribution by adding RandomReal
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
        {
            // Compute noise for this time step
            Raw(i,j,k) += real(RandomReal[k + Nz*(j + Ny*i)]);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    // Calculate the global maximum of Raw
    double Max = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,reduction(max:Max))
    {
        Max = max(Max, abs(Raw(i,j,k)));
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        Raw(i,j,k) /= Max;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}


void Noise::AddToDrivingForce(int tStep, Temperature& Temp, DrivingForce& DF)
{
    UpdateRaw(tStep);

    // Interpolate between old and new noise distribution by adding RandomReal
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        // Add noise to driving force
        for (auto it = DF.Raw(i,j,k).cbegin(); it < DF.Raw(i,j,k).cend(); ++it)
        {
            double Ampl = sqrt( 2 * PhysicalConstants::k_Boltzmann * Temp(i,j,k)/gamma);
            DF.Raw(i,j,k).add_asym1(it->indexA, it->indexB, Ampl * Raw(i,j,k));
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Noise::AddToDrivingForce(int tStep, DrivingForce& DF)
{
    UpdateRaw(tStep);

    // Interpolate between old and new noise distribution by adding RandomReal
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        // Add noise to driving force
        for (auto it = DF.Raw(i,j,k).cbegin(); it < DF.Raw(i,j,k).cend(); ++it)
        {
            DF.Raw(i,j,k).add_asym1(it->indexA, it->indexB, FixAmpl * Raw(i,j,k));
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Noise::WriteVTK(int tStep, const Settings& locSettings, const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = UserInterface::MakeFileName(locSettings.VTKDir, thisclassname+"_", tStep, ".vts");

    ListOfFields.push_back((VTK::Field_t){"Noise",  [this](int i,int j,int k){return Raw(i,j,k);}});
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

}// end of name space
