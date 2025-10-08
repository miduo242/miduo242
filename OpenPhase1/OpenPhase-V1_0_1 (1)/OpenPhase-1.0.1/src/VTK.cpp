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
 *   Main contributors :   Philipp Engels, Raphael Schiedung, Muhammad Adil Ali,
 *                         Oleg Shchyglo
 *
 */

#include "VTK.h"
#include "Settings.h"
#include "Base/UserInterface.h"
#include "Info.h"

namespace openphase
{

using namespace std;

vector<string> voigtcon {"xx", "yy", "zz", "yz", "xz", "xy"};
vector<string> matrixcon  {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"};


void VTK::Write(
    const std::string Filename,
    const Settings& locSettings,
    std::vector<Field_t> ListOfFields,
    const int precision,
    const int resolution)
{
    const int Nx = (resolution == 1) ? locSettings.Nx : (locSettings.dNx + 1)*locSettings.Nx;
    const int Ny = (resolution == 1) ? locSettings.Ny : (locSettings.dNy + 1)*locSettings.Ny;
    const int Nz = (resolution == 1) ? locSettings.Nz : (locSettings.dNz + 1)*locSettings.Nz;

#ifndef MPI_PARALLEL
    ofstream vtk_file(Filename.c_str());
    VTK::WriteHeader(vtk_file, Nx, Ny, Nz);
    {
        WritePointData(vtk_file,ListOfFields,Nx,Ny,Nz,precision);
    }
    VTK::WriteEndPointData(vtk_file);
    VTK::WriteCoordinates(vtk_file, locSettings, resolution);
    VTK::CloseFile(vtk_file);
#else
    const int TotalNx = (resolution == 1) ? locSettings.TotalNx : (locSettings.dNx + 1)*locSettings.TotalNx;
    const int TotalNy = (resolution == 1) ? locSettings.TotalNy : (locSettings.dNy + 1)*locSettings.TotalNy;
    const int TotalNz = (resolution == 1) ? locSettings.TotalNz : (locSettings.dNz + 1)*locSettings.TotalNz;

    const int OffsetX = (resolution == 1) ? locSettings.OffsetX : (locSettings.dNx + 1)*locSettings.OffsetX;
    const int OffsetY = (resolution == 1) ? locSettings.OffsetY : (locSettings.dNy + 1)*locSettings.OffsetY;
    const int OffsetZ = (resolution == 1) ? locSettings.OffsetZ : (locSettings.dNz + 1)*locSettings.OffsetZ;

    MPI_File fh;
    MPI_Status status;
    std::stringstream buffer;
    std::stringstream hbuffer;
    std::stringstream tbuffer;
    MPI_File_delete(Filename.c_str(), MPI_INFO_NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_open(MPI_COMM_WORLD, Filename.c_str(),
        MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_Offset myoffset = 0;
    MPI_File_set_view(fh, myoffset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);

    buffer << "<Piece Extent=\""
           << OffsetX << " " << Nx-1 + OffsetX << " "
           << OffsetY << " " << Ny-1 + OffsetY << " "
           << OffsetZ << " " << Nz-1 + OffsetZ << "\">\n";
    {
        WritePointData(buffer,ListOfFields,Nx,Ny,Nz,precision);
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, locSettings, resolution);
    buffer << "</Piece>\n";

    MPI_Barrier(MPI_COMM_WORLD);
    hbuffer << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";
    hbuffer << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    hbuffer << "<StructuredGrid WholeExtent=\""
            << 0 << " " << TotalNx-1 << " "
            << 0 << " " << TotalNy-1 << " "
            << 0 << " " << TotalNz-1 << "\"> \n";
    tbuffer << "</StructuredGrid> \n";
    tbuffer << "</VTKFile> \n";
    size_t buffersize = buffer.str().size();
    size_t headersize = hbuffer.str().size();
    if (MPI_RANK == 0)
    {
        const std::string tmp = hbuffer.str();
        const char* headdata = tmp.c_str();
        MPI_File_write_at(fh, myoffset, headdata, tmp.size(), MPI_CHAR, &status);
    }
    MPI_Bcast(&headersize, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    std::vector<size_t> sizes;
    sizes.resize(MPI_SIZE);
    MPI_Allgather(&buffersize, 1, MPI_LONG_LONG,
                  sizes.data(), 1, MPI_LONG_LONG,
                  MPI_COMM_WORLD);
    myoffset = headersize;
    for (int i = 0; i < MPI_RANK; ++i)
    {
        myoffset += sizes[i];
    }
    const std::string tmp = buffer.str();
    const char* bufferdata = tmp.c_str();
    MPI_File_write_at_all(fh, myoffset, bufferdata, tmp.size(), MPI_CHAR, &status);
    if (MPI_RANK == 0)
    {
        myoffset = headersize;
        for (int i = 0; i < MPI_SIZE; ++i)
        {
            myoffset += sizes[i];
        }
        const std::string tmp = tbuffer.str();
        const char* taildata = tmp.c_str();
        MPI_File_write_at(fh, myoffset, taildata, tmp.size(), MPI_CHAR, &status);
    }
    MPI_File_close(&fh);
#endif
}

void VTK::WriteDistorted(
    const std::string Filename,
    const Settings& locSettings,
    const ElasticProperties& EP,
    std::vector<Field_t> ListOfFields,
    const int precision,
    const int resolution)
{
    const int Nx = (resolution == 1) ? locSettings.Nx : (locSettings.dNx + 1)*locSettings.Nx;
    const int Ny = (resolution == 1) ? locSettings.Ny : (locSettings.dNy + 1)*locSettings.Ny;
    const int Nz = (resolution == 1) ? locSettings.Nz : (locSettings.dNz + 1)*locSettings.Nz;

#ifndef MPI_PARALLEL
    ofstream vtk_file(Filename.c_str());
    VTK::WriteHeader(vtk_file, Nx, Ny, Nz);
    {
        WritePointData(vtk_file,ListOfFields,Nx,Ny,Nz,precision);
    }
    VTK::WriteEndPointData(vtk_file);
    VTK::WriteCoordinatesDistorted(vtk_file,EP,locSettings,resolution);
    VTK::CloseFile(vtk_file);
#else
    const int TotalNx = (resolution == 1) ? locSettings.TotalNx : (locSettings.dNx + 1)*locSettings.TotalNx;
    const int TotalNy = (resolution == 1) ? locSettings.TotalNy : (locSettings.dNy + 1)*locSettings.TotalNy;
    const int TotalNz = (resolution == 1) ? locSettings.TotalNz : (locSettings.dNz + 1)*locSettings.TotalNz;

    const int OffsetX = (resolution == 1) ? locSettings.OffsetX : (locSettings.dNx + 1)*locSettings.OffsetX;
    const int OffsetY = (resolution == 1) ? locSettings.OffsetY : (locSettings.dNy + 1)*locSettings.OffsetY;
    const int OffsetZ = (resolution == 1) ? locSettings.OffsetZ : (locSettings.dNz + 1)*locSettings.OffsetZ;

    MPI_File fh;
    MPI_Status status;
    std::stringstream buffer;
    std::stringstream hbuffer;
    std::stringstream tbuffer;
    MPI_File_delete(Filename.c_str(), MPI_INFO_NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_open(MPI_COMM_WORLD, Filename.c_str(),
        MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_Offset myoffset = 0;
    MPI_File_set_view(fh, myoffset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);

    buffer << "<Piece Extent=\""
           << OffsetX << " " << Nx-1 + OffsetX << " "
           << OffsetY << " " << Ny-1 + OffsetY << " "
           << OffsetZ << " " << Nz-1 + OffsetZ << "\">\n";
    {
        WritePointData(buffer,ListOfFields,Nx,Ny,Nz,precision);
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinatesDistorted(buffer,EP,locSettings,resolution);
    buffer << "</Piece>\n";

    MPI_Barrier(MPI_COMM_WORLD);
    hbuffer << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";
    hbuffer << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    hbuffer << "<StructuredGrid WholeExtent=\""
            << 0 << " " << TotalNx-1 << " "
            << 0 << " " << TotalNy-1 << " "
            << 0 << " " << TotalNz-1 << "\"> \n";
    tbuffer << "</StructuredGrid> \n";
    tbuffer << "</VTKFile> \n";
    size_t buffersize = buffer.str().size();
    size_t headersize = hbuffer.str().size();
    if (MPI_RANK == 0)
    {
        const std::string tmp = hbuffer.str();
        const char* headdata = tmp.c_str();
        MPI_File_write_at(fh, myoffset, headdata, tmp.size(), MPI_CHAR, &status);
    }
    MPI_Bcast(&headersize, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    std::vector<size_t> sizes;
    sizes.resize(MPI_SIZE);
    MPI_Allgather(&buffersize, 1, MPI_LONG_LONG,
                  sizes.data(), 1, MPI_LONG_LONG,
                  MPI_COMM_WORLD);
    myoffset = headersize;
    for (int i = 0; i < MPI_RANK; ++i)
    {
        myoffset += sizes[i];
    }
    const std::string tmp = buffer.str();
    const char* bufferdata = tmp.c_str();
    MPI_File_write_at_all(fh, myoffset, bufferdata, tmp.size(), MPI_CHAR, &status);
    if (MPI_RANK == 0)
    {
        myoffset = headersize;
        for (int i = 0; i < MPI_SIZE; ++i)
        {
            myoffset += sizes[i];
        }
        const std::string tmp = tbuffer.str();
        const char* taildata = tmp.c_str();
        MPI_File_write_at(fh, myoffset, taildata, tmp.size(), MPI_CHAR, &status);
    }
    MPI_File_close(&fh);
#endif
}

}// namespace openphase
