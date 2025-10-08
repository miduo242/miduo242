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

#ifndef VTK_H
#define VTK_H

#include "Base/Includes.h"
#include "Settings.h"
#include "Mechanics/ElasticProperties.h"
#include "Tools.h"

namespace openphase
{

enum class VTKDataTypes                                                         ///< VTK data types
{
    PDScalars,                                                                  ///< Scalar
    PDVectors,                                                                  ///< Vector
    PDTensors                                                                   ///< Matrix/Tensor
};

typedef struct
{
    uint32_t nblocks;
    uint32_t blocksize;
    uint32_t lastblocksize;
    uint32_t compressedsize;
} header_singleblock;

struct OP_EXPORTS VTK                                                           ///< static class to write VTK data in xml format
{
 public:
    struct Field_t                                                              ///< Short hand for a field's name and function which will be written to file
    {
        std::string Name;
        std::function<std::any(const int, const int, const int)> Function;

        Field_t(const std::string& inp_Name,
                const std::function<std::any(const int, const int, const int)>& inp_Function) :
            Name(inp_Name), Function(inp_Function) {};
    };

    template <typename T, class function_t>
    static void ForEach (T& buffer,
            const long int Nx, const  long int Ny, const long int Nz,
            function_t Function)
    {
        for(long int k = 0; k < Nz; ++k)
        for(long int j = 0; j < Ny; ++j)
        for(long int i = 0; i < Nx; ++i)
        {
            buffer << Function(i,j,k) << "\n";
        }
        buffer << "</DataArray>\n";
    }

    template <typename T>
    static void WriteFieldHeader(
            T& buffer,
            std::string Name,
            std::string Type,
            size_t NComponents = 1,
            std::string format = "ascii")
    {
        buffer << "<DataArray type = \"" << Type
               << "\" Name = \"" << Name
               << "\" NumberOfComponents=\"" << NComponents
               << "\" format=\""<< format << "\">\n";
    }

    template <typename T, class container_t>
    static void WritePointData(
            T& buffer,
            container_t ListOfFields,
            const long int Nx, const  long int Ny, const long int Nz,
            const int precision = 16)
    {
        // Use type with highest dimensions to write the data to file
        bool check = true;
        for (auto Field : ListOfFields)
        if (Field.Function(0,0,0).type()==typeid(dMatrix3x3) or
            Field.Function(0,0,0).type()==typeid(dMatrix6x6))
        {
            buffer << "<PointData Tensors= \"TensorData\">\n";
            check = false;
            break;
        }
        if(check)
        for (auto Field : ListOfFields)
        if (Field.Function(0,0,0).type()==typeid(dVector3) or
            Field.Function(0,0,0).type()==typeid(dVector6) or
            Field.Function(0,0,0).type()==typeid(vStress) or
            Field.Function(0,0,0).type()==typeid(vStrain))
        {
            buffer << "<PointData Vectors= \"VectorData\">\n";
            check = false;
            break;
        }
        if(check)
        {
            buffer << "<PointData Scalars= \"ScalarData\">\n";
        }
        for (auto Field : ListOfFields)
        {
            if (Field.Function(0,0,0).type() == typeid(int))
            {
                buffer << std::fixed;
                buffer << std::setprecision(0);
                WriteFieldHeader(buffer, Field.Name, "Int32");
                ForEach(buffer,Nx,Ny,Nz,[&Field](int i,int j,int k){return std::any_cast<int>(Field.Function(i,j,k));});
            }
            else if (Field.Function(0,0,0).type() == typeid(size_t))
            {
                buffer << std::fixed;
                buffer << std::setprecision(0);
                WriteFieldHeader(buffer, Field.Name, "UInt64");
                ForEach(buffer,Nx,Ny,Nz,[&Field](int i,int j,int k){return std::any_cast<size_t>(Field.Function(i,j,k));});
            }
            else if (Field.Function(0,0,0).type() == typeid(double))
            {
                //out << std::scientific; //NOTE: this results in unnecessarily large files
                buffer << std::defaultfloat;
                buffer << std::setprecision(precision);
                WriteFieldHeader(buffer, Field.Name, "Float64");
                ForEach(buffer,Nx,Ny,Nz,[&Field](int i,int j,int k){return std::any_cast<double>(Field.Function(i,j,k));});
            }
            else if (Field.Function(0,0,0).type() == typeid(dVector3))
            {
                WriteFieldHeader(buffer, Field.Name, "Float64", 3);
                ForEach(buffer,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<dVector3>(Field.Function(i,j,k)).write(precision);});
            }
            else if (Field.Function(0,0,0).type() == typeid(dVector6))
            {
                WriteFieldHeader(buffer, Field.Name, "Float64", 6);
                ForEach(buffer,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<dVector6>(Field.Function(i,j,k)).write(precision);});
            }
            else if (Field.Function(0,0,0).type() == typeid(vStrain))
            {
                WriteFieldHeader(buffer, Field.Name, "Float64", 6);
                ForEach(buffer,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<vStrain>(Field.Function(i,j,k)).write(precision);});
            }
            else if (Field.Function(0,0,0).type() == typeid(vStress))
            {
                WriteFieldHeader(buffer, Field.Name, "Float64", 6);
                ForEach(buffer,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<vStress>(Field.Function(i,j,k)).write(precision);});
            }
            else if (Field.Function(0,0,0).type() == typeid(dMatrix3x3))
            {
                WriteFieldHeader(buffer, Field.Name, "Float64", 9);
                ForEach(buffer,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<dMatrix3x3>(Field.Function(i,j,k)).write(precision);});
            }
            else if (Field.Function(0,0,0).type() == typeid(dMatrix6x6))
            {
                WriteFieldHeader(buffer, Field.Name, "Float64", 36);
                ForEach(buffer,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<dMatrix6x6>(Field.Function(i,j,k)).write(precision);});
            }
        }
    }

    /// Writes a list of fields which returns a double into a VTK file
    static void Write(
            const std::string Filename,
            const Settings& locSettings,
            std::vector<Field_t> ListOfFields,
            const int precision = 16,
            const int resolution = 1);

    static void WriteDistorted(
            const std::string Filename,
            const Settings& locSettings,
            const ElasticProperties& EP,
            std::vector<Field_t> ListOfFields,
            const int precision = 16,
            const int resolution = 1);
    /// Short for a list of undefined data
    //typedef std::vector<std::function<void(std::stringstream& buffer)>> ListOfData_t;

    ///// Writes VTK data for file WriteVTKData function is provided
    //static void WriteData(
    //        const std::string Filename,
    //        const Settings& locSettings,
    //        ListOfData_t ListOfData);


    template <typename T>
    static void WriteHeader(T& buffer, const int Nx, const int Ny, const int Nz, int resolution = 1)
    {
        buffer << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";
        buffer << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        buffer << "<StructuredGrid WholeExtent=\""
                 << 0 << " " << resolution * Nx - 1 << " "
                 << 0 << " " << resolution * Ny - 1 << " "
                 << 0 << " " << resolution * Nz - 1 << "\">\n";
    }

    template <typename T>
    static void WriteBeginPointData(T& buffer, const std::vector<VTKDataTypes> PointDataTypes)
    {
        // Method requires vector that contains the types of the given point data.
        // Each type has to be given once.

        // Example vector initialization for an output with scalar and tensor data:
        // ' std::vector<int> DataTypes {PDScalars, PDTensors}; '

        buffer << "<PointData ";

        for(auto it = PointDataTypes.cbegin(); it != PointDataTypes.cend(); ++it)
        {
            switch(*it)
            {
                case VTKDataTypes::PDScalars:
                {
                    buffer << " Scalars= \"ScalarData\"";
                    break;
                }
                case VTKDataTypes::PDVectors:
                {
                    buffer << " Vectors= \"VectorData\"";
                    break;
                }
                case VTKDataTypes::PDTensors:
                {
                    buffer << " Tensors= \"TensorData\"";
                    break;
                }
                default:
                {
                    break;
                }
            }
        }
        buffer << ">\n";
    }

    template <typename T>
    static void WriteEndPointData(T& buffer)
    {
        buffer << "</PointData>\n";
    }

    template <typename T>
    static void WriteCoordinates(T& buffer, const Settings& locSettings, int resolution = 1)
    {
        int TotalNx = locSettings.TotalNx;
        int TotalNy = locSettings.TotalNy;
        int TotalNz = locSettings.TotalNz;

        int Nx = (resolution == 1) ? locSettings.Nx : (locSettings.dNx + 1)*locSettings.Nx;
        int Ny = (resolution == 1) ? locSettings.Ny : (locSettings.dNy + 1)*locSettings.Ny;
        int Nz = (resolution == 1) ? locSettings.Nz : (locSettings.dNz + 1)*locSettings.Nz;

        int OffsetX = (resolution == 1) ? locSettings.OffsetX : (locSettings.dNx + 1)*locSettings.OffsetX;
        int OffsetY = (resolution == 1) ? locSettings.OffsetY : (locSettings.dNy + 1)*locSettings.OffsetY;
        int OffsetZ = (resolution == 1) ? locSettings.OffsetZ : (locSettings.dNz + 1)*locSettings.OffsetZ;

        buffer << "<Points>\n";
        buffer << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        buffer << std::fixed;
        double a = 1.0;
        double b = 1.0;
        double c = 1.0;
        if(resolution == 2)
        {
            a = (TotalNx-1) ? (TotalNx*0.5 - 0.25)/(TotalNx) : 0;
            b = (TotalNy-1) ? (TotalNy*0.5 - 0.25)/(TotalNy) : 0;
            c = (TotalNz-1) ? (TotalNz*0.5 - 0.25)/(TotalNz) : 0;
            buffer << std::setprecision(1);
        }
        else
        {
            buffer << std::setprecision(0);
        }
        for(int k = 0; k < Nz; ++k)
        for(int j = 0; j < Ny; ++j)
        for(int i = 0; i < Nx; ++i)
        {
            buffer << (i + OffsetX)*a << " " << (j + OffsetY)*b << " " << (k + OffsetZ)*c << "\n";
        }
        buffer << "</DataArray>\n";
        buffer << "</Points>\n";
    }

    template <typename T>
    static void WriteCoordinatesDistorted(T& buffer, const  ElasticProperties& EP, const Settings& locSettings, int resolution = 1)
    {
        dMatrix3x3 locDefGrad = EP.AverageDeformationGradient;
        //dMatrix3x3 locDefGrad = Tools::AlignBaseAxes(EP.AverageDeformationGradient);
        int TotalNx = locSettings.TotalNx;
        int TotalNy = locSettings.TotalNy;
        int TotalNz = locSettings.TotalNz;

        int Nx = (resolution == 1) ? locSettings.Nx : (locSettings.dNx + 1)*locSettings.Nx;
        int Ny = (resolution == 1) ? locSettings.Ny : (locSettings.dNy + 1)*locSettings.Ny;
        int Nz = (resolution == 1) ? locSettings.Nz : (locSettings.dNz + 1)*locSettings.Nz;

        int OffsetX = (resolution == 1) ? locSettings.OffsetX : (locSettings.dNx + 1)*locSettings.OffsetX;
        int OffsetY = (resolution == 1) ? locSettings.OffsetY : (locSettings.dNy + 1)*locSettings.OffsetY;
        int OffsetZ = (resolution == 1) ? locSettings.OffsetZ : (locSettings.dNz + 1)*locSettings.OffsetZ;

        buffer << "<Points>\n";
        buffer << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        buffer << std::fixed;

        double a = 1.0;
        double b = 1.0;
        double c = 1.0;
        if(resolution == 2)
        {
            a = (TotalNx-1) ? (TotalNx*0.5 - 0.25)/(TotalNx) : 0;
            b = (TotalNy-1) ? (TotalNy*0.5 - 0.25)/(TotalNy) : 0;
            c = (TotalNz-1) ? (TotalNz*0.5 - 0.25)/(TotalNz) : 0;
        }

        for(int k = 0; k < Nz; ++k)
        for(int j = 0; j < Ny; ++j)
        for(int i = 0; i < Nx; ++i)
        {
            double x = i * a;
            double y = j * b;
            double z = k * c;
            dVector3 coordinates {(-0.5*TotalNx + x + a*OffsetX),
                                  (-0.5*TotalNy + y + b*OffsetY),
                                  (-0.5*TotalNz + z + c*OffsetZ)};
            int l = int(round(x));
            int m = int(round(y));
            int n = int(round(z));
            coordinates = locDefGrad*coordinates + EP.Displacements(l,m,n)*(1.0/EP.dx);
            buffer << std::setprecision(4);
            buffer << 0.5*TotalNx + coordinates[0] << " "
                   << 0.5*TotalNy + coordinates[1] << " "
                   << 0.5*TotalNz + coordinates[2] << "\n";
        }
        buffer << "</DataArray>\n";
        buffer << "</Points>\n";
    }

    template <typename T>
    static void WriteToFile(T& buffer, const std::string Filename)
    {
        buffer << "</StructuredGrid>\n";
        buffer << "</VTKFile>\n";

        std::ofstream vtk_file(Filename.c_str());
        vtk_file << buffer.rdbuf();
        vtk_file.close();
    }
    template <typename T>
    static void CloseFile(T& buffer)
    {
        buffer << "</StructuredGrid>\n";
        buffer << "</VTKFile>\n";
        buffer.close();
    }
};
}// namespace openphase
#endif
