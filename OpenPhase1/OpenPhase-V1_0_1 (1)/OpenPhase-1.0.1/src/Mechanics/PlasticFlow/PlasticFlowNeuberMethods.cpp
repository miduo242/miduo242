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
 *   File created :   2018
 *   Main contributors :   Johannes Goerler; Oleg Shchyglo
 *
 */

#include "Mechanics/PlasticFlow/PlasticFlowNeuberMethods.h"
#include "Mechanics/ElasticProperties.h"
#include "VTK.h"
#include "Base/UserInterface.h"

namespace openphase
{
    using namespace std;

    void PlasticFlowNeuberMethods::getNeuberDataRO(vStress& StressIn, vStrain& StrainIn, vStress& StressOut, vStrain& StrainOut, double YoungModulus, double YieldStress, double alpha)
    {
        for(int n = 0; n < 6; n++)
        {
            double strain_sign = 1;
            if(StrainIn[n] < 0.0) strain_sign = -1;
            double stress_sign = 1;
            if(StressIn[n] < 0.0) stress_sign = -1;
            double energy = fabs(StressIn[n] * StrainIn[n]);
            double newStress = sqrt(fabs(-0.5*YieldStress*YieldStress/alpha + 0.5*sqrt((4.0*YoungModulus*energy*alpha*YieldStress*YieldStress + pow(YieldStress,4))/pow(alpha,2))));
            double newStrain = newStress/YoungModulus*(1.0 + alpha*pow(newStress/YieldStress,2));
            StressOut[n] = newStress*stress_sign;
            StrainOut[n] = newStrain*strain_sign;
        }
    }
    vStrain PlasticFlowNeuberMethods::getNeuberStrains(vStrain StrainIn)
    {
        return StrainIn;            //Put Neuber equations here
    }

    vStress PlasticFlowNeuberMethods::getNeuberStresses(vStress StressIn)
    {
        return StressIn;            //Put Neuber equations here
    }

    void PlasticFlowNeuberMethods::writeNeuberStressVTK(const Settings& locSettings, ElasticProperties& EP, const int tStep)
    {
        stringstream buffer;
        vector<VTKDataTypes> DataTypes{VTKDataTypes::PDScalars};

        VTK::WriteHeader(buffer, EP.Nx, EP.Ny, EP.Nz);
        VTK::WriteBeginPointData(buffer, DataTypes);
        {
            PlasticFlowNeuberMethods::WriteStressesVTKData(EP, buffer);
        }
        VTK::WriteEndPointData(buffer);
        VTK::WriteCoordinates(buffer, locSettings);
        string FileName = UserInterface::MakeFileName(DefaultVTKDir, "NeuberStresses_", tStep, ".vts");
        VTK::WriteToFile(buffer, FileName);
    }

    void PlasticFlowNeuberMethods::writeNeuberStrainVTK(const Settings& locSettings, ElasticProperties& EP, const int tStep)
    {
        stringstream buffer;
        vector<VTKDataTypes> DataTypes{VTKDataTypes::PDScalars};
        VTK::WriteHeader(buffer, EP.Nx, EP.Ny, EP.Nz);
        VTK::WriteBeginPointData(buffer, DataTypes);
        {
            PlasticFlowNeuberMethods::WriteStrainsVTKData(EP, buffer);
        }
        VTK::WriteEndPointData(buffer);
        VTK::WriteCoordinates(buffer, locSettings);
        string FileName = UserInterface::MakeFileName(DefaultVTKDir, "NeuberStrains_", tStep, ".vts");
        VTK::WriteToFile(buffer, FileName);
    }

    void PlasticFlowNeuberMethods::WriteStressesVTKData(ElasticProperties& EP, stringstream& buffer)
    {
        vector<long int> dimV{EP.Stresses.sizeX(), EP.Stresses.sizeY(), EP.Stresses.sizeZ()};
        vector<int> compV{ 0, 1, 2, 3, 4, 5 };
        vector<string> compNameV{ "xx", "yy", "zz", "yz", "xz", "xy" };

        for (auto it = compV.cbegin(); it != compV.cend(); ++it)
        {
            string compname = "\"NeuberStress_" + compNameV[*it] + "\" ";
            buffer << "<DataArray type = \"Float64\" Name = " << compname <<
                "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
            for (int k = 0; k < EP.Stresses.sizeZ(); ++k)
            for (int j = 0; j < EP.Stresses.sizeY(); ++j)
            for (int i = 0; i < EP.Stresses.sizeX(); ++i)
            {
                buffer << getNeuberStresses(EP.Stresses(i, j, k))[*it] << "\n";
            }
            buffer << "</DataArray>" << endl;
        }

        buffer << "<DataArray type = \"Float64\" Name = \"" << "Pressure" <<
            "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for (int k = 0; k < EP.Stresses.sizeZ(); ++k)
        for (int j = 0; j < EP.Stresses.sizeY(); ++j)
        for (int i = 0; i < EP.Stresses.sizeX(); ++i)
        {
            buffer << getNeuberStresses(EP.Stresses(i, j, k)).Pressure() << endl;
        }
        buffer << "</DataArray>" << endl;

        buffer << "<DataArray type = \"Float64\" Name = \"" << "vMises" <<
            "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for (int k = 0; k < EP.Stresses.sizeZ(); ++k)
        for (int j = 0; j < EP.Stresses.sizeY(); ++j)
        for (int i = 0; i < EP.Stresses.sizeX(); ++i)
        {
            buffer << getNeuberStresses(EP.Stresses(i, j, k)).Mises() << endl;
        }
        buffer << "</DataArray>" << endl;
    }

    void PlasticFlowNeuberMethods::WriteStrainsVTKData(ElasticProperties& EP, stringstream& buffer)
    {
        vector<int> compV{ 0, 1, 2, 3, 4, 5 };
        vector<string> compNameV{ "1", "2", "3", "4", "5", "6" };

        for (auto it = compV.cbegin(); it != compV.cend(); ++it)
        {
            buffer << "<DataArray type = \"Float64\" Name = \""
                << "E_" << compNameV[*it]
                << "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
            for (int k = 0; k < EP.Nz; ++k)
            for (int j = 0; j < EP.Ny; ++j)
            for (int i = 0; i < EP.Nx; ++i)
            {
                buffer << getNeuberStrains(EP.TotalStrains(i, j, k))[*it] << endl;
            }
            buffer << "</DataArray>" << endl;
        }
    }
}
