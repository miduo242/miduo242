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
 *   Main contributors :   Oleg Shchyglo; Hesham Salama
 *
 */

#include "Tools/MicrostructureAnalysis.h"
#include "Info.h"
#include "PhaseField.h"
#include "Mechanics/SymmetryVariants.h"
#include <map>
#include "Tools.h"
#include "Base/UserInterface.h"



namespace openphase
{
using namespace std;

double MicrostructureAnalysis::GrainBoundaryStatistics(PhaseField& Phase, std::vector<dVector3> Facets, double DegreeTolerance)
{
    int totalInterface = 0;
    int totalFacet = 0;
    for(size_t n = 0; n < Facets.size(); n++)
    {
        Facets[n].normalize();
    }

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    {
        if(Phase.Interface(i,j,k))
        {
            NodeV3 locNormals = Phase.Normals(i,j,k);
            if(locNormals.size() == 1)
            for(auto alpha = locNormals.begin();
                     alpha != locNormals.end(); alpha++)
            {
                totalInterface++;

                for(size_t n = 0; n < Facets.size(); n++)
                {
                    dVector3 abNorm;
                    abNorm = alpha->vector3;
                    abNorm.normalize();
                    double cosTheta = Facets[n]*abNorm * 0.9999999;
                    double angle = acos(cosTheta) * 180.0/Pi;
                    if(fabs(angle) < DegreeTolerance and abNorm.abs() > 1.0 - 0.0000001 and abNorm.abs() < 1.0 + 0.0000001)
                    {
                        totalFacet++;
                    }
                }
            }
        }
    }
    STORAGE_LOOP_END
    if (totalInterface > 0)
    {
        cout << " (" << totalFacet << "/" << totalInterface << ") = ";
        return (double)totalFacet/(double)totalInterface;
    }
    else
    {
        return 0;
    }
}

void MicrostructureAnalysis::WriteEBSDDataQuaternions(const PhaseField& Phase, const int tStep, double scale)
{
    stringstream outbuffer;
    outbuffer << "Index\t"
              << "Phase\t"
              << "X\t"
              << "Y\t"
              << "Z\t"
              << "Quat_real\t"
              << "Quat_i\t"
              << "Quat_j\t"
              << "Quat_k\t" << endl;

    int index = 1;
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    {
        Quaternion  tempQuat;

        size_t pIndex = 0;
        // selecting the quaternion of the majority phase field.
        double value = 0.0;
        for(auto alpha = Phase.Fields(i, j, k).cbegin();
                 alpha != Phase.Fields(i, j, k).cend(); ++alpha)
        if(alpha->value > value)
        {
            tempQuat = Phase.FieldsStatistics[alpha->index].Orientation;
            value = alpha->value;
            pIndex = Phase.FieldsStatistics[alpha->index].Phase;
        }

        tempQuat.normalize();
        outbuffer << index       << "\t"
                  << pIndex + 1  << "\t"
                  << i*scale     << "\t"
                  << j*scale     << "\t"
                  << k*scale     << "\t"
                  << tempQuat[0] << "\t"
                  << tempQuat[1] << "\t"
                  << tempQuat[2] << "\t"
                  << tempQuat[3] << endl;
        index++;
    }
    STORAGE_LOOP_END

    string FileName = UserInterface::MakeFileName(DefaultRawDataDir, "EBSD_", tStep, ".dat");

    ofstream ebsd_file(FileName.c_str());
    ebsd_file << outbuffer.rdbuf();
    ebsd_file.close();
}

void MicrostructureAnalysis::WriteEBSDDataQuaternionsSlice(const PhaseField& Phase, const int tStep,
                            const char Axis, const int Position, double scale)
{
    stringstream outbuffer;
    outbuffer << "Index\t"
              << "Phase\t";
    if(Axis != 'X')
    {
        outbuffer << "X\t";
    }
    if(Axis != 'Y')
    {
        outbuffer << "Y\t";
    }
    if(Axis != 'Z')
    {
        outbuffer << "Z\t";
    }

    outbuffer << "Quat_real\t"
              << "Quat_i\t"
              << "Quat_j\t"
              << "Quat_k\t" << endl;

    int index = 1;
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    {
        if((Axis == 'X' and i == Position) or
           (Axis == 'Y' and j == Position) or
           (Axis == 'Z' and k == Position))
        {
            Quaternion  tempQuat;
            size_t pIndex = 0;

            // selecting the quaternion of the majority phase field.
            double value = 0.0;
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha != Phase.Fields(i, j, k).cend(); ++alpha)
            if(alpha->value > value)
            {
                pIndex = Phase.FieldsStatistics[alpha->index].Phase;
                tempQuat = Phase.FieldsStatistics[alpha->index].Orientation;
                value = alpha->value;
            }

            tempQuat.normalize();

            outbuffer << index       << "\t"
                      << pIndex + 1  << "\t";
            if(Axis != 'X')
            {
                outbuffer << i*scale << "\t";
            }
            if(Axis != 'Y')
            {
                outbuffer << j*scale << "\t";
            }
            if(Axis != 'Z')
            {
                outbuffer << k*scale << "\t";
            }
            outbuffer << tempQuat[0] << "\t"
                      << tempQuat[1] << "\t"
                      << tempQuat[2] << "\t"
                      << tempQuat[3] << endl;
            index++;
        }
    }
    STORAGE_LOOP_END
    stringstream sliceInd;
    sliceInd << Axis << "-" << Position << "_";
    string FileName = UserInterface::MakeFileName(DefaultRawDataDir, "EBSD_" + sliceInd.str(), tStep, ".dat");

    ofstream ebsd_file(FileName.c_str());
    ebsd_file << outbuffer.rdbuf();
    ebsd_file.close();
}

void MicrostructureAnalysis::WriteEBSDDataQuaternions(const PhaseField& Phase, const SymmetryVariants& SV, const int tStep, double scale)
{
    stringstream outbuffer;
    outbuffer << "Index\t"
              << "Phase\t"
              << "X\t"
              << "Y\t"
              << "Z\t"
              << "Quat_real\t"
              << "Quat_i\t"
              << "Quat_j\t"
              << "Quat_k\t" << endl;

    int index = 1;
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    {
        Quaternion  tempQuat;
        int pIndex = 0;

        // selecting the quaternion of the majority phase field.
        double value = 0.0;
        for(auto alpha = Phase.Fields(i, j, k).cbegin();
                 alpha != Phase.Fields(i, j, k).cend(); ++alpha)
        if(alpha->value > value)
        {
            pIndex = Phase.FieldsStatistics[alpha->index].Phase;
            int variant = Phase.FieldsStatistics[alpha->index].Variant;
            Quaternion locQ;
            dMatrix3x3 locRM = SV(pIndex,variant);
            locQ.set(locRM);
            tempQuat = Phase.FieldsStatistics[alpha->index].Orientation + locQ;
            value = alpha->value;
        }

        tempQuat.normalize();

        outbuffer << index       << "\t"
                  << pIndex + 1  << "\t"
                  << i*scale     << "\t"
                  << j*scale     << "\t"
                  << k*scale     << "\t"
                  << tempQuat[0] << "\t"
                  << tempQuat[1] << "\t"
                  << tempQuat[2] << "\t"
                  << tempQuat[3] << endl;
        index++;
    }
    STORAGE_LOOP_END

    string FileName = UserInterface::MakeFileName(DefaultRawDataDir, "EBSD_", tStep, ".dat");

    ofstream ebsd_file(FileName.c_str());
    ebsd_file << outbuffer.rdbuf();
    ebsd_file.close();
}

void MicrostructureAnalysis::WriteEBSDDataQuaternionsSlice(const PhaseField& Phase, const SymmetryVariants& SV, const int tStep,
                            const char Axis, const int Position, double scale)
{
    stringstream outbuffer;
    outbuffer << "Index\t"
              << "Phase\t";
    if(Axis != 'X')
    {
        outbuffer << "X\t";
    }
    if(Axis != 'Y')
    {
        outbuffer << "Y\t";
    }
    if(Axis != 'Z')
    {
        outbuffer << "Z\t";
    }

    outbuffer << "Quat_real\t"
              << "Quat_i\t"
              << "Quat_j\t"
              << "Quat_k\t" << endl;

    int index = 1;
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    {
        if((Axis == 'X' and i == Position) or
           (Axis == 'Y' and j == Position) or
           (Axis == 'Z' and k == Position))
        {
            Quaternion  tempQuat;
            size_t pIndex = 0;

            // selecting the quaternion of the majority phase field.
            double value = 0.0;
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha != Phase.Fields(i, j, k).cend(); ++alpha)
            if(alpha->value > value)
            {
                pIndex = Phase.FieldsStatistics[alpha->index].Phase;
                size_t variant = Phase.FieldsStatistics[alpha->index].Variant;
                Quaternion locQ;
                dMatrix3x3 locRM = SV(pIndex,variant);
                locQ.set(locRM);
                tempQuat = Phase.FieldsStatistics[alpha->index].Orientation + locQ;
                value = alpha->value;
            }

            tempQuat.normalize();

            outbuffer << index       << "\t"
                      << pIndex + 1  << "\t";
            if(Axis != 'X')
            {
                outbuffer << i*scale << "\t";
            }
            if(Axis != 'Y')
            {
                outbuffer << j*scale << "\t";
            }
            if(Axis != 'Z')
            {
                outbuffer << k*scale << "\t";
            }
            outbuffer << tempQuat[0] << "\t"
                      << tempQuat[1] << "\t"
                      << tempQuat[2] << "\t"
                      << tempQuat[3] << endl;
            index++;
        }
    }
    STORAGE_LOOP_END
    stringstream sliceInd;
    sliceInd << Axis << "-" << Position << "_";
    string FileName = UserInterface::MakeFileName(DefaultRawDataDir, "EBSD_" + sliceInd.str(), tStep, ".dat");

    ofstream ebsd_file(FileName.c_str());
    ebsd_file << outbuffer.rdbuf();
    ebsd_file.close();
}


// Reference sample direction (sd) input as an integer between 1 and 3. Options are [100], [010], and [001], respectively.
void MicrostructureAnalysis::WriteEBSDVTK(PhaseField& Phase, Crystallography& CR, const Settings& locSettings, EulerConvention locConvention, const int sd, const int tStep)
{
	Storage3D<dVector3, 0> tempRGB;
	tempRGB.Allocate(Phase.Nx, Phase.Ny, Phase.Nz, Phase.dNx, Phase.dNy, Phase.dNz, 1);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        EulerAngles tempAng;
        Quaternion tempQuat;

        size_t locPF = 0;
        double locVal = 0.0;
        for (auto beta = Phase.Fields(i,j,k).cbegin();
        		beta != Phase.Fields(i,j,k).cend();  ++beta)
        {
        	if(beta->value > locVal)
        	{
        		locVal = beta->value;
        		locPF = beta->index;
        	}
        }
        tempQuat =  Phase.FieldsStatistics[locPF].Orientation;
        tempAng.set(tempQuat, locConvention, false);

        dVector3 RGB;
        RGB = Tools::IPFColor(sd, tempAng, CR);
        tempRGB(i,j,k) = RGB;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
	WriteVTKRGB(tempRGB, locSettings, "sRGB", tStep);
}

void MicrostructureAnalysis::WriteSymmetryVariantsStatistics(size_t pIndex, const PhaseField& Phase, const SymmetryVariants& SV, const int tStep)
{
    string FileName = "SymmetryVariantsStatistics.dat";

    if(tStep == 0)
    {
        stringstream header;
        header << std::left << std::setw(12) << "tStep";
        for(size_t n = 0; n < SV.Nvariants(pIndex); n++)
        {
            stringstream variantN;
            variantN << "V" << n;
            header << std::right << std::setw(6) << variantN.str();
        }
        header << std::right << std::setw(8) << "Total" << endl;

        fstream out_file(FileName.c_str(), ios::out);
        out_file << header.rdbuf();
        out_file.close();
    }

    vector<int> Npfs(SV.Nvariants(pIndex),0);
    vector<double> Volumes(SV.Nvariants(pIndex),0.0);

    for(size_t n = 0; n < Phase.FieldsStatistics.size(); n++)
    if(Phase.FieldsStatistics[n].Exist and Phase.FieldsStatistics[n].Phase == pIndex)
    {
        size_t variant = Phase.FieldsStatistics[n].Variant;
        Npfs[variant] += 1.0;
        Volumes[variant] += Phase.FieldsStatistics[n].Volume;
    }
    stringstream outbuffer;

    outbuffer << std::left << std::setw(12) << tStep;
    int counter = 0;
    for(size_t n = 0; n < Npfs.size(); n++)
    {
        counter += Npfs[n];
        outbuffer << std::right << std::setw(6) << Npfs[n] /*<< " " << Volumes[n]*/;
    }
    outbuffer << std::right << std::setw(8) << counter << endl;

    fstream out_file(FileName.c_str(),ios::app);
    out_file << outbuffer.rdbuf();
    out_file.close();
}

void MicrostructureAnalysis::WriteGrainsStatistics(const PhaseField& Phase, const int tStep)
{
/// ++++++++++++++++++++++++++++++++++++++++++++
    size_t nPFs = Phase.FieldsStatistics.size();

    map<size_t, size_t> AllGrains;
    map <size_t,size_t> pairs;

    pairs.clear();
    AllGrains.clear();

    for (int i = 1; i < Phase.Nx+1; ++i)
    for (int j = 1; j < Phase.Ny+1; ++j)
    for (int k = 1; k < Phase.Nz+1; ++k)
    {
      if (Phase.Interface(i, j, k) && Phase.Fields(i, j, k).size() == 2)
      {
        int idx1 = Phase.Fields(i, j, k).cbegin()->index;
        int idx2 = (Phase.Fields(i, j, k).cbegin()+1)->index;
        pairs[nPFs*idx1 + idx2] = 1;
        pairs[nPFs*idx2 + idx1] = 1;
      }
      for(auto n = Phase.Fields(i, j, k).cbegin(); n < Phase.Fields(i, j, k).cend();++n)
      {
          AllGrains[n->index] += n->value;
      }
    }

    double AveSize = double(Phase.Nx*Phase.Ny*Phase.Nz)/AllGrains.size();
/// ++++++++++++++++++++++++++++++++++++++++++++ ///

    ofstream Fl1;
    ofstream Fl2;
    ofstream Fl3;
    if(tStep)
    {
        Fl1.open("SizeAveInfo.dat", ios::app);
        Fl2.open("SizeDetails.dat", ios::app);
        Fl3.open("NeighboInfo.dat", ios::app);
    }
    else
    {
        Fl1.open("SizeAveInfo.dat", ios::out);
        Fl2.open("SizeDetails.dat", ios::out);
        Fl3.open("NeighboInfo.dat", ios::out);
    }
    Fl1.precision(10);
    Fl2.precision(10);
    /// ++++++++++++++++++++++++++++++++++++++++++++
    Fl1 << tStep << " " << AllGrains.size() << " " << AveSize << " " << AveSize*Phase.dx*Phase.dx*Phase.dx << endl;
    Fl1.close();
    /// ++++++++++++++++++++++++++++++++++++++++++++
    Fl2 << tStep << " " << Phase.FieldsStatistics.size() << " ";
    for (size_t i = 0; i < nPFs; i++)
    {
        Fl2 << Phase.FieldsStatistics[i].Volume << " ";
    }
    Fl2 << endl;
    Fl2.close();
    /// ++++++++++++++++++++++++++++++++++++++++++++
    Fl3 << tStep << " " << Phase.FieldsStatistics.size() << " ";
    vector<size_t> sumpairs(nPFs);
    sumpairs.assign(nPFs, 0);
    for(size_t i = 0; i < nPFs; ++i)
    {
        for(size_t j = 0; j < nPFs; ++j)
        {
            if ((i != j) && (pairs.find(i*nPFs+j) != pairs.end())) sumpairs[i] += pairs.find(i*nPFs+j)->second;
        }
        Fl3 << sumpairs[i] << " ";
    }
    Fl3 << endl;
    Fl3.close();
}

dVector3 MicrostructureAnalysis::FindValuePosition(PhaseField& Phi, size_t index, double value, dVector3 start_position, dVector3 direction, double tolerance)
{
    direction.normalize();
    dVector3 current_position = start_position;
    double current_value = Phi.Fields.at(start_position[0],start_position[1],start_position[2]).get_value(index);
    double step = 1.0;

    double local_difference = 1.0;

    while(current_position[0] >= 0 and current_position[0] < Phi.Nx and
          current_position[1] >= 0 and current_position[1] < Phi.Ny and
          current_position[2] >= 0 and current_position[2] < Phi.Nz and
          local_difference > tolerance)
    {
        double local_value = current_value;
        current_position += direction*step;
        current_value = Phi.Fields.at(current_position[0],current_position[1],current_position[2]).get_value(index);
        if((current_value < value and local_value > value) or
           (current_value > value and local_value < value))
        {
            current_position -= direction*step;
            step *= 0.5;
            local_difference = fabs(local_value - current_value);
            current_value = local_value;
        }
    }
    return current_position;
};

}// namespace openphase
