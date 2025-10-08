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
 *   Main contributors :   Matthias Stratmann; Johannes Goerler
 *
 */

#include "Info.h"
#include "Base/UserInterface.h"
#include "Settings.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"
#include "Mechanics/ElasticProperties.h"
#include "PhaseField.h"
#include "Temperature.h"
#include "Nucleation.h"
#include "TextOutput.h"
#include "Composition.h"

namespace openphase
{
using namespace std;

bool TextOutput::FileExists(string filename)
{
    /** This function checks the availability of the specified text file. If the
    file is nonexistent, locked or other problems prevent this function from
    opening the file, it will return the value 'false'.*/

    bool result = false;

    ifstream my_file(filename);

    if(my_file.good())
    result = true;

    return result;
}

void TextOutput::PhasePercent(Composition& Cx, PhaseField& Phi,
                                Settings& OPSettings, string filename,
                                double time)
{
    /** This function will create tabulated data on the volume percent of each
    thermodynamic phase defined in the chemical property input. Each time this
    function is called, a new row will be written in the specified file name
    for the current time step. If the file is not present it will be created,
    if the phase already exists, the new data will be appended!*/
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    {
        if (!FileExists(filename))
        {
            ofstream file(filename, ios::out);                                      //Create file and write header (or overwrite)
            file << "# time ";
            for(size_t idx = 0; idx < Cx.Nphases; idx++)
            {
                file << Cx.PhaseNames[idx] << "  ";
            }
            file << endl;
            file.close();
        }
    }
    vector<double> tPhaseFrac(Cx.Nphases);
#ifdef MPI_PARALLEL
    double totalVolume = double(Phi.TotalNx * Phi.Ny * Phi.Nz);
#else
    double totalVolume = double(Phi.Nx * Phi.Ny * Phi.Nz);
#endif
    for(size_t idx = 0; idx < Phi.FieldsStatistics.size(); idx++)
    {
        tPhaseFrac[Phi.FieldsStatistics[idx].Phase]
                += Phi.FieldsStatistics[idx].Volume;
    }
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    {
        ofstream file(filename, ios::app);                                          //Write content of the file
        file << fixed << time << "   ";

        for(size_t idx = 0; idx < Cx.Nphases; idx++)
        {
            file << scientific  <<  100.0*tPhaseFrac[idx]/totalVolume << "    ";
        }
    file << endl;
    file.close();                                                               //Close file
    }
}

void TextOutput::AverageTemp(Temperature& Tx, string filename, double time)
{
    /** This function writes the average temperature of the system, each time it
    is called in a file.*/

    if (!FileExists(filename))
    {
        ofstream file(filename, ios::out);                                      //Create file and write header (or overwrite)
        file << "# Time Temperature" << endl;
        file.close();
    }

    double average = 0.0;
    double points  = 0.0;
    STORAGE_LOOP_BEGIN(i,j,k,Tx.Tx,0)
    {
        points++;
        average += Tx.Tx(i,j,k);
    }
    STORAGE_LOOP_END
#ifdef MPI_PARALLEL
    double  tempPoints = points;
    MPI_Reduce(&tempPoints,&points,1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
    double  tempaverage = average;
    MPI_Reduce(&tempaverage,&tempaverage,1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
    if(MPI_RANK == 0)
#endif
    {
        if(points > 0.0)
            average /= points;
        ofstream file(filename, ios::app);                                          //Write content of the file
        file.precision(10);
        file << time << "    " << average << endl;
        file.close();
    }
}

void TextOutput::GrainVolumes(PhaseField& Phi, string filename, double time)
{
    /** This function will write a list of all grain volumes in a file, each
    time this function is called in a new row. The first column is the time,
    the second row beginning with a "#" marks the thermodynamic index for each
    grain.*/

    size_t nPFs = Phi.FieldsStatistics.size();

    if (!FileExists(filename))
    {
        ofstream file(filename, ios::out);                                      //Create file and write header (or overwrite)
        file << "# Volume of each grain in the order of their field-index! "
             << "First row contains the simulation time. "
             << "The next line contains the thermodynamic phase field index "
             << "of each grain for identification" << endl;
        file << "#";
        for(size_t n = 0; n < nPFs; n++)
        {
            file << " " << Phi.FieldsStatistics[n].Phase;
        }
        file << endl;
        file.close();
    }

    map<int, double> AllGrains;

    STORAGE_LOOP_BEGIN(i,j,k,Phi.Fields,0)
    {
        for(auto n = Phi.Fields(i, j, k).cbegin();
                 n < Phi.Fields(i, j, k).cend(); ++n)
        {
            AllGrains[n->index] += n->value;
        }
    }
    STORAGE_LOOP_END

    vector<double> AllGrainsSize(nPFs);
    AllGrainsSize.assign(nPFs, 0.0);
    for (map<int, double>::iterator it = AllGrains.begin();
            it != AllGrains.end(); ++it)
    {
        AllGrainsSize[it->first] = it->second;
    }

    ofstream file(filename, ios::app);                                          //Write content of the file
    file.precision(10);
    file << time << " ";

    for (size_t i = 0; i < nPFs; i++)
    {
        file << AllGrainsSize[i] << " ";
    }
    file << endl;
    file.close();

    //This function was originally taken from PhaseField::WriteGrainsStatistics
}

void TextOutput::WriteValue(double Value,string filename, double time)
{
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
#endif
    {
        if (!FileExists(filename))
        {
            ofstream file(filename, ios::out);
            file << setw(12) << left << "time"
                 << setw(20) << right << "Value"
                 << endl;
            file.close();
        }

        ofstream file(filename, ios::app);

        file << scientific << setw(12) << left << time
             << setw(20) << right << Value
             << endl;
        file.close();
    }
}

void TextOutput::WriteMultipleValues(std::vector<string> Names, std::vector<double> value, string filename, double time)
{
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
#endif
    {
        if (!FileExists(filename))
        {
            ofstream file(filename, ios::out);
            file << setw(12) << left << "time";
            for (auto& N : Names)
            {
                file << setw(20) << right << N;
            }
            file << endl;
            file.close();
        }

        ofstream file(filename, ios::app);

        file << setw(12) << left << time;
        for (auto& V : value)
        {
          file << scientific << setw(20) << right << V;
        }
        file << endl;
        file.close();
    }
}

void TextOutput::AverageStress(ElasticProperties& EP, string filename,
                                double timeOrStrain)
{
    /** This function will create tabulated data on the average stress in the
    system. Each time this function is called, a new row will be written in
    the specified file name for the current time step. If the file is not
    present it will be created, if the phase already exists, the new data
    will be appended!*/

    vStress avgStress = EP.AverageStress;

#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
#endif
    {
        if (!FileExists(filename))
        {
            ofstream file(filename, ios::out);
            file << setw(12) <<  left  << "# strain"
                    << setw(20) <<  right << "sigma_1"
                    << setw(20) <<  right << "sigma_2"
                    << setw(20) <<  right << "sigma_3"
                    << setw(20) <<  right << "sigma_4"
                    << setw(20) <<  right << "sigma_5"
                    << setw(20) <<  right << "sigma_6"
                    << setw(20) <<  right << "Pressure"
                    << setw(20) <<  right << "Mises"
                    << setw(20) <<  right << "sigma_norm";;
            file << endl;
            file.close();
        }

        ofstream file(filename, ios::app);

        file << setw(12) <<  left  << timeOrStrain
                << scientific << setw(20) <<  right << avgStress[0]
                << setw(20) <<  right << avgStress[1]
                << setw(20) <<  right << avgStress[2]
                << setw(20) <<  right << avgStress[3]
                << setw(20) <<  right << avgStress[4]
                << setw(20) <<  right << avgStress[5]
                << setw(20) <<  right << avgStress.Pressure()
                << setw(20) <<  right << avgStress.Mises()
                << setw(20) <<  right << avgStress.norm();
        file << endl;
        file.close();
    }
}

void TextOutput::AverageStrain(ElasticProperties& EP, string filename,
                                double timeOrStrain)
{
    /** This function will create tabulated data on the average strain in the
    system. Each time this function is called, a new row will be written in
    the specified file name for the current time step. If the file is not
    present it will be created, if the phase already exists, the new data
    will be appended!*/

    vStrain avgStrain = EP.AverageStrain;

#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
#endif
    {
    if (!FileExists(filename))
    {
        ofstream file(filename, ios::out);
        file << setw(12) <<  left  << "# time"
             << setw(20) <<  right << "Epsilon_0"
             << setw(20) <<  right << "Epsilon_1"
             << setw(20) <<  right << "Epsilon_2"
             << setw(20) <<  right << "Epsilon_3"
             << setw(20) <<  right << "Epsilon_4"
             << setw(20) <<  right << "Epsilon_5";
        file << endl;
        file.close();
    }

    ofstream file(filename, ios::app);

        file << setw(12) <<  left  << timeOrStrain
                << scientific << setw(20) <<  right << avgStrain[0]
                << setw(20) <<  right << avgStrain[1]
                << setw(20) <<  right << avgStrain[2]
                << setw(20) <<  right << avgStrain[3]
                << setw(20) <<  right << avgStrain[4]
                << setw(20) <<  right << avgStrain[5];
    file << endl;
    file.close();
    }
}

void TextOutput::AverageDMatrix3x3(Storage3D<dMatrix3x3, 0>& Matrix, Settings& OP, string filename,
    double timeOrStrain)
{
    /** This function will create tabulated data on the average stress in the
    system. Each time this function is called, a new row will be written in
    the specified file name for the current time step. If the file is not
    present it will be created, if the phase already exists, the new data
    will be appended!*/

    int Nx = Matrix.sizeX();
    int Ny = Matrix.sizeY();
    int Nz = Matrix.sizeZ();

    dMatrix3x3 avgMatrix;
    avgMatrix.set_to_zero();

    for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
    for (int k = 0; k < Nz; k++)
    {
        avgMatrix += Matrix(i, j, k);
    }

#ifdef MPI_PARALLEL
    dMatrix3x3 temp = avgMatrix;
    MPI_Reduce(&temp(0,0),&avgMatrix(0,0),1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
    MPI_Reduce(&temp(0,1),&avgMatrix(0,1),1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
    MPI_Reduce(&temp(0,2),&avgMatrix(0,2),1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
    MPI_Reduce(&temp(1,0),&avgMatrix(1,0),1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
    MPI_Reduce(&temp(1,1),&avgMatrix(1,1),1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
    MPI_Reduce(&temp(1,2),&avgMatrix(1,2),1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
    MPI_Reduce(&temp(2,0),&avgMatrix(2,0),1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
    MPI_Reduce(&temp(2,1),&avgMatrix(2,1),1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
    MPI_Reduce(&temp(2,2),&avgMatrix(2,2),1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
    avgMatrix /= (OP.TotalNx*Matrix.sizeY()*Matrix.sizeZ());
#else
    avgMatrix /= (Matrix.sizeX()*Matrix.sizeY()*Matrix.sizeZ());
#endif

#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
#endif
    {
        if (!FileExists(filename))
        {
            ofstream file(filename, ios::out);
            file << setw(12) <<  left  << "# time"
                    << setw(20) <<  right << "00"
                    << setw(20) <<  right << "01"
                    << setw(20) <<  right << "02"
                    << setw(20) <<  right << "10"
                    << setw(20) <<  right << "11"
                    << setw(20) <<  right << "12"
                    << setw(20) <<  right << "20"
                    << setw(20) <<  right << "21"
                    << setw(20) <<  right << "22";
            file << endl;
            file.close();
        }

        ofstream file(filename, ios::app);

        file << setw(12) <<  left  << timeOrStrain
                << scientific << setw(20) <<  right << avgMatrix(0, 0)
                << setw(20) <<  right << avgMatrix(0, 1)
                << setw(20) <<  right << avgMatrix(0, 2)
                << setw(20) <<  right << avgMatrix(1, 0)
                << setw(20) <<  right << avgMatrix(1, 1)
                << setw(20) <<  right << avgMatrix(1, 2)
                << setw(20) <<  right << avgMatrix(2, 0)
                << setw(20) <<  right << avgMatrix(2, 1)
                << setw(20) <<  right << avgMatrix(2, 2);
        file << endl;
        file.close();
    }
}

void TextOutput::AverageDouble(Storage3D<double,0>& value, Settings& OP, std::string filename,
                          double timeOrStrain)
{
    /** This function will create tabulated data on the average value in the
    system. Each time this function is called, a new row will be written in
    the specified file name for the current time step. If the file is not
    present it will be created, if the phase already exists, the new data
    will be appended!*/

    double avgValue = 0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, value, 0,reduction(+:avgValue) )   // OMP start
    {
        avgValue += value(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
#ifdef MPI_PARALLEL
    double temp = avgValue;
    MPI_Reduce(&temp,&avgValue, 1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
    avgValue /=(OP.TotalNx*value.sizeY()*value.sizeZ());
#else
    avgValue /=( value.sizeX()*value.sizeY()*value.sizeZ());
#endif
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
#endif
    {
    if (!FileExists(filename))
    {
        ofstream file(filename, ios::out);
        file << "# time  "
             << "value";
        file << endl;
        file.close();
    }

    ofstream file(filename, ios::app);

    file << scientific << timeOrStrain << "   "
         << avgValue;
    file << endl;
    file.close();
    }
}

void TextOutput::maxElasticRotation(ElasticProperties& EP, PhaseField& Phase,Orientations& OR, int tStep, string Filename)
{
    double maxRotation = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Phase.Fields,0, reduction(max:maxRotation))
    {
        double locangle = 0;
        dVector3 locaxis;
        Tools::getAxisAngle(EP.DeformationGradientsTotal(i,j,k), locaxis, locangle);
        for (auto it = Phase.Fields(i, j, k).cbegin();
        it != Phase.Fields(i, j, k).cend(); ++it)
        {
            locangle *= 180 / M_PI;
            maxRotation = std::max(maxRotation,locangle);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
#ifdef MPI_PARALLEL
double tempmaxRotation = maxRotation;
MPI_Allreduce(&tempmaxRotation,&maxRotation, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    WriteValue(maxRotation, Filename, tStep);
}

void TextOutput::LineConcentration(Composition& Cx, PhaseField& Phi,
                                   string filename, double timestep,
                                   string type, string axis,
                                   int x, int y, int z)
{
    /** This function will create a separate file each time it is called, with
    a different file name according to the current time step. In this file, the
    tabulated composition data will be written down along a single line along a
    given axis with type = "X", "Y" or "Z". The line can be positioned with a
    point on it given with "x", "y" and "z". The output can be weight fraction
    with type = "WF", weight percent "WP", mole fraction "MF" and mole percent
    with "MP".*/

    int mode = -1;
    int dir = -1;
    int tablength = 16;
    std::stringstream ss;

    if(type == "WF")
    {
        mode = 0;
    }
    else if(type == "WP")
    {
        mode = 1;
    }
    else if(type == "MF")
    {
        mode = 2;
    }
    else if(type == "MP")
    {
        mode = 3;
    }
    else
    {
        cout << "Undefined type " << type
             << " in TextOutput::LineConcentration()! Chose \"WF\", \"WP\", "
             << "\"MF\" or \"MP\". Now exiting!" << endl;
        exit(1);
    }

    if(axis == "X")
    {
        dir = 0;
    }
    else if(axis == "Y")
    {
        dir = 1;
    }
    else if(axis == "Z")
    {
        dir = 2;
    }
    else
    {
        cout << "Undefined axis direction " << axis
             << " in TextOutput::LineConcentration()! Chose \"X\", \"Y\" or"
             << "\"Z\". Now exiting!" << endl;
        exit(1);
    }

    ss << fixed << std::setw(9) << std::setfill('0') << int(timestep);
    std::string s = ss.str();
    string name = filename + s + ".opd";

    ofstream file(name, ios::out);                                              //Create file and write header (or overwrite)

    file << std::setw(5) << std::setfill(' ') << "# ";
    switch (dir)
    {
        case 0:
        {
            file << "x";
            break;
        }
        case 1:
        {
            file << "y";
            break;
        }
        case 2:
        {
            file << "z";
            break;
        }
    }

    for(size_t comp = 0; comp < Cx.Ncomp; comp++)
    {
        string temp;
        string name(Cx.ElementNames[comp]);

        if(mode == 0)
        {
            temp = "wf." + name;
        }
        else if(mode == 1)
        {
            temp = "wp." + name;
        }
        else if(mode == 2)
        {
            temp = "mf." + name;
        }
        else if(mode == 3)
        {
            temp = "mp." + name;
        }

        file << std::setw(tablength) << std::setfill(' ') << temp;
    }

    file << endl;

    switch(dir)                                                                 // Write data to the file
    {
        case 0:
        {
            for(int x = 0; x < Cx.Nx; x++)
            {
                file << fixed << std::setw(5) << std::setfill(' ') << x;
                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                {
                    double val = 0.0;
                    switch(mode)
                    {
                        case 0:
                        {//WF
                            val = Cx.WeightFractionsTotal(x,y,z)({comp});
                            break;
                        }
                        case 1:
                        {//WP
                            val = Cx.WeightFractionsTotal(x,y,z)({comp})*100.0;
                            break;
                        }
                        case 2:
                        {//MF
                            val = Cx.MoleFractionsTotal(x,y,z)({comp});
                            break;
                        }
                        case 3:
                        {//MP
                            val = Cx.MoleFractionsTotal(x,y,z)({comp})*100.0;
                            break;
                        }
                    }
                    file << scientific << std::setw(tablength) << std::setfill(' ')
                         << val;
                }
                file << endl;
            }
            break;
        }
        case 1:
        {
            for(int y = 0; y < Cx.Ny; y++)
            {
                file << fixed << std::setw(5) << std::setfill(' ') << y;
                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                {
                    double val = 0.0;
                    switch(mode)
                    {
                        case 0:
                        {//WF
                            val = Cx.WeightFractionsTotal(x,y,z)({comp});
                            break;
                        }
                        case 1:
                        {//WP
                            val = Cx.WeightFractionsTotal(x,y,z)({comp})*100.0;
                            break;
                        }
                        case 2:
                        {//MF
                            val = Cx.MoleFractionsTotal(x,y,z)({comp});
                            break;
                        }
                        case 3:
                        {//MP
                            val = Cx.MoleFractionsTotal(x,y,z)({comp})*100.0;
                            break;
                        }
                    }
                    file << scientific << std::setw(tablength) << std::setfill(' ')
                         << val;
                }
                file << endl;
            }
            break;
        }
        case 2:
        {
            for(int z = 0; z < Cx.Nz; z++)
            {
                file << fixed << std::setw(5) << std::setfill(' ') << z;
                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                {
                    double val = 0.0;
                    switch(mode)
                    {
                        case 0:
                        {//WF
                            val = Cx.WeightFractionsTotal(x,y,z)({comp});
                            break;
                        }
                        case 1:
                        {//WP
                            val = Cx.WeightFractionsTotal(x,y,z)({comp})*100.0;
                            break;
                        }
                        case 2:
                        {//MF
                            val = Cx.MoleFractionsTotal(x,y,z)({comp});
                            break;
                        }
                        case 3:
                        {//MP
                            val = Cx.MoleFractionsTotal(x,y,z)({comp})*100.0;
                            break;
                        }
                    }
                    file << scientific << std::setw(tablength) << std::setfill(' ')
                         << val;
                }
                file << endl;
            }
            break;
        }
    }
    file.close();                                                               //Close file
}

void TextOutput::LocalPhaseComposition(Composition& Cx,
                                       PhaseField& Phi, string filename,
                                       double time, int x, int y, int z)
{
    /** This function will create tabulated data on the phase composition of
    each thermodynamic phase defined in the chemical property input. Each time
    this function is called, a new row will be written in the specified file name
    for the current time step. If the file is not present it will be created,
    if the phase already exists, the new data will be appended!*/

    if (!FileExists(filename))
    {
        ofstream file(filename, ios::out);                                      //Create file and write header (or overwrite)
        file << "# time ";
        for(size_t alpha = 0; alpha < Cx.Nphases; ++alpha)
        for(size_t comp = 0; comp < Cx.Ncomp; comp++)
        {
            file << Cx.PhaseNames[alpha] << "_"
                 << Cx.ElementNames[comp] << "    ";
        }
        file << endl;
        file.close();
    }

    ofstream file(filename, ios::app);                                          //Write content of the file
    file << scientific << time << "   ";

    for(size_t alpha = 0; alpha < Cx.Nphases; ++alpha)
    {
        if(Phi.Fractions(x,y,z)({alpha}) > 0.0)
        {
            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
            {
                file << Cx.MoleFractions(x,y,z)({alpha,comp}) << "   ";
            }
        }
        else
        {
            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
            file << 0.0 << "   ";
        }
    }

    file << endl;
    file.close();
}

}// namespace openphase
