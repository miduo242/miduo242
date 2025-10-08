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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev;
 *                         Philipp Engels
 *
 */

#ifndef VSTRAIN_H
#define VSTRAIN_H

#include <cfloat>
#include <cmath>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Base/dMatrix3x3.h"

namespace openphase
{

class vStrain
{
    /*
     *  This is a special version of the six component Voigt vector:
     *  the off diagonal elements of the strain matrix are multiplied by 2(!)
     *  while transforming the 3x3 strain matrix to a 6-component vector,
     *  and have to be divided by two if they are transformed to a 3x3 matrix
     *  again. Therefore the transformation functions and some other functions
     *  have to be specified in different way then for the stress-type
     *  Voigt vector. The differences are specified in the functions.
     */
 public:
    vStrain()
    {
        memset(data(), 0, 6*sizeof(double));
    };
    vStrain(const vStrain& rhs)
    {
        memmove(data(), rhs.const_data(), 6*sizeof(double));
    };
    vStrain& operator=(const vStrain& rhs)
    {
        memmove(data(), rhs.const_data(), 6*sizeof(double));
        return *this;
    };

    double& operator[](const int i)
    {
#ifdef DEBUG
        if(i > 5)
        {
            std::stringstream message;
            message << "Error in vStrain::operator[]\n"
                    << "Access beyond storage range. i = "
                    << i << " > 5"
                    << "\nTerminating!!!\n";
            throw std::logic_error(message.str());
        }
#endif
        return storage[i];
    };
    double const& operator[](const int i) const
    {
#ifdef DEBUG
        if(i > 5)
        {
            std::stringstream message;
            message << "Error in vStrain::operator[]\n"
                    << "Access beyond storage range. i = "
                    << i << " > 5"
                    << "\nTerminating!!!\n";
            throw std::logic_error(message.str());
        }
#endif
        return storage[i];
    };
    vStrain& set_to_zero(void)
    {
        memset(storage, 0, 6*sizeof(double));
        return *this;
    };
    void pack(std::vector<double>& buffer)
    {
        for (int i = 0; i < 6; ++i)
        {
            buffer.push_back(storage[i]);
        }
    };
    void unpack(std::vector<double>& buffer, size_t& it)
    {
        for (int i = 0; i < 6; ++i)
        {
            storage[i] = buffer[it]; ++it;
        }
    };
    vStrain& set_to_unity(void)
    {
        storage[0] = 1.0;
        storage[1] = 1.0;
        storage[2] = 1.0;
        storage[3] = 0.0;
        storage[4] = 0.0;
        storage[5] = 0.0;
        return *this;
    };
    bool operator==(const vStrain& rhs)
    {
        for (int i = 0; i < 6; i++)
        {
            if (fabs(storage[i] - rhs[i]) > FLT_EPSILON)
            {
                return false;
            };
        }
        return true;
    };
    vStrain operator*(const double m) const
    {
        vStrain tmp;
        tmp[0] = storage[0]*m;
        tmp[1] = storage[1]*m;
        tmp[2] = storage[2]*m;
        tmp[3] = storage[3]*m;
        tmp[4] = storage[4]*m;
        tmp[5] = storage[5]*m;
        return tmp;
    };
    vStrain& operator*=(const double m)
    {
        storage[0] *= m;
        storage[1] *= m;
        storage[2] *= m;
        storage[3] *= m;
        storage[4] *= m;
        storage[5] *= m;
        return *this;
    };
    vStrain operator/(const double m) const
    {
        vStrain tmp;
        double num = 1.0/m;
        tmp[0] = storage[0]*num;
        tmp[1] = storage[1]*num;
        tmp[2] = storage[2]*num;
        tmp[3] = storage[3]*num;
        tmp[4] = storage[4]*num;
        tmp[5] = storage[5]*num;
        return tmp;
    };
    vStrain& operator/=(const double m)
    {
        double num = 1.0/m;
        storage[0] *= num;
        storage[1] *= num;
        storage[2] *= num;
        storage[3] *= num;
        storage[4] *= num;
        storage[5] *= num;
        return *this;
    };
    vStrain operator+(const double m) const
    {
        vStrain tmp;
        tmp[0] = storage[0]+m;
        tmp[1] = storage[1]+m;
        tmp[2] = storage[2]+m;
        tmp[3] = storage[3]+m;
        tmp[4] = storage[4]+m;
        tmp[5] = storage[5]+m;
        return tmp;
    };
    vStrain operator-(const double m) const
    {
        vStrain tmp;
        tmp[0] = storage[0]-m;
        tmp[1] = storage[1]-m;
        tmp[2] = storage[2]-m;
        tmp[3] = storage[3]-m;
        tmp[4] = storage[4]-m;
        tmp[5] = storage[5]-m;
        return tmp;
    };
    vStrain operator-(const vStrain& rhs) const
    {
        vStrain tmp;
        tmp[0] = storage[0] - rhs[0];
        tmp[1] = storage[1] - rhs[1];
        tmp[2] = storage[2] - rhs[2];
        tmp[3] = storage[3] - rhs[3];
        tmp[4] = storage[4] - rhs[4];
        tmp[5] = storage[5] - rhs[5];
        return tmp;
    };
    vStrain& operator-=(const vStrain& rhs)
    {
        storage[0] -= rhs[0];
        storage[1] -= rhs[1];
        storage[2] -= rhs[2];
        storage[3] -= rhs[3];
        storage[4] -= rhs[4];
        storage[5] -= rhs[5];
        return *this;
    };
    vStrain operator+(const vStrain& rhs) const
    {
        vStrain tmp;
        tmp[0] = storage[0] + rhs[0];
        tmp[1] = storage[1] + rhs[1];
        tmp[2] = storage[2] + rhs[2];
        tmp[3] = storage[3] + rhs[3];
        tmp[4] = storage[4] + rhs[4];
        tmp[5] = storage[5] + rhs[5];
        return tmp;
    };
    vStrain& operator+=(const vStrain& rhs)
    {
        storage[0] += rhs[0];
        storage[1] += rhs[1];
        storage[2] += rhs[2];
        storage[3] += rhs[3];
        storage[4] += rhs[4];
        storage[5] += rhs[5];
        return *this;
    };
    double norm(void) const  /// Frobenius norm
    {
        /*
         * multiplication by 0.5 of the off diagonal elements before taking
         * square and doubling the result due to double appearance of the off
         * diagonal elements => 0.5^2 * 2 = 0.5
         */
        double tmp = storage[0] * storage[0]
                   + storage[1] * storage[1]
                   + storage[2] * storage[2]
                   + storage[3] * storage[3] * 0.5
                   + storage[4] * storage[4] * 0.5
                   + storage[5] * storage[5] * 0.5;
        return sqrt(tmp);
    };
    double double_contract(const vStrain& Bstrain) const // "double-dot product"
    {
        double tmp =
        storage[0] * Bstrain[0] +
        storage[1] * Bstrain[1] +
        storage[2] * Bstrain[2] +
        2.0*(0.5*storage[3]) * (0.5*Bstrain[3]) +
        2.0*(0.5*storage[4]) * (0.5*Bstrain[4]) +
        2.0*(0.5*storage[5]) * (0.5*Bstrain[5]);
        return tmp;
    };
    double max_abs(void) const
    {
        // Returns maximum absolute value
        double tempmax = 0.0;
        for (int i = 0; i < 6; i++)
        {
            if (std::abs(storage[i]) > tempmax)
            {
                tempmax = std::abs(storage[i]);
            }
        }
        return tempmax;
    };
    vStrain rotated(const dMatrix3x3& RotationMatrix)
    {
        dMatrix3x3 locStrainTensor = tensor().rotated(RotationMatrix);
        vStrain locStrain;
        locStrain[0] = locStrainTensor(0,0);
        locStrain[1] = locStrainTensor(1,1);
        locStrain[2] = locStrainTensor(2,2);
        locStrain[3] = locStrainTensor(1,2)*2.0;
        locStrain[4] = locStrainTensor(0,2)*2.0;
        locStrain[5] = locStrainTensor(0,1)*2.0;
        return locStrain;
    };
    vStrain& rotate(const dMatrix3x3& RotationMatrix)
    {
        dMatrix3x3 locStrainTensor = tensor().rotated(RotationMatrix);

        storage[0] = locStrainTensor(0,0);
        storage[1] = locStrainTensor(1,1);
        storage[2] = locStrainTensor(2,2);
        storage[3] = locStrainTensor(1,2)*2.0;
        storage[4] = locStrainTensor(0,2)*2.0;
        storage[5] = locStrainTensor(0,1)*2.0;

        return *this;
    };
    double trace() const
    {
        return storage[0] + storage[1] + storage[2];
    };
    double get_tensor(const int i, const int j) const
    {
        /*
         * multiplication by 0.5 of the off diagonal elements during translating
         * to a 3x3 matrix
         */
        return (storage[(i==j)?(i):(6-(i+j))] * ((i==j)?(1.0):(0.5)));
    };
    dMatrix3x3 tensor(void) const
    {
        /*
         * multiplication by 0.5 of the off diagonal elements during translating
         * to a 3x3 matrix
         */
        dMatrix3x3 tmp;
        tmp(0,0) = storage[0];
        tmp(0,1) = storage[5]*0.5;
        tmp(0,2) = storage[4]*0.5;
        tmp(1,0) = storage[5]*0.5;
        tmp(1,1) = storage[1];
        tmp(1,2) = storage[3]*0.5;
        tmp(2,0) = storage[4]*0.5;
        tmp(2,1) = storage[3]*0.5;
        tmp(2,2) = storage[2];
        return tmp;
    };
    std::string print(void) const
    {
        std::stringstream out;
        out << "< | ";
        for(int i = 0; i < 6; i++)
        {
            out << storage[i]<< " " << " | ";
        }
        out << " >";
        return out.str();
    };
    std::string write(const int precision = 16, const char sep = ' ') const
    {
        std::stringstream out;
        out << std::setprecision(precision) << std::defaultfloat;
        out << storage[0] << sep;
        out << storage[1] << sep;
        out << storage[2] << sep;
        out << storage[5] << sep;
        out << storage[3] << sep;
        out << storage[4] << sep;
        return out.str();
    };
    std::vector<float> writeCompressed() const
    {
        std::vector<float> out;
        out.push_back((float)storage[0]);
        out.push_back((float)storage[1]);
        out.push_back((float)storage[2]);
        out.push_back((float)storage[5]);
        out.push_back((float)storage[3]);
        out.push_back((float)storage[4]);
        return out;
    };
    std::vector<double> writeBinary() const
    {
        std::vector<double> out;
        out.push_back((double)storage[0]);
        out.push_back((double)storage[1]);
        out.push_back((double)storage[2]);
        out.push_back((double)storage[5]);
        out.push_back((double)storage[3]);
        out.push_back((double)storage[4]);
        return out;
    };
    double* data(void)
    {
        return &storage[0];
    };
    const double* const_data(void) const
    {
        return &storage[0];
    };
 protected:
 private:
    double storage[6];
};

}// namespace openphase
#endif
