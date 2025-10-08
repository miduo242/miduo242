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

#ifndef DMATRIX3X3_H
#define DMATRIX3X3_H

#include <cassert>
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

namespace openphase
{

class dMatrix3x3
{
 public:

    dMatrix3x3()
    {
        memset(storage, 0, 9*sizeof(double));
    };

    dMatrix3x3(const dMatrix3x3& rhs)
    {
        memmove(storage, rhs.const_data(), 9*sizeof(double));
    };

    dMatrix3x3(std::initializer_list<double> matinit)
    {
        assert(matinit.size() == 9 && "Initialization list size is not equal to storage range.");

        int ii = 0;
        int jj = 0;
        for (auto it = matinit.begin(); it != matinit.end(); it++)
        {
            if (ii == 3)
            {
                ii =  0;
                jj += 1;
            }
            storage[jj][ii] = *it;
            ii += 1;
        }
    }

    double& operator()(const size_t i, const size_t j)
    {
        assert(i < 3 && "Access beyond storage range");
        assert(j < 3 && "Access beyond storage range");

        return storage[i][j];
    };
    double const& operator()(const size_t i, const size_t j) const
    {
        assert(i < 3 && "Access beyond storage range");
        assert(j < 3 && "Access beyond storage range");

        return storage[i][j];
    };
    dMatrix3x3& set_to_zero(void)
    {
        memset(storage, 0, 9*sizeof(double));
        return *this;
    };
    dMatrix3x3& set_to_unity(void)
    {
        memset(storage, 0, 9*sizeof(double));
        storage[0][0] = 1.0;
        storage[1][1] = 1.0;
        storage[2][2] = 1.0;
        return *this;
    };
    void pack(std::vector<double>& buffer)
    {
        for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            buffer.push_back(storage[i][j]);
        }
    }
    void unpack(std::vector<double>& buffer, size_t& it)
    {
        for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            storage[i][j] = buffer[it]; ++it;
        }
    }

    void set(double r00, double r01, double r02,
             double r10, double r11, double r12,
             double r20, double r21, double r22)
    {
        storage[0][0] = r00;
        storage[0][1] = r01;
        storage[0][2] = r02;
        storage[1][0] = r10;
        storage[1][1] = r11;
        storage[1][2] = r12;
        storage[2][0] = r20;
        storage[2][1] = r21;
        storage[2][2] = r22;
    }

    dMatrix3x3& operator=(const dMatrix3x3& rhs)
    {
        memmove(storage, rhs.const_data(), 9*sizeof(double));
        return *this;
    };

    bool operator==(const dMatrix3x3& rhs)
    {
        for(int i=0; i < 3; i++)
        for(int j=0; j < 3; j++)
        {
            if(fabs(storage[i][j]-rhs(i,j)) > DBL_EPSILON)
            {
                return false;
            }
        }
        return true;
    };

    bool operator!=(const dMatrix3x3& rhs)
    {
        for(int i=0; i < 3; i++)
        for(int j=0; j < 3; j++)
        {
            if(fabs(storage[i][j]-rhs(i,j)) > DBL_EPSILON)
            {
                return true;
            }
        }
        return false;
    };

    dMatrix3x3 operator*(const double m) const
    {
        dMatrix3x3 tmp;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp(i,j) = storage[i][j]*m;
        }
        return tmp;
    };
    dMatrix3x3 operator/(const double m) const
    {
        dMatrix3x3 tmp;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp(i,j) = storage[i][j]/m;
        }
        return tmp;
    };

    dMatrix3x3 operator*(const dMatrix3x3& rhs) const
    {
        dMatrix3x3 tmp;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        for(int k = 0; k < 3; k++)
        {
            tmp(i,j) += storage[i][k]*rhs(k,j);
        }
        return tmp;
    };
    dMatrix3x3 operator+(const dMatrix3x3& rhs) const
    {
        dMatrix3x3 tmp;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp(i,j) = storage[i][j] + rhs(i,j);
        }
        return tmp;
    };
    dMatrix3x3 operator-(const dMatrix3x3& rhs) const
    {
        dMatrix3x3 tmp;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp(i,j) = storage[i][j] - rhs(i,j);
        }
        return tmp;
    };
    dMatrix3x3& operator+=(const dMatrix3x3& rhs)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            storage[i][j] += rhs(i,j);
        }
        return *this;
    };
    dMatrix3x3& operator-=(const dMatrix3x3& rhs)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            storage[i][j] -= rhs(i,j);
        }
        return *this;
    };
    dMatrix3x3& operator*=(const double m)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            storage[i][j] *= m;
        }
        return *this;
    };
    dMatrix3x3& operator/=(const double m)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            storage[i][j] /= m;
        }
        return *this;
    };
    double& operator[](const size_t i)
    {
        assert(i < 9 && "Access beyond storage range");

        int column = i % 3;
        int row = int((i - column)/3);
        return storage[row][column];
    };
    dMatrix3x3 H_product(const dMatrix3x3& rhs) const
    {
        dMatrix3x3 tmp;

        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp(i,j) = storage[i][j]*rhs(i,j);
        }
        return tmp;
    };
    inline double determinant(void) const
    {
        return (storage[0][0]*storage[1][1]*storage[2][2] +
                storage[0][1]*storage[1][2]*storage[2][0] +
                storage[0][2]*storage[1][0]*storage[2][1] -
                storage[0][2]*storage[1][1]*storage[2][0] -
                storage[0][1]*storage[1][0]*storage[2][2] -
                storage[0][0]*storage[1][2]*storage[2][1]);
    };
    dMatrix3x3& invert(void)
    {
        dMatrix3x3 tmp;

        double detInv = determinant();

        if(detInv != 0.0)
        {
            detInv = 1.0/detInv;
        }
        else
        {
            std::cerr << "dMatrix3x3: Can Not Compute Inverse Matrix.\n"
                      << this->print() << "Matrix is Singular !!!\n";
            std::terminate();
        }

        tmp(0,0) = (storage[1][1]*storage[2][2] -
                    storage[1][2]*storage[2][1])*detInv;
        tmp(1,0) =-(storage[1][0]*storage[2][2] -
                    storage[1][2]*storage[2][0])*detInv;
        tmp(2,0) = (storage[1][0]*storage[2][1] -
                    storage[1][1]*storage[2][0])*detInv;
        tmp(0,1) =-(storage[0][1]*storage[2][2] -
                    storage[0][2]*storage[2][1])*detInv;
        tmp(1,1) = (storage[0][0]*storage[2][2] -
                    storage[0][2]*storage[2][0])*detInv;
        tmp(2,1) =-(storage[0][0]*storage[2][1] -
                    storage[0][1]*storage[2][0])*detInv;
        tmp(0,2) = (storage[0][1]*storage[1][2] -
                    storage[1][1]*storage[0][2])*detInv;
        tmp(1,2) =-(storage[0][0]*storage[1][2] -
                    storage[0][2]*storage[1][0])*detInv;
        tmp(2,2) = (storage[0][0]*storage[1][1] -
                    storage[0][1]*storage[1][0])*detInv;

        memmove(storage, tmp.data(), 9*sizeof(double));
        return *this;
    };
    dMatrix3x3 inverted(void) const
    {
        dMatrix3x3 tmp;

        double detInv = determinant();

        if(detInv != 0.0)
        {
            detInv = 1.0/detInv;
        }
        else
        {
            std::cerr << "dMatrix3x3: Can Not Compute Inverse Matrix.\n"
                      << this->print() << "Matrix is Singular !!!\n";
            std::terminate();
        }

        tmp(0,0) = (storage[1][1]*storage[2][2] -
                    storage[1][2]*storage[2][1])*detInv;
        tmp(1,0) =-(storage[1][0]*storage[2][2] -
                    storage[1][2]*storage[2][0])*detInv;
        tmp(2,0) = (storage[1][0]*storage[2][1] -
                    storage[1][1]*storage[2][0])*detInv;
        tmp(0,1) =-(storage[0][1]*storage[2][2] -
                    storage[0][2]*storage[2][1])*detInv;
        tmp(1,1) = (storage[0][0]*storage[2][2] -
                    storage[0][2]*storage[2][0])*detInv;
        tmp(2,1) =-(storage[0][0]*storage[2][1] -
                    storage[0][1]*storage[2][0])*detInv;
        tmp(0,2) = (storage[0][1]*storage[1][2] -
                    storage[1][1]*storage[0][2])*detInv;
        tmp(1,2) =-(storage[0][0]*storage[1][2] -
                    storage[0][2]*storage[1][0])*detInv;
        tmp(2,2) = (storage[0][0]*storage[1][1] -
                    storage[0][1]*storage[1][0])*detInv;

        return tmp;
    };
    dMatrix3x3& transpose(void)
    {
        for(int i = 0; i < 2; i++)
        for(int j = i+1; j < 3; j++)
        {
              double tmp = storage[i][j];
              storage[i][j] = storage[j][i];
              storage[j][i] = tmp;
        }
        return *this;
    };
    dMatrix3x3 transposed(void) const
    {
        dMatrix3x3 TempMat;

        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            TempMat(i,j) = storage[j][i];
        }
        return TempMat;
    };
    dMatrix3x3& rotate(const dMatrix3x3& RotationMatrix)
    {
        dMatrix3x3 TempMat;

        memmove(TempMat.data(), storage, 9*sizeof(double));

        TempMat = RotationMatrix * (TempMat * RotationMatrix.transposed());

        memmove(storage, TempMat.data(), 9*sizeof(double));
        return *this;
    };
    dMatrix3x3 rotated(const dMatrix3x3& RotationMatrix) const
    {
        dMatrix3x3 Out;
        memmove(Out.data(), storage, 9*sizeof(double));

        Out = RotationMatrix * (Out * RotationMatrix.transposed());

        return Out;
    };

    dMatrix3x3& rotateU(const dMatrix3x3& RotationMatrix)
    {
        dMatrix3x3 TempMat;

        memmove(TempMat.data(), storage, 9*sizeof(double));

        TempMat = RotationMatrix * TempMat;

        memmove(storage, TempMat.data(), 9*sizeof(double));
        return *this;
    };
    dMatrix3x3 rotatedU(const dMatrix3x3& RotationMatrix) const
    {
        dMatrix3x3 Out;
        memmove(Out.data(), storage, 9*sizeof(double));

        Out = RotationMatrix * Out;

        return Out;
    };

    double double_contract(const dMatrix3x3& rHS) const
    {
        double tmp = 0.0;
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            tmp += storage[i][j]*rHS(i,j);
        }
        return tmp;
    };
    double trace(void) const
    {
        return storage[0][0] + storage[1][1] + storage[2][2];
    };
    double max(void) const
    {
        double value = -DBL_MAX;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            value = std::max(value,storage[i][j]);
        }
        return value;
    }
    double min(void) const
    {
        double value = DBL_MAX;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            value = std::min(value,storage[i][j]);
        }
        return value;
    }
    dMatrix3x3 get_sym(void) const
    {
        dMatrix3x3 TempMat;

        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            TempMat(i,j) = (storage[i][j] + storage[j][i]) * 0.5;
        }

        return TempMat;
    };
    dMatrix3x3 get_skew(void) const
    {
        dMatrix3x3 TempMat;

        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            TempMat(i,j) = (storage[i][j] - storage[j][i]) * 0.5;
        }

        return TempMat;
    };
    double norm(void) const                                                     ///< Frobenius norm
    {
        double tmp = 0.0;
        for(int i = 0; i < 3; i++)
        {
            tmp += storage[i][0] * storage[i][0]
                 + storage[i][1] * storage[i][1]
                 + storage[i][2] * storage[i][2];
        }

        return sqrt(tmp);
    };
    std::string print(void) const
    {
        std::stringstream out;
        for(int i = 0; i < 3; i++)
        {
            out << "||" << std::setprecision(6) << std::right
                        << std::setw(8) << storage[i][0] << " "
                        << std::setw(8) << storage[i][1] << " "
                        << std::setw(8) << storage[i][2] << "||\n";
        }
        return out.str();
    };
    std::string write(const int precision = 16, const char sep = ' ') const
    {
        std::stringstream out;
        out << std::setprecision(precision) << std::defaultfloat;
        for(int i = 0; i < 3; i++)
        {
           for(int j = 0; j < 3; j++)
           {
               out << storage[i][j] << sep;
           }
           out << "\n";
        }
        return out.str();
    };
    void read(std::fstream& inp)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            inp >> storage[i][j];
        }
    };
    std::vector<double> writeBinary() const
    {
        std::vector<double> out;
        for(int i = 0; i < 3; i++)
        {
           for(int j = 0; j < 3; j++)
           {
               out.push_back((double)storage[i][j]);
           }
        }
        return out;
    };
    std::vector<float> writeCompressed() const
    {
        std::vector<float> out;
        for(int i = 0; i < 3; i++)
        {
           for(int j = 0; j < 3; j++)
           {
               out.push_back((float)storage[i][j]);
           }
        }
        return out;
    };
    void read(std::stringstream& inp)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            inp >> storage[i][j];
        }
    };
    double * data(void)
    {
        return &storage[0][0];
    };
    double const * const_data(void) const
    {
        return &storage[0][0];
    };

    static dMatrix3x3 UnitTensor(void)
    {
        dMatrix3x3 unity;
        return unity.set_to_unity();
    };
    static dMatrix3x3 ZeroTensor()
    {
        dMatrix3x3 myZeroTensor;
        return myZeroTensor.set_to_zero();
    };

 protected:
 private:
    double storage[3][3];
};

}// namespace openphase
#endif //DMATRIX3X3_H
