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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev; Philipp Engels
 *
 */

#ifndef IVECTOR3_H
#define IVECTOR3_H

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

class iVector3
{
 public:

    iVector3()
    {
        memset(storage.data(), 0, 3*sizeof(int));
    };
    iVector3(const iVector3& vecinit)
    {
        storage = vecinit.storage;
    };
    iVector3(std::initializer_list<int> vecinit)
    {
        assert(vecinit.size() == 3 && "Initialization list size is not equal to storage range");
        int ii = 0;
        for (auto it = vecinit.begin(); it != vecinit.end(); it++)
        {
            storage[ii] = *it;
            ii += 1;
        }
    }
    int& operator[](const size_t i)
    {
        assert(i < 3 && "Access beyond storage range");
        return storage[i];
    };
    int const& operator[](const size_t i) const
    {
        assert(i < 3 && "Access beyond storage range");
        return storage[i];
    };
    bool operator<(iVector3 rhs) const
    {
        return length() < rhs.length();
    }
    bool operator>(iVector3 rhs) const
    {
        return length() > rhs.length();
    }
    bool operator==(iVector3 rhs) const
    {
        return length_sqr() == rhs.length_sqr();
    }
    int getX(void) const
    {
        return storage[0];
    };

    void setX(const int newX)
    {
        storage[0] = newX;
    };

    int getY(void) const
    {
        return storage[1];
    };

    void setY(const int newY)
    {
        storage[1] = newY;
    };

    int getZ(void) const
    {
        return storage[2];
    };

    void setZ(const int newX)
    {
        storage[2] = newX;
    };
    void pack(std::vector<int>& buffer)
    {
        for (int i = 0; i < 3; ++i)
        {
            buffer.push_back(storage[i]);
        }
    }
    void unpack(std::vector<int>& buffer, size_t& it)
    {
        for (int i = 0; i < 3; ++i)
        {
            storage[i] = buffer[it]; ++it;
        }
    }
    void set_to_zero(void)
    {
        memset(storage.data(), 0, 3*sizeof(int));
    };
    void set_to_unitX(void)
    {
        storage[0] = 1;
        storage[1] = 0;
        storage[2] = 0;
    };
    void set_to_unitY(void)
    {
        storage[0] = 0;
        storage[1] = 1;
        storage[2] = 0;
    };
    void set_to_unitZ(void)
    {
        storage[0] = 0;
        storage[1] = 0;
        storage[2] = 1;
    };
    iVector3 operator*(const int m) const
    {
        iVector3 tmp;
        tmp[0] = storage[0]*m;
        tmp[1] = storage[1]*m;
        tmp[2] = storage[2]*m;
        return tmp;
    };
    iVector3 operator/(const int m) const
    {
        iVector3 tmp;
        tmp[0] = storage[0]/m;
        tmp[1] = storage[1]/m;
        tmp[2] = storage[2]/m;
        return tmp;
    };
    int operator*(const iVector3& rhs) const
    {
        return storage[0]*rhs[0] + storage[1]*rhs[1] + storage[2]*rhs[2];
    };
    double abs() const
    {
        return sqrt(storage[0]*storage[0] +
                    storage[1]*storage[1] +
                    storage[2]*storage[2]);
    };
    double length(void) const
    {
        return sqrt(storage[0]*storage[0] +
                    storage[1]*storage[1] +
                    storage[2]*storage[2]);
    };
    int length_sqr(void) const
    {
        return (storage[0]*storage[0] +
                storage[1]*storage[1] +
                storage[2]*storage[2]);
    };
    iVector3 cross(const iVector3& rhs) const
    {
        iVector3 tmp;
        tmp[0] = storage[1] * rhs[2] - storage[2] * rhs[1];
        tmp[1] = storage[2] * rhs[0] - storage[0] * rhs[2];
        tmp[2] = storage[0] * rhs[1] - storage[1] * rhs[0];
        return tmp;
    };
    iVector3 operator+(const iVector3& rhs) const
    {
        iVector3 tmp;
        tmp[0] = storage[0] + rhs[0];
        tmp[1] = storage[1] + rhs[1];
        tmp[2] = storage[2] + rhs[2];
        return tmp;
    };
    iVector3 operator-(const iVector3& rhs) const
    {
        iVector3 tmp;
        tmp[0] = storage[0] - rhs[0];
        tmp[1] = storage[1] - rhs[1];
        tmp[2] = storage[2] - rhs[2];
        return tmp;
    };
    iVector3& operator*=(const int m)
    {
        storage[0] *= m;
        storage[1] *= m;
        storage[2] *= m;
        return *this;
    };
    iVector3& operator/=(const int m)
    {
        storage[0] /= m;
        storage[1] /= m;
        storage[2] /= m;
        return *this;
    };
    iVector3& operator-=(const iVector3& rhs)
    {
        storage[0] -= rhs[0];
        storage[1] -= rhs[1];
        storage[2] -= rhs[2];
        return *this;
    };
    iVector3& operator+=(const iVector3& rhs)
    {
        storage[0] += rhs[0];
        storage[1] += rhs[1];
        storage[2] += rhs[2];
        return *this;
    };
    iVector3& operator=(const iVector3& rhs)
    {
        storage = rhs.storage;
        return *this;
    };
    iVector3& operator=(const int rhs[3])
    {
        memmove(storage.data(), rhs, 3*sizeof(int));
        return *this;
    };
    dVector3 normalized(void) const
    {
        double norm = length();
        dVector3 tmp;
        if(norm != 0.0)
        {
            double norm_inv = 1.0/norm;
            tmp[0] = storage[0]*norm_inv;
            tmp[1] = storage[1]*norm_inv;
            tmp[2] = storage[2]*norm_inv;
        }
        else
        {
            tmp.set_to_zero();
        }
        return tmp;
    };
    iVector3 Xreflected(void) const
    {
        iVector3 Out(*this);
        Out[0] *= -1.0;
        return Out;
    };
    iVector3 Yreflected(void) const
    {
        iVector3 Out(*this);
        Out[1] *= -1.0;
        return Out;
    };
    iVector3 Zreflected(void) const
    {
        iVector3 Out(*this);
        Out[2] *= -1.0;
        return Out;
    };
    std::string print(void) const
    {
        std::stringstream out;

        out << "(" << storage[0] << ", "
                   << storage[1] << ", "
                   << storage[2] << ")";
        return out.str();
    };

    std::string write(const int precision = 16, const char sep = ' ') const
    {
        std::stringstream out;
        out << std::setprecision(precision) << std::scientific;
        for(int i = 0; i < 3; i++)
        {
            out << storage[i] << sep;
        }
        return out.str();
    };
/*
    int* data(void)
    {
        return storage.data();
    };
    const int* const_data(void) const
    {
        return storage.data();
    };*/
 protected:
 private:
    std::array<int, 3> storage;
};

}// namespace openphase
#endif
