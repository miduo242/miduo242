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
 *   Main contributors :   Muhammad Adil Ali; Oleg Shchyglo
 *
 */

#ifndef DVECTORN_H
#define DVECTORN_H

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace openphase
{

class dVectorN
{
 public:
    dVectorN()
    {
    };
    dVectorN(const size_t size): storage(size, 0.0)
    {
    };
    dVectorN(const size_t size, const double value): storage(size, value)
    {
    };
    dVectorN(const dVectorN& other): storage(other.storage)
    {
    };
    dVectorN(std::initializer_list<double> vecinit)
    {
        storage.resize(vecinit.size());
        int ii = 0;
        for (auto it = vecinit.begin(); it != vecinit.end(); it++)
        {
            storage[ii] = *it;
            ii += 1;
        }
    }
    void set_to_zero()
    {
        std::fill(storage.begin(), storage.end(), 0.0);
    };
    void set_to_value(const double value)
    {
        std::fill(storage.begin(), storage.end(), value);
    };
    void Allocate(const size_t size)
    {
        if(storage.size() != size)
        {
            storage.resize(size);
            set_to_zero();
        }
    };
    double norm()
    {
    double dnorm = 0.;
    for(size_t i = 0; i < storage.size(); i++)
    {
           dnorm += storage[i]*storage[i];
        }
    return sqrt(dnorm);
    }
    void normalize()
    {
    double dnorm = norm();
    if (dnorm != 0)
    for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] /= dnorm;
        }
    }
    double operator*(const dVectorN& rhs) const
    {
    double dsum = 0.;
    for(size_t i = 0; i < storage.size(); i++)
        {
            dsum += storage[i]*rhs[i];
        }
         return dsum;
    }
    double& operator[](const size_t i)
    {
        assert(i < storage.size() && "Access beyond storage range");
        return storage[i];
    };
    double const& operator[](const size_t i) const
    {
        assert(i < storage.size() && "Access beyond storage range");
        return storage[i];
    };
    dVectorN& operator=(const dVectorN& rhs)
    {
        assert((rhs.size() == storage.size() or storage.size() == 0) && "Sizes of the vectors are not equal");
        storage = rhs.storage;
        return *this;
    };
    dVectorN operator*(const double m) const
    {
        dVectorN tmp(storage.size(),0.0);
        for(size_t i = 0; i < storage.size(); i++)
        {
            tmp[i] = storage[i]*m;
        }
        return tmp;
    };
    dVectorN operator/(const double m) const
    {
        if(m == 0.0)
        {
            std::stringstream message;
            message << "Error in dVectorN::operator/ :\n"
                    << "Division by zero!"
                    << "\nTerminating!!!\n";
            std::cerr << message.str();
            std::terminate();
        }
        dVectorN tmp(storage.size(),0.0);
        for(size_t i = 0; i < storage.size(); i++)
        {
            tmp[i] = storage[i]/m;
        }
        return tmp;
    };
    dVectorN& operator*=(const double m)
    {
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] *= m;
        }
        return *this;
    };
    dVectorN& operator/=(const double m)
    {
        if(m == 0.0)
        {
            std::stringstream message;
            message << "Error in dVectorN::operator/= :\n"
                    << "Division by zero!"
                    << "\nTerminating!!!\n";
            std::cerr << message.str();
            std::terminate();
        }
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] /= m;
        }
        return *this;
    };
    dVectorN operator+(dVectorN& rhs) const
    {
        assert(rhs.size() == storage.size() && "Sizes of the vectors are not equal");
        dVectorN tmp(storage.size(),0.0);
        for(size_t i = 0; i < storage.size(); i++)
        {
            tmp[i] = storage[i] + rhs[i];
        }
        return tmp;
    };
    dVectorN& operator+=(dVectorN& rhs)
    {
        assert(rhs.size() == storage.size() && "Sizes of the vectors are not equal");
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] = storage[i] + rhs[i];
        }
        return *this;
    };
    dVectorN operator-(dVectorN& rhs) const
    {
        assert(rhs.size() == storage.size() && "Sizes of the vectors are not equal");
        dVectorN tmp(storage.size(),0.0);
        for(size_t i = 0; i < storage.size(); i++)
        {
            tmp[i] = storage[i] - rhs[i];
        }
        return tmp;
    };
    dVectorN& operator-=(dVectorN& rhs)
    {
        assert(rhs.size() == storage.size() && "Sizes of the vectors are not equal");
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] = storage[i] - rhs[i];
        }
        return *this;
    };
    std::string print(void) const
    {
        std::stringstream out;
        out << "< | ";
        for(size_t i = 0; i < storage.size(); i++)
        {
            out << storage[i]<< " " << " | ";
        }
        out << " >";
        return out.str();
    };
    size_t size(void) const
    {
        return storage.size();
    }
    double* data(void)
    {
        return storage.data();
    }
    const double* data(void) const
    {
        return storage.data();
    }
    void Read(std::istream& inp)
    {
        size_t size = 0;
        inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));
        storage.resize(size);
        for(size_t n = 0; n < size; n++)
        {
            inp.read(reinterpret_cast<char*>(&storage[n]), sizeof(double));
        }
    }
    void Write(std::ostream& outp)
    {
        size_t size = storage.size();
        outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));
        for(size_t n = 0; n < size; n++)
        {
            outp.write(reinterpret_cast<const char*>(&storage[n]), sizeof(double));
        }
    }
 protected:
 private:
    std::vector<double> storage;
};

} // namespace openphase
#endif
