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
 *                         Raphael Schiedung
 *
 */

#ifndef TENSOR_H
#define TENSOR_H

#include <cfloat>
#include <cmath>
#include <cstring>
#include <fstream>
#include <array>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace openphase
{
template <class T, int Rank>
class Tensor                                                                    /// Tensor template class. It can handle any type of values but POD is preferred
{
 public:
    ~Tensor()
    {
        if(allocated)
        {
            delete[] locData;
        }
    }
    Tensor()
    {
        locData = nullptr;
        Dimensions.fill(0);
        allocated = false;
        totSize = 0;
    }
    Tensor(std::array<size_t, Rank> locDimensions)
    {
        Dimensions = locDimensions;
        totSize = 1;
        for(size_t n = 0; n < Rank; n++)
        {
            totSize *= Dimensions[n];
        }
        if(totSize)
        {
            locData = new T[totSize] ();
            allocated = true;
        }
        else
        {
            locData = nullptr;
            allocated = false;
        }
    }
    Tensor(T* data_ptr, std::array<size_t, Rank> locDimensions)
    {
        locData = data_ptr;
        allocated = false;
        Dimensions = locDimensions;
        totSize = 1;
        for(size_t n = 0; n < Rank; n++)
        {
            totSize *= Dimensions[n];
        }
    }
    Tensor(const Tensor<T, Rank>& locTensor)
    {
        Dimensions = locTensor.Dimensions;
        totSize = locTensor.totSize;
        if(totSize)
        {
            locData = new T[totSize] ();
            for(size_t n = 0; n < totSize; n++)
            {
                locData[n] = locTensor[n];
            }
            allocated = true;
        }
        else
        {
            locData = nullptr;
            allocated = false;
        }
    }
    T& operator()(const std::array<size_t,Rank> locPosition)
    {
#ifdef DEBUG
        for(size_t n = 0; n < Rank; n++)
        if(locPosition[n] >= Dimensions[n])
        {
            std::stringstream message;
            message << "Error in Tensor<T," << Rank << ">::operator()\n"
                    << "Access beyond the size of the tensor.\n"
                    << "Requested position [" << n << "] = " << locPosition[n]
                    << " > allowed size of " << Dimensions[n] << "\n";
            throw std::logic_error(message.str());
        }
#endif
        return locData[Index(locPosition)];
    }
    T const& operator()(const std::array<size_t,Rank> locPosition) const
    {
#ifdef DEBUG
        for(size_t n = 0; n < Rank; n++)
        if(locPosition[n] >= Dimensions[n])
        {
            std::stringstream message;
            message << "Error in Tensor<T," << Rank << ">::operator()\n"
                    << "Access beyond the size of the tensor.\n"
                    << "Requested position [" << n << "] = " << locPosition[n]
                    << " >= than allowed size of " << Dimensions[n]  << "\n";
            throw std::logic_error(message.str());
        }
#endif
        return locData[Index(locPosition)];
    }
    T& operator[](const size_t position)
    {
#ifdef DEBUG
        if(position >= totSize)
        {
            std::stringstream message;
            message << "Error in Tensor<T," << Rank << ">::operator[]\n"
                    << "Access beyond the total size of the tensor.\n"
                    << "Requested position " << position
                    << " >= tensor total size of " << totSize << "\n";
            throw std::logic_error(message.str());
        }
#endif
        return locData[position];
    }
    T const& operator[](const size_t position) const
    {
#ifdef DEBUG
        if(position >= totSize)
        {
            std::stringstream message;
            message << "Error in Tensor<T," << Rank << ">::operator[]\n"
                    << "Access beyond the total size of the tensor.\n"
                    << "Requested position "
                    << position << " >= tensor total size of " << totSize
                    << "\n" << std::endl;
            throw std::logic_error(message.str());
        }
#endif
        return locData[position];
    }
    void Assign(T* data_ptr, std::array<size_t, Rank> locDimensions)
    {
        locData = data_ptr;
        allocated = false;
        Dimensions = locDimensions;
        totSize = 1;
        for(size_t n = 0; n < Rank; n++)
        {
            totSize *= Dimensions[n];
        }
    }
    size_t Allocate(const std::array<size_t, Rank> locDimensions)
    {
        if(!allocated)
        {
            Dimensions = locDimensions;
            totSize = 1;
            for(size_t n = 0; n < Rank; n++)
            {
                totSize *= Dimensions[n];
            }
            if(totSize)
            {
                locData = new T[totSize] ();
                allocated = true;
                return sizeof(T)*totSize;
            }
            else
            {
                locData = nullptr;
                allocated = false;
            }
        }
        else
        {
            std::stringstream message;
            message << "ERROR: Tensor<T," << Rank
                    << ">: Allocation attempt of already allocated Tensor!\n"
                    << "Reallocate() should be used instead!\n";
#ifdef DEBUG
            throw std::logic_error(message.str());
#else
            std::cerr << message.str();
            std::exit(EXIT_FAILURE);
#endif
        }
        return 0;
    }
    void Reallocate(const std::array<size_t, Rank> locDimensions)
    {
        if(allocated)
        {
            delete[] locData;
        }

        Dimensions = locDimensions;
        totSize = 1;
        for(size_t n = 0; n < Rank; n++)
        {
            totSize *= Dimensions[n];
        }

        if(totSize)
        {
            locData = new T[totSize] ();
            allocated = true;
        }
        else
        {
            locData = nullptr;
            allocated = false;
        }
    }
    void set_to_value(T value)
    {
        if(allocated)
        for (size_t i = 0; i < totSize; i++)
        {
            locData[i] = value;
        }
    }
    void set_to_zero(void)
    {
        if(allocated)
        for (size_t i = 0; i < totSize; i++)
        {
            locData[i] = T();
        }
    }
    size_t size(size_t n) const
    {
        if(n < Rank)
        {
            return Dimensions[n];
        }
        else
        {
            return 0;
        }
    }
    size_t size(void) const
    {
        return totSize;
    }
    size_t rank(void) const
    {
        return Rank;
    }
    const T* data(void) const
    {
        return locData;
    }
    T* data(void)
    {
        return locData;
    }
    Tensor<T, Rank>& operator=(const Tensor<T, Rank>& locTensor)
    {
        if (this == &locTensor)
        {
            return *this;
        }
        if (!locTensor.allocated)
        {
            if(allocated)
            {
                delete[] locData;
            }
            allocated = false;
            return *this;
        }
        if (!allocated)
        {
            Dimensions = locTensor.Dimensions;
            totSize = 1;
            for(size_t n = 0; n < Rank; n++)
            {
                totSize *= Dimensions[n];
            }
            if(totSize)
            {
                locData = new T[totSize] ();
                allocated = true;
            }
        }
        if (locTensor.Dimensions == Dimensions)
        {
            for(size_t n = 0; n < totSize; n++)
            {
                locData[n] = locTensor[n];
            }
        }
        else
        {
            std::stringstream message;
            message << "ERROR: Tensor<T," << Rank
                    << ">: Different tensor size in assignment operator!\n";
#ifdef DEBUG
            throw std::logic_error(message.str());
#else
            std::cerr << message.str();
            std::exit(EXIT_FAILURE);
#endif
        }
            return *this;
    }
    template<typename T2>
    Tensor<T, Rank>& operator+=(const Tensor<T2, Rank>& locTensor)
    {
        if ((not allocated) and (locTensor.allocated))
        {
            Dimensions = locTensor.Dimensions;
            totSize = 1;
            for(size_t n = 0; n < Rank; n++)
            {
                totSize *= Dimensions[n];
            }
            if(totSize)
            {
                locData = new T[totSize] ();
                allocated = true;
            }
        }
        if (locTensor.Dimensions == Dimensions)
        {
            for(size_t n = 0; n < totSize; n++)
            {
                locData[n] += (T)locTensor[n];
            }
        }
        else
        {
            std::stringstream message;
            message << "ERROR: Tensor<T," << Rank
                << ">: Different tensor size in += operator!" << std::endl;

            message << "Dimensions lside: ";
                for (size_t i = 0; i < Rank; i++)
                    message << Dimensions[i] << " ";
            message << std::endl;
            message << "Dimensions rside: ";
                for (size_t i = 0; i < Rank; i++)
                    message << locTensor.Dimensions[i] << " \n";
#ifdef DEBUG
            throw std::logic_error(message.str());
#else
            std::cerr << message.str();
            std::exit(EXIT_FAILURE);
#endif
        }
        return *this;
    }
    template<typename T2>
    Tensor<T, Rank>& operator-=(const Tensor<T2, Rank>& locTensor)
    {
        if ((not allocated) and (locTensor.allocated))
        {
            Dimensions = locTensor.Dimensions;
            totSize = 1;
            for(size_t n = 0; n < Rank; n++)
            {
                totSize *= Dimensions[n];
            }
            if(totSize)
            {
                locData = new T[totSize] ();
                allocated = true;
            }
        }
        if (locTensor.Dimensions == Dimensions)
        {
            for(size_t n = 0; n < totSize; n++)
            {
                locData[n] += (T)locTensor[n];
            }
        }
        else
        {
            std::stringstream message;
            message << "ERROR: Tensor<T," << Rank
                << ">: Different tensor size in -= operator!"
                << std::endl;

            message << "Dimensions lside: ";
                for (size_t i = 0; i < Rank; i++)
                    message << Dimensions[i] << " ";
            message << std::endl;
            message << "Dimensions rside: ";
                for (size_t i = 0; i < Rank; i++)
                    message << locTensor.Dimensions[i] << " \n";
#ifdef DEBUG
            throw std::logic_error(message.str());
#else
            std::cerr << message.str();
            std::exit(EXIT_FAILURE);
#endif
        }
        return *this;
    }
    template<typename T2>
    Tensor<T, Rank>& operator/=(const T2 number)
    {
        if (allocated)
        for(size_t n = 0; n < totSize; n++)
        {
            locData[n] /= (T)number;
        }
        return *this;
    }
    template<typename T2>
    Tensor<T, Rank>& operator*=(const T2 number)
    {
        if (allocated)
        for(size_t n = 0; n < totSize; n++)
        {
            locData[n] *= (T)number;
        }
        return *this;
    }
    Tensor<T, Rank> operator+(const Tensor<T, Rank>& locTensor) const
    {
        if (locTensor.Dimensions == Dimensions)
        {
            Tensor<T,Rank> myReturn(locTensor);

            for(size_t n = 0; n < totSize; n++)
            {
                myReturn[n] += locData[n];
            }
            return myReturn;
        }
        else
        {
            std::stringstream message;
            message << "ERROR: Tensor<T," << Rank
                    << ">: Different tensor size in summation operator!\n";
#ifdef DEBUG
            throw std::logic_error(message.str());
#else
            std::cerr << message.str();
            std::exit(EXIT_FAILURE);
#endif
        }
        return *this;
    }
    Tensor<T, Rank> operator-(const Tensor<T, Rank>& locTensor) const
    {
        if (locTensor.Dimensions == Dimensions)
        {
            Tensor<T,Rank> myReturn(locTensor);

            for(size_t n = 0; n < totSize; n++)
            {
                myReturn[n] -= locData[n];
            }
            return myReturn;
        }
        else
        {
            std::stringstream message;
            message << "ERROR: Tensor<T," << Rank
                    << ">: Different tensor size in subtraction operator!\n";
#ifdef DEBUG
            throw std::logic_error(message.str());
#else
            std::cerr << message.str();
            std::exit(EXIT_FAILURE);
#endif
        }
        return *this;
    }
    Tensor<T, Rank> operator*(double val) const
    {
        Tensor<T, Rank> myReturn(*this);
        for(size_t n = 0; n < totSize; n++)
        {
            myReturn[n] *= val;
        }
        return myReturn;
    }
    Tensor<T, Rank> operator/(double val) const
    {
        Tensor<T, Rank> myReturn(*this);
        for(size_t n = 0; n < totSize; n++)
        {
            myReturn[n] /= val;
        }
        return myReturn;
    }
    void normalize()
    {
        double sum = 0.0;
        for (size_t n = 0; n < totSize; n++)
        {
            sum += locData[n];
        }
        for (size_t n = 0; n < totSize; n++)
        {
            locData[n] /= sum;
        }
    }
    bool IsAllocated() const
    {
        return allocated;
    }

 protected:
    std::array<size_t, Rank> Dimensions;
    size_t totSize;
    bool allocated;
    T*  locData;

 private:
    size_t Index(const std::array<size_t,Rank> position) const
    {
        size_t locIndex = *(position.begin());
        for(size_t n = 1; n < Rank; n++)
        {
            locIndex *= Dimensions[n];
            locIndex += *(position.begin() + n);
        }
        return locIndex;
    }
};

}// namespace openphase
#endif
