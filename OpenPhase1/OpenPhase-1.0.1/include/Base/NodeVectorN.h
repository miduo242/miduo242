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
 *   Main contributors :    Muhammad Adil Ali; Oleg Shchyglo
 *
 */

#ifndef NodeVectorN_H
#define NodeVectorN_H

#include "Base/Includes.h"
namespace openphase
{

class VecEntry                                                              	///< Individual vector entry. Used in the NodeVectorN class as a storage unit.
{
public:
    VecEntry(): index(0)
    {
    };
    VecEntry(size_t size, size_t idx): value(size)
    {
        index = idx;
    };
    VecEntry(const VecEntry& other): value(other.value)
    {
        index = other.index;
    };
    dVectorN value;
    size_t index;
};

/***********************************************************************************/

class NodeVectorN                                                           	///< Basic class to store the vector valued quantities for each phase field.
{
public:

    NodeVectorN()
    {
        SIZE_X = 0;
    };
    NodeVectorN(const size_t size)
    {
        SIZE_X = size;
    };
    NodeVectorN(const NodeVectorN& other): VectorFields(other.VectorFields)
    {
        SIZE_X = other.SIZE_X;
    };
    void Allocate(const size_t size)
    {
        SIZE_X = size;
    }
    size_t size_of_array()
    {
        return SIZE_X;
    }
    void set_to_value(const size_t idx, const int value);                       ///< Set all entries to value for phase field index idx.
    void set_to_zero(const size_t idx);                                         ///< Set all entries to zero for phase field index idx.
    void set(const size_t idx, const size_t ii, const double value);            ///< Set component ii for phase field index idx to value.
    void set(const size_t idx, dVectorN value);                                 ///< Set dVectorN for phase field index idx.
    double   get(const size_t idx, const size_t ii) const;                      ///< Return value of component ii for phase field index idx.
    dVectorN get(const size_t idx) const;                                       ///< Return values for phase field index idx.
    void         add(const size_t idx, dVectorN value);                         ///< Increment value using two indeces.
    NodeVectorN  add(const NodeVectorN& n) const;                               ///< Add two nodes.
    void         add(const size_t idx, const size_t ii, const double value);    ///< Add in component ii for phase field index idx to value.

    NodeVectorN& operator=(const NodeVectorN& n);
    NodeVectorN& operator+=( NodeVectorN n);
    NodeVectorN& operator-=( NodeVectorN n);
    NodeVectorN& operator*=(const double n);
    NodeVectorN  operator+(const  NodeVectorN& n) const;                        ///< Plus operator. Takes as input another Node type entry.
    NodeVectorN  operator-(const  NodeVectorN& n) const;                        ///< Minus operator. Takes as input another Node type entry.
    NodeVectorN  operator*(const double n) const;                               ///< Multiply all fields by a number.

    void   clear(){VectorFields.clear();};                                      ///< Empties the vector fields.
    size_t  size() const {return VectorFields.size();};                         ///< Returns the size of vector fields.

    typedef typename std::vector<VecEntry>::iterator iterator;                  ///< Iterator over the vector fields
    typedef typename std::vector<VecEntry>::const_iterator citerator;           ///< Constant iterator over the vector fields
    iterator begin() {return VectorFields.begin();};                            ///< Iterator to the begin of vector fields
    iterator end()   {return VectorFields.end();};                              ///< Iterator to the end of vector fields
    citerator cbegin() const {return VectorFields.begin();};                    ///< Constant iterator to the begin of vector fields
    citerator cend()   const {return VectorFields.end();};                      ///< Constant iterator to the end of vector fields

    void pack(std::vector<double>& buffer)
    {
        buffer.push_back(SIZE_X);
        buffer.push_back(VectorFields.size());
        for (size_t i = 0; i < VectorFields.size(); ++i)
        {
            buffer.push_back(VectorFields[i].index);
            for (size_t n = 0; n < SIZE_X; ++n)
            {
                buffer.push_back(VectorFields[i].value[n]);
            }
        }
    }
    void unpack(std::vector<double>& buffer, size_t& it)
    {
        SIZE_X = buffer[it]; ++it;
        VectorFields.resize(buffer[it]); ++it;
        for (size_t i = 0; i < VectorFields.size(); ++i)
        {
            VectorFields[i].index = buffer[it]; ++it;
            VectorFields[i].value.Allocate(SIZE_X);
            for (size_t n = 0; n < SIZE_X; ++n)
            {
                VectorFields[i].value[n] = buffer[it]; ++it;
            }
        }
    }
    void Write(std::ostream& outp)
    {
        size_t size = VectorFields.size();
        outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&SIZE_X), sizeof(size_t));

        for(auto &Field : VectorFields)
        {
            outp.write(reinterpret_cast<const char*>(&Field.index), sizeof(size_t));
            Field.value.Write(outp);
        }
    }
    void Read(std::istream& inp)
    {
        size_t size = 0;
        inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&SIZE_X), sizeof(size_t));
        VectorFields.resize(size);

        for(auto &Field : VectorFields)
        {
            inp.read(reinterpret_cast<char*>(&Field.index), sizeof(size_t));
            Field.value.Read(inp);
        }

    }
protected:
private:
    std::vector<VecEntry> VectorFields;
    size_t SIZE_X;
};

inline void NodeVectorN::set_to_value(const size_t idx, const int value)
{
    for (iterator i = begin(); i < end(); ++i)
    {
        if(i->index == idx)
        {
            i->value.set_to_value(value);
            return;
        }
    }
    VecEntry NewEntry(SIZE_X, idx);
    NewEntry.value.set_to_value(value);
    VectorFields.push_back(NewEntry);
}

inline void NodeVectorN::set_to_zero(const size_t idx)
{
    for (iterator i = begin(); i < end(); ++i)
    {
        if(i->index == idx)
        {
            for(unsigned int ii = 0; ii < SIZE_X; ii++)
            {
                i->value[ii] = 0.0;
            }
            return;
        }
    }
    VecEntry NewEntry(SIZE_X,idx);
    NewEntry.value.set_to_zero();
    VectorFields.push_back(NewEntry);
}

inline void NodeVectorN::set(const size_t idx, const size_t ii, const double value)
{
    for (iterator i = begin(); i < end(); ++i)
    {
        if(i->index == idx)
        {
            i->value[ii] = value;
            return;
        }
    }
    VecEntry NewEntry(SIZE_X,idx);
    NewEntry.value.set_to_zero();
    NewEntry.value[ii] = value;
    VectorFields.push_back(NewEntry);
}

inline void NodeVectorN::set(const size_t idx, dVectorN value)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == idx)
    {
        if(SIZE_X != value.size())
        {
            std::cerr << "Error in NodeVectorN::set\n"
                    << "Sizes of the vectors are not equal. rhs.size() = "
                    <<  SIZE_X << " and lhs.size() = " << value.size()
                    << "\nTerminating!!!\n";;
            std::terminate();
        }
        else
        {
            i->value = value;
            return;
        }
    }
    VecEntry NewEntry(SIZE_X,idx);
    NewEntry.value = value;
    VectorFields.push_back(NewEntry);
}

inline double NodeVectorN::get(const size_t idx, const size_t ii) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->index == idx)
    {
        return i->value[ii];
    }
    return 0.0;
}

inline dVectorN NodeVectorN::get(const size_t idx) const
{
    dVectorN returndV(SIZE_X); returndV.set_to_zero();
    for (citerator i = cbegin(); i < cend(); ++i)
    {
        if (i->index == idx)
        {
            returndV = i->value;
            break;
        }
    }
    return returndV;
}

inline void NodeVectorN::add(const size_t idx, dVectorN value)
{
    for (auto i = begin(); i < end(); ++i)
    if (i->index == idx)
    {
        for(size_t ss = 0; ss < SIZE_X; ++ss)
        {
            i->value[ss] += value[ss];
        }
        return;
    }

    VecEntry NewEntry(SIZE_X,idx);
    NewEntry.value = value;
    VectorFields.push_back(NewEntry);
}

inline NodeVectorN NodeVectorN::add(const NodeVectorN& n) const
{
    NodeVectorN result = n;
    for (auto i = cbegin(); i < cend(); ++i)
    {
        result.add(i->index, i->value);
    }
    return result;
}

inline void NodeVectorN::add(const size_t idx, const size_t ii, const double value)
{
    for (iterator i = begin(); i < end(); ++i)
    {
        if(i->index == idx)
        {
                i->value[ii] += value;
                return;
        }
    }
    VecEntry NewEntry(SIZE_X,idx);
    NewEntry.value.set_to_zero();
    NewEntry.value[ii] += value;
    VectorFields.push_back(NewEntry);
}


inline NodeVectorN& NodeVectorN::operator=(const NodeVectorN& n)
{
    VectorFields = n.VectorFields;
    return *this;
}

inline NodeVectorN& NodeVectorN::operator+=( NodeVectorN n)
{
    for (auto beta : n.VectorFields)
    {
        bool inside = false;
        for (auto alpha : VectorFields)
        {
            if (alpha.index == beta.index)
            {
                alpha.value += beta.value;
                inside = true;
                break;
            }
        }
        if(!inside)
        {
            VecEntry NewEntry(beta.value.size(),beta.index);
            NewEntry.value = beta.value;
            VectorFields.push_back(NewEntry);
        }
    }
    return *this;
}

inline NodeVectorN& NodeVectorN::operator-=(NodeVectorN n)
{
    for (auto beta : n.VectorFields)
    {
        bool inside = false;
        for (auto alpha : VectorFields)
        {
            if (alpha.index == beta.index)
            {
                alpha.value -= beta.value;
                inside = true;
                break;
            }
        }
        if(!inside)
        {
            VecEntry NewEntry(beta.value.size(),beta.index);
            NewEntry.value = beta.value;
            VectorFields.push_back(NewEntry);
        }
    }
    return *this;
}

inline NodeVectorN& NodeVectorN::operator*=(const double n)
{
    for (iterator i = begin(); i < end(); ++i)
    {
        i->value *= n;
    }
    return *this;
}

inline NodeVectorN NodeVectorN::operator+(const NodeVectorN& n) const
{
    NodeVectorN result = n;
    for (auto i = cbegin(); i < cend(); ++i)
    {
        result.add(i->index, i->value);
    }
    return result;
}

inline NodeVectorN NodeVectorN::operator-(const NodeVectorN& n) const
{
    NodeVectorN result = n;
    for (auto i = result.begin(); i < result.end(); ++i)
    for(size_t ss = 0; ss < SIZE_X; ++ss)
    {
        i->value[ss] = -i->value[ss];
    }

    for (auto i = cbegin(); i < cend(); ++i)
    {
        result.add(i->index, i->value);
    }
    return result;
}

inline NodeVectorN NodeVectorN::operator*(const double n) const
{
    NodeVectorN result = *this;
    for(auto i = result.begin(); i < result.end(); ++i)
    {
        i->value *= n;
    }
    return result;
}

} //namespace openphase
#endif
