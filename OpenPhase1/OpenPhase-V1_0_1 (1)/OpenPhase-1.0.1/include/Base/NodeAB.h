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
 *   File created :   2009
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */

#ifndef NODEAB_H
#define NODEAB_H

#include "Base/Includes.h"

namespace openphase
{
/**********************************************************/

struct DoubleIndexFieldEntry                                                    ///< Structure for storing the field entry with two indices. Used in the NodeAB class as a storage unit.
{
    size_t indexA;                                                              ///< First index.
    size_t indexB;                                                              ///< Second index.
    double value1;                                                              ///< First value.
    double value2;                                                              ///< Second value.

    DoubleIndexFieldEntry()                                                     ///< Default constructor.
    {
        indexA = 0;
        indexB = 0;
        value1 = 0.0;
        value2 = 0.0;
    }
    DoubleIndexFieldEntry(const DoubleIndexFieldEntry& entry)                   ///< Copy constructor.
    {
        indexA = entry.indexA;
        indexB = entry.indexB;
        value1 = entry.value1;
        value2 = entry.value2;
    }
    DoubleIndexFieldEntry& operator=(const DoubleIndexFieldEntry& entry)        ///< Copy operator.
    {
        indexA = entry.indexA;
        indexB = entry.indexB;
        value1 = entry.value1;
        value2 = entry.value2;
        return *this;
    }
};

/**********************************************************/
class NodeAB                                                                    ///< Stores the fields with two indices and two values at a grid point. Provides access and manipulation methods for the stored field entries.
{
 public:
    NodeAB()                                                                    ///< Constructor used to allocate space for at least 3 fields.
    {
        Fields.reserve(3);
    }
    NodeAB(const NodeAB& n)                                                     ///< Copy constructor.
    {
        Fields.reserve(3);
        Fields = n.Fields;
    }
    NodeAB  operator+(const NodeAB& n) const;                                   ///< Plus operator. Takes as input another Node type entry.
    NodeAB  operator-(const NodeAB& n) const;                                   ///< Minus operator. Takes as input another Node type entry.
    NodeAB  operator*(const double  n) const;                                   ///< Multiplies all fields by a number.
    NodeAB& operator=(const NodeAB& n);                                         ///< Assignment operator

    NodeAB& operator+=(const NodeAB& n);                                        ///< Plus-equal operator. Takes as input another Node type entry.
    NodeAB& operator-=(const NodeAB& n);                                        ///< Minus-equal operator. Takes as input another Node type entry.
    NodeAB& operator*=(const double  n);                                        ///< Multiplies all fields by a number.

    void    set1(const size_t n, const size_t m, const double value);           ///< Sets value1 following strict index order.
    void    set2(const size_t n, const size_t m, const double value);           ///< Sets value2 following strict index order.
    double  get1(const size_t n, const size_t m) const;                         ///< Returns value1 following strict index order.
    double  get2(const size_t n, const size_t m) const;                         ///< Returns value2 following strict index order.
    void    add1(const size_t n, const size_t m,  const double value);          ///< Increments value1 following strict index order.
    void    add2(const size_t n, const size_t m,  const double value);          ///< Increments value2 following strict index order.

    void    set_sym1(const size_t n, const size_t m, const double value);       ///< Sets value1 in symmetric case: f(n,m) = f(m,n).
    void    set_sym2(const size_t n, const size_t m, const double value);       ///< Sets value2 in symmetric case: f(n,m) = f(m,n).
    double  get_sym1(const size_t n, const size_t m) const;                     ///< Returns value1 in symmetric case: f(n,m) = f(m,n).
    double  get_sym2(const size_t n, const size_t m) const;                     ///< Returns value2 in symmetric case: f(n,m) = f(m,n).
    void    add_sym1(const size_t n, const size_t m,  const double value);      ///< Increments value1 in symmetric case: f(n,m) = f(m,n).
    void    add_sym2(const size_t n, const size_t m,  const double value);      ///< Increments value2 in symmetric case: f(n,m) = f(m,n).

    void    set_asym1(const size_t n, const size_t m, const double value);      ///< Sets value1 in antisymmetric case: f(n,m) = -f(m,n).
    void    set_asym2(const size_t n, const size_t m, const double value);      ///< Sets value2 in antisymmetric case: f(n,m) = -f(m,n).
    double  get_asym1(const size_t n, const size_t m) const;                    ///< Returns value1 in antisymmetric case: f(n,m) = -f(m,n).
    double  get_asym2(const size_t n, const size_t m) const;                    ///< Returns value2 in antisymmetric case: f(n,m) = -f(m,n).
    void    add_asym1(const size_t n, const size_t m,  const double value);     ///< Increments value1 in antisymmetric case: f(n,m) = -f(m,n).
    void    add_asym2(const size_t n, const size_t m,  const double value);     ///< Increments value2 in antisymmetric case: f(n,m) = -f(m,n).

    void    set_pair(const size_t n, const size_t m, const double value1, const double value2);///< Sets value1 and value2 following strict index order.
    void    add_pair(const size_t n, const size_t m, const double value1, const double value2);///< Increments value1 and value2 following strict index order.
    std::pair<double,double> get_pair(const size_t n, const size_t m) const;                   ///< Returns pair<value1,value2> following strict index order.

    void    set_sym_pair(const size_t n, const size_t m, const double value1, const double value2);///< Sets value1 and value2 in symmetric case: f(n,m) = f(m,n).
    void    add_sym_pair(const size_t n, const size_t m, const double value1, const double value2);///< Increments value1 and value2 in symmetric case: f(n,m) = f(m,n).
    std::pair<double,double> get_sym_pair(const size_t n, const size_t m) const;                   ///< Returns pair<value1,value2>  in symmetric case: f(n,m) = f(m,n).

    void    set_asym_pair(const size_t n, const size_t m, const double value1, const double value2);///< Sets value1 and value2 in antisymmetric case: f(n,m) = -f(m,n).
    void    add_asym_pair(const size_t n, const size_t m, const double value1, const double value2);///< Increments value1 and value2 in antisymmetric case: f(n,m) = -f(m,n).
    std::pair<double,double> get_asym_pair(const size_t n, const size_t m) const;                   ///< Returns pair<value1,value2>  in antisymmetric case: f(n,m) = -f(m,n).

    void    add1(const NodeAB& value);                                          ///< Add value1 of two nodes following strict index order.
    void    add2(const NodeAB& value);                                          ///< Add value2 of two nodes following strict index order.

    void    add_sym1(const NodeAB& value);                                      ///< Add value1 of two nodes in symmetric case: f(n,m) = f(m,n).
    void    add_sym2(const NodeAB& value);                                      ///< Add value2 of two nodes in symmetric case: f(n,m) = f(m,n).

    void    add_asym1(const NodeAB& value);                                     ///< Add value1 of two nodes in antisymmetric case: f(n,m) = -f(m,n).
    void    add_asym2(const NodeAB& value);                                     ///< Add value2 of two nodes in antisymmetric case: f(n,m) = -f(m,n).

    void    add_sym_pairs(const NodeAB& value);                                 ///< Add value1 and value2 of two nodes in symmetric case: f(n,m) = f(m,n).
    void    add_asym_pairs(const NodeAB& value);                                ///< Add value1 and value2 of two nodes in antisymmetric case: f(n,m) = -f(m,n).

    void    add_sym1_exist(const NodeAB& value);                                ///< Add only value1 existing in two nodes simultaneously in symmetric case: f(n,m) = f(m,n).
    void    add_sym2_exist(const NodeAB& value);                                ///< Add only value2 existing in two nodes simultaneously in symmetric case: f(n,m) = f(m,n).

    void    add_asym1_exist(const NodeAB& value);                               ///< Add only value1 existing in two nodes simultaneously in antisymmetric case: f(n,m) = -f(m,n).
    void    add_asym2_exist(const NodeAB& value);                               ///< Add only value2 existing in two nodes simultaneously in antisymmetric case: f(n,m) = -f(m,n).

    void    clear() {Fields.clear();};                                          ///< Empties the field storage. Sets flag to 0.
    size_t  size() const {return Fields.size();};                               ///< Returns the size of storage.
    typedef std::vector<DoubleIndexFieldEntry>::iterator iterator;              ///< Iterator over storage vector
    typedef std::vector<DoubleIndexFieldEntry>::const_iterator citerator;       ///< Constant iterator over storage vector
    iterator  begin() {return Fields.begin();};                                 ///< Iterator to the begin of storage vector
    iterator  end()   {return Fields.end();};                                   ///< Iterator to the end of storage vector
    citerator cbegin() const {return Fields.cbegin();};                         ///< Constant iterator to the begin of storage vector
    citerator cend()   const {return Fields.cend();};                           ///< Constant iterator to the end of storage vector
    iterator erase(iterator it) {return Fields.erase(it);};                     ///< Erase a single record pointed by iterator it
    DoubleIndexFieldEntry& front(void) {return Fields.front();};                ///< Reference to the first FieldEntry.
    const DoubleIndexFieldEntry& front(void) const {return Fields.front();};    ///< Constant reference to the first FieldEntry.

    void pack(std::vector<double>& buffer);
    void unpack(std::vector<double>& buffer, size_t& it);
    void Read(std::istream& inp);
    void Write(std::ostream& outp) const;

 protected:
 private:
    std::vector<DoubleIndexFieldEntry> Fields;                                  ///< Fields storage vector.
};

/***************************************************************/

inline NodeAB NodeAB::operator+(const NodeAB& n) const
{
    NodeAB result = n;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add_pair(i->indexA, i->indexB, i->value1, i->value2);
    }
    return result;
}

inline NodeAB NodeAB::operator-(const NodeAB& n) const
{
    NodeAB result = n;

    for (auto i = result.begin(); i < result.end(); ++i)
    {
        i->value1 *= -1.0;
        i->value2 *= -1.0;
    }

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add_pair(i->indexA, i->indexB, i->value1, i->value2);
    }
    return result;
}

inline NodeAB& NodeAB::operator+=(const NodeAB& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_pair(i->indexA, i->indexB, i->value1, i->value2);
    }
    return *this;
}

inline NodeAB& NodeAB::operator-=(const NodeAB& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_pair(i->indexA, i->indexB, -i->value1, -i->value2);
    }
    return *this;
}

inline NodeAB NodeAB::operator*(const double n) const
{
    NodeAB result = *this;

    for(auto i = result.begin(); i < result.end(); ++i)
    {
        i->value1 *= n;
        i->value2 *= n;
    }
    return result;
}

inline NodeAB& NodeAB::operator*=(const double n)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    {
        i->value1 *= n;
        i->value2 *= n;
    }
    return *this;
}

inline NodeAB& NodeAB::operator=(const NodeAB& n)
{
    Fields = n.Fields;
    return *this;
}

//*****************************************************************************

inline void NodeAB::set1(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m))
    {
        i->value1 = value;
        return;
    }

    DoubleIndexFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value;
    NewEntry.value2 = 0.0;

    Fields.push_back(NewEntry);
}

inline void NodeAB::set2(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m))
    {
        i->value2 = value;
        return;
    }

    DoubleIndexFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = 0.0;
    NewEntry.value2 = value;

    Fields.push_back(NewEntry);
}

inline double NodeAB::get1(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m))
    {
        return i->value1;
    }
    return 0.0;
}

inline double NodeAB::get2(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m))
    {
        return i->value2;
    }
    return 0.0;
}

inline void NodeAB::set_sym1(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->value1 = value;
        return;
    }

    DoubleIndexFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value;
    NewEntry.value2 = 0.0;

    Fields.push_back(NewEntry);
}

inline void NodeAB::set_sym2(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->value2 = value;
        return;
    }

    DoubleIndexFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = 0.0;
    NewEntry.value2 = value;

    Fields.push_back(NewEntry);
}

inline void NodeAB::add_sym1(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->value1 += value;
        return;
    }

    DoubleIndexFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value;
    NewEntry.value2 = 0.0;

    Fields.push_back(NewEntry);
}

inline void NodeAB::add_sym2(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->value2 += value;
        return;
    }

    DoubleIndexFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = 0.0;
    NewEntry.value2 = value;

    Fields.push_back(NewEntry);
}

/***************************************************************/
inline double NodeAB::get_sym1(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        return i->value1;
    }
    return 0.0;
}

inline double NodeAB::get_sym2(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        return i->value2;
    }
    return 0.0;
}

inline void NodeAB::set_asym1(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        i->value1 = value;
        return;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        i->value1 = -value;
        return;
    }

    DoubleIndexFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value;
    NewEntry.value2 = 0.0;

    Fields.push_back(NewEntry);
}

inline void NodeAB::set_asym2(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        i->value2 = value;
        return;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        i->value2 = -value;
        return;
    }

    DoubleIndexFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = 0.0;
    NewEntry.value2 = value;

    Fields.push_back(NewEntry);
}

inline void NodeAB::add_asym1(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        i->value1 += value;
        return;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        i->value1 -= value;
        return;
    }

    DoubleIndexFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value;
    NewEntry.value2 = 0.0;

    Fields.push_back(NewEntry);
}

inline void NodeAB::add_asym2(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        i->value2 += value;
        return;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        i->value2 -= value;
        return;
    }

    DoubleIndexFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = 0.0;
    NewEntry.value2 = value;

    Fields.push_back(NewEntry);
}

inline double NodeAB::get_asym1(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        return i->value1;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        return -i->value1;
    }
    return 0.0;
}

inline double NodeAB::get_asym2(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        return i->value2;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        return -i->value2;
    }
    return 0.0;
}

inline void NodeAB::set_pair(const size_t n, const size_t m, const double value1, const double value2)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m))
    {
        i->value1 = value1;
        i->value2 = value2;
        return;
    }

    DoubleIndexFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value1;
    NewEntry.value2 = value2;

    Fields.push_back(NewEntry);
}

inline void NodeAB::add_pair(const size_t n, const size_t m, const double value1, const double value2)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m))
    {
        i->value1 += value1;
        i->value2 += value2;
        return;
    }

    DoubleIndexFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value1;
    NewEntry.value2 = value2;

    Fields.push_back(NewEntry);
}

inline std::pair<double,double> NodeAB::get_pair(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m))
    {
        return std::pair<double,double>(i->value1, i->value2);
    }
    return std::pair<double,double>(0.0, 0.0);
}

inline void NodeAB::set_sym_pair(const size_t n, const size_t m, const double value1, const double value2)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->value1 = value1;
        i->value2 = value2;
        return;
    }

    DoubleIndexFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value1;
    NewEntry.value2 = value2;

    Fields.push_back(NewEntry);
}

inline void NodeAB::add_sym_pair(const size_t n, const size_t m, const double value1, const double value2)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->value1 += value1;
        i->value2 += value2;
        return;
    }

    DoubleIndexFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value1;
    NewEntry.value2 = value2;

    Fields.push_back(NewEntry);
}

inline std::pair<double,double> NodeAB::get_sym_pair(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        return std::pair<double,double>(i->value1, i->value2);
    }
    return std::pair<double,double>(0.0, 0.0);
}

inline void NodeAB::set_asym_pair(const size_t n, const size_t m, const double value1, const double value2)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        i->value1 = value1;
        i->value2 = value2;
        return;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        i->value1 = -value1;
        i->value2 = -value2;
        return;
    }

    DoubleIndexFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value1;
    NewEntry.value2 = value2;

    Fields.push_back(NewEntry);
}

inline void NodeAB::add_asym_pair(const size_t n, const size_t m, const double value1, const double value2)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        i->value1 += value1;
        i->value2 += value2;
        return;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        i->value1 -= value1;
        i->value2 -= value2;
        return;
    }

    DoubleIndexFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value1;
    NewEntry.value2 = value2;

    Fields.push_back(NewEntry);
}

inline std::pair<double,double> NodeAB::get_asym_pair(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        return std::pair<double,double>(i->value1, i->value2);
    }
    else if(i->indexA == m and i->indexB == n)
    {
        return std::pair<double,double>(-i->value1, -i->value2);
    }
    return std::pair<double,double>(0.0, 0.0);
}

//*****************************************************************************

inline void NodeAB::add_sym1(const NodeAB& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_sym1(i->indexA, i->indexB, i->value1);
    }
}

inline void NodeAB::add_sym2(const NodeAB& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_sym2(i->indexA, i->indexB, i->value2);
    }
}

inline void NodeAB::add_asym1(const NodeAB& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_asym1(i->indexA, i->indexB, i->value1);
    }
}

inline void NodeAB::add_asym2(const NodeAB& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_asym2(i->indexA, i->indexB, i->value2);
    }
}

inline void NodeAB::add_sym_pairs(const NodeAB& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_sym_pair(i->indexA, i->indexB, i->value1, i->value2);
    }
}

inline void NodeAB::add_asym_pairs(const NodeAB& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_asym_pair(i->indexA, i->indexB, i->value1, i->value2);
    }
}

inline void NodeAB::add_sym1_exist(const NodeAB& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        double locVal = get_sym1(i->indexA, i->indexB);
        if(locVal != 0.0 && i->value1 != 0.0)
        {
            set_sym1(i->indexA, i->indexB, i->value1 + locVal);
        }
    }
}

inline void NodeAB::add_asym1_exist(const NodeAB& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        double locVal = get_asym1(i->indexA, i->indexB);
        if(locVal != 0.0 && i->value1 != 0.0)
        {
            set_asym1(i->indexA, i->indexB, i->value1 + locVal);
        }
    }
}

inline void NodeAB::add_sym2_exist(const NodeAB& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        double locVal = get_sym2(i->indexA, i->indexB);
        if(locVal != 0.0 && i->value2 != 0.0)
        {
            set_sym2(i->indexA, i->indexB, i->value2 + locVal);
        }
    }
}

inline void NodeAB::add_asym2_exist(const NodeAB& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        double locVal = get_asym2(i->indexA, i->indexB);
        if(locVal != 0.0 && i->value2 != 0.0)
        {
            set_asym2(i->indexA, i->indexB, i->value2 + locVal);
        }
    }
}

inline void NodeAB::pack(std::vector<double>& buffer)
{
    buffer.push_back(Fields.size());
    for(auto it = Fields.begin(); it != Fields.end();++it)
    {
        buffer.push_back(it->indexA);
        buffer.push_back(it->indexB);
        buffer.push_back(it->value1);
        buffer.push_back(it->value2);
    }
}

inline void NodeAB::unpack(std::vector<double>& buffer, size_t& it)
{
    clear();
    size_t size = buffer[it]; ++it;
    Fields.resize(size);
    for(size_t i = 0; i < size; ++i)
    {
        Fields[i].indexA = buffer[it]; ++it;
        Fields[i].indexB = buffer[it]; ++it;
        Fields[i].value1 = buffer[it]; ++it;
        Fields[i].value2 = buffer[it]; ++it;
    }
}

inline void NodeAB::Read(std::istream& inp)
{
    size_t size = 0;
    inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));
    Fields.resize(size);

    for(auto &Field : Fields)
    {
        inp.read(reinterpret_cast<char*>(&Field.indexA), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&Field.indexB), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&Field.value1), sizeof(double));
        inp.read(reinterpret_cast<char*>(&Field.value2), sizeof(double));
    }
}

inline void NodeAB::Write(std::ostream& outp) const
{
    size_t size = Fields.size();
    outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

    for(auto &Field : Fields)
    {
        outp.write(reinterpret_cast<const char*>(&Field.indexA), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&Field.indexB), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&Field.value1), sizeof(double));
        outp.write(reinterpret_cast<const char*>(&Field.value2), sizeof(double));
    }
}

}//namespace openphase
#endif
