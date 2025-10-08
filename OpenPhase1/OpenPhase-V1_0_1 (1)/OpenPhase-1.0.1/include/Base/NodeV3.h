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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */


#ifndef NODEV3_H
#define NODEV3_H

//#include "Base/Includes.h"
#include "Base/dVector3.h"
namespace openphase
{

struct Vector3Entry                                                             /// Individual vector entry. Used in the NodeV class as a storage unit.
{
    dVector3 vector3;
    double& operator[](size_t index)
    {
        return vector3[index];
    }
    const double& operator[](size_t index) const
    {
        return vector3[index];
    }
    double& X(void)
    {
        return vector3[0];
    }
    const double& X(void) const
    {
        return vector3[0];
    }

    double& Y(void)
    {
        return vector3[1];
    }
    const double& Y(void) const
    {
        return vector3[1];
    }

    double& Z(void)
    {
        return vector3[2];
    }
    const double& Z(void) const
    {
        return vector3[2];
    }
    union                                                                       /// First index
    {
        size_t index;
        size_t indexA;
    };
    size_t     indexB;                                                          /// Second index
};

/***************************************************************/

class NodeV3                                                                    /// Basic class to store the vector valued quantities for each phase field. Provide access and manipulation methods for the entrys.
{
 public:
    NodeV3  operator+(const  NodeV3& n) const;                                  /// Plus operator.
    NodeV3  operator-(const  NodeV3& n) const;                                  /// Minus operator.
    NodeV3  operator*(const  double n) const;                                   /// Multiply operator.

    NodeV3  Xreflected() const;                                                 /// Returns the NodeV with the sign of X components of all vectors changed
    NodeV3  Yreflected() const;                                                 /// Returns the NodeV with the sign of Y components of all vectors changed
    NodeV3  Zreflected() const;                                                 /// Returns the NodeV with the sign of Z components of all vectors changed

    void     set     (const size_t n, const dVector3);                          /// Set all components using dVector3 type
    dVector3 get     (const size_t n) const;                                    /// Returns all components as dVector3 type
    void     add     (const size_t n, const dVector3);                          /// Add all components using dVector3 type

    void     set_asym(const size_t n, const size_t m, const dVector3);          /// Set entry assuming f(n,m) = -f(m,n)
    void     set_sym (const size_t n, const size_t m, const dVector3);          /// Set entry assuming f(n,m) = f(m,n)
    void     set     (const size_t n, const size_t m, const dVector3);          /// Set entry assuming f(n,m) != f(m,n)

    dVector3 get_asym(const size_t n, const size_t m) const;                    /// Return entry assuming f(n,m) = -f(m,n)
    dVector3 get_sym (const size_t n, const size_t m) const;                    /// Return entry assuming f(n,m) = f(m,n)
    dVector3 get     (const size_t n, const size_t m) const;                    /// Return entry assuming f(n,m) != f(m,n)

    void     add_asym(const size_t n, const size_t m, const dVector3);          /// Add to an entry assuming f(n,m) = -f(m,n)
    void     add_sym (const size_t n, const size_t m, const dVector3);          /// Add to an entry assuming f(n,m) = f(m,n)
    void     add     (const size_t n, const size_t m, const dVector3);          /// Add to an entry assuming f(n,m) != f(m,n)

    void    clear(){Vector3Fields.clear();};                                    /// Empty the vector fields.
    size_t  size() const {return Vector3Fields.size();};                        /// Return the size of vector fields.
    typedef std::vector<Vector3Entry>::iterator iterator;                       /// Iterator over the vector fields
    typedef std::vector<Vector3Entry>::const_iterator citerator;                /// Constant iterator over the vector fields
    iterator begin() {return Vector3Fields.begin();};                           /// Iterator to the begin of vector fields
    iterator end()   {return Vector3Fields.end();};                             /// Iterator to the end of vector fields
    citerator cbegin() const {return Vector3Fields.cbegin();};                  /// Constant iterator to the begin of vector fields
    citerator cend() const   {return Vector3Fields.cend();};                    /// Constant iterator to the end of vector fields

 protected:
 private:
    std::vector<Vector3Entry> Vector3Fields;                                    /// List of nonvanishing vector fields.
};

inline NodeV3 NodeV3::operator+(const  NodeV3& n) const
{
    NodeV3 result = n;
    for (citerator i = Vector3Fields.cbegin(); i < Vector3Fields.cend(); ++i)
    {
        result.add(i->indexA, i->indexB, i->vector3);
    }
    return result;
}

inline NodeV3 NodeV3::operator-(const NodeV3& n) const
{
    NodeV3 result = n*(-1.0);

    for (citerator i = Vector3Fields.cbegin(); i < Vector3Fields.cend(); ++i)
    {
        result.add(i->indexA, i->indexB, i->vector3);
    }
    return result;
}

inline NodeV3 NodeV3::operator*(const double n) const
{
    NodeV3 result = *this;
    for (iterator i = result.begin(); i < result.end(); ++i)
    {
        i->vector3 *= n;
    }
    return result;
}

inline NodeV3 NodeV3::Xreflected() const
{
    NodeV3 result = *this;

    for (iterator i = result.begin(); i < result.end(); ++i)
    {
        i->vector3[0] *= -1.0;
    }
    return result;
}

inline NodeV3 NodeV3::Yreflected() const
{
    NodeV3 result = *this;

    for (iterator i = result.begin(); i < result.end(); ++i)
    {
        i->vector3[1] *= -1.0;
    }
    return result;
}

inline NodeV3 NodeV3::Zreflected() const
{
    NodeV3 result = *this;

    for (iterator i = result.begin(); i < result.end(); ++i)
    {
        i->vector3[2] *= -1.0;
    }
    return result;
}

inline void NodeV3::set(const size_t n, const dVector3 value)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == n)
    {
        i->vector3 = value;
        return;
    }

    Vector3Entry NewEntry;
    NewEntry.index  = n;
    NewEntry.indexB = 0;
    NewEntry.vector3 = value;

    Vector3Fields.push_back(NewEntry);
}

inline dVector3 NodeV3::get(const size_t n) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->index == n)
    {
        return i->vector3;
    }
    return dVector3::ZeroVector();
}

inline void NodeV3::add(const size_t n, const dVector3 value)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == n)
    {
        i->vector3 += value;
        return;
    }

    Vector3Entry NewEntry;
    NewEntry.index  = n;
    NewEntry.indexB = 0;
    NewEntry.vector3 = value;

    Vector3Fields.push_back(NewEntry);
}

inline void NodeV3::set_asym(const size_t n, const size_t m, const dVector3 value)
{
    for (iterator i = begin(); i < end(); ++i)
    {
        if (i->indexA == n && i->indexB == m)
        {
            i->vector3 = value;
            return;
        }
        if (i->indexA == m && i->indexB == n)
        {
            i->vector3 = value*(-1.0);
            return;
        }
    }

    Vector3Entry NewEntry;
    NewEntry.index  = n;
    NewEntry.indexB = m;
    NewEntry.vector3 = value;

    Vector3Fields.push_back(NewEntry);
}

inline void NodeV3::set_sym(const size_t n, const size_t m, const dVector3 value)
{
    for (iterator i = begin(); i < end(); ++i)
    if ((i->indexA == n && i->indexB == m) or
        (i->indexA == m && i->indexB == n))
    {
        i->vector3 = value;
        return;
    }

    Vector3Entry NewEntry;
    NewEntry.index  = n;
    NewEntry.indexB = m;
    NewEntry.vector3 = value;

    Vector3Fields.push_back(NewEntry);
}

inline void NodeV3::set(const size_t n, const size_t m, const dVector3 value)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        i->vector3 = value;
        return;
    }

    Vector3Entry NewEntry;
    NewEntry.index  = n;
    NewEntry.indexB = m;
    NewEntry.vector3 = value;

    Vector3Fields.push_back(NewEntry);
}

inline dVector3 NodeV3::get_asym(const size_t n, const size_t m) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    {
        if (i->indexA == n && i->indexB == m)
        {
            return i->vector3;
        }
        if (i->indexA == m && i->indexB == n)
        {
            return i->vector3*(-1.0);
        }
    }
    return dVector3::ZeroVector();
}

inline dVector3 NodeV3::get_sym(const size_t n, const size_t m) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if ((i->indexA == n && i->indexB == m) or
        (i->indexA == m && i->indexB == n))
    {
        return i->vector3;
    }
    return dVector3::ZeroVector();
}

inline dVector3 NodeV3::get(const size_t n, const size_t m) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        return i->vector3;
    }
    return dVector3::ZeroVector();
}

inline void NodeV3::add_asym(const size_t n, const size_t m, const dVector3 value)
{
    for (iterator i = begin(); i < end(); ++i)
    {
        if (i->indexA == n && i->indexB == m)
        {
            i->vector3 += value;
            return;
        }
        if (i->indexA == m && i->indexB == n)
        {
            i->vector3 -= value;
            return;
        }
    }

    Vector3Entry NewEntry;
    NewEntry.index  = n;
    NewEntry.indexB = m;
    NewEntry.vector3 = value;

    Vector3Fields.push_back(NewEntry);
}

inline void NodeV3::add_sym(const size_t n, const size_t m, const dVector3 value)
{
    for (iterator i = begin(); i < end(); ++i)
    if ((i->indexA == n && i->indexB == m) or
        (i->indexA == m && i->indexB == n))
    {
        i->vector3 += value;
        return;
    }

    Vector3Entry NewEntry;
    NewEntry.index  = n;
    NewEntry.indexB = m;
    NewEntry.vector3 = value;

    Vector3Fields.push_back(NewEntry);
}

inline void NodeV3::add(const size_t n, const size_t m, const dVector3 value)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        i->vector3 += value;
        return;
    }

    Vector3Entry NewEntry;
    NewEntry.index  = n;
    NewEntry.indexB = m;
    NewEntry.vector3 = value;

    Vector3Fields.push_back(NewEntry);
}
} //namespace openphase
#endif
