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
 *   Main contributors :   Oleg Shchyglo, Philipp Engels, Raphael Schiedung
 *
 */

#ifndef NODEMATRIX3X3_H
#define NODEMATRIX3X3_H

#include "Base/Includes.h"
#include "ElasticityTensors.h"
namespace openphase
{

struct Matrix3x3Entry                                                           /// Individual entries. Used in the NodeMatrix3x3 class as a storage unit.
{
    dMatrix3x3 PhaseMatrix;
    union                                                                       /// First index
    {
        size_t index;
        size_t indexA;
    };
    size_t     indexB;                                                          /// Second index
};

/***************************************************************/
class NodeMatrix3x3                                                             /// Basic class to store the matrices for each phase field. Provide access and manipulation methods for the entrys.
{
 public:
    void       add (const size_t idx, const dMatrix3x3 inMatrix);               /// Adds Matrix3X3 with index idx (author: Raphael Schiedung)
    dMatrix3x3 get (const size_t idx) const;                                    /// Return value of component idx for phase field index n.
    void       set (const size_t idx, const dMatrix3x3 inMatrix);               /// Set components idx for phase field index n.

    void       add (const size_t idx, const size_t idy, const dMatrix3x3 inMatrix);/// Adds Matrix3X3 with index idx and idy (author: Raphael Schiedung)
    dMatrix3x3 get (const size_t idx, const size_t idy) const;                     /// Returns Matrix3X3 with index idx and idy (author: Raphael Schiedung)
    void       set (const size_t idx, const size_t idy, const dMatrix3x3 inMatrix);/// Stores Matrix3X3 with index idx and idy (author: Raphael Schiedung)

    void   clear(){MatrixFields.clear();};                                      /// Empties the vector fields.
    size_t  size() const {return MatrixFields.size();};                         /// Returns the size of vector fields.
    typedef std::vector<Matrix3x3Entry>::iterator iterator;                     /// Iterator over the vector fields
    typedef std::vector<Matrix3x3Entry>::const_iterator citerator;              /// Constant iterator over the vector fields
    iterator begin() {return MatrixFields.begin();};                            /// Iterator to the begin of vector fields
    iterator end()   {return MatrixFields.end();};                              /// Iterator to the end of vector fields
    citerator cbegin() const {return MatrixFields.begin();};                    /// Constant iterator to the begin of vector fields
    citerator cend()   const {return MatrixFields.end();};                      /// Constant iterator to the end of vector fields

 protected:
 private:
    std::vector<Matrix3x3Entry> MatrixFields;                                   /// List of nonvanishing vector fields.
};

inline void NodeMatrix3x3::add(const size_t idx, const dMatrix3x3 inMatrix)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == idx)
    {
        i->PhaseMatrix += inMatrix;
        return;
    }

    Matrix3x3Entry NewEntry;
    NewEntry.index  = idx;
    NewEntry.indexB = 0;
    NewEntry.PhaseMatrix = inMatrix;

    MatrixFields.push_back(NewEntry);
}

inline dMatrix3x3 NodeMatrix3x3::get(const size_t idx) const
{
    dMatrix3x3 ZeroMatrix;
    ZeroMatrix.set_to_zero();
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->index == idx)
    {
        return i->PhaseMatrix;
    }
    return ZeroMatrix;
}

inline void NodeMatrix3x3::set(const size_t idx, const dMatrix3x3 inMatrix)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == idx)
    {
        i->PhaseMatrix = inMatrix;
        return;
    }

    Matrix3x3Entry NewEntry;
    NewEntry.index  = idx;
    NewEntry.indexB = 0;
    NewEntry.PhaseMatrix = inMatrix;

    MatrixFields.push_back(NewEntry);
}

inline void NodeMatrix3x3::add(const size_t idx, const size_t idy,
        const dMatrix3x3 inMatrix)
{
    for (iterator i = begin(); i < end(); ++i)
    if ((i->indexA == idx) and (i->indexB == idy))
    {
        i->PhaseMatrix += inMatrix;
        return;
    }

    Matrix3x3Entry NewEntry;
    NewEntry.indexA = idx;
    NewEntry.indexB = idy;
    NewEntry.PhaseMatrix = inMatrix;

    MatrixFields.push_back(NewEntry);
}

inline dMatrix3x3 NodeMatrix3x3::get(const size_t idx, const size_t idy) const

{
    dMatrix3x3 ZeroMatrix;
    ZeroMatrix.set_to_zero();
    for (citerator i = cbegin(); i < cend(); ++i)
    if ((i->indexA == idx) and (i->indexB == idy))
    {
        return i->PhaseMatrix;
    }
    return ZeroMatrix;
}

inline void NodeMatrix3x3::set(const size_t idx, const size_t idy,
        const dMatrix3x3 inMatrix)
{
    for (iterator i = begin(); i < end(); ++i)
    if ((i->indexA == idx) and (i->indexB == idy))
    {
        i->PhaseMatrix = inMatrix;
        return;
    }

    Matrix3x3Entry NewEntry;
    NewEntry.indexA = idx;
    NewEntry.indexB = idy;
    NewEntry.PhaseMatrix = inMatrix;

    MatrixFields.push_back(NewEntry);
}
} //namespace openphase
#endif
