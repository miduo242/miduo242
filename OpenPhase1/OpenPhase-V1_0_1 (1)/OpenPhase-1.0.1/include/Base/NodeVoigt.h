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
 *   Main contributors :   Oleg Shchyglo, Philipp Engels
 *
 */

#ifndef NODEVOIGT_H
#define NODEVOIGT_H

#include "Base/Includes.h"
namespace openphase
{

struct VoigtEntry                                                               /// Individual entries. Used in the Node12 class as a storage unit.
{
    dVector6 PhaseVoigt;
    size_t index;
};

/***************************************************************/
class NodeVoigt                                                                 /// Basic class to store the vector valued quantities for each phase field. Provide access and manipulation methods for the entrys.
{
 public:
    void   set      (const size_t idx, const dVector6 inVoigt);                    /// Set components idx for phase field index n.
    dVector6 get      (const size_t idx) const;                                    /// Return value of component idx for phase field index n.

    void   clear(){VoigtFields.clear();};                                       /// Empties the vector fields.
    size_t  size() const {return VoigtFields.size();};                           /// Returns the size of vector fields.
    typedef std::vector<VoigtEntry>::iterator iterator;                         /// Iterator over the vector fields
    typedef std::vector<VoigtEntry>::const_iterator citerator;                  /// Constant iterator over the vector fields
    iterator begin() {return VoigtFields.begin();};                             /// Iterator to the begin of vector fields
    iterator end()   {return VoigtFields.end();};                               /// Iterator to the end of vector fields
    citerator cbegin() const {return VoigtFields.begin();};                     /// Constant iterator to the begin of vector fields
    citerator cend()   const {return VoigtFields.end();};                       /// Constant iterator to the end of vector fields

 protected:
 private:
    std::vector<VoigtEntry> VoigtFields;                                        /// List of nonvanishing vector fields.
};

inline void NodeVoigt::set(const size_t idx, const dVector6 inVoigt)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == idx)
    {
        i->PhaseVoigt = inVoigt;
        return;
    }

    VoigtEntry NewEntry;
    NewEntry.index = idx;
    NewEntry.PhaseVoigt = inVoigt;

    VoigtFields.push_back(NewEntry);
}

inline dVector6 NodeVoigt::get(const size_t idx) const
{
    dVector6 emptyvStrain;
    emptyvStrain.set_to_zero();
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->index == idx)
    {
        return i->PhaseVoigt;
    }
    return emptyvStrain;
}

} //namespace openphase
#endif
