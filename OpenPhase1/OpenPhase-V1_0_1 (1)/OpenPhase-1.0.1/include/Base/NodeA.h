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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich;
 *                         Reza Darvishi Kamachali; Dmitry Medvedev
 *
 */

#ifndef NODEA_H
#define NODEA_H

#include "Base/Includes.h"

namespace openphase
{
/********************************* Declaration *******************************/

struct SingleIndexFieldEntry                                                    ///< Structure for storing the single index entries. Used in the NodeA class as a storage unit.
{
    size_t index;                                                               ///< Field index.
    double value;                                                               ///< Stored value.
    SingleIndexFieldEntry()                                                     ///< Default constructor.
    {
        index = 0;
        value = 0.0;
    }
    SingleIndexFieldEntry(const SingleIndexFieldEntry& entry)                   ///< Copy constructor.
    {
        index = entry.index;
        value = entry.value;
    }
    SingleIndexFieldEntry& operator=(const SingleIndexFieldEntry& entry)        ///< Copy operator.
    {
        index = entry.index;
        value = entry.value;
        return *this;
    }
};

class NodeA                                                                     ///< Stores the single-valued fields at a grid point. Provides access and manipulation methods for the entries.
{
 public:
    NodeA()                                                                     ///< Constructor used to allocate space for at least 3 fields
    {
        Fields.reserve(3);
    }

    NodeA(const NodeA& n)                                                       ///< Constructor used to allocate space for at least 3 fields
    {
        Fields.reserve(3);
        Fields = n.Fields;
    }

    void clear()                                                                ///< Empties the field storage.
    {
        Fields.clear();
    };

    double  operator[](const size_t n) const;                                   ///< Index operator for accessing the n's field value
    NodeA   operator+(const NodeA& n) const;                                    ///< Plus operator. Takes as input another NodeA type entry.
    NodeA   operator-(const NodeA& n) const;                                    ///< Minus operator. Takes as input another NodeA type entry.
    NodeA   operator*(const double n) const;                                    ///< Multiplies all fields by a number.
    NodeA&  operator=(const NodeA& n);                                          ///< Assignment operator.

    NodeA&  operator+=(const NodeA& n);                                         ///< Plus-equal operator. Takes as input another NodeA type entry.
    NodeA&  operator-=(const NodeA& n);                                         ///< Minus-equal operator. Takes as input another NodeA type entry.
    NodeA&  operator*=(const double n);                                         ///< Multiply all fields by a number.

    bool    present(const size_t idx) const;                                    ///< Returns true if the field with a given index is present in the node, false otherwise.

    void    set_value(const size_t n, const double value);                      ///< Sets value.
    void    add_value(const size_t n, const double value);                      ///< Increments value.
    double  get_value(const size_t n) const;                                    ///< Returns value.

    void    add_values         (const NodeA& value);                            ///< Adds values of two nodes.
    void    add_existing_values(const NodeA& value);                            ///< Adds only values existing in two nodes simultaneously.

    SingleIndexFieldEntry get_max(void) const;                                  ///< Returns FieldEntry with max value

    typedef std::vector<SingleIndexFieldEntry>::iterator iterator;              ///< Iterator over storage vector.
    typedef std::vector<SingleIndexFieldEntry>::const_iterator citerator;       ///< Constant iterator over storage vector.
    iterator  begin() {return Fields.begin();};                                 ///< Iterator to the begin of storage vector.
    iterator  end()   {return Fields.end();};                                   ///< Iterator to the end of storage vector.
    citerator cbegin() const {return Fields.cbegin();};                         ///< Constant iterator to the begin of storage vector.
    citerator cend()   const {return Fields.cend();};                           ///< Constant iterator to the end of storage vector.
    size_t    size() const {return Fields.size();};                             ///< Returns the size of storage.
    iterator  erase(iterator it) {return Fields.erase(it);};                    ///< Erase a single record pointed by iterator it.
    SingleIndexFieldEntry& front(void) {return Fields.front();};                ///< Reference to the first FieldEntry.
    const SingleIndexFieldEntry& front(void) const {return Fields.front();};    ///< Constant reference to the first FieldEntry.

    void pack(std::vector<double>& buffer);
    void unpack(std::vector<double>& buffer, size_t& it);
    void Read(std::istream& inp);                                               ///< Reads NodeA content from the input stream.
    void Write(std::ostream& outp) const;                                       ///< Writes NodeA content to the output stream.

 protected:
 private:
    std::vector<SingleIndexFieldEntry> Fields;                                  ///< Storage vector.
};

/******************************* Implementation ******************************/

inline double NodeA::operator[](const size_t n) const
{
    for(auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if(i->index == n) return i->value;
    }
    return 0.0;
}

inline bool NodeA::present(const size_t index) const
{
    for(auto i = Fields.cbegin(); i != Fields.cend(); ++i)
    {
        if(i->index == index) return true;
    }
    return false;
}

inline void NodeA::set_value(const size_t n, const double value)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->value = value;
        return;
    }
    SingleIndexFieldEntry NewEntry;
    NewEntry.index      = n;
    NewEntry.value      = value;
    Fields.push_back(NewEntry);
}

inline void NodeA::add_value(const size_t n, const double value)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->value += value;
        return;
    }
    SingleIndexFieldEntry NewEntry;
    NewEntry.index      = n;
    NewEntry.value      = value;
    Fields.push_back(NewEntry);
}

inline double NodeA::get_value(const size_t n) const
{
    for(auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if(i->index == n) return i->value;
    }
    return 0.0;
}

inline NodeA NodeA::operator+(const NodeA& n) const
{
    NodeA result = n;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add_value(i->index, i->value);
    }
    return result;
}

inline NodeA NodeA::operator-(const NodeA& n) const
{
    NodeA result = n;

    for (auto i = result.begin(); i < result.end(); ++i)
    {
        i->value *= -1.0;
    }

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add_value(i->index, i->value);
    }
    return result;
}

inline NodeA NodeA::operator*(const double n) const
{
    NodeA result = *this;

    for(auto i = result.begin(); i < result.end(); ++i)
    {
        i->value *= n;
    }
    return result;
}

inline NodeA& NodeA::operator+=(const NodeA& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_value(i->index, i->value);
    }
    return *this;
}

inline NodeA& NodeA::operator-=(const NodeA& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_value(i->index, -i->value);
    }
    return *this;
}

inline NodeA& NodeA::operator*=(const double n)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    {
        i->value *= n;
    }
    return *this;
}

inline NodeA& NodeA::operator=(const NodeA& n)
{
    Fields = n.Fields;
    return *this;
}

inline void NodeA::add_existing_values(const NodeA& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        double locVal = get_value(i->index);
        if(locVal != 0.0 and i->value != 0.0)
        {
            set_value(i->index, i->value + locVal);
        }
    }
}

inline SingleIndexFieldEntry NodeA::get_max(void) const
{
    SingleIndexFieldEntry returnFieldEntry;

    returnFieldEntry.value = 0.0;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if (fabs(i->value) > returnFieldEntry.value)
    {
        returnFieldEntry = *i;
    }
    return returnFieldEntry;
}

inline void NodeA::pack(std::vector<double>& buffer)
{
    buffer.push_back(Fields.size());
    for(auto it = Fields.begin(); it != Fields.end(); ++it)
    {
        buffer.push_back(it->index);
        buffer.push_back(it->value);
    }
}

inline void NodeA::unpack(std::vector<double>& buffer, size_t& it)
{
    clear();
    size_t size = buffer[it]; ++it;
    Fields.resize(size);
    for(size_t i = 0; i < size; ++i)
    {
        Fields[i].index = buffer[it]; ++it;
        Fields[i].value = buffer[it]; ++it;
    }
}

inline void NodeA::Read(std::istream& inp)
{
    size_t size = 0;
    inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));
    Fields.resize(size);
    for(auto &Field : Fields)
    {
        inp.read(reinterpret_cast<char*>(&Field.index), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&Field.value), sizeof(double));
    }
}

inline void NodeA::Write(std::ostream& outp) const
{
    size_t size = Fields.size();
    outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

    for(auto &Field : Fields)
    {
        outp.write(reinterpret_cast<const char*>(&Field.index), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&Field.value), sizeof(double));
    }
}

}//namespace openphase
#endif
