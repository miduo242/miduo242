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
 *   Main contributors :   Philipp Engels; Efim Borukhovich; Hesham Salama
 *
 */

#ifndef QUATERNION_H
#define QUATERNION_H

#include <cassert>

#include "Base/Includes.h"
#include "Base/dMatrix3x3.h"
#include "Base/EulerAngles.h"

namespace openphase
{

class OP_EXPORTS Quaternion                                                                ///< Stores and manages Quaternions
{
 public:

    Quaternion()                                                                ///< Standard constructor
    {
        s = 1.0;
        v.set_to_zero();
        RotationMatrix.set_to_unity();
    }
    Quaternion(const Quaternion& rhS)
    {
        s = rhS.s;
        v = rhS.v;
        setRotationMatrix();
    }
    Quaternion(std::initializer_list<double> vecinit)
    {
        assert(vecinit.size() == 4 && "Initialization list size is not equal to storage range");

        int ii = 0;
        s = *vecinit.begin();
        for (auto it = vecinit.begin()+1; it != vecinit.end(); it++)
        {
            v[ii] = *it;
            ii += 1;
        }
        setRotationMatrix();
    }

    double& operator[](const size_t index)
    {
        assert(index < 4 && "Access beyond storage range");

        if (index == 0)
        {
            return s;
        };
        return v[index-1];
    }

    double const& operator[](const size_t index) const
    {
        assert(index < 4 && "Access beyond storage range");

        if (index == 0){return s;};
        return v[index-1];
    }
    void set(const double sIn, const double x, const double y, const double z)
    {
        s = sIn;
        v[0] = x;
        v[1] = y;
        v[2] = z;
        setRotationMatrix();
    }

    void set(const double sIn, const dVector3 vIn)
    {
        s = sIn;
        v = vIn;
        setRotationMatrix();
    }

    void set(const dMatrix3x3& RotMatrix)
    {
        assert(RotMatrix.determinant() >= 1.0 - 3.0*DBL_EPSILON && "matrix is not a proper rotation matrix");
        assert(RotMatrix.determinant() <= 1.0 + 3.0*DBL_EPSILON && "matrix is not a proper rotation matrix");

        // Found at http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
        RotationMatrix = RotMatrix;

        double trace = RotMatrix.trace();

        if(1.0 + trace > 0)
        {
            double ss = sqrt(trace + 1.0)*2.0;
            set(0.25*ss,
               (RotMatrix(2,1) - RotMatrix(1,2))/ss,
               (RotMatrix(0,2) - RotMatrix(2,0))/ss,
               (RotMatrix(1,0) - RotMatrix(0,1))/ss);
        }
        else
        {
            if (RotMatrix(0,0) > RotMatrix(1,1) && RotMatrix(0,0) > RotMatrix(2,2))
            {
                double ss = 2.0 * sqrt( 1.0 + RotMatrix(0,0) - RotMatrix(1,1) - RotMatrix(2,2));
                set((RotMatrix(2,1) - RotMatrix(1,2))/ss,
                     0.25 * ss,
                    (RotMatrix(0,1) + RotMatrix(1,0))/ss,
                    (RotMatrix(0,2) + RotMatrix(2,0))/ss);
            }
            else if (RotMatrix(1,1) > RotMatrix(2,2))
            {
                double ss = 2.0 * sqrt(1.0 + RotMatrix(1,1) - RotMatrix(0,0) - RotMatrix(2,2));
                set((RotMatrix(0,2) - RotMatrix(2,0))/ss,
                    (RotMatrix(0,1) + RotMatrix(1,0))/ss,
                     0.25 * ss,
                    (RotMatrix(1,2) + RotMatrix(2,1))/ss);
            }
            else
            {
                double ss = 2.0 * sqrt(1.0 + RotMatrix(2,2) - RotMatrix(0,0) - RotMatrix(1,1));
                set((RotMatrix(1,0) - RotMatrix(0,1))/ss,
                    (RotMatrix(0,2) + RotMatrix(2,0))/ss,
                    (RotMatrix(1,2) + RotMatrix(2,1))/ss,
                     0.25 * ss);
            }
        }
    }

    void set_entry(const int idx, const double val)
    {
        if (idx == 0)
        {
            s = val;
        }
        else
        {
            v[idx] = val;
        }
        setRotationMatrix();
    }
    void set(dVector3 Axis, const double Angle)
    {
        Axis.normalize();
        // https://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/index.htm
        s = cos(Angle/2);
        v[0] = Axis[0] * sin(Angle/2);
        v[1] = Axis[1] * sin(Angle/2);
        v[2] = Axis[2] * sin(Angle/2);
        normalize();
        setRotationMatrix();
    }
    void set_to_zero()
    {
        s    = 1.0;
        v[0] = 0.0;
        v[1] = 0.0;
        v[2] = 0.0;
        RotationMatrix.set_to_unity();
    }
    Quaternion& operator=(const Quaternion rhS)
    {
        s = rhS.s;
        v = rhS.v;
        setRotationMatrix();
        return *this;
    }
    Quaternion  operator+(const Quaternion rhS) const
    {
        Quaternion result;

        result.s = s + rhS.s;
        result.v = v + rhS.v;
        //result.setRotationMatrix();                                           // will be set in assignment operator
        return result;
    }
    Quaternion& operator+=(const Quaternion rhS)
    {
        s = s + rhS.s;
        v = v + rhS.v;
        setRotationMatrix();
        return *this;
    }
    Quaternion  operator-(const Quaternion rhS) const
    {
        Quaternion result;
        result.s = s - rhS.s;
        result.v = v - rhS.v;
        //result.setRotationMatrix();                                           // will be set in assignment operator
        return result;
    }
    Quaternion& operator-=(const Quaternion rhS)
    {
        s = s - rhS.s;
        v = v - rhS.v;
        setRotationMatrix();
        return *this;
    }
    Quaternion  operator*(const double scalar) const
    {
        Quaternion result;
        result.s = s * scalar;
        result.v = v * scalar;
        //result.setRotationMatrix();                                           // will be set in assignment operator
        return result;
    }
    Quaternion& operator*=(const double scalar)
    {
        s *= scalar;
        v *= scalar;
        setRotationMatrix();
        return *this;
    }
    Quaternion  operator*(const Quaternion rhS) const
    {
        Quaternion result;
        result.v[0] =  v[0]*rhS.s    + v[1]*rhS.v[2] - v[2]*rhS.v[1] + s*rhS.v[0];
        result.v[1] = -v[0]*rhS.v[2] + v[1]*rhS.s    + v[2]*rhS.v[0] + s*rhS.v[1];
        result.v[2] =  v[0]*rhS.v[1] - v[1]*rhS.v[0] + v[2]*rhS.s    + s*rhS.v[2];
        result.s    = -v[0]*rhS.v[0] - v[1]*rhS.v[1] - v[2]*rhS.v[2] + s*rhS.s;
        //result.setRotationMatrix();                                           // will be set in assignment operator
        return result;
    }
    Quaternion& operator*=(const Quaternion rhS)
    {
        s = s*rhS.s - v*rhS.v;
        v = v.cross(rhS.v) + rhS.v*s + v*rhS.s;
        setRotationMatrix();
        return *this;
    }
    Quaternion operator/(const double divisor)
    {
        Quaternion result;
        result.s = s/divisor;
        result.v = v/divisor;
        //result.setRotationMatrix();                                           // will be set in assignment operator
        return result;
    }
    Quaternion& operator/=(const double divisor)
    {
        s = s/divisor;
        v = v/divisor;
        setRotationMatrix();
        return *this;
    }
    double length() const
    {
        return sqrt(s*s + v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }
    Quaternion& normalize()
    {
        *this /= sqrt(s*s + v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        return *this;
    }
    Quaternion  normalized() const
    {
        Quaternion result;
        double norm = length();
        result.s = s/norm;
        result.v = v/norm;
        return result;
    }
    Quaternion& conjugate()
    {
        v = v*(-1.0);
        setRotationMatrix();
        return *this;
    }
    Quaternion  conjugated() const
    {
        Quaternion result;
        result.s = s;
        result.v = v*(-1.0);
        return result;
    }
    Quaternion& invert()
    {
        // if |q|=1 then q inverse = q conjugate
        v = v*(-1.0);
        *this /= (s*s + v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        return *this;
    }
    Quaternion  inverted() const
    {
        Quaternion result = *this;
        return result = result.conjugate()/(s*s+v*v);
    }
    void pack(std::vector<double>& buffer)
    {
        buffer.push_back(s);
        v.pack(buffer);
    }
    void unpack(std::vector<double>& buffer, size_t& it)
    {
        s = buffer[it]; ++it;
        v.unpack(buffer, it);
    }
    void setRotationMatrix(void);
    dMatrix3x3 getRotationMatrix(const bool Active = true);
    EulerAngles getEulerAngles(bool Passive = false);

    //EulerAngles getEulerAngles(const EulerConvention EConvention) const;

    static Quaternion lerp(const Quaternion& rhSQ1, const Quaternion& rhSQ2,
                                                            const double t);
    static Quaternion slerp(const Quaternion& a, const Quaternion& b,
                                                            const double t);
    std::string print(void) const;
    std::string print_matrix(void) const;

    dMatrix3x3 RotationMatrix;
 protected:
 private:

    double s;                                                                   /// Real component
    dVector3 v;                                                                 /// Imaginary components
};
}
#endif
