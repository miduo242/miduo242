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
 *                         Philipp Engels
 *
 */

#ifndef MATRIXVECTOROPERATORS_H
#define MATRIXVECTOROPERATORS_H

#include "Base/dMatrix3x3.h"
#include "Base/dMatrix6x6.h"
#include "Base/dVector3.h"
#include "Base/vStrain.h"
#include "Base/vStress.h"

namespace openphase
{

extern inline dVector3 operator*(const dMatrix3x3& lhs, const dVector3& rhs)
{
    dVector3 tmp;
    for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
    {
        tmp[i] += lhs(i,j)*rhs[j];
    }
    return tmp;
};

extern inline double operator*(const iVector3& lhs, const dVector3& rhs)
{
    double tmp = 0.0;
    for(int i = 0; i < 3; i++)
    {
        tmp += lhs[i]*rhs[i];
    }
    return tmp;
};

extern inline double operator*(const dVector3& lhs, const iVector3& rhs)
{
    double tmp = 0.0;
    for(int i = 0; i < 3; i++)
    {
        tmp += lhs[i]*rhs[i];
    }
    return tmp;
};

extern inline dVector3 operator+(const iVector3& lhs, const dVector3& rhs)
{
    dVector3 tmp;
    for(int i = 0; i < 3; i++)
    {
        tmp[i] = lhs[i]+rhs[i];
    }
    return tmp;
};

extern inline dVector3 operator+(const dVector3& lhs, const iVector3& rhs)
{
    dVector3 tmp;
    for(int i = 0; i < 3; i++)
    {
        tmp[i] = lhs[i]+rhs[i];
    }
    return tmp;
};

extern inline dVector3 operator-(const iVector3& lhs, const dVector3& rhs)
{
    dVector3 tmp;
    for(int i = 0; i < 3; i++)
    {
        tmp[i] = lhs[i]-rhs[i];
    }
    return tmp;
};

extern inline dVector3 operator-(const dVector3& lhs, const iVector3& rhs)
{
    dVector3 tmp;
    for(int i = 0; i < 3; i++)
    {
        tmp[i] = lhs[i]-rhs[i];
    }
    return tmp;
};

extern inline dMatrix6x6 outer(const dMatrix3x3& lhs, const dMatrix3x3& rhs)
{
    dMatrix6x6 tmp;

    tmp(0,0) += lhs(0,0)*rhs(0,0);
    tmp(0,5) += lhs(0,0)*rhs(0,1);
    tmp(0,4) += lhs(0,0)*rhs(0,2);
    tmp(0,5) += lhs(0,0)*rhs(1,0);
    tmp(0,1) += lhs(0,0)*rhs(1,1);
    tmp(0,3) += lhs(0,0)*rhs(1,2);
    tmp(0,4) += lhs(0,0)*rhs(2,0);
    tmp(0,3) += lhs(0,0)*rhs(2,1);
    tmp(0,2) += lhs(0,0)*rhs(2,2);
    tmp(5,0) += lhs(0,1)*rhs(0,0);
    tmp(5,5) += lhs(0,1)*rhs(0,1);
    tmp(5,4) += lhs(0,1)*rhs(0,2);
    tmp(5,5) += lhs(0,1)*rhs(1,0);
    tmp(5,1) += lhs(0,1)*rhs(1,1);
    tmp(5,3) += lhs(0,1)*rhs(1,2);
    tmp(5,4) += lhs(0,1)*rhs(2,0);
    tmp(5,3) += lhs(0,1)*rhs(2,1);
    tmp(5,2) += lhs(0,1)*rhs(2,2);
    tmp(4,0) += lhs(0,2)*rhs(0,0);
    tmp(4,5) += lhs(0,2)*rhs(0,1);
    tmp(4,4) += lhs(0,2)*rhs(0,2);
    tmp(4,5) += lhs(0,2)*rhs(1,0);
    tmp(4,1) += lhs(0,2)*rhs(1,1);
    tmp(4,3) += lhs(0,2)*rhs(1,2);
    tmp(4,4) += lhs(0,2)*rhs(2,0);
    tmp(4,3) += lhs(0,2)*rhs(2,1);
    tmp(4,2) += lhs(0,2)*rhs(2,2);
    tmp(5,0) += lhs(1,0)*rhs(0,0);
    tmp(5,5) += lhs(1,0)*rhs(0,1);
    tmp(5,4) += lhs(1,0)*rhs(0,2);
    tmp(5,5) += lhs(1,0)*rhs(1,0);
    tmp(5,1) += lhs(1,0)*rhs(1,1);
    tmp(5,3) += lhs(1,0)*rhs(1,2);
    tmp(5,4) += lhs(1,0)*rhs(2,0);
    tmp(5,3) += lhs(1,0)*rhs(2,1);
    tmp(5,2) += lhs(1,0)*rhs(2,2);
    tmp(1,0) += lhs(1,1)*rhs(0,0);
    tmp(1,5) += lhs(1,1)*rhs(0,1);
    tmp(1,4) += lhs(1,1)*rhs(0,2);
    tmp(1,5) += lhs(1,1)*rhs(1,0);
    tmp(1,1) += lhs(1,1)*rhs(1,1);
    tmp(1,3) += lhs(1,1)*rhs(1,2);
    tmp(1,4) += lhs(1,1)*rhs(2,0);
    tmp(1,3) += lhs(1,1)*rhs(2,1);
    tmp(1,2) += lhs(1,1)*rhs(2,2);
    tmp(3,0) += lhs(1,2)*rhs(0,0);
    tmp(3,5) += lhs(1,2)*rhs(0,1);
    tmp(3,4) += lhs(1,2)*rhs(0,2);
    tmp(3,5) += lhs(1,2)*rhs(1,0);
    tmp(3,1) += lhs(1,2)*rhs(1,1);
    tmp(3,3) += lhs(1,2)*rhs(1,2);
    tmp(3,4) += lhs(1,2)*rhs(2,0);
    tmp(3,3) += lhs(1,2)*rhs(2,1);
    tmp(3,2) += lhs(1,2)*rhs(2,2);
    tmp(4,0) += lhs(2,0)*rhs(0,0);
    tmp(4,5) += lhs(2,0)*rhs(0,1);
    tmp(4,4) += lhs(2,0)*rhs(0,2);
    tmp(4,5) += lhs(2,0)*rhs(1,0);
    tmp(4,1) += lhs(2,0)*rhs(1,1);
    tmp(4,3) += lhs(2,0)*rhs(1,2);
    tmp(4,4) += lhs(2,0)*rhs(2,0);
    tmp(4,3) += lhs(2,0)*rhs(2,1);
    tmp(4,2) += lhs(2,0)*rhs(2,2);
    tmp(3,0) += lhs(2,1)*rhs(0,0);
    tmp(3,5) += lhs(2,1)*rhs(0,1);
    tmp(3,4) += lhs(2,1)*rhs(0,2);
    tmp(3,5) += lhs(2,1)*rhs(1,0);
    tmp(3,1) += lhs(2,1)*rhs(1,1);
    tmp(3,3) += lhs(2,1)*rhs(1,2);
    tmp(3,4) += lhs(2,1)*rhs(2,0);
    tmp(3,3) += lhs(2,1)*rhs(2,1);
    tmp(3,2) += lhs(2,1)*rhs(2,2);
    tmp(2,0) += lhs(2,2)*rhs(0,0);
    tmp(2,5) += lhs(2,2)*rhs(0,1);
    tmp(2,4) += lhs(2,2)*rhs(0,2);
    tmp(2,5) += lhs(2,2)*rhs(1,0);
    tmp(2,1) += lhs(2,2)*rhs(1,1);
    tmp(2,3) += lhs(2,2)*rhs(1,2);
    tmp(2,4) += lhs(2,2)*rhs(2,0);
    tmp(2,3) += lhs(2,2)*rhs(2,1),
    tmp(2,2) += lhs(2,2)*rhs(2,2);

    return tmp;
}

extern inline vStrain VoigtStrain(const dMatrix3x3& locStrainTensor)
{
    vStrain locStrain;
    locStrain[0] = locStrainTensor(0,0);
    locStrain[1] = locStrainTensor(1,1);
    locStrain[2] = locStrainTensor(2,2);
    locStrain[3] = locStrainTensor(1,2)*2.0;
    locStrain[4] = locStrainTensor(0,2)*2.0;
    locStrain[5] = locStrainTensor(0,1)*2.0;
    return locStrain;
}

extern inline vStress VoigtStress(const dMatrix3x3& locStressTensor)
{
    vStress locStress;
    locStress[0] = locStressTensor(0,0);
    locStress[1] = locStressTensor(1,1);
    locStress[2] = locStressTensor(2,2);
    locStress[3] = locStressTensor(1,2);
    locStress[4] = locStressTensor(0,2);
    locStress[5] = locStressTensor(0,1);
    return locStress;
}

// Return columns of a matrix
extern inline std::vector<dVector3> Col(const dMatrix3x3& Mat)
{
    std::vector<dVector3> Col;
    for(int i = 0; i < 3; i++)
    {
        dVector3 tmp;
        for(int j = 0; j < 3; j++)
        {
            tmp[j] = Mat(j,i);
        }
        Col.push_back(tmp);
    }
    return Col;
}

extern inline vStrain operator*(const dMatrix6x6& locCompliance, const vStress& locStress)
{
    vStrain locStrain;
    for(int i = 0; i < 6; i++)
    for(int j = 0; j < 6; j++)
    {
        locStrain[i] += locCompliance(i,j)*locStress[j];
    }
    return locStrain;
}

extern inline vStress operator*(const dMatrix6x6& locStiffness, const vStrain& locStrain)
{
    vStress locStress;
    for(int i = 0; i < 6; i++)
    for(int j = 0; j < 6; j++)
    {
        locStress[i] += locStiffness(i,j)*locStrain[j];
    }
    return locStress;
};

extern inline double operator*(const vStrain& locStrain, const vStress& locStress)
{
    double locEnergy = locStrain[0] * locStress[0]
                     + locStrain[1] * locStress[1]
                     + locStrain[2] * locStress[2]
                     + locStrain[3] * locStress[3]
                     + locStrain[4] * locStress[4]
                     + locStrain[5] * locStress[5];
    return locEnergy;
};

extern inline double operator*(const vStress& locStress, const vStrain& locStrain)
{
    double locEnergy = locStrain[0] * locStress[0]
                     + locStrain[1] * locStress[1]
                     + locStrain[2] * locStress[2]
                     + locStrain[3] * locStress[3]
                     + locStrain[4] * locStress[4]
                     + locStrain[5] * locStress[5];
    return locEnergy;
};

}// namespace openphase
#endif
