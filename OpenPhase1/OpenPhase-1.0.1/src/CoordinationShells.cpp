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
 *   File created :   2015
 *   Main contributors :   Oleg Shchyglo;
 *
 */

#include "CoordinationShells.h"
#include "Settings.h"

using namespace std;

namespace openphase
{

void CoordinationShells::Initialize(const Settings& locSettings)
{
    int MAXsize = std::max(std::max(locSettings.TotalNx, locSettings.Ny), locSettings.Nz);
    for(int x = 0; x < MAXsize; x++)
    for(int y = x; y < MAXsize; y++)
    for(int z = y; z < MAXsize; z++)
    {
        CShells.push_back(iVector3({x,y,z}));
    }

    std::sort(CShells.begin(),CShells.end());
}
vector<iVector3> CoordinationShells::operator[](size_t n)
{
    vector<iVector3> myReturn;

    int x_val = CShells[n][0];
    int y_val = CShells[n][1];
    int z_val = CShells[n][2];

    if(x_val == 0)
    {
        if(y_val == 0) // x == 0; y == 0;
        {
            if(z_val == 0) // x == 0; y == 0; z == 0; (1)
            {
                myReturn.push_back({ 0, 0, 0});
            }
            else // x == 0; y == 0; z != 0; (6)
            {
                myReturn.push_back({     0,     0, z_val});
                myReturn.push_back({     0, z_val,     0});
                myReturn.push_back({ z_val,     0,     0});

                myReturn.push_back({     0,     0,-z_val});
                myReturn.push_back({     0,-z_val,     0});
                myReturn.push_back({-z_val,     0,     0});
            }
        }
        else if(y_val == z_val) // x == 0; y == z; (12)
        {
            myReturn.push_back({     0, y_val, y_val});
            myReturn.push_back({ y_val,     0, y_val});
            myReturn.push_back({ y_val, y_val,     0});

            myReturn.push_back({     0, y_val,-y_val});
            myReturn.push_back({ y_val,     0,-y_val});
            myReturn.push_back({ y_val,-y_val,     0});

            myReturn.push_back({     0,-y_val, y_val});
            myReturn.push_back({-y_val,     0, y_val});
            myReturn.push_back({-y_val, y_val,     0});

            myReturn.push_back({     0,-y_val,-y_val});
            myReturn.push_back({-y_val,     0,-y_val});
            myReturn.push_back({-y_val,-y_val,     0});
        }
        else // x == 0; y != z; (24)
        {
            myReturn.push_back({     0, y_val, z_val});
            myReturn.push_back({ y_val,     0, z_val});
            myReturn.push_back({ y_val, z_val,     0});

            myReturn.push_back({     0, y_val,-z_val});
            myReturn.push_back({ y_val,     0,-z_val});
            myReturn.push_back({ y_val,-z_val,     0});

            myReturn.push_back({     0,-y_val, z_val});
            myReturn.push_back({-y_val,     0, z_val});
            myReturn.push_back({-y_val, z_val,     0});

            myReturn.push_back({     0,-y_val,-z_val});
            myReturn.push_back({-y_val,     0,-z_val});
            myReturn.push_back({-y_val,-z_val,     0});

            myReturn.push_back({     0, z_val, y_val});
            myReturn.push_back({ z_val,     0, y_val});
            myReturn.push_back({ z_val, y_val,     0});

            myReturn.push_back({     0, z_val,-y_val});
            myReturn.push_back({ z_val,     0,-y_val});
            myReturn.push_back({ z_val,-y_val,     0});

            myReturn.push_back({     0,-z_val, y_val});
            myReturn.push_back({-z_val,     0, y_val});
            myReturn.push_back({-z_val, y_val,     0});

            myReturn.push_back({     0,-z_val,-y_val});
            myReturn.push_back({-z_val,     0,-y_val});
            myReturn.push_back({-z_val,-y_val,     0});
        }
    }
    else if (x_val == y_val and y_val == z_val) // x == y == z; (8)
    {
        myReturn.push_back({ x_val, x_val, x_val});

        myReturn.push_back({ x_val, x_val,-x_val});
        myReturn.push_back({ x_val,-x_val, x_val});
        myReturn.push_back({-x_val, x_val, x_val});

        myReturn.push_back({ x_val,-x_val,-x_val});
        myReturn.push_back({-x_val, x_val,-x_val});
        myReturn.push_back({-x_val,-x_val, x_val});

        myReturn.push_back({-x_val,-x_val,-x_val});
    }
    else if (x_val == y_val) // x == y != z; (24)
    {
        myReturn.push_back({ x_val, x_val, z_val});
        myReturn.push_back({ x_val, z_val, x_val});
        myReturn.push_back({ z_val, x_val, x_val});

        myReturn.push_back({ x_val, x_val,-z_val});
        myReturn.push_back({ x_val,-z_val, x_val});
        myReturn.push_back({-z_val, x_val, x_val});

        myReturn.push_back({-x_val, x_val, z_val});
        myReturn.push_back({-x_val, z_val, x_val});
        myReturn.push_back({ z_val,-x_val, x_val});

        myReturn.push_back({-x_val, x_val,-z_val});
        myReturn.push_back({-x_val,-z_val, x_val});
        myReturn.push_back({-z_val,-x_val, x_val});

        myReturn.push_back({ x_val,-x_val, z_val});
        myReturn.push_back({ x_val, z_val,-x_val});
        myReturn.push_back({ z_val, x_val,-x_val});

        myReturn.push_back({ x_val,-x_val,-z_val});
        myReturn.push_back({ x_val,-z_val,-x_val});
        myReturn.push_back({-z_val, x_val,-x_val});

        myReturn.push_back({-x_val,-x_val, z_val});
        myReturn.push_back({-x_val, z_val,-x_val});
        myReturn.push_back({ z_val,-x_val,-x_val});

        myReturn.push_back({-x_val,-x_val,-z_val});
        myReturn.push_back({-x_val,-z_val,-x_val});
        myReturn.push_back({-z_val,-x_val,-x_val});
    }
    else if (y_val == z_val) // x != y == z; (24)
    {
        myReturn.push_back({ x_val, y_val, y_val});
        myReturn.push_back({ y_val, x_val, y_val});
        myReturn.push_back({ y_val, y_val, x_val});

        myReturn.push_back({ x_val, y_val,-y_val});
        myReturn.push_back({ y_val, x_val,-y_val});
        myReturn.push_back({ y_val,-y_val, x_val});

        myReturn.push_back({ x_val,-y_val, y_val});
        myReturn.push_back({-y_val, x_val, y_val});
        myReturn.push_back({-y_val, y_val, x_val});

        myReturn.push_back({ x_val,-y_val,-y_val});
        myReturn.push_back({-y_val, x_val,-y_val});
        myReturn.push_back({-y_val,-y_val, x_val});

        myReturn.push_back({-x_val, y_val, y_val});
        myReturn.push_back({ y_val,-x_val, y_val});
        myReturn.push_back({ y_val, y_val,-x_val});

        myReturn.push_back({-x_val, y_val,-y_val});
        myReturn.push_back({ y_val,-x_val,-y_val});
        myReturn.push_back({ y_val,-y_val,-x_val});

        myReturn.push_back({-x_val,-y_val, y_val});
        myReturn.push_back({-y_val,-x_val, y_val});
        myReturn.push_back({-y_val, y_val,-x_val});

        myReturn.push_back({-x_val,-y_val,-y_val});
        myReturn.push_back({-y_val,-x_val,-y_val});
        myReturn.push_back({-y_val,-y_val,-x_val});
    }
    else if (x_val == z_val) // x == z; x != y; (24)
    {
        myReturn.push_back({ x_val, y_val, x_val});
        myReturn.push_back({ y_val, x_val, x_val});
        myReturn.push_back({ x_val, x_val, y_val});

        myReturn.push_back({ x_val, y_val,-x_val});
        myReturn.push_back({ y_val,-x_val, x_val});
        myReturn.push_back({-x_val, x_val, y_val});

        myReturn.push_back({-x_val, y_val, x_val});
        myReturn.push_back({-y_val, x_val, x_val});
        myReturn.push_back({ x_val,-x_val, y_val});

        myReturn.push_back({-x_val, y_val,-x_val});
        myReturn.push_back({-y_val,-x_val, x_val});
        myReturn.push_back({-x_val,-x_val, y_val});

        myReturn.push_back({ x_val,-y_val, x_val});
        myReturn.push_back({ y_val, x_val,-x_val});
        myReturn.push_back({ x_val, x_val,-y_val});

        myReturn.push_back({ x_val,-y_val,-x_val});
        myReturn.push_back({ y_val,-x_val,-x_val});
        myReturn.push_back({-x_val, x_val,-y_val});

        myReturn.push_back({-x_val,-y_val, x_val});
        myReturn.push_back({-y_val, x_val,-x_val});
        myReturn.push_back({ x_val,-x_val,-y_val});

        myReturn.push_back({-x_val,-y_val,-x_val});
        myReturn.push_back({-y_val,-x_val,-x_val});
        myReturn.push_back({-x_val,-x_val,-y_val});
    }
    else // x != y != z; (48)
    {
        myReturn.push_back({ x_val, y_val, z_val});
        myReturn.push_back({ y_val, x_val, z_val});
        myReturn.push_back({ y_val, z_val, x_val});

        myReturn.push_back({ x_val, y_val,-z_val});
        myReturn.push_back({ y_val, x_val,-z_val});
        myReturn.push_back({ y_val,-z_val, x_val});

        myReturn.push_back({ x_val,-y_val, z_val});
        myReturn.push_back({-y_val, x_val, z_val});
        myReturn.push_back({-y_val, z_val, x_val});

        myReturn.push_back({ x_val,-y_val,-z_val});
        myReturn.push_back({-y_val, x_val,-z_val});
        myReturn.push_back({-y_val,-z_val, x_val});

        myReturn.push_back({ x_val, z_val, y_val});
        myReturn.push_back({ z_val, x_val, y_val});
        myReturn.push_back({ z_val, y_val, x_val});

        myReturn.push_back({ x_val, z_val,-y_val});
        myReturn.push_back({ z_val, x_val,-y_val});
        myReturn.push_back({ z_val,-y_val, x_val});

        myReturn.push_back({ x_val,-z_val, y_val});
        myReturn.push_back({-z_val, x_val, y_val});
        myReturn.push_back({-z_val, y_val, x_val});

        myReturn.push_back({ x_val,-z_val,-y_val});
        myReturn.push_back({-z_val, x_val,-y_val});
        myReturn.push_back({-z_val,-y_val, x_val});

        myReturn.push_back({-x_val, y_val, z_val});
        myReturn.push_back({ y_val,-x_val, z_val});
        myReturn.push_back({ y_val, z_val,-x_val});

        myReturn.push_back({-x_val, y_val,-z_val});
        myReturn.push_back({ y_val,-x_val,-z_val});
        myReturn.push_back({ y_val,-z_val,-x_val});

        myReturn.push_back({-x_val,-y_val, z_val});
        myReturn.push_back({-y_val,-x_val, z_val});
        myReturn.push_back({-y_val, z_val,-x_val});

        myReturn.push_back({-x_val,-y_val,-z_val});
        myReturn.push_back({-y_val,-x_val,-z_val});
        myReturn.push_back({-y_val,-z_val,-x_val});

        myReturn.push_back({-x_val, z_val, y_val});
        myReturn.push_back({ z_val,-x_val, y_val});
        myReturn.push_back({ z_val, y_val,-x_val});

        myReturn.push_back({-x_val, z_val,-y_val});
        myReturn.push_back({ z_val,-x_val,-y_val});
        myReturn.push_back({ z_val,-y_val,-x_val});

        myReturn.push_back({-x_val,-z_val, y_val});
        myReturn.push_back({-z_val,-x_val, y_val});
        myReturn.push_back({-z_val, y_val,-x_val});

        myReturn.push_back({-x_val,-z_val,-y_val});
        myReturn.push_back({-z_val,-x_val,-y_val});
        myReturn.push_back({-z_val,-y_val,-x_val});
    }
    return myReturn;
}

void sortIndex(double& a, double& b, double& c)
{
    if (a > b)
    {
        std::swap(a, b);
    }
    if (b > c)
    {
        std::swap(b, c);
    }
    if (a > b)
    {
        std::swap(a, b);
    }
}

std::vector<dVector3> CoordinationShells::findPermutations(std::vector<dVector3>& locFacets, size_t n)
{
    std::vector<dVector3> myReturn;

    double x_val = locFacets[n].getX();
    double y_val = locFacets[n].getY();
    double z_val = locFacets[n].getZ();

    sortIndex(x_val,y_val,z_val);

    if(x_val == 0)
    {
        if(y_val == 0) // x == 0; y == 0;
        {
            if(z_val == 0) // x == 0; y == 0; z == 0; (1)
            {
                myReturn.push_back({ 0, 0, 0});
            }
            else // x == 0; y == 0; z != 0; (6)
            {
                myReturn.push_back({     0,     0, z_val});
                myReturn.push_back({     0, z_val,     0});
                myReturn.push_back({ z_val,     0,     0});

                myReturn.push_back({     0,     0,-z_val});
                myReturn.push_back({     0,-z_val,     0});
                myReturn.push_back({-z_val,     0,     0});
            }
        }
        else if(y_val == z_val) // x == 0; y == z; (12)
        {
            myReturn.push_back({     0, y_val, y_val});
            myReturn.push_back({ y_val,     0, y_val});
            myReturn.push_back({ y_val, y_val,     0});

            myReturn.push_back({     0, y_val,-y_val});
            myReturn.push_back({ y_val,     0,-y_val});
            myReturn.push_back({ y_val,-y_val,     0});

            myReturn.push_back({     0,-y_val, y_val});
            myReturn.push_back({-y_val,     0, y_val});
            myReturn.push_back({-y_val, y_val,     0});

            myReturn.push_back({     0,-y_val,-y_val});
            myReturn.push_back({-y_val,     0,-y_val});
            myReturn.push_back({-y_val,-y_val,     0});
        }
        else // x == 0; y != z; (24)
        {
            myReturn.push_back({     0, y_val, z_val});
            myReturn.push_back({ y_val,     0, z_val});
            myReturn.push_back({ y_val, z_val,     0});

            myReturn.push_back({     0, y_val,-z_val});
            myReturn.push_back({ y_val,     0,-z_val});
            myReturn.push_back({ y_val,-z_val,     0});

            myReturn.push_back({     0,-y_val, z_val});
            myReturn.push_back({-y_val,     0, z_val});
            myReturn.push_back({-y_val, z_val,     0});

            myReturn.push_back({     0,-y_val,-z_val});
            myReturn.push_back({-y_val,     0,-z_val});
            myReturn.push_back({-y_val,-z_val,     0});

            myReturn.push_back({     0, z_val, y_val});
            myReturn.push_back({ z_val,     0, y_val});
            myReturn.push_back({ z_val, y_val,     0});

            myReturn.push_back({     0, z_val,-y_val});
            myReturn.push_back({ z_val,     0,-y_val});
            myReturn.push_back({ z_val,-y_val,     0});

            myReturn.push_back({     0,-z_val, y_val});
            myReturn.push_back({-z_val,     0, y_val});
            myReturn.push_back({-z_val, y_val,     0});

            myReturn.push_back({     0,-z_val,-y_val});
            myReturn.push_back({-z_val,     0,-y_val});
            myReturn.push_back({-z_val,-y_val,     0});
        }
    }
    else if (x_val == y_val and y_val == z_val) // x == y == z; (8)
    {
        myReturn.push_back({ x_val, x_val, x_val});

        myReturn.push_back({ x_val, x_val,-x_val});
        myReturn.push_back({ x_val,-x_val, x_val});
        myReturn.push_back({-x_val, x_val, x_val});

        myReturn.push_back({ x_val,-x_val,-x_val});
        myReturn.push_back({-x_val, x_val,-x_val});
        myReturn.push_back({-x_val,-x_val, x_val});

        myReturn.push_back({-x_val,-x_val,-x_val});
    }
    else if (x_val == y_val) // x == y != z; (24)
    {
        myReturn.push_back({ x_val, x_val, z_val});
        myReturn.push_back({ x_val, z_val, x_val});
        myReturn.push_back({ z_val, x_val, x_val});

        myReturn.push_back({ x_val, x_val,-z_val});
        myReturn.push_back({ x_val,-z_val, x_val});
        myReturn.push_back({-z_val, x_val, x_val});

        myReturn.push_back({-x_val, x_val, z_val});
        myReturn.push_back({-x_val, z_val, x_val});
        myReturn.push_back({ z_val,-x_val, x_val});

        myReturn.push_back({-x_val, x_val,-z_val});
        myReturn.push_back({-x_val,-z_val, x_val});
        myReturn.push_back({-z_val,-x_val, x_val});

        myReturn.push_back({ x_val,-x_val, z_val});
        myReturn.push_back({ x_val, z_val,-x_val});
        myReturn.push_back({ z_val, x_val,-x_val});

        myReturn.push_back({ x_val,-x_val,-z_val});
        myReturn.push_back({ x_val,-z_val,-x_val});
        myReturn.push_back({-z_val, x_val,-x_val});

        myReturn.push_back({-x_val,-x_val, z_val});
        myReturn.push_back({-x_val, z_val,-x_val});
        myReturn.push_back({ z_val,-x_val,-x_val});

        myReturn.push_back({-x_val,-x_val,-z_val});
        myReturn.push_back({-x_val,-z_val,-x_val});
        myReturn.push_back({-z_val,-x_val,-x_val});
    }
    else if (y_val == z_val) // x != y == z; (24)
    {
        myReturn.push_back({ x_val, y_val, y_val});
        myReturn.push_back({ y_val, x_val, y_val});
        myReturn.push_back({ y_val, y_val, x_val});

        myReturn.push_back({ x_val, y_val,-y_val});
        myReturn.push_back({ y_val, x_val,-y_val});
        myReturn.push_back({ y_val,-y_val, x_val});

        myReturn.push_back({ x_val,-y_val, y_val});
        myReturn.push_back({-y_val, x_val, y_val});
        myReturn.push_back({-y_val, y_val, x_val});

        myReturn.push_back({ x_val,-y_val,-y_val});
        myReturn.push_back({-y_val, x_val,-y_val});
        myReturn.push_back({-y_val,-y_val, x_val});

        myReturn.push_back({-x_val, y_val, y_val});
        myReturn.push_back({ y_val,-x_val, y_val});
        myReturn.push_back({ y_val, y_val,-x_val});

        myReturn.push_back({-x_val, y_val,-y_val});
        myReturn.push_back({ y_val,-x_val,-y_val});
        myReturn.push_back({ y_val,-y_val,-x_val});

        myReturn.push_back({-x_val,-y_val, y_val});
        myReturn.push_back({-y_val,-x_val, y_val});
        myReturn.push_back({-y_val, y_val,-x_val});

        myReturn.push_back({-x_val,-y_val,-y_val});
        myReturn.push_back({-y_val,-x_val,-y_val});
        myReturn.push_back({-y_val,-y_val,-x_val});
    }
    else if (x_val == z_val) // x == z; x != y; (24)
    {
        myReturn.push_back({ x_val, y_val, x_val});
        myReturn.push_back({ y_val, x_val, x_val});
        myReturn.push_back({ x_val, x_val, y_val});

        myReturn.push_back({ x_val, y_val,-x_val});
        myReturn.push_back({ y_val,-x_val, x_val});
        myReturn.push_back({-x_val, x_val, y_val});

        myReturn.push_back({-x_val, y_val, x_val});
        myReturn.push_back({-y_val, x_val, x_val});
        myReturn.push_back({ x_val,-x_val, y_val});

        myReturn.push_back({-x_val, y_val,-x_val});
        myReturn.push_back({-y_val,-x_val, x_val});
        myReturn.push_back({-x_val,-x_val, y_val});

        myReturn.push_back({ x_val,-y_val, x_val});
        myReturn.push_back({ y_val, x_val,-x_val});
        myReturn.push_back({ x_val, x_val,-y_val});

        myReturn.push_back({ x_val,-y_val,-x_val});
        myReturn.push_back({ y_val,-x_val,-x_val});
        myReturn.push_back({-x_val, x_val,-y_val});

        myReturn.push_back({-x_val,-y_val, x_val});
        myReturn.push_back({-y_val, x_val,-x_val});
        myReturn.push_back({ x_val,-x_val,-y_val});

        myReturn.push_back({-x_val,-y_val,-x_val});
        myReturn.push_back({-y_val,-x_val,-x_val});
        myReturn.push_back({-x_val,-x_val,-y_val});
    }
    else // x != y != z; (48)
    {
        myReturn.push_back({ x_val, y_val, z_val});
        myReturn.push_back({ y_val, x_val, z_val});
        myReturn.push_back({ y_val, z_val, x_val});

        myReturn.push_back({ x_val, y_val,-z_val});
        myReturn.push_back({ y_val, x_val,-z_val});
        myReturn.push_back({ y_val,-z_val, x_val});

        myReturn.push_back({ x_val,-y_val, z_val});
        myReturn.push_back({-y_val, x_val, z_val});
        myReturn.push_back({-y_val, z_val, x_val});

        myReturn.push_back({ x_val,-y_val,-z_val});
        myReturn.push_back({-y_val, x_val,-z_val});
        myReturn.push_back({-y_val,-z_val, x_val});

        myReturn.push_back({ x_val, z_val, y_val});
        myReturn.push_back({ z_val, x_val, y_val});
        myReturn.push_back({ z_val, y_val, x_val});

        myReturn.push_back({ x_val, z_val,-y_val});
        myReturn.push_back({ z_val, x_val,-y_val});
        myReturn.push_back({ z_val,-y_val, x_val});

        myReturn.push_back({ x_val,-z_val, y_val});
        myReturn.push_back({-z_val, x_val, y_val});
        myReturn.push_back({-z_val, y_val, x_val});

        myReturn.push_back({ x_val,-z_val,-y_val});
        myReturn.push_back({-z_val, x_val,-y_val});
        myReturn.push_back({-z_val,-y_val, x_val});

        myReturn.push_back({-x_val, y_val, z_val});
        myReturn.push_back({ y_val,-x_val, z_val});
        myReturn.push_back({ y_val, z_val,-x_val});

        myReturn.push_back({-x_val, y_val,-z_val});
        myReturn.push_back({ y_val,-x_val,-z_val});
        myReturn.push_back({ y_val,-z_val,-x_val});

        myReturn.push_back({-x_val,-y_val, z_val});
        myReturn.push_back({-y_val,-x_val, z_val});
        myReturn.push_back({-y_val, z_val,-x_val});

        myReturn.push_back({-x_val,-y_val,-z_val});
        myReturn.push_back({-y_val,-x_val,-z_val});
        myReturn.push_back({-y_val,-z_val,-x_val});

        myReturn.push_back({-x_val, z_val, y_val});
        myReturn.push_back({ z_val,-x_val, y_val});
        myReturn.push_back({ z_val, y_val,-x_val});

        myReturn.push_back({-x_val, z_val,-y_val});
        myReturn.push_back({ z_val,-x_val,-y_val});
        myReturn.push_back({ z_val,-y_val,-x_val});

        myReturn.push_back({-x_val,-z_val, y_val});
        myReturn.push_back({-z_val,-x_val, y_val});
        myReturn.push_back({-z_val, y_val,-x_val});

        myReturn.push_back({-x_val,-z_val,-y_val});
        myReturn.push_back({-z_val,-x_val,-y_val});
        myReturn.push_back({-z_val,-y_val,-x_val});
    }
    return myReturn;
}


} // namespace openphase
