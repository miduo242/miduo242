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
 *   File created :   2021
 *   Main contributors :   Oleg Shchyglo; Raphael Schiedung; Marvin Tegeler;
 *                         Matthias Stratmann
 *
 */

#ifndef TEMPERATURE1DEXTENSION_H
#define TEMPERATURE1DEXTENSION_H

#include "Base/Includes.h"
#include "Info.h"

namespace openphase
{

class Temperature;
class HeatDiffusion;

class Temperature1Dextension                                                    ///< 1D temperature field extension storage
{
 public:
    void Initialize(size_t size, iVector3 direction)
    {
        if(size <= 0) // Checks for sufficient extension size
        {
            std::stringstream message;
            message << "1D temperature extension size = " << size << " is incorrect!\n"
                    << "It should be in the range [1, inf) for correct operation\n";
            Info::WriteExit(message.str(), "Temperature1Dextension", "Initialize()");
            exit(1);
        }
        Data.resize(size+2,0.0);
        DataOld.resize(size+2,0.0);
        Direction = direction;
        Qdot = 0.0;
    }
    Temperature1Dextension& operator=(const Temperature1Dextension& RHS)
    {
        Qdot    = RHS.Qdot;
        Data    = RHS.Data;
        DataOld = RHS.DataOld;
        Direction = RHS.Direction;
        Xbounds = RHS.Xbounds;
        Ybounds = RHS.Ybounds;
        Zbounds = RHS.Zbounds;

        return *this;
    }
    void setBC()                                                                ///< Sets adiabatic boundary condition at the far end of the extension
    {
        Data[size()-1] = Data[size()-2];
    }
    void setBC(double value)                                                    ///< Sets temperature at the far end of the extension to the specified value
    {
        Data[size()-1] = value;
    }
    bool isActive()                                                             ///< Returns true if extension was activated, false otherwise
    {
        return Data.size() > 0;
    };
    size_t size() const
    {
        return Data.size();
    }
    void store_temporary(void)
    {
        DataOld = Data;
    }
    void read(std::ifstream& out)
    {
        out.read(reinterpret_cast<char*>(Data.data()),Data.size()*sizeof(double));
    }
    void write(std::ofstream& out)
    {
        out.write(reinterpret_cast<char*>(Data.data()),Data.size()*sizeof(double));
    }
    void SetInitial(const Temperature& Tx);
    void PerformImplicitIteration(Temperature& Tx, HeatDiffusion& HD, double& residual, double dt);

    double Qdot;                                                                ///< Heat source at the far end of the extension

    std::vector<double> Data;                                                   ///< Data storage array
    std::vector<double> DataOld;                                                ///< Temporary data storage for iterative heat diffusion solver
 protected:
 private:

    iVector3 Direction;                                                         ///< Selects the extension's direction: 0 -> direction inactive, 1 -> upper boundary extension, -1 -> lower boundary extension

    std::pair<int, int> Xbounds;                                                ///< X axis loop bounds
    std::pair<int, int> Ybounds;                                                ///< Y axis loop bounds
    std::pair<int, int> Zbounds;                                                ///< Z axis loop bounds
};
} // namespace openphase
#endif
