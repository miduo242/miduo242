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
 *   Main contributors :   Oleg Shchyglo
 *
 */

#include "UserDrivingForce.h"
#include "Info.h"
#include "Settings.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "Temperature.h"

namespace openphase
{
using namespace std;

UserDrivingForce::UserDrivingForce(Settings& locSettings,
                                               const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void UserDrivingForce::Initialize(Settings& locSettings)
{
    thisclassname = "UserDrivingForce";

    Nphases = locSettings.Nphases;

    Mode.Allocate(Nphases,Nphases);
    Value.Allocate(Nphases,Nphases);
    Teq.Allocate(Nphases,Nphases);
    LatentHeat.Allocate(Nphases,Nphases);

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void UserDrivingForce::ReadInput(std::string InputFileName)
{
    Info::WriteLineInsert("UserDrivingForce input");
    Info::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };

    std::stringstream data;
	data << inp.rdbuf();
	ReadInput(data);

    inp.close();
}

void UserDrivingForce::ReadInput(std::stringstream& inp)
{
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    // Reading driving force modes and parameters for each phase pair
    for (size_t n = 0; n < Nphases; n++)
    for (size_t m = n; m < Nphases; m++)
    {
        stringstream idx;
        idx << "_" << n << "_" << m;
        if(n == m)
        {
            Mode(n, m) = UserDrivingForceModes::None;
        }
        else
        {
            string locMode = UserInterface::ReadParameterK(inp, moduleLocation, "UDF_Mode" + idx.str(), false, "NONE");

            if(locMode == "NONE")    Mode(n, m) = Mode(m, n) = UserDrivingForceModes::None;
            if(locMode == "VALUE")   Mode(n, m) = Mode(m, n) = UserDrivingForceModes::Value;
            if(locMode == "FORMULA") Mode(n, m) = Mode(m, n) = UserDrivingForceModes::Formula;
        }

        switch(Mode(n, m))
        {
            case UserDrivingForceModes::Value:
            {
                Value(n,m) = UserInterface::ReadParameterD(inp, moduleLocation, "UDF_Value" + idx.str(), true, 0.0);
                Value(m,n) = -Value(n,m);

                LatentHeat(n,m) = LatentHeat(m,n) = 0.0;
                Teq(n,m)        = Teq(m,n)        = 1.0;                        /// Set to 1 to avoid accidental division by zero
                break;
            }
            case UserDrivingForceModes::Formula:
            {
                Value(n,m) = 0.0;
                Value(m,n) = -Value(n,m);

                LatentHeat(n,m) = UserInterface::ReadParameterD(inp, moduleLocation, "UDF_LatentHeat" + idx.str(), true, 0.0);
                LatentHeat(m,n) = -LatentHeat(n,m);
                Teq(n,m) = UserInterface::ReadParameterD(inp, moduleLocation, "UDF_Teq" + idx.str(), true, 1.0);
                Teq(m,n) = Teq(n,m);
                break;
            }
            case UserDrivingForceModes::None:
            default:
            {
                Value(n,m)       = Value(m,n)      = 0.0;
                LatentHeat(n,m)  = LatentHeat(m,n) = 0.0;
                Teq(n,m)         = Teq(m,n)        = 1.0;                       /// Set to 1 to avoid accidental division by zero
                break;
            }
        }
    }

    Info::WriteLine();
    Info::WriteBlankLine();
}

void UserDrivingForce::SetDrivingForce(PhaseField& Phase,
                                       DrivingForce& dGab)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend()-1; ++alpha)
            for(auto beta = alpha + 1;
                     beta != Phase.Fields(i,j,k).cend(); ++beta)
            {
                int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

                if(Mode(pIndexA,pIndexB) == UserDrivingForceModes::Value)
                {
                    dGab.Raw(i,j,k).add_asym1(alpha->index, beta->index, Value(pIndexA,pIndexB));
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void UserDrivingForce::SetDrivingForce(PhaseField& Phase,
                                       DrivingForce& dGab,
                                       Temperature& Tx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend()-1; ++alpha)
            for(auto beta = alpha + 1;
                     beta != Phase.Fields(i,j,k).cend(); ++beta)
            {
                int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

                switch(Mode(pIndexA,pIndexB))
                {
                    case UserDrivingForceModes::None:
                    {
                        break;
                    }
                    case UserDrivingForceModes::Formula:
                    {
                        double dG_AB = LatentHeat(pIndexA,pIndexB)*
                                       (Tx(i,j,k) - Teq(pIndexA,pIndexB))/
                                        Teq(pIndexA,pIndexB);

                        dGab.Raw(i,j,k).add_asym1(alpha->index, beta->index, dG_AB);
                        break;
                    }
                    case UserDrivingForceModes::Value:
                    {
                        dGab.Raw(i,j,k).add_asym1(alpha->index, beta->index, Value(pIndexA,pIndexB));
                        break;
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void UserDrivingForce::SetDrivingForce(PhaseField& Phase, DrivingForce& dGab,
                                       int indexA, int indexB, double dGvalue)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Interface(i,j,k))
        {
            if (Phase.Fields(i,j,k).present(indexA) and
                Phase.Fields(i,j,k).present(indexB))
            {
                dGab.Raw(i,j,k).add_asym1(indexA, indexB, dGvalue);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

}// namespace openphase
