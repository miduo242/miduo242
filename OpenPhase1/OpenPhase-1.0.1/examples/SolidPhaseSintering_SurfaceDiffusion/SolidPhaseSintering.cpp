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
 */

#include "BoundaryConditions.h"
#include "DoubleObstacle.h"
#include "Initializations.h"
#include "InterfaceDiffusion.h"
#include "InterfaceProperties.h"
#include "Mechanics/ElasticProperties.h"
#include "Mechanics/ElasticitySolverSpectral.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include <iostream>

namespace op = openphase;

int main()
{
    // Detects floating point errors on Linux machines
    //feenableexcept(FE_INEXACT);  // NOTE: THIS DETECTS ROUNDING ERRORS
    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID);    // domain error occurred
    feenableexcept(FE_OVERFLOW);   // the result was too large
    feenableexcept(FE_UNDERFLOW);  // the result was too small

    op::Settings OPSettings(op::DefaultInputFileName);

    op::BoundaryConditions       BC    (OPSettings);
    op::DoubleObstacle           DO    (OPSettings);
    op::ElasticProperties        EP    (OPSettings);
    op::ElasticitySolverSpectral ES    (OPSettings);
    op::InterfaceDiffusion       ID    (OPSettings);
    op::InterfaceProperties      IP    (OPSettings);
    op::PhaseField               Phase (OPSettings);
    op::RunTimeControl           RTC   (OPSettings);


    if (RTC.Restart == true)
    {
        Phase.Read(BC, RTC.tStart);
    }
    else
    {
        const double R  = (OPSettings.Nx+1)/5;
        const double x0 = (OPSettings.Nx+1)/2;
        const double y0 = (OPSettings.Ny+1)/2;
        const double z0 = (OPSettings.Nz+1)/2;

        op::Initializations::Single(Phase, 0, BC, OPSettings);
        op::Initializations::Sphere(Phase, 1, R, x0-R, y0, z0, BC, OPSettings);
        op::Initializations::Sphere(Phase, 1, R, x0+R, y0, z0, BC, OPSettings);
    }
    Phase.SetBoundaryConditions(BC);

    // Initialize elastic properties
    op::Info::WriteSimple("Starting simulation...");
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        if (RTC.WriteToScreen())
        {
            //Solve for Mechanical Equilibrium
            EP.SetGrainsProperties(Phase);
            EP.SetEffectiveElasticConstants(Phase);
            EP.SetEffectiveTransformationStretches(Phase, IP);
            ES.Solve(EP, BC, RTC.dt);

            Phase.WriteVTK            (RTC.tStep, OPSettings);
            EP.WriteStressesVTK       (RTC.tStep, OPSettings);
            EP.WriteElasticStrainsVTK (RTC.tStep, OPSettings);
            EP.WriteTotalStrainsVTK   (RTC.tStep, OPSettings);
        }

        if(RTC.WriteRawData())
        {
            Phase.Write(RTC.tStep);
        }

        if(RTC.WriteVTK())
        {
            op::Info::WriteTimeStep(RTC.tStep, RTC.nSteps);
            op::Info::Write("Simulation Time  [s]", RTC.SimulationTime);
            op::Info::Write("Interface Energy [J]", DO.Energy(Phase,IP));
        }

        IP.Set(Phase);
        ID.CalculatePhaseFieldIncrements(Phase, IP);
        Phase.MergeIncrements(BC, RTC.dt, false);

    }
    return 0;
}
