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

using namespace openphase;

int main(int argc, char** argv)
 {

    Settings OPSettings(DefaultInputFileName);

    BoundaryConditions       BC    (OPSettings);
    DoubleObstacle           DO    (OPSettings);
    ElasticProperties        EP    (OPSettings);
    ElasticitySolverSpectral ES    (OPSettings);
    InterfaceDiffusion       ID    (OPSettings);
    InterfaceProperties      IP    (OPSettings);
    PhaseField               Phase (OPSettings);
    RunTimeControl           RTC   (OPSettings);

    // Phase field Indies
    int Cu  = 0;
    int ZnO = 0;

    if (RTC.Restart)
    {
        Phase.Read(BC, RTC.tStart);
    }
    else
    {
        const int Nx = OPSettings.Nx;
        const int Ny = OPSettings.Ny;
        const int Nz = OPSettings.Nz;
        Initializations::Single(Phase, 0, BC, OPSettings);

        // Initialize Nanocluster
        const double CWidth  = Nx*4/6.;
        const double CHeight = Nz/4;
        const double CRadius = (4*pow(CHeight,2) + pow(CWidth,2))/(8*CHeight);
        const double s0      = 0.4;   // Relative surface position
        const double x0      = Nx/2.0;
        const double y0      = Ny/2.0;
        const double z0      = s0 * Nz - (CRadius - CHeight);
        Cu = Initializations::Sphere(Phase, 1, CRadius, x0, y0, z0, BC, OPSettings, false);

        // Initialize Substrate
        const double Lx = Nx + 2*OPSettings.iWidth;
        const double Ly = Ny + 2*OPSettings.iWidth;
        const double Lz = Nz/2.0;
        ZnO = Initializations::Rectangular(Phase, 2, Lx, Ly, Lz, Nx/2., Ny/2., Nz/4., BC, OPSettings, false);
        Phase.Finalize(BC);
    }

    // Initialize properties
    EP.SetGrainsProperties(Phase);
    IP.Set(Phase);

    std::string FileName = OPSettings.TextDir + "TimeLog.csv";
    std::fstream log(FileName, std::ios::trunc | std::ios::out);
    Info::Write("Starting simulation...");
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        if(RTC.WriteToScreen())
        {
            const double AvInterfaceEnergie = DO.Energy(Phase,IP);
            const double dV        = pow(OPSettings.dx, 3);
            const double VolumeCu  = dV*Phase.FieldsStatistics[Cu].Volume;
            const double VolumeZnO = dV*Phase.FieldsStatistics[ZnO].Volume;

            std::array<std::stringstream,2> line;
            line[1] << std::scientific << std::setprecision(16);
            Info::WriteTimeStep(RTC.tStep, RTC.nSteps);
            Info::WriteWithLog(line, RTC.tStep, "Simulation Time  [  s]", RTC.SimulationTime);
            Info::WriteWithLog(line, RTC.tStep, "Interface energy [  J]", AvInterfaceEnergie);
            Info::WriteWithLog(line, RTC.tStep, "Volume Cu        [m^3]", VolumeCu);
            Info::WriteWithLog(line, RTC.tStep, "Volume ZnO       [m^3]", VolumeZnO);
            Info::WriteLineToLogfile(log, line, RTC.tStep);
        }

        if(RTC.WriteRawData())
        {
            Phase.Write(RTC.tStep);
        }

        if (RTC.WriteVTK())
        {
            Phase.WriteVTK(RTC.tStep, OPSettings);

            //Solve for Mechanical Equilibrium
            EP.SetGrainsProperties(Phase);
            EP.SetEffectiveElasticConstants(Phase);
            EP.SetEffectiveTransformationStretches(Phase, IP);
            EP.WriteEigenStrainsVTK (RTC.tStep, OPSettings);
            ES.Solve(EP, BC, RTC.dt);

            EP.WriteStressesVTK       (RTC.tStep, OPSettings);
            EP.WriteElasticStrainsVTK (RTC.tStep, OPSettings);
            EP.WriteTotalStrainsVTK   (RTC.tStep, OPSettings);
        }

        IP.Set(Phase);
        ID.CalculatePhaseFieldIncrements(Phase, IP);
        Phase.MergeIncrements(BC, RTC.dt, false);
    }
    return 0;
}
