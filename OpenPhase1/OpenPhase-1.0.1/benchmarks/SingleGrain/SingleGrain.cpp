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

#include "Settings.h"
#include "RunTimeControl.h"
#include "InterfaceProperties.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "Initializations.h"
#include "BoundaryConditions.h"
#include "Tools/CSVParser.h"
#include "Tools/TimeInfo.h"
#include "Info.h"

using namespace std;
using namespace openphase;
int main()
{
    Settings                    OPSettings;
    OPSettings.ReadInput();

    RunTimeControl              RTC(OPSettings);
    PhaseField                  Phi(OPSettings);
    DoubleObstacle              DO(OPSettings);
    InterfaceProperties         IP(OPSettings);
    BoundaryConditions          BC(OPSettings);

    // Initialize phase-fields
    Initializations::Single(Phi, 0, BC, OPSettings);
    int locIndex = Initializations::Sphere(Phi, 1, OPSettings.Nx*0.4,
    (OPSettings.Nx)/2, (OPSettings.Ny)/2, (OPSettings.Nz)/2, BC, OPSettings);

    //-------------------------------------------------//

    double Vol0 = Phi.FieldsStatistics[locIndex].Volume;
    double R0 = exp(log(0.75/Pi*Vol0)/3.0);
    double R02 = R0*R0;

    // Create output file
    std::vector<std::string> headernames {"tstep", "R^2(simulation)", "R^2(analytic)", "Ratio"};
    std::vector<double> dataInit {0, R02, R02, 1};
    CSVParser::WriteHeader("R_2_graph.dat", headernames);
    CSVParser::WriteData("R_2_graph.dat", dataInit);

    cout << "Entering the Time Loop!!!" << endl;

    for(RTC.tStep = 0; RTC.tStep < RTC.nSteps+1; RTC.IncrementTimeStep())
    {
        IP.Set(Phi);
        DO.CalculatePhaseFieldIncrements(Phi, IP);
        Phi.NormalizeIncrements(BC, RTC.dt);
        Phi.MergeIncrements(BC, RTC.dt);

        if (RTC.WriteVTK())
        {
            /// Write data in VTK format
            Phi.WriteVTK(RTC.tStep,OPSettings);
            /// Output of the grain radius
            double Vol = Phi.FieldsStatistics[locIndex].Volume;
            double Radius = exp(log(0.75/Pi*Vol)/3.0);

            std::vector<double> datat;

            datat.push_back(double(RTC.tStep+1));
            datat.push_back(Radius*Radius);
            double RadiusTheor = R02 - 4.0*IP.InterfaceMobility(0,1).MaxMobility*
                                           IP.InterfaceEnergy(0,1).MaxEnergy*(RTC.tStep+1)*
                                           RTC.dt/OPSettings.dx/OPSettings.dx;
            if (RadiusTheor >= 0) datat.push_back(RadiusTheor);
            else datat.push_back(0.0);
            if (Radius > 0) {datat.push_back((Radius*Radius)/(R02 - 4.0*IP.InterfaceMobility(0,1).MaxMobility*
                                                                        IP.InterfaceEnergy(0,1).MaxEnergy*(RTC.tStep+1)*
                                                                        RTC.dt/OPSettings.dx/OPSettings.dx));};
            CSVParser::WriteData("R_2_graph.dat", datat);
        }

        if (RTC.WriteToScreen())
        {
            int I_En = DO.AverageEnergyDensity(Phi, IP);
            std::string message; message = Info::GetStandard("InterfaceEnergy", to_string(I_En));
            Info::WriteTimeStep(RTC.tStep, RTC.nSteps, message);
        }
    }
    return 0;
}
