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

#include "Info.h"
#include "Base/Includes.h"
#include "Settings.h"
#include "RunTimeControl.h"
#include "PhaseField.h"
#include "Initializations.h"
#include "BoundaryConditions.h"
#include "Mechanics/ElasticProperties.h"
#include "Mechanics/ElasticitySolverSpectral.h"

using namespace std;
using namespace openphase;

int main(int argc, char *argv[])
{
    Settings                    OPSettings;
    OPSettings.ReadInput();

    RunTimeControl              RTC(OPSettings);
    BoundaryConditions          BC(OPSettings);
    PhaseField                  Phi(OPSettings);
    ElasticitySolverSpectral    ES(OPSettings);
    ElasticProperties           EP(OPSettings);

    double centerx = (OPSettings.Nx+1)/2.0;
    double centery = (OPSettings.Ny+1)/2.0;
    double centerz = (OPSettings.Nz+1)/2.0;

    int iRadius = 10;

    Initializations::Single(Phi, 0, BC, OPSettings);
    Initializations::Sphere(Phi, 1, iRadius, centerx, centery, centerz, BC, OPSettings);
    Info::WriteStandard("Initial Microstructure", "Set");

    EP.SetGrainsProperties(Phi);
    Info::WriteStandard("Grains Properties", "Set");

    EP.SetEffectiveElasticConstants(Phi);
    EP.SetEffectiveTransformationStretches(Phi);

    Info::WriteStandard("Calculation", "Started");

    int niterations = ES.Solve(EP, BC, RTC.dt);

    Info::WriteStandard("Solver iterations", std::to_string(niterations));
    Info::WriteStandard("Elastic energy density", std::to_string(EP.Energy()));

    Phi.WriteDistortedVTK(0,OPSettings,EP,true);
    EP.WriteStressesVTK(0,OPSettings);
    EP.WriteTotalStrainsVTK(0,OPSettings);

    fstream out("LineStresses.dat", ios::out);
    string dlm = ",";

    out << std::fixed;
    out << "X" << dlm << "Sigma_11" << dlm << "Sigma_33" << dlm << "Eshelby_Sigma_11" << dlm << "Eshelby_Sigma_33" << endl;

    double C11 = EP.PhaseElasticConstants[0](0,0);
    double C12 = EP.PhaseElasticConstants[0](0,1);
    double nu = C12/(C12 + C11);

    vStrain locEigenStrain = EP.EigenStrains(centerx,centery,centerz);

    double epsilon = locEigenStrain[0];
    double sigma_0 = 2.0/3.0*epsilon*(C11 + 2.0*C12)*(1.0 - 2.0*nu)/(1.0 - nu);
    for (int i = centerx; i < centerx + iRadius; i++)
    {
        dMatrix3x3 locStress = EP.Stresses(i, centery, centerz).tensor();

        double sigma_11 = -sigma_0;
        double sigma_33 = -sigma_0;
        out << double(i - centerx)/(OPSettings.Nx-1)   << dlm
            << locStress(0,0)                          << dlm
            << locStress(2,2)                          << dlm
            << sigma_11                                << dlm
            << sigma_33                                << endl;
    }
    for (int i = centerx + iRadius; i < OPSettings.Nx; i++)
    {
        dMatrix3x3 locStress = EP.Stresses(i, centery, centerz).tensor();

        double sigma_11 = -sigma_0*pow(iRadius/double(i-centerx), 3);
        double sigma_33 = 0.5*sigma_0*pow(iRadius/double(i-centerx), 3);

        out << double((i - centerx))/(OPSettings.Nx-1) << dlm
            << locStress(0,0)                          << dlm
            << locStress(2,2)                          << dlm
            << sigma_11                                << dlm
            << sigma_33                                << endl;
    }
    out.close();

    Info::WriteLine();
    simulation_end();
    return 0;
}
