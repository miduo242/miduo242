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
#include "Info.h"
#include "Initializations.h"
#include "Magnetism/LinearMagneticSolver.h"
#include "Mechanics/ElasticitySolverSpectral.h"
#include "Mechanics/ElasticProperties.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "Tools/TimeInfo.h"

namespace op = openphase;

int main(int argc, char *argv[])
{
    std::string InputFileName = op::DefaultInputFileName;
    if (argc > 1) InputFileName = argv[1];

    // Detects floating point errors on Linux machines
    //feenableexcept(FE_INEXACT);  // NOTE: THIS DETECTS ROUNDING ERRORS
    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID);    // domain error occurred
    //feenableexcept(FE_OVERFLOW);   // the result was too large
    //feenableexcept(FE_UNDERFLOW);  // the result was too small

    op::Settings                 OPSettings (InputFileName);

    op::BoundaryConditions       BC    (OPSettings);
    op::ElasticProperties        EP    (OPSettings);
    op::ElasticitySolverSpectral ES    (OPSettings, BC);
    //op::ElasticitySolverMatrix   ESM   (OPSettings);
    op::LinearMagneticSolver     MS    (OPSettings);
    op::PhaseField               Phi   (OPSettings);
    op::RunTimeControl           RTC   (OPSettings);
    op::TimeInfo                 Timer (OPSettings, "Execution Time Statistics");

    std::fstream inpF(InputFileName, std::ios::in);
    std::stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();

    const double dx = OPSettings.dx;
    const int  moduleLocation = op::UserInterface::FindModuleLocation(inp, "MagnetoactiveElastomerLinear");
    const double R1           = op::UserInterface::ReadParameterD(inp, moduleLocation, std::string( "dR1"))/dx;
    const double R2           = op::UserInterface::ReadParameterD(inp, moduleLocation, std::string( "dR2"))/dx;
    const double R3           = op::UserInterface::ReadParameterD(inp, moduleLocation, std::string( "dR3"))/dx;
    const double distance     = op::UserInterface::ReadParameterD(inp, moduleLocation, std::string("dDis"))/dx;
    const double H0           = op::UserInterface::ReadParameterD(inp, moduleLocation, std::string("dH0"));
    const double Phi0         = op::UserInterface::ReadParameterD(inp, moduleLocation, std::string("dPhi0"));
    const double DPhi         = op::UserInterface::ReadParameterD(inp, moduleLocation, std::string("dDPhi"));
    const double accuracy     = op::UserInterface::ReadParameterD(inp, moduleLocation, std::string("dAcc"));
    const size_t MinItr       = op::UserInterface::ReadParameterD(inp, moduleLocation, std::string("iMinItr"));

    const size_t Nx = OPSettings.Nx;
    const size_t Ny = OPSettings.Ny;
    const size_t Nz = OPSettings.Nz;

    const double x0 = (Nx>1) ? Nx/2.0 : 0;
    const double y0 = (Ny>1) ? Ny/2.0 : 0;
    const double z0 = (Nz>1) ? Nz/2.0 : 0;

    const double x1 = Nx/2.0 - distance/2;
    const double x2 = Nx/2.0 + distance/2;

    op::Info::Write("Set initial conditions");
    if(RTC.Restart)
    {
        // Read restart input
        EP.Read  (BC, RTC.tStart);
        Phi.Read (BC, RTC.tStart);
    }
    else
    {
       RTC.tStart = 0;
       op::Initializations::Single(Phi, 0, BC, OPSettings);
       op::Initializations::CylinderSimple(Phi, 1, R3, Nz+20, 2, x0, y0, z0, BC,OPSettings);
       op::Initializations::Sphere(Phi, 2, R1, x1, y0, z0, BC, OPSettings);
       op::Initializations::Sphere(Phi, 2, R2, x2, y0, z0, BC, OPSettings);
    }
    MS.SetEffectiveSusceptibility(Phi, BC);

    BC.SetX(MS.chi);
    BC.SetY(MS.chi);
    BC.SetZ(MS.chi);

    EP.SetGrainsProperties(Phi);
    EP.SetEffectiveElasticConstants(Phi);
    EP.SetEffectiveTransformationStretches(Phi);

    op::Info::Write("Start Calculation");
    std::string FileName = OPSettings.TextDir + "TimeLog.csv";
    std::fstream log(FileName, std::ios::trunc | std::ios::out);
    log << std::scientific << std::setprecision(16);
    for(int tStep = RTC.tStart; tStep <= RTC.nSteps; tStep++)
    {


        // Solve for magnetic field and force density
        MS.H0x = H0*std::cos(M_PI*(Phi0+tStep*DPhi)/180.0);
        MS.H0y = H0*std::sin(M_PI*(Phi0+tStep*DPhi)/180.0);
        MS.Solve(BC,0.0001);
        MS.CalcForceDensity(EP, BC);

        // Solve for elastic deformation
        size_t Iteration = 0;
        double DeltaDistance    = 0.0;
        double DeltaDistanceOld = 0.0;
        const int j0 = std::round(y0);
        const int k0 = std::round(z0);
        std::stringstream FName;
        FName << OPSettings.TextDir << "DeltaDistance" << tStep << ".csv";
        std::fstream logDist(FName.str(), std::ios::trunc | std::ios::out);
        auto EndSolve = [&tStep, &OPSettings, &EP, &j0, &k0, &x1, &x2, &dx, &Iteration, &DeltaDistance, &DeltaDistanceOld, &logDist, &MinItr, &accuracy]()
        {
            DeltaDistance = 0.0;
            for (int i = std::round(x1); i < std::round(x2); i++)
            {
                DeltaDistance += EP.ElasticStrains(i,j0,k0)[0]*dx;
            }
            double diff = std::abs((DeltaDistanceOld - DeltaDistance)/DeltaDistance);
            std::cout << "Elastic Solver Interation: "<<  Iteration << "," << "\tParticle Displacement: " << DeltaDistance << "," << "\t Relative Particle Displacement Change: " << diff << "\n";
            logDist   << Iteration << "," << DeltaDistance << "," << diff << "\n";
            //EP.WriteElasticStrainsVTK (tStep, OPSettings);

            DeltaDistanceOld = DeltaDistance;
            Iteration++;
            if (diff < accuracy and Iteration > MinItr) return true;
            else return false;
        };
        ES.Solve(EP, BC, RTC.dt, EndSolve);

        //ESM.Solve(EP,BC,1.0e-7);

        if (!(tStep%RTC.tScreenWrite))
        {
            op::dVector3 H0 = MS.H0(0,0,0); // Resulting external H-field
            op::dVector3 H  = MS.H(0,0,0);  // Resulting external H-field
            op::dVector3 B  = MS.B(0,0,0);  // Resulting external B-field
            op::dVector3 M  = MS.M(std::round(x1), j0, k0); // Resulting particle magnetisation

            DeltaDistance = 0.0;
            for (int i = std::round(x1); i < std::round(x2); i++)
            {
                DeltaDistance += EP.ElasticStrains(i,j0,k0)[0]*dx;
            }

            std::array<std::stringstream,2> line;
            line[1] << std::scientific << std::setprecision(16);
            op::Info::WriteTimeStep(tStep, RTC.nSteps);
            op::Info::WriteWithLog(line, tStep, "Phi", Phi0+tStep*DPhi);
            op::Info::WriteWithLog(line, tStep, "Applied external H0", H0.length());
            op::Info::WriteWithLog(line, tStep, "Resulting external H",  H.length());
            op::Info::WriteWithLog(line, tStep, "Resulting external B",  B.length());
            op::Info::WriteWithLog(line, tStep, "Magnetisation of particle", M.length());
            op::Info::WriteWithLog(line, tStep, "DeltaDistance", DeltaDistance);
            op::Info::WriteWithLog(line, tStep, "Elastic Energy", EP.Energy());
            op::Info::WriteLineToLogfile(log, line, tStep);
        }

        //  Output to VTK file
        if (RTC.WriteVTK())
        {
            EP.WriteElasticStrainsVTK (tStep, OPSettings);
            //EP.WriteStressesVTK       (tStep, OPSettings);
            //EP.WriteForceDensityVTK   (tStep, OPSettings);
            //ES.WriteVTK               (tStep, OPSettings);
            MS.WriteVTK               (tStep, OPSettings);
            //Phi.WriteDistortedVTK     (tStep, OPSettings, ES, 1000000);
            Phi.WriteVTK              (tStep, OPSettings);
        }  //  Output to file

        //  Output raw Data to file
        if (RTC.WriteRawData())
        {
            EP.Write(tStep);
            Phi.Write(tStep);
        }
    } //end time loop
    return 0;
}
