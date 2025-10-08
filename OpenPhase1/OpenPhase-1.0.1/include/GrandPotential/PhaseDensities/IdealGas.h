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
 *   Main contributors :   Raphael Schiedung
 *
 */

#ifndef GRANDPOTENTIALIDEALGAS_H
#define GRANDPOTENTIALIDEALGAS_H

#include "GrandPotential/PhaseDensities/PhaseDensity.h"

namespace openphase
{
class Settings;
namespace GrandPotential
{struct IdealGas: PhaseDensity                                                  ///< EOS for ideal gas
{
    IdealGas(size_t PhaseIdxInp): PhaseDensity(PhaseIdxInp){};

    void Initialize(Settings& locSettings) override;                            ///<  Initializes global settings
    void ReadInput(std::stringstream& InputFile, int moduleLocation) override;  ///<  Reads input parameters from a file

    double PhasePressure (double Temperature, const Tensor<double,1>& ChemicalPotential) const override
    {
        using namespace PhysicalConstants;
        using namespace std;
        assert(Temperature > 0);
        double locPressure = 0.0;
        for (size_t comp = 0; comp < Ncomp; ++comp)
        {
            const double& M  = MolarMass[comp];
            const double& T  = Temperature;
            const double& mu = ChemicalPotential({comp});
            locPressure += 2*sqrt(2.0*M_PI*M*R*T)*M_PI*exp(mu/R/T)*M*R*R*T*T/N_A/N_A/N_A/N_A/h/h/h;
        }
        return locPressure;
    }
    double PhaseConcentration (double Temperature, double ChemicalPotential, size_t comp) const override
    {
        using namespace PhysicalConstants;
        using namespace std;
        assert(Temperature > 0);
        assert(comp < Ncomp);
        const double& M  = MolarMass[comp];
        const double& T  = Temperature;
        const double& mu = ChemicalPotential;
        return 2*sqrt(2.0*M_PI*M*R*T)*M_PI*exp(mu/R/T)*M*R*T/N_A/N_A/N_A/N_A/h/h/h;
    }
    double PhaseSusceptibility (double Temperature, double ChemicalPotential, size_t comp) const override
    {
        using namespace PhysicalConstants;
        using namespace std;
        assert(Temperature > 0);
        assert(comp < Ncomp);
        const double& M  = MolarMass[comp];
        const double& T  = Temperature;
        const double& mu = ChemicalPotential;
        return 2*sqrt(2.0*M_PI*M*R*T)*M_PI*exp(mu/R/T)*M/N_A/N_A/N_A/N_A/h/h/h;
    }

    double PhasePressure (double height, double Temperature, const Tensor<double,1>& ChemicalPotential) const override
    {
        using namespace PhysicalConstants;
        using namespace std;
        assert(Temperature > 0);
        double locPressure = 0.0;
        for (size_t comp = 0; comp < Ncomp; ++comp)
        {
            const double& M  = MolarMass[comp];
            const double& T  = Temperature;
            const double& g  = GravityAcceleration;
            const double& mu = ChemicalPotential({comp});
            locPressure += 2*sqrt(2.0*M_PI*M*R*T)*M_PI*exp((mu-M*g*height)/R/T)*M*R*R*T*T/N_A/N_A/N_A/N_A/h/h/h;
        }
        return locPressure;
    }
    double PhaseConcentration  (double height, double Temperature, double ChemicalPotential, size_t comp) const override
    {
        using namespace PhysicalConstants;
        using namespace std;
        assert(Temperature > 0);
        assert(comp < Ncomp);
        const double& M  = MolarMass[comp];
        const double& T  = Temperature;
        const double& g  = GravityAcceleration;
        const double& mu = ChemicalPotential;
        return 2*sqrt(2.0*M_PI*M*R*T)*M_PI*exp((mu-M*g*height)/R/T)*M*R*T/N_A/N_A/N_A/N_A/h/h/h;
    }
    double PhaseSusceptibility  (double height, double Temperature, double ChemicalPotential, size_t comp) const override
    {
        using namespace PhysicalConstants;
        using namespace std;
        assert(Temperature > 0);
        assert(comp < Ncomp);
        const double& M  = MolarMass[comp];
        const double& T  = Temperature;
        const double& g  = GravityAcceleration;
        const double& mu = ChemicalPotential;
        return 2*sqrt(2.0*M_PI*M*R*T)*M_PI*exp((mu-M*g*height)/R/T)*M/N_A/N_A/N_A/N_A/h/h/h;
    }

    static constexpr auto thisclassname = "GrandPotential::PhaseDensities::IdealGas";

 protected:
    std::vector<double> MolarMass;                                              ///< Mass of Species
    double GravityAcceleration;                                                 ///< Gravitational acceleration [m/s^2]
};
}// namespace openphase::GrandPotential
}// namespace openphase
#endif
