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
 *   File created :   2019
 *   Main contributors :   Oleg Shchyglo; Alexander Monas; Efim Borukhovich
 *
 */

#ifndef INTERFACEMOBILITYMODEL_H
#define INTERFACEMOBILITYMODEL_H

#include "Base/Includes.h"
#include "Base/UserInterface.h"
#include "CoordinationShells.h"
namespace openphase
{
enum class InterfaceMobilityModels{Ext, Iso, Cubic, HexBoettger, HexSun, HexYang, Faceted};

class InterfaceMobilityModel                                                    ///< Interface mobility models implementation class
{
 public:
    InterfaceMobilityModels Model;

    InterfaceMobilityModel()
    {
        Model = InterfaceMobilityModels::Iso;
        Mobility = 0.0;
        MaxMobility = 0.0;
        ActivationEnergy = 0.0;
        Epsilon1 = 0.0;
        Epsilon2 = 0.0;
        Epsilon3 = 0.0;
        Epsilon4 = 0.0;
        FamilyFacets = 0;
    };
    void ReadInputIso(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the isotropic model
    {
        Mobility = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Mu") + counter);
        ActivationEnergy = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Q") + counter, false, 0.0);

        MaxMobility = Mobility;
        Model = InterfaceMobilityModels::Iso;
    };
    void ReadInputCubic(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the cubic model
    {
        Mobility = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Mu") + counter);
        Epsilon1 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM") + counter);
        ActivationEnergy = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Q") + counter, false, 0.0);

        MaxMobility = Mobility;//*std::min(1.0, (1.0 - Epsilon1));
        Model = InterfaceMobilityModels::Cubic;
    };
    void ReadInputHexBoettger(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the hexagonal model by Boettger et al
    {
        Mobility = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Mu") + counter);
        Epsilon1 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM") + counter);
        ActivationEnergy = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Q") + counter, false, 0.0);

        MaxMobility = Mobility;//*std::min(1.0, (1.0 - Epsilon1));
        Model = InterfaceMobilityModels::HexBoettger;
    };
    void ReadInputHexSun(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the hexagonal model by Sun et al
    {
        Mobility = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Mu") + counter);
        Epsilon1 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM20") + counter, false, -0.026);
        Epsilon2 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM40") + counter, false,  0.0);
        Epsilon3 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM60") + counter, false,  0.0);
        Epsilon4 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM66") + counter, false,  0.003);

        ActivationEnergy = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Q") + counter, false, 0.0);

        MaxMobility = Mobility;
        Model = InterfaceMobilityModels::HexSun;
    };
    void ReadInputHexYang(std::stringstream& inp, int moduleLocation, std::string counter) ///< Read parameters of the hexagonal model by Yang et al
    {
        Mobility = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Mu") + counter);
        ActivationEnergy = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Q") + counter, false, 0.0);

        Epsilon1 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM1") + counter, false,-0.02);
        Epsilon2 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM2") + counter, false, 0.15);
        Epsilon3 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM3") + counter, false, 0.15);

        MaxMobility = Mobility;
        Model = InterfaceMobilityModels::HexYang;
    };
    void ReadInputFaceted(std::stringstream& inp, int moduleLocation, std::string counter) ///< Read parameters of the Faceted model
    {
        FamilyFacets = UserInterface::ReadParameterI(inp, moduleLocation, std::string("FamilyOfFacets") + counter);
        locFacets.resize(FamilyFacets);
        FacetVector.resize(FamilyFacets);
        MobilityFacet.resize(FamilyFacets,0);
        EpsilonFacet.resize(FamilyFacets,0);

        for(size_t NoFacets = 0 ; NoFacets < FamilyFacets; NoFacets++)
        {
            std::stringstream converter;
            converter << "_" << NoFacets;
            std::string counter1 = converter.str();

            locFacets[NoFacets] = UserInterface::ReadParameterV3(inp, moduleLocation, std::string("F") + counter + counter1 , true);
            MobilityFacet[NoFacets] = UserInterface::ReadParameterD(inp, moduleLocation, std::string("MuF") + counter + counter1);
            EpsilonFacet[NoFacets] = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonMF") + counter + counter1);
            ActivationEnergy = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Q") + counter + counter1, false, 0.0);

        }
        for(size_t i = 0; i < locFacets.size() ; i++)
        {
            FacetVector[i] = CoordinationShells::findPermutations(locFacets, i);
        }

        for(size_t i = 0; i < locFacets.size(); i++)
        {
            for(size_t j = 0 ; j < FacetVector[i].size() ; j++)
            {
                FacetVector[i][j].normalize();
            }
        }

        double maxEpsilonFacet = 0.0;
        for(size_t i = 0; i < EpsilonFacet.size() ; i++)
        {
            if ( EpsilonFacet[i] > maxEpsilonFacet)
            {
                maxEpsilonFacet = EpsilonFacet[i];
            }
        }
        double maxlocMobilityFacet = 0.0;
        for(size_t i = 0; i < MobilityFacet.size() ; i++)
        {
            if ( MobilityFacet[i] > maxlocMobilityFacet)
            {
                maxlocMobilityFacet = MobilityFacet[i];
            }
        }
        MaxMobility = maxlocMobilityFacet;
        Model = InterfaceMobilityModels::Faceted;
    };

    double Calculate(dVector3& locNormal)                                       ///< Calculate anisotropic interface mobility for cubic symmetry grains
    {
        double locMobility = 0.0;
        switch(Model)
        {
            case InterfaceMobilityModels::Iso:
            {
                locMobility = Mobility;
                break;
            }
            case InterfaceMobilityModels::Cubic:
            {
                locMobility = Mobility*(1.0 - Epsilon1 * (1.5 - 2.5*(pow(locNormal[0], 4) +
                                                                     pow(locNormal[1], 4) +
                                                                     pow(locNormal[2], 4))));
                break;
            }
            case InterfaceMobilityModels::HexBoettger:
            {
                locMobility = Mobility*(1.0 + Epsilon1 *
                            (      std::pow(locNormal[0], 6) -
                                   std::pow(locNormal[1], 6) -
                            15.0 * std::pow(locNormal[0], 4) * locNormal[1]*locNormal[1] +
                            15.0 * std::pow(locNormal[1], 4) * locNormal[0]*locNormal[0] +
                            (5.0 * std::pow(locNormal[2], 4) -
                             5.0 * std::pow(locNormal[2], 2) +
                                   std::pow(locNormal[2], 6))));
                break;
            }
            case InterfaceMobilityModels::HexSun:
            {
                locMobility = Mobility*(1.0
                        -Epsilon1*sqrt(5.0/16.0/Pi)*(3.0*locNormal[2]*locNormal[2] - 1.0)
                        -Epsilon2*3.0/16.0/sqrt(Pi)*(35.0*pow(locNormal[2],4) -
                                                     30.0*locNormal[2]*locNormal[2] + 3.0)
                        -Epsilon3*sqrt(13.0/Pi)/32.0*(231.0*pow(locNormal[2],6) -
                                                      315.0*pow(locNormal[2],4) +
                                                      105.0*locNormal[2]*locNormal[2] - 5.0)
                        -Epsilon4*sqrt(6006.0/Pi)/64.0*(pow(locNormal[0],6) -
                                                        15.0*pow(locNormal[0],4)*locNormal[1]*locNormal[1] +
                                                        15.0*locNormal[0]*locNormal[0]*pow(locNormal[1],4) -
                                                        pow(locNormal[1],6)));
                break;
            }
            case InterfaceMobilityModels::HexYang:
            {
                locMobility = Mobility*(1.0
                            - Epsilon1 * pow((3.0*locNormal[2]*locNormal[2] - 1.0), 2)
                            - Epsilon2 * pow((locNormal[0]*locNormal[0]*locNormal[0] - 3.0*locNormal[0]*locNormal[1]*locNormal[1]),2)
                            * pow((9.0*locNormal[2]*locNormal[2] - 1.0 + Epsilon3),2));
                break;
            }
            case InterfaceMobilityModels::Faceted:
            {
                double Inclination = 0.0;
                double MobilityA = 0.0;
                double EpsilonMA = 0.0;
                CalculateInclination(locNormal, Inclination, MobilityA, EpsilonMA);
                locMobility = MobilityA * (EpsilonMA + (1.0 - EpsilonMA) * fabs(tan(Inclination)) * tanh(pow(fabs(tan(Inclination)),-1)));
                break;
            }
            case InterfaceMobilityModels::Ext:
            {
                break;
            }
        }
        return locMobility;
    };

    ///< Calculate inclination angle between facet normal and interface normal
    void CalculateInclination(dVector3& locGrad, double& Inclination, double& MobilityA, double& EpsilonMA)
    {
        Inclination = M_PI;
        dVector3 Facet;
        for(size_t i = 0; i < locFacets.size(); i++)
        {
            for(size_t j = 0 ; j < FacetVector[i].size() ; j++)
            {
                Facet =  FacetVector[i][j];
                double dotVec = acos((locGrad * Facet) / (locGrad.abs() * Facet.abs()));
                if(fabs(dotVec) <= fabs(Inclination))
                {
                    Inclination = dotVec;
                    MobilityA = MobilityFacet[i];
                    EpsilonMA = EpsilonFacet[i];
                }
            }
        }
    }

    size_t FamilyFacets;                                                        ///< Number of type of facets, e.g. <100>, <111>, <110>, etc...
    std::vector<dVector3> locFacets;                                            ///< Facet type as a vector [100], [111], [110], etc
    std::vector<std::vector<dVector3>> FacetVector;                             ///< Set of planes for different facet type ...
    std::vector<double> MobilityFacet;                                          ///< Interface mobility prefactor per Facet
    std::vector<double> EpsilonFacet;                                           ///< Interface mobility anisotropy parameter per Facet

    double Mobility;                                                            ///< Interface mobility prefactor
    double MaxMobility;                                                         ///< Minimum interface mobility for a phase pair
    double Epsilon1;                                                            ///< Interface mobility anisotropy parameter
    double Epsilon2;                                                            ///< Interface mobility anisotropy parameter
    double Epsilon3;                                                            ///< Interface mobility anisotropy parameter
    double Epsilon4;                                                            ///< Interface mobility anisotropy parameter
    double ActivationEnergy;                                                    ///< Interface mobility activation energy
};
}// namespace openphase
#endif
