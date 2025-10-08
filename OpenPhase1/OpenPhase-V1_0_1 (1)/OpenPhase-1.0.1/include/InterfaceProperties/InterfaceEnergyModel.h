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

#ifndef INTERFACEENERGYMODEL_H
#define INTERFACEENERGYMODEL_H

#include "Base/Includes.h"
#include "Base/UserInterface.h"
#include "CoordinationShells.h"
namespace openphase
{
class ProbabilityDistribuitons;

enum class InterfaceEnergyModels{Ext, Iso, Cubic, HexBoettger, HexSun, HexYang, Faceted};

class InterfaceEnergyModel                                                      ///< Interface mobility models implementation class
{
 public:
    InterfaceEnergyModels Model;

    InterfaceEnergyModel() // @suppress("Class members should be properly initialized")
    {
        Model = InterfaceEnergyModels::Iso;

        Energy = 0.0;
        MaxEnergy = 0.0;
        Epsilon1 = 0.0;
        Epsilon2 = 0.0;
        Epsilon3 = 0.0;
        Epsilon4 = 0.0;
        FamilyFacets = 0;
        SecondFacetFlag = false;
    };
    void ReadInputIso(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the isotropic model
    {
        Energy = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Sigma") + counter);

        MaxEnergy = Energy;
        Model = InterfaceEnergyModels::Iso;
    };
    void ReadInputCubic(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the cubic model
    {
        Energy   = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Sigma") + counter);
        Epsilon1 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE") + counter);

        MaxEnergy = Energy;
        Model = InterfaceEnergyModels::Cubic;
    };
    void ReadInputHexBoettger(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the hexagonal model by Boettger et al
    {
        Energy   = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Sigma") + counter);
        Epsilon1 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE") + counter);

        MaxEnergy = Energy;
        Model = InterfaceEnergyModels::HexBoettger;
    };
    void ReadInputHexSun(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the hexagonal model by Sun et al
    {
        Energy   = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Sigma") + counter);
        Epsilon1 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE20") + counter, false, -0.026);
        Epsilon2 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE40") + counter, false,  0.0);
        Epsilon3 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE60") + counter, false,  0.0);
        Epsilon4 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE66") + counter, false,  0.003);

        MaxEnergy = Energy;
        Model = InterfaceEnergyModels::HexSun;
    };
    void ReadInputHexYang(std::stringstream& inp, int moduleLocation, std::string counter) ///< Read parameters of the hexagonal model by Yang et al
    {
        Energy   = UserInterface::ReadParameterD(inp, moduleLocation, std::string("Sigma") + counter);
        Epsilon1 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE1") + counter, false,-0.02);
        Epsilon2 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE2") + counter, false, 0.15);
        Epsilon3 = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE3") + counter, false, 0.15);

        MaxEnergy = Energy;
        Model = InterfaceEnergyModels::HexYang;
    };
    void ReadInputFaceted(std::stringstream& inp, int moduleLocation, std::string counter) ///< Read parameters of the Faceted model
    {
        FamilyFacets = UserInterface::ReadParameterI(inp, moduleLocation, std::string("FamilyOfFacets") + counter);
        SecondFacetFlag = UserInterface::ReadParameterB(inp, moduleLocation, std::string("SingleFacetFlag") + counter, false, false);

        locFacets.resize(FamilyFacets);
        FacetVector.resize(FamilyFacets);
        EnergyFacet.resize(FamilyFacets,0);
        EpsilonFacet.resize(FamilyFacets,0);

        for(size_t NoFacets = 0 ; NoFacets < FamilyFacets; NoFacets++)
        {
            std::stringstream converter;
            converter << "_" << NoFacets;
            std::string counter1 = converter.str();

            locFacets[NoFacets] = UserInterface::ReadParameterV3(inp, moduleLocation, std::string("F") + counter + counter1 , true);
            EnergyFacet[NoFacets] = UserInterface::ReadParameterD(inp, moduleLocation, std::string("SigmaF") + counter + counter1);
            EpsilonFacet[NoFacets] = UserInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonF") + counter + counter1);
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

        double maxlocEnergyFacet = 0.0;
        for(size_t i = 0; i < EnergyFacet.size() ; i++)
        {
            if ( EnergyFacet[i] > maxlocEnergyFacet)
            {
                maxlocEnergyFacet = EnergyFacet[i];
            }
        }
        MaxEnergy = maxlocEnergyFacet;
        Model = InterfaceEnergyModels::Faceted;
    };
    double Calculate(dVector3& locNormal, bool faceted_b = false)                                       ///< Calculate anisotropic interface mobility for cubic symmetry grains
    {
        double locEnergy = 0.0;
        switch(Model)
        {
            case InterfaceEnergyModels::Iso:
            {
                locEnergy = Energy;
                break;
            }
            case InterfaceEnergyModels::Cubic:
            {
                locEnergy = Energy*(1.0 + Epsilon1 * (1.5 - 2.5*(pow(locNormal[0], 4) +
                                                                 pow(locNormal[1], 4) +
                                                                 pow(locNormal[2], 4))));
                break;
            }
            case InterfaceEnergyModels::HexBoettger:
            {
                locEnergy = Energy*(1.0 - Epsilon1 *
                            (      std::pow(locNormal[0], 6) -
                                   std::pow(locNormal[1], 6) -
                            15.0 * std::pow(locNormal[0], 4) * locNormal[1]*locNormal[1] +
                            15.0 * std::pow(locNormal[1], 4) * locNormal[0]*locNormal[0] +
                            (5.0 * std::pow(locNormal[2], 4) -
                             5.0 * std::pow(locNormal[2], 2) +
                                   std::pow(locNormal[2], 6))));
                break;
            }
            case InterfaceEnergyModels::HexSun:
            {
                locEnergy = Energy*(1.0
                        + Epsilon1*sqrt(5.0/16.0/Pi)*(3.0*locNormal[2]*locNormal[2] - 1.0)
                        + Epsilon2*3.0/16.0/sqrt(Pi)*(35.0*pow(locNormal[2],4) -
                                                     30.0*locNormal[2]*locNormal[2] + 3.0)
                        + Epsilon3*sqrt(13.0/Pi)/32.0*(231.0*pow(locNormal[2],6) -
                                                      315.0*pow(locNormal[2],4) +
                                                      105.0*locNormal[2]*locNormal[2] - 5.0)
                        + Epsilon4*sqrt(6006.0/Pi)/64.0*(pow(locNormal[0],6) -
                                                        15.0*pow(locNormal[0],4)*locNormal[1]*locNormal[1] +
                                                        15.0*locNormal[0]*locNormal[0]*pow(locNormal[1],4) -
                                                        pow(locNormal[1],6)));
                break;
            }
            case InterfaceEnergyModels::HexYang:
            {
                locEnergy = Energy*(1.0
                            + Epsilon1 * pow((3.0*locNormal[2]*locNormal[2] - 1.0), 2)
                            + Epsilon2 * pow((locNormal[0]*locNormal[0]*locNormal[0] - 3.0*locNormal[0]*locNormal[1]*locNormal[1]),2)
                            * pow((9.0*locNormal[2]*locNormal[2] - 1.0 + Epsilon3),2));
                break;
            }
            case InterfaceEnergyModels::Faceted:
            {
                double Inclination = 0.0;
                double EnergyA = 0.0;
                double EpsilonA = 0.0;
                int FacetN = 0;
                CalculateInclination(locNormal, Inclination, EnergyA, EpsilonA, FacetN);
                locEnergy = EnergyA * EpsilonA * EpsilonA * pow(sqrt((EpsilonA*EpsilonA*pow(cos(Inclination),2.0)) + pow(sin(Inclination),2.0)),-3.0);
                if(SecondFacetFlag)
                {
                    locEnergy = EnergyA * pow(EpsilonA + (1.0 - EpsilonA) * fabs(tan(Inclination)) * tanh(pow(fabs(tan(Inclination)),-1)),-1);
                }
                break;
            }
            case InterfaceEnergyModels::Ext:
            {
                break;
            }
        }
        return locEnergy;
    };

    ///< Calculate inclination angle between facet normal and interface normal
    void CalculateInclination(dVector3& locGrad, double& Inclination, double& EnergyA, double& EpsilonA, int FacetN)
    {
        Inclination = Pi;
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
                    EnergyA = EnergyFacet[i];
                    EpsilonA = EpsilonFacet[i];
                    FacetN = i;
                }
            }
        }
    }


    size_t FamilyFacets;                                                        ///< Number of type of facets, e.g. <100>, <111>, <110>, etc...
    bool SecondFacetFlag;
    std::vector<dVector3> locFacets;                                            ///< Facet type as a vector [100], [111], [110], etc
    std::vector<std::vector<dVector3>> FacetVector;                             ///< Set of planes for different facet type ...
    std::vector<double> EnergyFacet;                                            ///< Interface energy prefactor per Facet
    std::vector<double> EpsilonFacet;                                           ///< Interface energy anisotropy parameter per Facet

    double Energy;                                                              ///< Interface energy prefactor
    double MaxEnergy;                                                           ///< Minimum interface energy for a phase pair
    double Epsilon1;                                                            ///< Interface energy anisotropy parameter
    double Epsilon2;                                                            ///< Interface energy anisotropy parameter
    double Epsilon3;                                                            ///< Interface energy anisotropy parameter
    double Epsilon4;                                                            ///< Interface energy anisotropy parameter
};
}// namespace openphase
#endif
