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
 *   File created :   2012
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Raphael Schiedung
 *
 */

#ifndef INTERFACEPROPERTIES_H
#define INTERFACEPROPERTIES_H

#include "Base/Includes.h"
#include "Base/NodeAB.h"
#include "InterfaceProperties/InterfaceEnergyModel.h"
#include "InterfaceProperties/InterfaceMobilityModel.h"

namespace openphase
{
class ElasticProperties;
class PhaseField;
class Settings;
class Temperature;
class Orientations;
class Nucleation;

class OP_EXPORTS InterfaceProperties : public OPObject                                     ///< Interface properties module.
{
 public:

    size_t Nphases;

    int Nx;                                                                     ///< X dimension
    int Ny;                                                                     ///< Y dimension
    int Nz;                                                                     ///< Z dimension
    int dNx;                                                                    ///< Active X dimension
    int dNy;                                                                    ///< Active Y dimension
    int dNz;                                                                    ///< Active Z dimension
    double dx;                                                                  ///< Grid spacing

    double TripleJunctionFactor;                                                ///< Parameter in front of the triple junction energy term
    double maxSigma;                                                            ///< Maximum interface energy
    double maxMu;                                                               ///< Maximum interface mobility

    double maxSigmaPhase;                                                       ///< Maximum interface energy from Input File of all the phases
    double R;                                                                   ///< Gas constant

    Resolutions Resolution;                                                     ///< Grid resolution for phase fields. Can be Single or Double

    Storage3D< NodeAB, 0 > IntProperties;                                       ///< 3D interface properties storage
    Storage3D< NodeAB, 0 > IntPropertiesDR;                                     ///< 3D interface properties storage
    Matrix<InterfaceEnergyModel> InterfaceEnergy;                               ///< Interface energy model storage for all pairs of phases
    Matrix<InterfaceMobilityModel> InterfaceMobility;                           ///< Interface mobility model storage for all pairs of phases
    Matrix< bool   > RespectParentBoundaries;                                   ///< Prevents growth of product phase beyond parent grain boundaries

    std::string VTKDir;                                                         ///< Directory-path added in front of VTK files
    std::string RawDataDir;                                                     ///< Directory-path added in front of restart files

    InterfaceProperties(){};
    InterfaceProperties(Settings& locSettings,
                        const std::string InputFileName = DefaultInputFileName)
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }
    void Initialize(Settings& locSettings) override;                            ///< Initializes all variables, allocates storages
    void ReInitialize(const PhaseField& Phase);                                 ///< Reallocates 3D storages using new system dimensions
    void ReadInput(const std::string InputFileName) override;
    void ReadInput(std::stringstream& inp) override;
    void Coarsen(const PhaseField& Phase);
    void Refine(const PhaseField& Phase);

    int Bcells(void) const
    {
        return IntProperties.Bcells();
    }

    int BcellsDR(void) const
    {
        return IntPropertiesDR.Bcells();
    }

    void SetFacetOrientation(const PhaseField&Phase);

    void set_energy(const int i, const int j, const int k, const int alpha, const int beta, const double value)
    {
        IntProperties(i,j,k).set_sym1(alpha, beta, value);
    };
    void set_energy_and_mobility(const int i, const int j, const int k, const int alpha, const int beta, const double Evalue, const double Mvalue)
    {
        IntProperties(i,j,k).set_sym_pair(alpha, beta, Evalue, Mvalue);
    };
    void add_energy(const int i, const int j, const int k, const int alpha, const int beta, const double value)
    {
        IntProperties(i,j,k).add_sym1(alpha, beta, value);
    };
    double get_energy(const int i, const int j, const int k, const int alpha, const int beta) const
    {
        return IntProperties(i,j,k).get_sym1(alpha, beta);
    };
    double get_energy(const ElasticProperties& EP, const int i, const int j, const int k, const int alpha, const int beta) const;

    void set_mobility(const int i, const int j, const int k, const int alpha, const int beta, const double value)
    {
        IntProperties(i,j,k).set_sym2(alpha, beta, value);
    };
    void add_mobility(const int i, const int j, const int k, const int alpha, const int beta, const double value)
    {
        IntProperties(i,j,k).add_sym2(alpha, beta, value);
    };
    double get_mobility(const int i, const int j, const int k, const int alpha, const int beta) const
    {
        return IntProperties(i,j,k).get_sym2(alpha, beta);
    };
    void clear(const int i, const int j, const int k)
    {
        IntProperties(i,j,k).clear();
    };

    /// Double resolution methods
    void set_energy_DR(const int i, const int j, const int k, const int alpha, const int beta, const double value)
    {
        IntPropertiesDR(i,j,k).set_sym1(alpha, beta, value);
    };
    void set_energy_and_mobility_DR(const int i, const int j, const int k, const int alpha, const int beta, const double Evalue, const double Mvalue)
    {
        IntPropertiesDR(i,j,k).set_sym_pair(alpha, beta, Evalue, Mvalue);
    };
    void add_energy_DR(const int i, const int j, const int k, const int alpha, const int beta, const double value)
    {
        IntPropertiesDR(i,j,k).add_sym1(alpha, beta, value);
    };
    double get_energy_DR(const int i, const int j, const int k, const int alpha, const int beta) const
    {
        return IntPropertiesDR(i,j,k).get_sym1(alpha, beta);
    };
    void set_mobility_DR(const int i, const int j, const int k, const int alpha, const int beta, const double value)
    {
        IntPropertiesDR(i,j,k).set_sym2(alpha, beta, value);
    };
    void add_mobility_DR(const int i, const int j, const int k, const int alpha, const int beta, const double value)
    {
        IntPropertiesDR(i,j,k).add_sym2(alpha, beta, value);
    };
    double get_mobility_DR(const int i, const int j, const int k, const int alpha, const int beta) const
    {
        return IntPropertiesDR(i,j,k).get_sym2(alpha, beta);
    };
    void clear_DR(const int i, const int j, const int k)
    {
        IntPropertiesDR(i,j,k).clear();
    };

    void SetSR(const PhaseField& Phase);                                          ///< Sets both, interface energy and mobility
    void SetMobilityThermalEffectSR(const PhaseField& Phase, const Temperature& Tx);///< Sets interface mobility in the entire simulation domain using user selected interface mobility models and considering temperature effect

    void SetDR(const PhaseField& Phase);                                        ///< Sets both, interface energy and mobility
    void SetMobilityThermalEffectDR(const PhaseField& Phase, const Temperature& Tx);///< Sets interface mobility in the entire simulation domain using user selected interface mobility models and considering temperature effect

    void SetwMD(const PhaseField& Phase, Orientations& OR);    // Sets both, interface energy and mobility but energy is anisotropic based on Disorientation angle

    void ReduceMobilityForNucleation(PhaseField& Phi, BoundaryConditions& BC,
                                     Nucleation& Nuc);

    void Set(const PhaseField& Phase)                                           ///< Sets both, interface energy and mobility
    {
        switch(Resolution)
        {
            case Resolutions::Single:
            {
                SetSR(Phase);
                break;
            }
            case Resolutions::Double:
            {
                SetDR(Phase);
                Coarsen(Phase);
                break;
            }
        }
    }

    void Set(const PhaseField& Phase, const Temperature& Tx)                    ///< Sets both, interface energy and mobility considering temperature effect
    {
        switch(Resolution)
        {
            case Resolutions::Single:
            {
                SetSR(Phase);
                SetMobilityThermalEffectSR(Phase, Tx);
                break;
            }
            case Resolutions::Double:
            {
                SetDR(Phase);
                SetMobilityThermalEffectDR(Phase, Tx);
                Coarsen(Phase);
                break;
            }
        }
    }
    double ReportMaximumTimeStep(void)                                          ///< Returns maximum time step for phase-field equation
    {
        return 0.25*dx*dx/(maxMu*maxSigma);
    }
    void WriteVTK(const PhaseField& Phase, const int tStep);                    ///< Writes effective interface properties in VTK format

    void WriteFacetArea(const PhaseField& Phase, double& tStep, double& DegreeTolerance, std::string OutFile);
 protected:
 private:
};

}// namespace openphase
#endif
