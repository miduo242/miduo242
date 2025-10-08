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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Johannes Goerler
 *
 */

#include "Mechanics/ElasticProperties.h"
#include "Mechanics/ElasticitySolverSpectral.h"
#include "Settings.h"
#include "Mechanics/SymmetryVariants.h"
#include "Mechanics/ElasticityModels/ElasticitySteinbach.h"
#include "Mechanics/ElasticityModels/ElasticityKhachaturyan.h"
#include "Mechanics/PlasticFlow/PlasticFlowNeuberMethods.h"
#include "Info.h"
#include "PhaseField.h"
#include "Composition.h"
#include "Orientations.h"
#include "DrivingForce.h"
#include "BoundaryConditions.h"
#include "Temperature.h"
#include "Base/UserInterface.h"
#include "Velocities.h"
#include "VTK.h"
#include "InterfaceProperties.h"
#include "Tools.h"
#include "AdvectionHR/AdvectionHR.h"

namespace openphase
{
using namespace std;

ElasticProperties::ElasticProperties(Settings& locSettings,
                                     const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void ElasticProperties::Initialize(Settings& locSettings)
{
    thisclassname = "ElasticProperties";

    KeepAspectRatio = false;
    PreventShear = false;

    TotalNx = locSettings.TotalNx;
    OffsetX = locSettings.OffsetX;
    TotalNy = locSettings.TotalNy;
    OffsetY = locSettings.OffsetY;
    TotalNz = locSettings.TotalNz;
    OffsetZ = locSettings.OffsetZ;

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;

    dNx = locSettings.dNx;
    dNy = locSettings.dNy;
    dNz = locSettings.dNz;

    dx = locSettings.dx;

    EModel = ElasticityModels::Khachaturyan;

    Nphases = locSettings.Nphases;
    Ncomp = locSettings.Ncomp;
    Names = locSettings.ElementNames;

    StrainAccuracy = 1.0e-6;
    MAXIterations = 100;

    AppliedStrainMask.set_to_zero();
    AppliedStressMask.set_to_zero();

    AppliedStrain.set_to_zero();
    AppliedStrainOLD.set_to_zero();

    AppliedStress.set_to_zero();
    AppliedStressOLD.set_to_zero();

    AverageElasticConstants.set_to_zero();
    MAXElasticConstants.set_to_zero();

    AverageStrain.set_to_zero();
    RemeshedStrain.set_to_zero();
    StrainToRemesh.set_to_zero();
    AverageDeformationGradient.set_to_unity();

    size_t Bcells = locSettings.Bcells;

    Stresses.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
    DeformationJumps.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
    DeformationGradientsTotal.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
    DeformationGradientsEigen.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
    EffectiveElasticConstants.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Bcells);
    Displacements.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, 0);

    PhaseElasticConstants.Allocate(Nphases);
    PhaseCompliences.Allocate(Nphases);
    PhaseTransformationStretches.Allocate(Nphases);

    PhaseAlpha.Allocate(Nphases);
    PhaseGamma.Allocate(Nphases);
    Tref.Allocate(Nphases);
    PoissonRatio.Allocate(Nphases);

    ElasticConstants.Allocate(Nphases);
    Compliances.Allocate(Nphases);
    TransformationStretches.Allocate(Nphases);

    Alpha.Allocate(Nphases);
    Gamma.Allocate(Nphases);

    if(Ncomp > 0)
    {
        PhaseKappa.Allocate({Nphases, Ncomp});
        PhaseLambda.Allocate({Nphases, Ncomp});
        Cref.Allocate({Nphases, Ncomp});

        Kappa.Allocate({Nphases, Ncomp});
        Lambda.Allocate({Nphases, Ncomp});
    }

    NeuberCorrection.resize(Nphases);

    for(size_t alpha = 0; alpha != Nphases; alpha++)
    {
        PhaseAlpha[alpha].set_to_zero();
        PhaseGamma[alpha].set_to_zero();

        Tref[alpha] = 0.0;

        PhaseTransformationStretches[alpha].set_to_unity();
        PhaseElasticConstants[alpha].set_to_zero();
        PhaseCompliences[alpha].set_to_zero();

        for(size_t comp = 0; comp != Ncomp; comp++)
        {
            Cref({alpha, comp}) = 0.0;
            PhaseLambda({alpha, comp}).set_to_zero();
            PhaseKappa({alpha, comp}).set_to_zero();
        }
    }

    MAXElasticConstants.set_to_zero();

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Stresses,Bcells,)
    {
        DeformationGradientsTotal(i,j,k).set_to_unity();                        //NOTE: zero strain
        DeformationGradientsEigen(i,j,k).set_to_unity();                        //NOTE: zero strain
        Stresses(i,j,k).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Variants.Initialize(locSettings);

    VTKDir = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;

    Info::Write(thisclassname, "Initialized");
}

void ElasticProperties::ReadInput(const string InputFileName)
{
    Info::WriteLineInsert("ElasticProperties input");
    Info::Write("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in);

    if (!inp)
    {
        Info::WriteExit("File " + InputFileName + " could not be opened",thisclassname, "ReadInput()");
        exit(1);
    };
    std::stringstream data;
    data << inp.rdbuf();

    ReadInput(data);
   
    inp.close();

    Variants.ReadInput(InputFileName);

    Info::WriteLine();
}

void ElasticProperties::ReadInput(stringstream& inp)
{
    for(size_t alpha = 0; alpha != Nphases; alpha++)
    {
        PhaseElasticConstants[alpha].set_to_zero();
    }

    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    StrainAccuracy = UserInterface::ReadParameterD(inp, moduleLocation, "StrainAccuracy", false, StrainAccuracy);
    MAXIterations  = UserInterface::ReadParameterD(inp, moduleLocation, "MAXIterations", false, MAXIterations);

    string tmp1 = UserInterface::ReadParameterK(inp, moduleLocation, "EModel", false, "KHACHATURYAN");

    if (tmp1 == "KHACHATURYAN")
    {
        EModel = ElasticityModels::Khachaturyan;
    }
    else if (tmp1 == "STEINBACH")
    {
        EModel = ElasticityModels::Steinbach;
    }
    else if (tmp1 == "VOIGT")
    {
        EModel = ElasticityModels::Voigt;
    }
    else if (tmp1 == "REUSS")
    {
        EModel = ElasticityModels::Reuss;
    }
    else
    {
        Info::WriteWarning("No or wrong elasticity model specified!\nThe default \"Khachaturyan\" model is used!", thisclassname, "ReadInput()");
    }

    /// Mechanical boundary conditions

    std::vector<std::string> BCnames {"BCX", "BCY","BCZ","BCYZ","BCXZ","BCXY"};
    std::vector<std::string> BCvalueNames {"BCValueX", "BCValueY","BCValueZ","BCValueYZ","BCValueXZ","BCValueXY"};

    for(int n = 0; n < 6; n++)
    {
        string tmp2 = UserInterface::ReadParameterK(inp, moduleLocation, BCnames[n], 0);

        if (tmp2 == "FREEBOUNDARIES")
        {
            AppliedStressMask[n] = 1.0;
            AppliedStress[n] = 0.0;
        }
        if (tmp2 == "APPLIEDSTRAIN")
        {
            AppliedStrain[n] = UserInterface::ReadParameterD(inp, moduleLocation, BCvalueNames[n]);
            AppliedStrainMask[n] = 1.0;
        }
        if (tmp2 == "APPLIEDSTRAINRATE")
        {
            AppliedStrainRate[n] = UserInterface::ReadParameterD(inp, moduleLocation, BCvalueNames[n]);
            AppliedStrainMask[n] = 1.0;
        }
        if (tmp2 == "APPLIEDSTRESS")
        {
            AppliedStressMask[n] = 1.0;
            AppliedStress[n] = UserInterface::ReadParameterD(inp, moduleLocation, BCvalueNames[n]);
        }
    }

    string tmp3 = UserInterface::ReadParameterK(inp, moduleLocation, "Restrict", 0);

    if (tmp3 == "ASPECTRATIO")
    {
        KeepAspectRatio = true;
    }

    if (tmp3 == "SHEAR")
    {
        PreventShear = true;
    }

    for(size_t alpha = 0; alpha != Nphases; alpha++)
    {
        PhaseElasticConstants[alpha].set_to_zero();
    }

    // Reading elastic constants
    for(size_t pIndex = 0; pIndex < Nphases; pIndex++)
    {
        // Read elastic moduli if provided
        std::stringstream KConv;
        std::stringstream EConv;
        std::stringstream LConv;
        std::stringstream GConv;
        std::stringstream NConv;
        std::stringstream MConv;

        KConv << "K_" << pIndex;
        EConv << "E_" << pIndex;
        LConv << "L_" << pIndex;
        GConv << "G_" << pIndex;
        NConv << "Nu_"<< pIndex;
        MConv << "M_" << pIndex;

        double K      = UserInterface::ReadParameterD(inp, moduleLocation, KConv.str(), false, 0); ///< Bulk modulus
        double E      = UserInterface::ReadParameterD(inp, moduleLocation, EConv.str(), false, 0); ///< Young's modulus
        double lambda = UserInterface::ReadParameterD(inp, moduleLocation, LConv.str(), false, 0); ///< First Lame
        double G      = UserInterface::ReadParameterD(inp, moduleLocation, GConv.str(), false, 0); ///< Shear modulus
        double nu     = UserInterface::ReadParameterD(inp, moduleLocation, NConv.str(), false, 0); ///< Poisson's ratio
        double M      = UserInterface::ReadParameterD(inp, moduleLocation, MConv.str(), false, 0); ///< P-wave modulus

        if      (K      != 0.0 and E      != 0.0) {lambda = 3.0*K*(3.0*K-E)/(9.0*K-E);    G = 3.0*K*E/(9.0*K-E);}
        else if (K      != 0.0 and lambda != 0.0) {                                       G = 3*(K-lambda)/2.0;}
        else if (K      != 0.0 and G      != 0.0) {lambda = K-2.0*G/3.0;}
        else if (K      != 0.0 and nu     != 0.0) {lambda = 3.0*K*nu/(1+nu);              G = 3*K*(1-2.0*nu)/(2.0*(1.0+nu));}
        else if (K      != 0.0 and M      != 0.0) {lambda = (3.0*K-M)/2;                  G = 3.0*(M-K)/4.0;}
        else if (E      != 0.0 and lambda != 0.0) {                                       G = (E-3.0*lambda+std::sqrt(E*E+9.0*lambda*lambda+2.0*E*lambda))/4.0;}
        else if (E      != 0.0 and G      != 0.0) {lambda = G*(E-2.0*G)/(3.0*G-E);}
        else if (E      != 0.0 and nu     != 0.0) {lambda = E*nu/((1.0+nu)*(1.0-2.0*nu)); G = E/(2.0*(1.0+nu));}
        else if (E      != 0.0 and M      != 0.0) {Info::WriteExit("Choice of elastic moduli is not unique", thisclassname, "ReadInput"); std::exit(EXIT_FAILURE);}
        else if (lambda != 0.0 and nu     != 0.0) {                                       G = lambda*(1.0-2.0*nu)/(2.0*nu);}
        else if (lambda != 0.0 and M      != 0.0) {                                       G = (M-lambda)/2.0;}
        else if (G      != 0.0 and nu     != 0.0) {lambda = 2.0*G*nu/(1.0-2.0*nu);}
        else if (G      != 0.0 and M      != 0.0) {lambda = M-2.0*G;}
        else if (nu     != 0.0 and M      != 0.0) {lambda = M*nu/(1.0-nu);                G = M*(1.0-2.0*nu)/(2.0*(1.0-nu));}

        if (std::fabs(lambda) > DBL_EPSILON and std::fabs(G) > DBL_EPSILON)
        {
            PhaseElasticConstants[pIndex](0,0) = 2.0*G + lambda;
            PhaseElasticConstants[pIndex](1,1) = 2.0*G + lambda;
            PhaseElasticConstants[pIndex](2,2) = 2.0*G + lambda;

            PhaseElasticConstants[pIndex](0,1) = lambda;
            PhaseElasticConstants[pIndex](0,2) = lambda;
            PhaseElasticConstants[pIndex](1,2) = lambda;

            PhaseElasticConstants[pIndex](1,0) = lambda;
            PhaseElasticConstants[pIndex](2,0) = lambda;
            PhaseElasticConstants[pIndex](2,1) = lambda;

            PhaseElasticConstants[pIndex](3,3) = G;
            PhaseElasticConstants[pIndex](4,4) = G;
            PhaseElasticConstants[pIndex](5,5) = G;

        }
        else
        {
            stringstream converter;
            converter << "C" << "_" << pIndex;
            PhaseElasticConstants[pIndex] = UserInterface::ReadParameterM6x6(inp, moduleLocation, converter.str(),false,dMatrix6x6::ZeroTensor());
            if(PhaseElasticConstants[pIndex] == dMatrix6x6::ZeroTensor())
            {
                for(int ii =  1; ii <= 6; ii++)
                for(int jj = ii; jj <= 6; jj++)
                {
                    stringstream Cij;
                    Cij << "C" << ii << jj << "_" << pIndex;

                    PhaseElasticConstants[pIndex](ii-1,jj-1) = UserInterface::ReadParameterD(inp, moduleLocation, Cij.str(),false,0);
                    if(ii != jj)
                    {
                        PhaseElasticConstants[pIndex](jj-1,ii-1) = PhaseElasticConstants[pIndex](ii-1,jj-1);
                    }
                }
            }
        }
        PoissonRatio[pIndex] = PhaseElasticConstants[pIndex](0,1) / (PhaseElasticConstants[pIndex](0,1) + PhaseElasticConstants[pIndex](0,0));

        // Check if stiffness constants are read.
        if (PhaseElasticConstants[pIndex].norm() <= DBL_EPSILON)
        {
            std::string message = "Elastic constants for phase " + to_string(pIndex) + " not set.";
            Info::WriteExit(message, thisclassname, "ReadInput()");
            exit(3);
        }
        else
        {
            PhaseCompliences[pIndex] = PhaseElasticConstants[pIndex].inverted();
        }
    }

    for(size_t pIndex = 0; pIndex < Nphases; pIndex++)
    {
        stringstream converter;
        converter << "U_" << pIndex;
        PhaseTransformationStretches[pIndex] = UserInterface::ReadParameterM3x3(inp, moduleLocation, converter.str(),true,dMatrix3x3::UnitTensor());
    }
    // Considering external forces
    ConsiderExternalForces = UserInterface::ReadParameterB(inp, moduleLocation, "ConsiderExternalForces", false, false);
    if(ConsiderExternalForces)
    {
        ForceDensity.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, 1);
    }

    // Reading chemo-mechanical coupling parameters if considered
    ChemoMechanicalCoupling = UserInterface::ReadParameterB(inp, moduleLocation, "ChemoMechanicalCoupling",false,false);
    if(ChemoMechanicalCoupling)
    {
        int counter = 0;
        for(size_t comp = 0; comp < Ncomp; comp++)
        for(size_t alpha = 0; alpha != Nphases; alpha++)
        {
            stringstream converter;
            converter << Names[comp] << "_" << alpha;
            string counter = converter.str();
            Cref({alpha, comp}) = UserInterface::ReadParameterD(inp, moduleLocation, string("Cref_") + counter, false, 0.0);
        }

        for (size_t pIndex = 0; pIndex < Nphases; pIndex++)
        {
            for (size_t comp = 0; comp < Ncomp; comp++)
            {
                stringstream converter;
                converter << "Kappa_" << pIndex << "_" << Names[comp];
                PhaseKappa({ pIndex, comp }) = UserInterface::ReadParameterM6x6(inp, moduleLocation, converter.str(),false,dMatrix6x6::ZeroTensor());

                if (PhaseKappa({pIndex, comp}).norm() != 0.0)
                {
                    counter++;
                }
            }
        }

        for (size_t pIndex = 0; pIndex < Nphases; pIndex++)
        for (size_t comp = 0; comp < Ncomp; comp++)
        {
            stringstream converter;
            converter << "Lambda_" << pIndex << "_" << Names[comp];
            PhaseLambda({ pIndex, comp }) = UserInterface::ReadParameterM3x3(inp, moduleLocation, converter.str(),false,dMatrix3x3::ZeroTensor());

            if (PhaseLambda({pIndex, comp}).norm() != 0.0)
            {
                counter++;
            }
        }

        if(counter == 0)
        {
            Info::WriteWarning("ChemoMechanicalCoupling is ON but no coupling parameters specified", thisclassname, "ReadInput()");
        }
    }
    // Reading thermo-mechanical coupling parameters if considered
    ThermoMechanicalCoupling = UserInterface::ReadParameterB(inp, moduleLocation, "ThermoMechanicalCoupling",false,false);
    if(ThermoMechanicalCoupling)
    {
        int counter = 0;
        for(size_t alpha = 0; alpha != Nphases; alpha++)
        {
            stringstream converter;
            converter << "Tref_" << alpha;
            Tref[alpha] = UserInterface::ReadParameterD(inp, moduleLocation, converter.str(), false, 0.0);
        }

        for (size_t pIndex = 0; pIndex < Nphases; pIndex++)
        {
            stringstream converter;
            converter << "Gamma_" << pIndex;

            double ScalarGamma = UserInterface::ReadParameterD(inp, moduleLocation, converter.str(), false, 0.0);
            std::cout << pIndex << std::endl;
            if (ScalarGamma != 0.0)
            {
                for (int ii = 0; ii < 6; ii++)
                for (int jj = 0; jj < 6; jj++)
                {
                    PhaseGamma[pIndex](ii,jj) = ScalarGamma;
                }
            }
            else
            {
                for (int ii = 1; ii <= 6; ii++)
                for (int jj = ii; jj <= 6; jj++)
                {
                    converter << "Gamma" << ii << jj << "_" << pIndex;

                    PhaseGamma[pIndex](ii - 1, jj - 1) =
                        UserInterface::ReadParameterD(inp, moduleLocation, converter.str(), false, 0.0);
                    if (ii != jj)
                    {
                        PhaseGamma[pIndex](jj - 1, ii - 1) = PhaseGamma[pIndex](ii - 1, jj - 1);
                    }
                }
            }

            if (PhaseGamma[pIndex].norm() != 0.0)
            {
                counter++;

                stringstream GammaNX;
                GammaNX << "Gamma_" << pIndex;
                Info::Write(GammaNX.str(),PhaseGamma[pIndex],6);
            }
        }

        for (size_t pIndex = 0; pIndex < Nphases; pIndex++)
        {
            stringstream converter;
            for (int ii = 1; ii <= 3; ii++)
            for (int jj = ii; jj <= 3; jj++)
            {
                converter << "Alpha" << ii << jj << "_" << pIndex;

                PhaseAlpha[pIndex](ii - 1, jj - 1) =
                    UserInterface::ReadParameterD(inp, moduleLocation, converter.str(), false, 0.0);
                if (ii != jj)
                {
                    PhaseAlpha[pIndex](jj - 1, ii - 1) = PhaseAlpha[pIndex](ii - 1, jj - 1);
                }
            }
            if(PhaseAlpha[pIndex].norm() != 0.0)
            {
                counter++;
                stringstream AlphaNX;
                AlphaNX << "Alpha_" << pIndex;
                Info::Write(AlphaNX.str(),PhaseAlpha[pIndex],6);
            }
        }

        if(counter == 0)
        {
            Info::WriteWarning("ThermoMechanicalCoupling is ON but no coupling parameters specified", thisclassname, "ReadInput()");
        }
    }
    // Reading Neuber parameters
    for(size_t n = 0; n < Nphases; n++)
    {
        stringstream converter;
        converter << n;
        NeuberCorrection[n].Active = UserInterface::ReadParameterB(inp, moduleLocation, "NeuberCorrection_" + converter.str(), false, false);
        if (NeuberCorrection[n].Active)
        {
            NeuberCorrection[n].YoungsModulus = UserInterface::ReadParameterD(inp, moduleLocation, "YoungsModulus_" + converter.str(), true, 1.0);
            NeuberCorrection[n].YieldStrength = UserInterface::ReadParameterD(inp, moduleLocation, "YieldStrength_" + converter.str(), true, 1.0);
            NeuberCorrection[n].Hardening     = UserInterface::ReadParameterD(inp, moduleLocation, "Hardening_" + converter.str(), true, 1.0);
        }
    }

    Info::WriteLine();
    Info::WriteBlankLine();
}

void ElasticProperties::Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)
{
    ///Changing box size
    dMatrix3x3 stretch;
    stretch.set_to_zero();
    stretch(0,0) = double(newNx)/double(Nx);
    stretch(1,1) = double(newNy)/double(Ny);
    stretch(2,2) = double(newNz)/double(Nz);

    RemeshedStrain[0] += stretch(0,0) - 1.0;
    RemeshedStrain[1] += stretch(1,1) - 1.0;
    RemeshedStrain[2] += stretch(2,2) - 1.0;

    StrainToRemesh[0] -= stretch(0,0) - 1.0;
    StrainToRemesh[1] -= stretch(1,1) - 1.0;
    StrainToRemesh[2] -= stretch(2,2) - 1.0;

    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    Stresses.Remesh(Nx, Ny, Nz);

    DeformationGradientsTotal.Remesh(Nx, Ny, Nz);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, DeformationGradientsTotal, 0,)
    {
        DeformationGradientsTotal(i,j,k) = stretch*DeformationGradientsTotal(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    DeformationGradientsEigen.Remesh(Nx, Ny, Nz);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, DeformationGradientsEigen, 0,)
    {
        DeformationGradientsEigen(i,j,k) = stretch*DeformationGradientsEigen(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Displacements.Reallocate(Nx, Ny, Nz);
    EffectiveElasticConstants.Remesh(Nx, Ny, Nz);

    SetBoundaryConditions(BC);

    //WriteRemeshingData(tStep, sim_time);
    Info::Write(thisclassname, "Remeshed");
}

void ElasticProperties::SetGrainsProperties(PhaseField& Phase)
{
    size_t size = Phase.FieldsStatistics.size();
    if(TransformationStretches.size() != size)
    {
        TransformationStretches.Reallocate(size);
        ElasticConstants.Reallocate(size);
        Compliances.Reallocate(size);
        if(Ncomp)
        {
            Lambda.Reallocate({size, Ncomp});
            Kappa.Reallocate({size, Ncomp});
        }
        Alpha.Reallocate(size);
        Gamma.Reallocate(size);
    }
    for(size_t alpha = 0; alpha != size; alpha++)
    if(Phase.FieldsStatistics[alpha].Exist)
    {
        size_t pIndex = Phase.FieldsStatistics[alpha].Phase;
        size_t vIndex = Phase.FieldsStatistics[alpha].Variant;

        TransformationStretches[alpha] = PhaseTransformationStretches[pIndex];
        ElasticConstants[alpha] = PhaseElasticConstants[pIndex];
        Alpha[alpha] = PhaseAlpha[pIndex];
        Gamma[alpha] = PhaseGamma[pIndex];

        if(Variants.set)
        {
            TransformationStretches[alpha].rotate(Variants(pIndex, vIndex));
            ElasticConstants[alpha].rotate(Variants(pIndex, vIndex));
            Alpha[alpha].rotate(Variants(pIndex, vIndex));
            Gamma[alpha].rotate(Variants(pIndex, vIndex));
        }

        TransformationStretches[alpha].rotate(Phase.FieldsStatistics[alpha].Orientation.RotationMatrix);
        ElasticConstants[alpha].rotate(Phase.FieldsStatistics[alpha].Orientation.RotationMatrix);
        Alpha[alpha].rotate(Phase.FieldsStatistics[alpha].Orientation.RotationMatrix);
        Gamma[alpha].rotate(Phase.FieldsStatistics[alpha].Orientation.RotationMatrix);

        Compliances[alpha] = ElasticConstants[alpha].inverted();

        for(size_t comp = 0; comp != Ncomp; comp++)
        {
            Lambda({alpha, comp}) = PhaseLambda({pIndex, comp});
            Kappa({alpha, comp})  = PhaseKappa({pIndex, comp});

            if(Variants.set)
            {
                Lambda({alpha, comp}).rotate(Variants(pIndex, vIndex));
                Kappa({alpha, comp}).rotate(Variants(pIndex, vIndex));
            }

            Lambda({alpha, comp}).rotate(Phase.FieldsStatistics[alpha].Orientation.RotationMatrix);
            Kappa({alpha, comp}).rotate(Phase.FieldsStatistics[alpha].Orientation.RotationMatrix);
        }
    }

    for(size_t alpha = 0; alpha != size; alpha++)
    if(Phase.FieldsStatistics[alpha].Exist)
    for(int n = 0; n < 6; n++)
    for(int m = 0; m < 6; m++)
    {
        MAXElasticConstants(n,m) = std::max(ElasticConstants[alpha](n,m), MAXElasticConstants(n,m));
    }
}

void ElasticProperties::SetBoundaryConditions(const BoundaryConditions& BC)
{
    // Only periodic BC are correct. For non periodic boundary conditions gradients should be treated differently.

    if(dNx) BC.SetX(DeformationGradientsTotal);
    if(dNy) BC.SetY(DeformationGradientsTotal);
    if(dNz) BC.SetZ(DeformationGradientsTotal);

    if(dNx) BC.SetX(DeformationGradientsEigen);
    if(dNy) BC.SetY(DeformationGradientsEigen);
    if(dNz) BC.SetZ(DeformationGradientsEigen);

    if(dNx) BC.SetX(Stresses);
    if(dNy) BC.SetY(Stresses);
    if(dNz) BC.SetZ(Stresses);

    if(dNx) BC.SetX(EffectiveElasticConstants);
    if(dNy) BC.SetY(EffectiveElasticConstants);
    if(dNz) BC.SetZ(EffectiveElasticConstants);
}

ElasticProperties& ElasticProperties::operator= (const ElasticProperties& rhs)
{
    // protect against self-assignment and copy of unitialized object
    if (this != &rhs and rhs.thisclassname == "ElasticProperties")
    {
        thisclassname = rhs.thisclassname;

        Nx = rhs.Nx;
        Ny = rhs.Ny;
        Nz = rhs.Nz;

        dNx = rhs.dNx;
        dNy = rhs.dNy;
        dNz = rhs.dNz;

        dx = rhs.dx;

        Nphases = rhs.Nphases;
        Ncomp   = rhs.Ncomp;
        Names = rhs.Names;
        EModel = rhs.EModel;

        if (DeformationGradientsTotal.IsNotAllocated())
        {
            DeformationGradientsTotal.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, rhs.DeformationGradientsTotal.Bcells());
            DeformationGradientsEigen.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, rhs.DeformationGradientsEigen.Bcells());
            Stresses.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, rhs.Stresses.Bcells());
            EffectiveElasticConstants.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, rhs.EffectiveElasticConstants.Bcells());
            Displacements.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, 0);
        }
        else if (not DeformationGradientsTotal.IsSize(rhs.Nx, rhs.Ny, rhs.Nz))
        {
            DeformationGradientsTotal.Reallocate(Nx, Ny, Nz);
            DeformationGradientsEigen.Reallocate(Nx, Ny, Nz);
            Stresses.Reallocate(Nx, Ny, Nz);
            EffectiveElasticConstants.Reallocate(Nx, Ny, Nz);
            Displacements.Reallocate(Nx, Ny, Nz);
        }

        if (PhaseElasticConstants.IsNotAllocated())
        {
            PhaseElasticConstants.Allocate(Nphases);
            PhaseTransformationStretches.Allocate(Nphases);

            PhaseCompliences.Allocate(Nphases);
            PhaseKappa.Allocate({Nphases, Ncomp});
            PhaseLambda.Allocate({Nphases, Ncomp});
            PhaseGamma.Allocate(Nphases);
            PhaseAlpha.Allocate(Nphases);
            Tref.Allocate(Nphases);
            Cref.Allocate({Nphases, Ncomp});
        }
        else if (PhaseElasticConstants.size() != rhs.Nphases)
        {
            PhaseElasticConstants.Reallocate(Nphases);
            PhaseTransformationStretches.Reallocate(Nphases);
            PhaseCompliences.Reallocate(Nphases);
            PhaseKappa.Reallocate({Nphases, Ncomp});
            PhaseLambda.Reallocate({Nphases, Ncomp});
            PhaseGamma.Reallocate(Nphases);
            PhaseAlpha.Reallocate(Nphases);
            Tref.Reallocate(Nphases);
            Cref.Reallocate({Nphases, Ncomp});
        }
        for(size_t alpha = 0; alpha != Nphases; alpha++)
        {
            PhaseAlpha[alpha] = rhs.PhaseAlpha[alpha];
            PhaseGamma[alpha] = rhs.PhaseGamma[alpha];

            Tref[alpha] = rhs.Tref[alpha];

            PhaseTransformationStretches[alpha] = rhs.PhaseTransformationStretches[alpha];
            PhaseElasticConstants[alpha] = rhs.PhaseElasticConstants[alpha];
            PhaseCompliences[alpha] = rhs.PhaseCompliences[alpha];

            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                PhaseLambda({alpha, comp}) = rhs.PhaseLambda({alpha, comp});
                PhaseKappa({alpha, comp}) = rhs.PhaseKappa({alpha, comp});
                Cref({alpha, comp}) = rhs.Cref({alpha, comp});
            }
        }
        ElasticConstants = rhs.ElasticConstants;
        TransformationStretches = rhs.TransformationStretches;
        Compliances = rhs.Compliances;
        Alpha = rhs.Alpha;
        Gamma = rhs.Gamma;
        Lambda = rhs.Lambda;
        Kappa = rhs.Kappa;

        RemeshedStrain = rhs.RemeshedStrain;
        StrainToRemesh = rhs.StrainToRemesh;
        AverageStrain = rhs.AverageStrain;
        AppliedStress = rhs.AppliedStress;
        AppliedStrain = rhs.AppliedStrain;
        AppliedStrainOLD = rhs.AppliedStrainOLD;
        MAXElasticConstants = rhs.MAXElasticConstants;

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DeformationGradientsTotal,DeformationGradientsTotal.Bcells(),)
        {
            DeformationGradientsTotal(i,j,k) = rhs.DeformationGradientsTotal(i,j,k);
            DeformationGradientsEigen(i,j,k) = rhs.DeformationGradientsEigen(i,j,k);
            Stresses(i,j,k) = rhs.Stresses(i,j,k);

            EffectiveElasticConstants(i,j,k) = rhs.EffectiveElasticConstants(i,j,k);
            Displacements(i,j,k) = rhs.Displacements(i,j,k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    return *this;
}

void ElasticProperties::SetEffectiveTransformationStretches(PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, DeformationGradientsEigen, DeformationGradientsEigen.Bcells(),)
    {
        DeformationGradientsEigen(i,j,k).set_to_zero();
        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        {
            DeformationGradientsEigen(i,j,k) += TransformationStretches[alpha->index]*alpha->value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::SetEffectiveTransformationStretches(PhaseField& Phase, InterfaceProperties& IP)
{
    SetEffectiveTransformationStretches(Phase);
    CalculateInterfaceStress(Phase, IP);
}

void ElasticProperties::SetEffectiveTransformationStretches(PhaseField& Phase, Composition& Cx)
{
    SetEffectiveTransformationStretches(Phase);
    CalculateVegardsExpansion(Phase, Cx);
}

void ElasticProperties::SetEffectiveTransformationStretches(PhaseField& Phase, Temperature& Tx)
{
    SetEffectiveTransformationStretches(Phase);
    CalculateThermalExpansion(Phase, Tx);
}

void ElasticProperties::SetEffectiveTransformationStretches(PhaseField& Phase, Composition& Cx, Temperature& Tx)
{
    SetEffectiveTransformationStretches(Phase);
    CalculateVegardsExpansion(Phase, Cx);
    CalculateThermalExpansion(Phase, Tx);
}
void ElasticProperties::CalculateVegardsExpansion(PhaseField& Phase, Composition& Cx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, DeformationGradientsEigen, DeformationGradientsEigen.Bcells(),)
    {
        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        {
            size_t pIndex = Phase.FieldsStatistics[alpha->index].Phase;

            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                double delta = (Cx.MoleFractions(i,j,k)({pIndex, comp}) - Cref({pIndex, comp}));
                DeformationGradientsEigen(i,j,k) += Lambda({alpha->index, comp})*delta*alpha->value;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::CalculateThermalExpansion(PhaseField& Phase, Temperature& Tx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DeformationGradientsEigen,0,)
    {
        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        {
            int pIndex = Phase.FieldsStatistics[alpha->index].Phase;
            double delta = (Tx(i,j,k) - Tref[pIndex]);
            DeformationGradientsEigen(i,j,k) += Alpha[alpha->index]*delta*alpha->value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::CalculateInterfaceStress(PhaseField& Phase, InterfaceProperties& IP)
{
    //Raphael Schiedung, Ingo Steinbach, and Fathollah Varnik.
    //"Multi-phase-field method for surface tension induced elasticity."
    //Physical Review B 97.3 (2018): 035410.
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DeformationGradientsEigen,0,)
    if (Phase.Interface(i,j,k))
    {
        const double pre1 = 4.0/Phase.Eta;
        const double pre2 = Phase.Eta * Phase.Eta/Pi;

        vStress locInterfaceStress;
        // Calculate phase gradients and normals
        const NodeV3 locNormals   = Phase.Normals(i,j,k);
        const NodeV3 locGradients = Phase.Fields(i,j,k).get_gradients();

        for (auto it = locNormals.cbegin(); it < locNormals.cend(); ++it)
        {
            // Calculate phase field part of the interface stress
            const double PhiA = abs(Phase.Fields(i,j,k).get_value(it->indexA));
            const double PhiB = abs(Phase.Fields(i,j,k).get_value(it->indexB));

            if (PhiA * PhiB < 0.02) continue; 

            const double interpol = pre1 * (pre2 *
                    (locGradients.get(it->indexA) *
                     locGradients.get(it->indexB)) + PhiA * PhiB);

            // Calculate local projection matrix
            dMatrix3x3 Projection;
            Projection.set_to_unity();
            Projection(0,0) -= it->X()*it->X();
            Projection(1,0) -= it->Y()*it->X();
            Projection(2,0) -= it->Z()*it->X();
            Projection(0,1) -= it->X()*it->Y();
            Projection(1,1) -= it->Y()*it->Y();
            Projection(2,1) -= it->Z()*it->Y();
            Projection(0,2) -= it->X()*it->Z();
            Projection(1,2) -= it->Y()*it->Z();
            Projection(2,2) -= it->Z()*it->Z();

            // Calculate local interface stress
            locInterfaceStress -= VoigtStress(Projection) * interpol *
                    IP.get_energy(i,j,k,it->indexA,it->indexB);
        }

        //TODO account for derivatives Interface energy with respect to the strain
        vStrain locEigenStrainInterface =
            EffectiveElasticConstants(i,j,k).inverted() * locInterfaceStress *(-1.);

        DeformationGradientsEigen(i,j,k) += locEigenStrainInterface.tensor();

    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::CalculateAverageElasticConstants(void)
{
    dMatrix6x6 locAverageElasticConstants;
    AverageElasticConstants.set_to_zero();

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EffectiveElasticConstants, 0, reduction(dMatrix6x6SUM: locAverageElasticConstants))
    {
        locAverageElasticConstants += EffectiveElasticConstants(i, j, k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    AverageElasticConstants = locAverageElasticConstants;

#ifdef MPI_PARALLEL
    MPI_Allreduce(locAverageElasticConstants.data(),AverageElasticConstants.data(),36,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif

    AverageElasticConstants /= double(TotalNx*TotalNy*TotalNz);
}

void ElasticProperties::SetEffectiveElasticConstants(PhaseField& Phase)
{
    switch(EModel)
    {
        case ElasticityModels::Khachaturyan:
        {
            ElasticityKhachaturyan::SetEffectiveElasticConstants(Phase, *this);
            break;
        }
        case ElasticityModels::Steinbach:
        {
            ElasticitySteinbach::SetEffectiveElasticConstants(Phase, *this);
            break;
        }
        default:
        {
            Info::WriteExit("Nonexistent elasticity model selected", thisclassname, "SetEffectiveElasticConstants()");
            exit(1);
        }
    }
    CalculateAverageElasticConstants();
}

void ElasticProperties::SetEffectiveElasticConstants(PhaseField& Phase, Composition& Cx)
{
    switch(EModel)
    {
        case ElasticityModels::Khachaturyan:
        {
            ElasticityKhachaturyan::SetEffectiveElasticConstants(Phase, *this, Cx);
            break;
        }
        case ElasticityModels::Steinbach:
        {
            ElasticitySteinbach::SetEffectiveElasticConstants(Phase, *this, Cx);
            break;
        }
        default:
        {
            Info::WriteExit("Nonexistent elasticity model selected", thisclassname, "SetEffectiveElasticConstants()");
            exit(1);
        }
    }
    CalculateAverageElasticConstants();
}

void ElasticProperties::SetEffectiveElasticConstants(PhaseField& Phase, Temperature& Tx)
{
    switch(EModel)
    {
        case ElasticityModels::Khachaturyan:
        {
            ElasticityKhachaturyan::SetEffectiveElasticConstants(Phase, *this, Tx);
            break;
        }
        case ElasticityModels::Steinbach:
        {
            ElasticitySteinbach::SetEffectiveElasticConstants(Phase, *this, Tx);
            break;
        }
        default:
        {
            Info::WriteExit("Nonexistent elasticity model selected", thisclassname, "SetEffectiveElasticConstants()");
            exit(1);
        }
    }
    CalculateAverageElasticConstants();
}

void ElasticProperties::CalculateDrivingForce(PhaseField& Phase, Composition& Cx, DrivingForce& dGab) const
{
    switch(EModel)
    {
        case ElasticityModels::Khachaturyan:
        {
            ElasticityKhachaturyan::CalculateDrivingForce(Phase, *this, Cx, dGab);
            break;
        }
        case ElasticityModels::Steinbach:
        {
            ElasticitySteinbach::CalculateDrivingForce(Phase, *this, Cx, dGab);
            break;
        }
        default:
        {
            Info::WriteExit("Nonexistent elasticity model selected", thisclassname, "CalculateDrivingForce()");
            exit(1);
        }
    }
}

void ElasticProperties::CalculateDrivingForce(PhaseField& Phase, Temperature& Tx, DrivingForce& dGab) const
{
    switch(EModel)
    {
        case ElasticityModels::Khachaturyan:
        {
            ElasticityKhachaturyan::CalculateDrivingForce(Phase, *this, Tx, dGab);
            break;
        }
        case ElasticityModels::Steinbach:
        {
            ElasticitySteinbach::CalculateDrivingForce(Phase, *this, Tx, dGab);
            break;
        }
        default:
        {
            Info::WriteExit("Nonexistent elasticity model selected", thisclassname, "CalculateDrivingForce()");
            exit(1);
        }
    }
}

void ElasticProperties::CalculateDrivingForce(PhaseField& Phase, DrivingForce& dGab) const
{
    switch(EModel)
    {
        case ElasticityModels::Khachaturyan:
        {
            ElasticityKhachaturyan::CalculateDrivingForce(Phase, *this, dGab);
            break;
        }
        case ElasticityModels::Steinbach:
        {
            ElasticitySteinbach::CalculateDrivingForce(Phase, *this, dGab);
            break;
        }
        default:
        {
            Info::WriteExit("Nonexistent elasticity model selected", thisclassname, "CalculateDrivingForce()");
            exit(1);
        }
    }
}

vStrain ElasticProperties::CalculateNeuberCorrection(
        const vStrain ElasticStrains, const dMatrix6x6 locCompliance,
        const int i, const int j, const int k,
        const size_t pIndexA, const size_t pIndexB) const
{
    vStrain ElasticStrainsCorr = ElasticStrains;
    if (NeuberCorrection[pIndexA].Active or NeuberCorrection[pIndexB].Active)
    {
        vStress ElasticStresses;
        vStress deviatoricStresses;
        vStrain deviatoricStrains;
        vStrain dummyStrain;
        double StressTrace = Stresses(i, j, k).trace();
        double StrainTrace = ElasticStrains.trace();
        for(int n = 0; n < 3; n++)
        {
            deviatoricStresses[n] = Stresses(i, j, k)[n] - (1.0/3.0)*StressTrace;
            deviatoricStrains[n] = ElasticStrains[n] - (1.0/3.0)*StrainTrace;
        }
        if(NeuberCorrection[pIndexA].Active and !NeuberCorrection[pIndexB].Active)
        {
            PlasticFlowNeuberMethods::getNeuberDataRO(
                    deviatoricStresses,
                    deviatoricStrains,
                    ElasticStresses,
                    dummyStrain,
                    NeuberCorrection[pIndexA].YoungsModulus,
                    NeuberCorrection[pIndexA].YieldStrength,
                    NeuberCorrection[pIndexA].Hardening);
        }
        if(!NeuberCorrection[pIndexA].Active and NeuberCorrection[pIndexB].Active)
        {
            PlasticFlowNeuberMethods::getNeuberDataRO(
                    deviatoricStresses,
                    deviatoricStrains,
                    ElasticStresses,
                    dummyStrain,
                    NeuberCorrection[pIndexB].YoungsModulus,
                    NeuberCorrection[pIndexB].YieldStrength,
                    NeuberCorrection[pIndexB].Hardening);
        }
        if(NeuberCorrection[pIndexA].Active and NeuberCorrection[pIndexB].Active)
        {
            if(NeuberCorrection[pIndexA].YieldStrength <= NeuberCorrection[pIndexB].YieldStrength)
            {
                PlasticFlowNeuberMethods::getNeuberDataRO(
                        deviatoricStresses,
                        deviatoricStrains,
                        ElasticStresses,
                        dummyStrain,
                        NeuberCorrection[pIndexA].YoungsModulus,
                        NeuberCorrection[pIndexA].YieldStrength,
                        NeuberCorrection[pIndexA].Hardening);
            }
            else
            {
                PlasticFlowNeuberMethods::getNeuberDataRO(
                        deviatoricStresses,
                        deviatoricStrains,
                        ElasticStresses,
                        dummyStrain,
                        NeuberCorrection[pIndexB].YoungsModulus,
                        NeuberCorrection[pIndexB].YieldStrength,
                        NeuberCorrection[pIndexB].Hardening);
            }
        }
        for(int n = 0; n < 3; n++)
        {
            ElasticStresses[n] += (1.0/3.0)*StressTrace;
        }
        ElasticStrainsCorr = locCompliance*ElasticStresses;
    }
    return ElasticStrainsCorr;
}

void ElasticProperties::CalculateDeformationJumps(const PhaseField& Phase)
{
    NodeV3 Jumps;
    double residual = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DeformationJumps,0,)
    {
        if(Phase.Interface(i,j,k))
        {
            do
            {
                NodeV3 locJumps;
                residual = 0.0;
                for(auto alpha = Phase.Fields(i, j, k).cbegin();
                         alpha != Phase.Fields(i, j, k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta != Phase.Fields(i, j, k).cend();  ++beta)
                if(alpha->value != 0.0 and beta->value != 0.0)
                {
                    dVector3 locNormalAB = Phase.Normals(i, j, k).get_asym(alpha->index, beta->index);

                    vStrain dStrainA = (TotalStrains(i,j,k) - StrainSmall(TransformationStretches[alpha->index]));
                    vStrain dStrainB = (TotalStrains(i,j,k) - StrainSmall(TransformationStretches[ beta->index]));

                    if(Phase.Fields(i, j, k).size() > 2)
                    for(auto gamma = Phase.Fields(i, j, k).cbegin();
                             gamma != Phase.Fields(i, j, k).cend(); ++gamma)
                    if(gamma->value != 0.0 and gamma != alpha and gamma != beta)
                    {
                        dVector3 locNormalAG = Phase.Normals(i, j, k).get_asym(alpha->index, gamma->index);
                        dVector3 locNormalBG = Phase.Normals(i, j, k).get_asym( beta->index, gamma->index);

                        dVector3 locJumpA = Jumps.get_sym(alpha->index, gamma->index);
                        dVector3 locJumpB = Jumps.get_sym( beta->index, gamma->index);

                        dMatrix3x3 locFjumpA = locJumpA.dyadic(locNormalAG);
                        dMatrix3x3 locFjumpB = locJumpB.dyadic(locNormalBG);

                        /*if(alpha->index + beta->index == 3)
                        {
                            dStrainA -= (locFjumpA + locFjumpA.transposed()).VoigtStrain()*gamma->value*0.5;
                            dStrainB += (locFjumpB + locFjumpB.transposed()).VoigtStrain()*gamma->value*0.5;
                        }
                        else*/
                        {
                            dStrainA -= VoigtStrain(locFjumpA + locFjumpA.transposed())*gamma->value*0.5;
                            dStrainB -= VoigtStrain(locFjumpB + locFjumpB.transposed())*gamma->value*0.5;
                        }
                    }

                    dMatrix6x6 Cij = ElasticConstants[alpha->index]*(1.0 - alpha->value) +
                                     ElasticConstants[ beta->index]*(1.0 - beta->value);

                    dVector3 locJump;
                    dVector3 locJumpOLD = Jumps.get_sym(alpha->index, beta->index);
                    if(locNormalAB.abs() > DBL_EPSILON)
                    {
                        dMatrix3x3 projCij = Cij.project(locNormalAB);
                        locJump = projCij.inverted()*(ElasticConstants[ beta->index]*dStrainB - ElasticConstants[alpha->index]*dStrainA).tensor()*locNormalAB*(-1.0);
                    }
                    locJumps.set_sym(alpha->index, beta->index, locJump);
                    residual = max((locJumpOLD - locJump).abs()*beta->value, residual);
                    if(Phase.Fields(i, j, k).size() > 2)
                    {
                        cout << alpha->index << " (" << alpha->value << ") " << " " << beta->index << " (" << beta->value << "); residual " << " (" << residual<< ") " << ": " << locJumps.get_sym(alpha->index, beta->index).print() << endl;
                        getchar();
                    }
                }
                DeformationJumps(i,j,k) = locJumps;
            }
            while(residual > 1.0e-6);
        }
        else
        {
            DeformationJumps(i,j,k).clear();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::Advect(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi, const BoundaryConditions& BC, const double dt, const double tStep)
{
    Adv.AdvectField(DeformationGradientsTotal, Vel, BC, dx, dt, tStep);
    Adv.AdvectField(DeformationGradientsEigen, Vel, BC, dx, dt, tStep);
    Adv.AdvectField(Stresses, Vel, BC, dx, dt, tStep, true);
}

void ElasticProperties::GetAverageDeformationGradient()
{
    AverageDeformationGradient = AverageStrain.tensor() + dMatrix3x3::UnitTensor();
}

void ElasticProperties::WriteTotalRotationsVTK(const int tStep, const Settings& locSettings, const int precision) const
{
    auto CalculateRotations = [&](int i, int j, int k)
    {
        double locAngle = 0.0;
        dVector3 locAxis;
        Tools::getAxisAngle(DeformationGradientsTotal(i,j,k),locAxis,locAngle);
        locAngle *= 180.0 / Pi;
        return std::make_pair(locAngle, locAxis);
    };

    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"Angle", [=](int i,int j,int k){return CalculateRotations(i,j,k).first;}});
    ListOfFields.push_back((VTK::Field_t) {"Axis",  [=](int i,int j,int k){return CalculateRotations(i,j,k).second;}});

    std::string Filename = UserInterface::MakeFileName(VTKDir, "Rotations_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields);
}

void ElasticProperties::Write(int tStep) const
{
#ifdef MPI_PARALLEL
    string FileName = UserInterface::MakeFileName(RawDataDir, thisclassname + "_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = UserInterface::MakeFileName(RawDataDir, thisclassname + "_", tStep, ".dat");
#endif
    fstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        Info::WriteExit("File " + FileName + " could not be created",thisclassname);
        exit(1);
    };

    out.write(reinterpret_cast<const char*>(&(AverageStrain)), sizeof(vStrain));

    STORAGE_LOOP_BEGIN(i,j,k,DeformationGradientsTotal,DeformationGradientsTotal.Bcells())
    {
        out.write(reinterpret_cast<const char*>(&(DeformationGradientsTotal(i,j,k))), sizeof(dMatrix3x3));
        out.write(reinterpret_cast<const char*>(&(DeformationGradientsEigen(i,j,k))), sizeof(dMatrix3x3));
    }
    STORAGE_LOOP_END

    out.close();
}

void ElasticProperties::WriteDeformationGradientsTotal(int tStep) const
{
#ifdef MPI_PARALLEL
    string FileName = UserInterface::MakeFileName(RawDataDir, "DeformationGradientsTotal_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = UserInterface::MakeFileName(RawDataDir, "DeformationGradientsTotal_", tStep, ".dat");
#endif
    fstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        Info::WriteExit("File " + FileName + " could not be created",thisclassname);
        exit(1);
    };

    STORAGE_LOOP_BEGIN(i,j,k,DeformationGradientsTotal,DeformationGradientsTotal.Bcells())
    {
        out.write(reinterpret_cast<const char*>(&(DeformationGradientsTotal(i,j,k))), sizeof(dMatrix3x3));
    }
    STORAGE_LOOP_END

    out.close();
}

void ElasticProperties::WriteDeformationGradientsEigen(const int tStep) const
{
#ifdef MPI_PARALLEL
    string FileName = UserInterface::MakeFileName(RawDataDir, "DeformationGradientsEigen_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = UserInterface::MakeFileName(RawDataDir, "DeformationGradientsEigen_", tStep, ".dat");
#endif

    fstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        Info::WriteExit("File " + FileName + " could not be created", thisclassname);
        exit(1);
    };

    STORAGE_LOOP_BEGIN(i,j,k,DeformationGradientsEigen,DeformationGradientsEigen.Bcells())
    {
        out.write(reinterpret_cast<const char*>(&(DeformationGradientsEigen(i,j,k))), sizeof(dMatrix3x3));
    }
    STORAGE_LOOP_END

    out.close();
}

void ElasticProperties::WriteStresses(int tStep) const
{
#ifdef MPI_PARALLEL
    string FileName = UserInterface::MakeFileName(RawDataDir, "Stresses_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = UserInterface::MakeFileName(RawDataDir, "Stresses_", tStep, ".dat");
#endif
    fstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        Info::WriteExit("File " + FileName + " could not be created",thisclassname);
        exit(1);
    };

    STORAGE_LOOP_BEGIN(i,j,k,Stresses,Stresses.Bcells())
    {
        out.write(reinterpret_cast<const char*>(&(Stresses(i,j,k))), sizeof(vStress));
    }
    STORAGE_LOOP_END
    out.close();
}

void ElasticProperties::Read(const BoundaryConditions& BC, const int tStep)
{
#ifdef MPI_PARALLEL
    string FileName = UserInterface::MakeFileName(RawDataDir, thisclassname + "_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = UserInterface::MakeFileName(RawDataDir, thisclassname + "_", tStep, ".dat");
#endif

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteExit("File " + FileName + " could not be opened",thisclassname);
        exit(1);
    };

    inp.read(reinterpret_cast<char*>(&(AverageStrain)), sizeof(vStrain));

    STORAGE_LOOP_BEGIN(i,j,k,DeformationGradientsTotal,DeformationGradientsTotal.Bcells())
    {
        inp.read(reinterpret_cast<char*>(&(DeformationGradientsTotal(i,j,k))), sizeof(dMatrix3x3));
        inp.read(reinterpret_cast<char*>(&(DeformationGradientsEigen(i,j,k))), sizeof(dMatrix3x3));
    }
    STORAGE_LOOP_END
    SetBoundaryConditions(BC);
    Info::Write(thisclassname, "Binary input loaded");
}

void ElasticProperties::ReadDeformationGradientsTotal(const int tStep)
{
#ifdef MPI_PARALLEL
    string FileName = UserInterface::MakeFileName(RawDataDir, "DeformationGradientsTotal_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = UserInterface::MakeFileName(RawDataDir, "DeformationGradientsTotal_", tStep, ".dat");
#endif

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteExit("File " + FileName + " could not be opened",thisclassname);
        exit(1);
    };

    STORAGE_LOOP_BEGIN(i,j,k,DeformationGradientsTotal,DeformationGradientsTotal.Bcells())
    {
        inp.read(reinterpret_cast<char*>(&(DeformationGradientsTotal(i,j,k))), sizeof(dMatrix3x3));
    }
    STORAGE_LOOP_END
    Info::Write("TotalStrains", "Binary input loaded");
}

void ElasticProperties::ReadDeformationGradientsEigen(const int tStep)
{
#ifdef MPI_PARALLEL
    string FileName = UserInterface::MakeFileName(RawDataDir, "DeformationGradientsEigen_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = UserInterface::MakeFileName(RawDataDir, "DeformationGradientsEigen_", tStep, ".dat");
#endif

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteExit("File " + FileName + " could not be opened",thisclassname);
        exit(1);
    };

    STORAGE_LOOP_BEGIN(i,j,k,DeformationGradientsEigen,DeformationGradientsEigen.Bcells())
    {
        inp.read(reinterpret_cast<char*>(&(DeformationGradientsEigen(i,j,k))), sizeof(dMatrix3x3));
    }
    STORAGE_LOOP_END
    Info::Write("TotalStrains", "Binary input read successfully");
}

void ElasticProperties::ReadStresses(int tStep)
{
#ifdef MPI_PARALLEL
    string FileName = UserInterface::MakeFileName(RawDataDir, "Stresses_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = UserInterface::MakeFileName(RawDataDir, "Stresses_", tStep, ".dat");
#endif
    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteExit("File " + FileName + " could not be opened",thisclassname);
        exit(1);
    };

    STORAGE_LOOP_BEGIN(i,j,k,Stresses,Stresses.Bcells())
    {
        inp.read(reinterpret_cast<char*>(&(Stresses(i,j,k))), sizeof(vStress));
    }
    STORAGE_LOOP_END
}

void ElasticProperties::WriteTotalStrainsVTK(const int tStep,
        const Settings& locSettings, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"TotalStrain", [this](int i,int j,int k){return TotalStrains(i, j, k);}});

    std::string Filename = UserInterface::MakeFileName(VTKDir, "TotalStrain_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void ElasticProperties::WriteStressesVTK(const int tStep,
        const Settings& locSettings, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"Stresses", [this](int i,int j,int k){return Stresses(i, j, k);}});
    ListOfFields.push_back((VTK::Field_t) {"von Mises", [this](int i,int j,int k){return Stresses(i, j, k).Mises();}});
    ListOfFields.push_back((VTK::Field_t) {"Pressure", [this](int i,int j,int k){return Stresses(i, j, k).Pressure();}});

    std::string Filename = UserInterface::MakeFileName(VTKDir, "Stresses_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void ElasticProperties::WriteElasticStrainsVTK(const int tStep,
        const Settings& locSettings, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"ElasticStrains", [this](int i,int j,int k){return ElasticStrains(i,j,k);}});

    std::string Filename = UserInterface::MakeFileName(VTKDir, "ElasticStrains_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void ElasticProperties::WriteEigenStrainsVTK(const int tStep,
        const Settings& locSettings, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"EigenStrains", [this](int i,int j,int k){return EigenStrains(i,j,k);}});

    std::string Filename = UserInterface::MakeFileName(VTKDir, "EigenStrains_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void ElasticProperties::WriteDeformationGradientsTotalVTK(const int tStep,
        const Settings& locSettings, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"DeformationGradientsTotal", [this](int i,int j,int k){return DeformationGradientsTotal(i,j,k);}});
    std::string Filename = UserInterface::MakeFileName(VTKDir, "DeformationGradientsTotal_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void ElasticProperties::WriteEffectiveElasticConstantsVTK(const int tStep,
        const Settings& locSettings, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"EffectiveElasticConstants", [this](int i,int j,int k){return EffectiveElasticConstants(i, j, k);}});

    std::string Filename = UserInterface::MakeFileName(VTKDir, "EffectiveElasticConstants_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void ElasticProperties::WriteForceDensityVTK(const int tStep,
        const Settings& locSettings, const int precision) const
{
    if(ConsiderExternalForces)
    {
        std::vector<VTK::Field_t> ListOfFields;
        ListOfFields.push_back((VTK::Field_t) {"ForceDensity", [this](int i,int j,int k){return ForceDensity(i, j, k);}});

        std::string Filename = UserInterface::MakeFileName(VTKDir, "ForceDensity_", tStep, ".vts");
        VTK::Write(Filename, locSettings, ListOfFields, precision);
    }
    else
    {
        Info::WriteExit("ForceDensity Storage not allocated!", "ElasticProperties", "WriteForceDensityVTK");
        exit(EXIT_FAILURE);
    }
}

void ElasticProperties::WriteCauchyStressesVTK(const int tStep,
        const Settings& locSettings, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"CauchyStresses", [this](int i,int j,int k){return CauchyStress(i,j,k);}});
    ListOfFields.push_back((VTK::Field_t) {"Pressure", [this](int i,int j,int k){return CauchyStress(i,j,k).Pressure();}});
    ListOfFields.push_back((VTK::Field_t) {"vonMises", [this](int i,int j,int k){return CauchyStress(i,j,k).Mises();}});

    std::string Filename = UserInterface::MakeFileName(VTKDir, "CauchyStresses_", tStep, ".vts");
    VTK::WriteDistorted(Filename, locSettings, *this, ListOfFields, precision);
}

void ElasticProperties::WriteDisplacementsVTK(const int tStep, const Settings& OPSettings,
                                 const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"U", [this](int i, int j, int k){return Displacements(i,j,k);}});

    std::string Filename = UserInterface::MakeFileName(OPSettings.VTKDir, "Displacements_", tStep, ".vts");
    VTK::Write(Filename, OPSettings, ListOfFields, precision);
}

void ElasticProperties::PrintPointStatistics(const int x, const int y, const int z) const
{
    Info::WriteCoordinate(x, y, z, dx);
    stringstream outstream;
    outstream << "Stresses    : " <<     Stresses(x, y, z).print() << endl;
    outstream << "TotalStrains: " << TotalStrains(x, y, z).print() << endl;
    outstream << "EigenStrains: " << EigenStrains(x, y, z).print() << endl;
    outstream << "ElasticConstants: " << endl
              << EffectiveElasticConstants(x, y, z).print() << endl;

    Info::WriteSimple(outstream.str());
}

vStress ElasticProperties::WriteStressStrainData(std::string filename, std::string LDflag) const
{
    vStress aveStress;
    aveStress.set_to_zero();
    vStrain aveStrain;
    aveStrain.set_to_zero();
    STORAGE_LOOP_BEGIN(i,j,k,Stresses,0)
    {
        aveStress += Stresses(i,j,k);
        aveStrain += TotalStrains(i,j,k);
    }
    STORAGE_LOOP_END
    aveStress /= (Nx*Ny*Nz);
    aveStrain /= (Nx*Ny*Nz);

    ofstream outputFile;
    outputFile.open(filename, ios::app|ios::out);
    outputFile << aveStrain[0]  << ", " << aveStress[0] << ", " <<
                  aveStrain[1]  << ", " << aveStress[1] << ", " <<
                  aveStrain[2]  << ", " << aveStress[2] << ", " <<
                  aveStrain[3]  << ", " << aveStress[3] << ", " <<
                  aveStrain[4]  << ", " << aveStress[4] << ", " <<
                  aveStrain[5]  << ", " << aveStress[5] << ", " <<
                                               aveStress.Mises() << endl;
    outputFile.close();
    return aveStress;
}

double ElasticProperties::Energy(void) const
{
    double Energy = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DeformationGradientsTotal, 0, reduction(+:Energy))
    {
        Energy += EnergyDensity(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    #ifdef MPI_PARALLEL
    double sEnergy = Energy;
    MPI_Allreduce(&sEnergy,&Energy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    #endif

    return Energy*dx*dx*dx;
}

double ElasticProperties::AverageEnergyDensity(void) const
{
    double Energy = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DeformationGradientsTotal, 0, reduction(+:Energy))
    {
        Energy += EnergyDensity(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    #ifdef MPI_PARALLEL
    double sEnergy = Energy;
    MPI_Allreduce(&sEnergy,&Energy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    #endif

    return Energy/double(TotalNx*TotalNy*TotalNz);
}

void ElasticProperties::CalculateInterfaceEnergyContribution(PhaseField& Phase, InterfaceProperties& IP) const
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DeformationGradientsTotal, 0,)
    {
        if(Phase.Fields(i,j,k).flag)
        {
            for(auto alpha  = Phase.Fields(i, j, k).cbegin();
                     alpha != Phase.Fields(i, j, k).cend() - 1; ++alpha)
            for(auto  beta  = alpha + 1;
                      beta != Phase.Fields(i, j, k).cend(); ++beta)
            {
                IP.add_energy(i,j,k,alpha->index, beta->index, EnergyDensity(i,j,k)*Phase.Eta);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

double ElasticProperties::EnergyDensity(int i, int j, int k) const
{
    vStrain locElasticStrains = ElasticStrains(i,j,k);
    return 0.5*(locElasticStrains*(EffectiveElasticConstants(i,j,k)*locElasticStrains));
}

void ElasticProperties::WriteEnergyDensityVTK(const int tStep, const Settings& locSettings, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"Elastic Energy Density", [this] (int i, int j, int k){return EnergyDensity(i,j,k);}});
    std::string Filename = UserInterface::MakeFileName(locSettings.VTKDir, "ElasticEnergyDensity_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
}// namespace openphase
