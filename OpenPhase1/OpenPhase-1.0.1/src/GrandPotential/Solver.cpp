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


#include "Base/FindRoot.h"
#include "Base/SolveLinearSystem.h"
#include "Base/SparseMatrix.h"
#include "Base/UserInterface.h"
#include "BoundaryConditions.h"
#include "GrandPotential/Density.h"
#include "GrandPotential/Solver.h"
#include "Info.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "PhysicalConstants.h"
#include "Settings.h"
#include "VTK.h"

namespace openphase::GrandPotential
{
Solver::Solver(const Settings& locSettings)
{
    Initialize(locSettings);
    ReadInput(DefaultInputFileName);
}
Solver::Solver(const Settings& locSettings, std::string filename)
{
    Initialize(locSettings);
    ReadInput(filename);
}
void Solver::Initialize(const Settings& locSettings)
{
    TotalNx = locSettings.TotalNx;

    Nx  = locSettings.Nx;
    Ny  = locSettings.Ny;
    Nz  = locSettings.Nz;

    dNx = locSettings.dNx;
    dNy = locSettings.dNy;
    dNz = locSettings.dNz;

    ElementNames = locSettings.ElementNames;
    Ncomp        = locSettings.Ncomp;
    Nphases      = locSettings.Nphases;
    PhaseNames   = locSettings.PhaseNames;
    RawDataDir   = locSettings.RawDataDir;
    TextDir      = locSettings.TextDir;
    VTKDir       = locSettings.VTKDir;
    dx           = locSettings.dx;

    CInitial           .Allocate({Nphases, Ncomp});
    PhaseMobilities    .Allocate({Nphases, Ncomp});
    InterfaceMobilities.Allocate({Nphases, Nphases, Ncomp});

    dPhaseMobilities_dConcentration.Allocate({Nphases, Ncomp, Ncomp});

    InitialChemicalPotential.assign(Ncomp, 0.0);
    ElementMasses           .assign(Ncomp, 0.0);
    TOC0                    .assign(Ncomp, 0.0);

    size_t Bcells = locSettings.Bcells;
    if (Bcells < 2)
    {
        Info::WriteExit("Too few boundary cell! Two are required!", thisclassname, "ReadInput");
        std::exit(EXIT_FAILURE);
    };

    Concentrations      .Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Ncomp}, Bcells);
    ChemicalPotential   .Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Ncomp}, Bcells);
    ChemicalPotentialDot.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Ncomp}, Bcells);
    DiffusionFlux       .Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Ncomp}, Bcells);
    Mobilities          .Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Ncomp}, Bcells);

    Info::WriteStandard(thisclassname, "Initialized");
}
void Solver::ReadInput(std::string filename)
{
    std::fstream inp(filename.c_str(), std::ios::in);
    if (!inp)
    {
        Info::WriteExit("File \"" + filename + "\" could not be opened", thisclassname, "ReadInput");
        std::exit(EXIT_FAILURE);
    };
    std::stringstream inp_data;
    inp_data << inp.rdbuf();
    inp.close();

    Info::WriteLine();
    Info::WriteLineInsert("Grand Potential Solver");
    Info::WriteStandard("Source", filename);

    int moduleLocation = UserInterface::FindModuleLocation(inp_data, thisclassname);
    // Read component masses
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        std::string converter = ""+ElementNames[comp];
        ElementMasses            [comp] = UserInterface::ReadParameterD(inp_data, moduleLocation, "MASS_"+converter, false, 0.0);
        InitialChemicalPotential [comp] = UserInterface::ReadParameterD(inp_data, moduleLocation, "ICP_"+converter, false, 0.0);
    }

    // Read initial concentrations
    UseInitialPressure = UserInterface::ReadParameterB(inp_data, moduleLocation, "UIP",false,false);
    if (UseInitialPressure)
    {
        InitialPressure = UserInterface::ReadParameterD(inp_data, moduleLocation, "IP");
    }
    else
    {
        for(size_t PhaseIdx = 0; PhaseIdx < Nphases; PhaseIdx++)
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            std::string converter = ""+PhaseNames[PhaseIdx]+"_"+ElementNames[comp];
            CInitial ({PhaseIdx, comp}) = UserInterface::ReadParameterD(inp_data, moduleLocation, "CI_"+converter);
        }
    }

    // Read bulk mobilities
    for(size_t PhaseIdx = 0; PhaseIdx < Nphases; PhaseIdx++)
    for(size_t comp     = 0; comp     < Ncomp;       comp++)
    {
        std::string converter = ""+PhaseNames[PhaseIdx]+"_"+ElementNames[comp];
        PhaseMobilities({PhaseIdx, comp}) = UserInterface::ReadParameterD(inp_data, moduleLocation, "M0_"+converter);
    }

    // Read bulk mobilities concentration dependency
    bool PhaseMobilityConcentrationCoupling = UserInterface::ReadParameterB(inp_data, moduleLocation, "dMdc");
    if (PhaseMobilityConcentrationCoupling)
    for(size_t PhaseIdx = 0; PhaseIdx < Nphases; PhaseIdx++)
    for(size_t compA    = 0; compA     < Ncomp;     compA++)
    for(size_t compB    = 0; compB     < Ncomp;     compB++)
    {
        std::string converter = ""+PhaseNames[PhaseIdx]+"_"+ElementNames[compA]+"_"+ElementNames[compB];
        double valueAB = UserInterface::ReadParameterD(inp_data, moduleLocation, "dMdc_"+converter);
        dPhaseMobilities_dConcentration({PhaseIdx, compA, compB}) = valueAB;
        dPhaseMobilities_dConcentration({PhaseIdx, compB, compA}) = valueAB;
    }

    // Read interface mobilities
    bool InterfaceMobilityConstants = UserInterface::ReadParameterB(inp_data, moduleLocation, "IM");
    if (InterfaceMobilityConstants)
    for(size_t alpha = 0;     alpha < Nphases; alpha++)
    for(size_t beta  = alpha; beta  < Nphases;  beta++)
    for(size_t comp  = 0;     comp  < Ncomp;    comp++)
    {
        std::string converterAB = PhaseNames[alpha] + "_" + PhaseNames[beta ] + "_" + ElementNames[comp];
        //std::string converterBA = PhaseNames[beta ] + "_" + PhaseNames[alpha] + "_" + ElementNames[comp];
        double valueAB = UserInterface::ReadParameterD(inp_data, moduleLocation, "IM_" + converterAB);
        InterfaceMobilities({alpha, beta, comp}) = valueAB;
        InterfaceMobilities({beta, alpha, comp}) = valueAB;
    }

    UseImplicitSolver = UserInterface::ReadParameterB(inp_data, moduleLocation, "Implicit");
    if (UseImplicitSolver)
    {
        ChemicalPotentialAccuracy = UserInterface::ReadParameterD(inp_data, moduleLocation, "ACC");
        MaxIterations             = UserInterface::ReadParameterI(inp_data, moduleLocation, "MAXI");

        ChemicalPotentialOld .Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Ncomp}, boundary);
        ChemicalPotentialDot2.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, {Ncomp}, boundary);
    }
    else
    {
        ChemicalPotentialAccuracy = 0.0;
        MaxIterations             = 0;
    }

    ConserveTOC = UserInterface::ReadParameterB(inp_data, moduleLocation, "TOC");
    if (ConserveTOC)
    {
        TOCAccuracy      = UserInterface::ReadParameterD(inp_data, moduleLocation, "TOCACC");
        TOCMaxIterations = UserInterface::ReadParameterI(inp_data, moduleLocation, "TOCMAXI");
    }
    else
    {
        TOCAccuracy      = 0.0;
        TOCMaxIterations = 0;
    }

    Info::WriteLine();
}
void Solver::Write(const char* filename) const
{
    std::ofstream out;
#ifdef MPI_PARALLEL
    for (int proc = 0; proc < MPI_SIZE; proc++)
    {
        if (MPI_RANK == proc)
        {
            if (MPI_RANK == 0)
            {
                out.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
            }
            else
            {
                out.open(filename, std::ios::out | std::ios::binary | std::ios::app);
            }
            if (MPI_RANK == 0)
            {
                out.write(reinterpret_cast<const char*>(&TotalNx), sizeof(long int));
                out.write(reinterpret_cast<const char*>(&Ny     ), sizeof(long int));
                out.write(reinterpret_cast<const char*>(&Nz     ), sizeof(long int));
                out.write(reinterpret_cast<const char*>(&Ncomp  ), sizeof(size_t));
                for(size_t comp = 0; comp < Ncomp; comp++)
                {
                    out.write(reinterpret_cast<const char*>(&TOC0[comp]), sizeof(size_t));
                }
            }
            for (long int i = 0; i < Nx;    i++)
            for (long int j = 0; j < Ny;    j++)
            for (long int k = 0; k < Nz;    k++)
            for (size_t   n = 0; n < Ncomp; n++)
            {
               out.write(reinterpret_cast<const char*>(&ChemicalPotential(i,j,k)({n})), sizeof(double));
            }
            out.close();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#else
    out.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
    out.write(reinterpret_cast<const char*>(&TotalNx), sizeof(long int));
    out.write(reinterpret_cast<const char*>(&Ny     ), sizeof(long int));
    out.write(reinterpret_cast<const char*>(&Nz     ), sizeof(long int));
    out.write(reinterpret_cast<const char*>(&Ncomp  ), sizeof(size_t));
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        out.write(reinterpret_cast<const char*>(&TOC0[comp]), sizeof(size_t));
    }
    for (long int i = 0; i < Nx;    i++)
    for (long int j = 0; j < Ny;    j++)
    for (long int k = 0; k < Nz;    k++)
    for (size_t   n = 0; n < Ncomp; n++)
    {
       out.write(reinterpret_cast<const char*>(&ChemicalPotential(i,j,k)({n})), sizeof(double));
    }
    out.close();
#endif
}
void Solver::Write(const std::string& filename) const
{
    Write(filename.c_str());
}
void Solver::Write(long tStep) const
{
    std::string filename = UserInterface::MakeFileName(RawDataDir, thisclassname+'_', tStep, ".dat");
    Write(filename);
}
void Solver::Write() const
{
    std::string filename = RawDataDir+thisclassname+".dat";
    Write(filename);
}
void Solver::Read(const BoundaryConditions& BC, const PhaseField& Phase, const Density& omega, const char* filename)
{
    std::ifstream inp;
#ifdef MPI_PARALLEL
    long pos = 0;
    for (int proc = 0; proc < MPI_SIZE; proc++)
    {
        if (MPI_RANK == proc)
        {
            inp.open(filename, std::ios::in | std::ios::binary);
            if (MPI_RANK == 0)
            {
                long int locTotalNx = 0;
                long int locNy      = 0;
                long int locNz      = 0;
                size_t   locNcomp   = 0;
                inp.read(reinterpret_cast<char*>(&locTotalNx), sizeof(long int));
                inp.read(reinterpret_cast<char*>(&locNy     ), sizeof(long int));
                inp.read(reinterpret_cast<char*>(&locNz     ), sizeof(long int));
                inp.read(reinterpret_cast<char*>(&locNcomp  ), sizeof(size_t));
                if (locTotalNx != TotalNx or locNy != Ny or locNz != Nz or locNcomp != Ncomp)
                {
                    std::stringstream message;
                    message << thisclassname << " Inconsistent system dimensions!\n"
                            << "Filename: " << filename << "\n"
                            << "Input data dimensions:    (" << locTotalNx << "," << locNy << "," << locNz << ","  << locNcomp << ").\n"
                            << "Required data dimensions: (" <<    TotalNx << "," <<    Ny << "," <<    Nz << ","  <<    Ncomp << ").\n";
                   throw std::runtime_error(message.str());
                }
                for(size_t comp = 0; comp < Ncomp; comp++)
                {
                    inp.read(reinterpret_cast<char*>(&TOC0[comp]), sizeof(size_t));
                    assert(TOC0[comp] != 0);
                }
            }
            else inp.seekg(pos);

            for (long int i = 0; i < Nx;    i++)
            for (long int j = 0; j < Ny;    j++)
            for (long int k = 0; k < Nz;    k++)
            for (size_t   n = 0; n < Ncomp; n++)
            {
                inp.read(reinterpret_cast<char*>(&ChemicalPotential(i,j,k)({n})), sizeof(double));
            }
            pos = inp.tellg();
            inp.close();
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&pos,1,MPI_LONG,proc,MPI_COMM_WORLD);
        if (proc == 0)
        {
            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                MPI_Bcast(&(TOC0[comp]),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
                assert(TOC0[comp] != 0);
            }
        }
    }
#else
    inp.open(filename, std::ios::in | std::ios::binary);
    long int locTotalNx = 0;
    long int locNy      = 0;
    long int locNz      = 0;
    size_t   locNcomp   = 0;
    inp.read(reinterpret_cast<char*>(&locTotalNx), sizeof(long int));
    inp.read(reinterpret_cast<char*>(&locNy     ), sizeof(long int));
    inp.read(reinterpret_cast<char*>(&locNz     ), sizeof(long int));
    inp.read(reinterpret_cast<char*>(&locNcomp  ), sizeof(size_t));
    if (locTotalNx != TotalNx or locNy != Ny or locNz != Nz or locNcomp != Ncomp)
    {
        std::stringstream message;
        message << thisclassname << " Inconsistent system dimensions!\n"
                << "Filename: " << filename << "\n"
                << "Input data dimensions:    (" << locTotalNx << "," << locNy << "," << locNz << ","  << locNcomp << ").\n"
                << "Required data dimensions: (" <<    TotalNx << "," <<    Ny << "," <<    Nz << ","  <<    Ncomp << ").\n";
       throw std::runtime_error(message.str());
    }
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        inp.read(reinterpret_cast<char*>(&TOC0[comp]), sizeof(size_t));
    }
    for (long int i = 0; i < Nx;    i++)
    for (long int j = 0; j < Ny;    j++)
    for (long int k = 0; k < Nz;    k++)
    for (size_t   n = 0; n < Ncomp; n++)
    {
        inp.read(reinterpret_cast<char*>(&ChemicalPotential(i,j,k)({n})), sizeof(double));
    }
    inp.close();
#endif

    BC.SetX(ChemicalPotential);
    BC.SetY(ChemicalPotential);
    BC.SetZ(ChemicalPotential);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,ChemicalPotential.Bcells()-1,)
    {
        CalculateLocalConcentrations (i,j,k,Phase,omega);
        CalculateLocalMobilities     (i,j,k,Phase);
        CaclulateLocalDiffusionFlux  (i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Info::WriteStandard(thisclassname, "Binary input loaded");
}
void Solver::Read(const BoundaryConditions& BC, const PhaseField& Phase, const Density& omega, const std::string& filename)
{
    Read(BC,Phase,omega,filename.c_str());
}
void Solver::Read(const BoundaryConditions& BC, const PhaseField& Phase, const Density& omega, const  long tStep)
{
    std::string filename = UserInterface::MakeFileName(RawDataDir, thisclassname+'_', tStep, ".dat");
    Read(BC,Phase,omega,filename);
}
void Solver::Read(const BoundaryConditions& BC, const PhaseField& Phase, const Density& omega)
{
    std::string filename = RawDataDir+thisclassname+".dat";
    Read(BC,Phase,omega,filename);
}
double Solver::GrandPotential(const PhaseField& Phase, const Density& omega) const
{
    return CalculateVolumeIntegral([this, &Phase, &omega](long i, long j, long k){return omega(i,j,k,Phase);},dx);
}
void Solver::SetInitialConcentration(const PhaseField& Phase, const Density& omega)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,0,)
    {
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            // Calculate initial concentration
            double FinalConcentration = 0;
            for(auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            {
                size_t PhaseIdx = Phase.FieldsStatistics[alpha->index].Phase;
                FinalConcentration += alpha->value*CInitial({PhaseIdx,comp});
            }

            ChemicalPotential(i,j,k)({comp}) = 0.0;

            // Because an initial chemical potential is needed as input
            // a chemical potential which maps into the desired initial
            // concentration is calculated using newtons method
            double precision = DBL_EPSILON;
            double f,df, dmu = 0;
            int n = 0;
            do
            {
                CalculateLocalConcentrations(i,j,k,Phase,omega);
                f  = Concentrations(i,j,k)({comp}) - FinalConcentration;
                df = Susceptibility(i,j,k,comp,Phase,omega);
                if (df == 0.0) break;
                dmu = -f/df/2;
                ChemicalPotential(i,j,k)({comp})+=dmu;
                n++;
                //std::cout << n << " " << f << "\n";
            }
            while (std::abs(dmu) != 0.0 and std::abs(dmu) > std::abs(ChemicalPotential(i,j,k)({comp}))*precision and n < 1000);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void Solver::SetInitialPressure(const PhaseField& Phase, const Density& omega)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,0,)
    {
        // Calculate initial concentration
        ChemicalPotential(i,j,k)({0}) = InitialChemicalPotential[0];

        // Because an initial chemical potential is needed as input
        // a chemical potential which maps into the desired initial
        // concentration is calculated using newtons method
        double precision = 0.1;
        auto PressureDelta = [this, &Phase, &omega, i, j, k](double mu)
        {
            //TODO mu should not be unused!!!!
            return (-omega(i,j,k,Phase)) - InitialPressure;
        };
        auto dPressureDelta_dChemicalPotential = [ &Phase, &omega, i, j, k]([[maybe_unused]] double mu)
        {
            double tmp = 0;
            for(auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            {
                size_t PhaseIdx = Phase.FieldsStatistics[alpha->index].Phase;
                tmp -= alpha->value*omega(PhaseIdx).dChemicalPotential(i,j,k,0);
            }
            return tmp;
        };

        try
        {
            FindRoot::Newton(PressureDelta, dPressureDelta_dChemicalPotential, ChemicalPotential(i,j,k)({0}), precision, 1000);
        }
        catch (std::runtime_error& ecep)
        {
            Info::WriteWarning(ecep.what(),thisclassname,"SetInitialPressure");
        }
        if (std::abs(-omega(i,j,k,Phase)-InitialPressure) > 1.0 )
        {
            std::cerr << i << " " << j << " " << k << " " << (-omega(i,j,k,Phase)) <<  "\n";
            std::abort();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Info::Write("Set Pressure sucessufully!");
}
double Solver::TotalAmountOfComponent(size_t comp) const
{
    return CalculateVolumeIntegral([this,comp](long i, long j, long k){return Concentrations(i,j,k)({comp});},dx);
}

void Solver::SetInitial(const PhaseField& Phase, const Density& omega, const BoundaryConditions& BC, double Temp)
{
    //TODO use either (mu0_0, mu0_1,..) or (rho0, c_0, c_1, ...) to determine initial state
    if   (UseInitialPressure) SetInitialPressure      (Phase, omega);
    else                      SetInitialConcentration (Phase, omega);

    BC.SetX(ChemicalPotential);
    BC.SetY(ChemicalPotential);
    BC.SetZ(ChemicalPotential);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,ChemicalPotential.Bcells()-1,)
    {
        CalculateLocalConcentrations(i,j,k,Phase,omega);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        TOC0[comp] = TotalAmountOfComponent(comp);
        if (std::abs(TOC0[comp]) < DBL_EPSILON)
        {
            std::stringstream message;
            message << "There is no  amount of " << ElementNames[comp] << " present!\n"
                    << "This may be unintentional and hence an error!";
            Info::WriteWarning(message.str(), thisclassname, "SetInitial");
        }
    }
}
void Solver::EnforceConservationOfTOC(const PhaseField& Phase, const Density& omega)
{
    for (std::size_t comp = 0; comp < Ncomp; comp++)
    {
        auto residual = [&] (double& delta)
        {
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
            {
                ChemicalPotential(i,j,k)({comp}) += delta;
                CalculateLocalConcentrations(i,j,k,comp,Phase,omega);
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            double r = (TotalAmountOfComponent(comp) - TOC0[comp])/TOC0[comp];

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,comp)
            {
                ChemicalPotential(i,j,k)({comp}) -= delta;
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            return r;
        };

        double delta = 0;
        try
        {
            FindRoot::Secant(residual, delta, TOCAccuracy, TOCMaxIterations);
        }
        catch (std::runtime_error& ecep)
        {
            Info::WriteWarning(ecep.what(),thisclassname,"EnforceConservation");
        }

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
        {
            ChemicalPotential(i,j,k)({comp}) += delta;
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}
void Solver::CaclulateLocalDiffusionFlux(long i, long j, long k)
{
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        // Calculate chemical potential gradient
        DiffusionFlux(i,j,k)({comp}) = {0.0,0.0,0.0};
        dVector3 locForce = {0.0,0.0,0.0};
        if (dNx) locForce[0] -= (ChemicalPotential(i+1,j  ,k   )({comp}) - ChemicalPotential(i-1,j  ,k  )({comp}))/dx/2.0;
        if (dNy) locForce[1] -= (ChemicalPotential(i  ,j+1,k   )({comp}) - ChemicalPotential(i  ,j-1,k  )({comp}))/dx/2.0;
        if (dNz) locForce[2] -= (ChemicalPotential(i  ,j  ,k +1)({comp}) - ChemicalPotential(i  ,j  ,k-1)({comp}))/dx/2.0;

        DiffusionFlux(i,j,k)({comp}) += locForce*Mobilities(i,j,k)({comp});
    }
}
double Solver::MassDensity(long i, long j, long k, const PhaseField& Phase, const Density& omega) const
{
    double locMassDensity = 0;
    for (size_t comp = 0; comp < Ncomp; comp++)
    {
        locMassDensity += Concentration(i,j,k,comp,Phase,omega)*ElementMasses[comp];
    }
    return locMassDensity;
}
double Solver::Susceptibility(long i, long j, long k, size_t comp, const PhaseField& Phase, const Density& omega) const
{
    // Calculate calculate local bulk susceptibility coefficient
    double locSusceptibility = 0.0;
    for (auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
    {
        size_t PhaseIdx = Phase.FieldsStatistics[alpha->index].Phase;
        double locPhaseSusceptibility = -omega(PhaseIdx).dChemicalPotential2(i,j,k,comp);
        locSusceptibility += alpha->value*locPhaseSusceptibility;
    }
    return locSusceptibility;
}
double Solver::Concentration(long i, long j, long k, size_t comp, const PhaseField& Phase, const Density& omega) const
{
    double locConcentration = 0;
    for (auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
    {
        size_t PhaseIdx = Phase.FieldsStatistics[alpha->index].Phase;
        locConcentration -= alpha->value*omega(PhaseIdx).dChemicalPotential(i,j,k,comp);
    }
    return locConcentration;
}
void Solver::CalculateLocalConcentrations(long i, long j, long k, size_t comp, const PhaseField& Phase, const Density& omega)
{
    Concentrations(i,j,k)({comp}) = 0.0;
    for (auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
    {
        size_t PhaseIdx = Phase.FieldsStatistics[alpha->index].Phase;
        Concentrations(i,j,k)({comp}) -= alpha->value*omega(PhaseIdx).dChemicalPotential(i,j,k,comp);
    }
}
void Solver::CalculateLocalConcentrations(long i, long j, long k, const PhaseField& Phase, const Density& omega)
{
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        CalculateLocalConcentrations(i,j,k,comp,Phase,omega);
    }
}
void Solver::CalculateLocalIncrements(long i, long j, long k, const PhaseField& Phase, const Density& omega, double dt)
{
    for (size_t comp = 0; comp < Ncomp; comp++)
    {
        // Calculate change of concentration due to diffusion
        double ConcentrationDot = 0.0;
        if (dNx) ConcentrationDot -= (DiffusionFlux(i+1,j  ,k   )({comp})[0]-DiffusionFlux(i-1,j  ,k   )({comp})[0])/dx/2.0;
        if (dNy) ConcentrationDot -= (DiffusionFlux(i  ,j+1,k   )({comp})[1]-DiffusionFlux(i  ,j-1,k   )({comp})[1])/dx/2.0;
        if (dNz) ConcentrationDot -= (DiffusionFlux(i  ,j  ,k +1)({comp})[2]-DiffusionFlux(i  ,j  ,k -1)({comp})[2])/dx/2.0;

        double locInverseSusceptibility = 1/Susceptibility(i,j,k,comp,Phase,omega);
        ChemicalPotentialDot(i,j,k)({comp}) += ConcentrationDot*locInverseSusceptibility;

        // Calculate change of concentration due to phase-transformation
        //NodeA locPhaseDot = Phase.Dot(i,j,k, dt);
        for (auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
        {
            size_t PhaseIdx = Phase.FieldsStatistics[alpha->index].Phase;
            double locPhaseConcentration = -omega(PhaseIdx).dChemicalPotential(i,j,k,comp);
            double locPhaseDot = Phase.Dot1(i,j,k, dt).get_value(alpha->index)+Phase.Dot2(i,j,k, dt).get_value(alpha->index);
            ChemicalPotentialDot(i,j,k)({comp}) -= locPhaseConcentration*locPhaseDot*locInverseSusceptibility;
        }
    }
}
void Solver::CalculateLocalIncrements1(long i, long j, long k, const PhaseField& Phase, const Density& omega, double dt)
{
    for (size_t comp = 0; comp < Ncomp; comp++)
    {
        // Calculate change of concentration due to diffusion
        double ConcentrationDot = 0.0;
        if (dNx) ConcentrationDot -= (DiffusionFlux(i+1,j  ,k   )({comp})[0]-DiffusionFlux(i-1,j  ,k   )({comp})[0])/dx/2.0;
        if (dNy) ConcentrationDot -= (DiffusionFlux(i  ,j+1,k   )({comp})[1]-DiffusionFlux(i  ,j-1,k   )({comp})[1])/dx/2.0;
        if (dNz) ConcentrationDot -= (DiffusionFlux(i  ,j  ,k +1)({comp})[2]-DiffusionFlux(i  ,j  ,k -1)({comp})[2])/dx/2.0;

        double locInverseSusceptibility = 1/Susceptibility(i,j,k,comp,Phase,omega);
        ChemicalPotentialDot(i,j,k)({comp}) += ConcentrationDot*locInverseSusceptibility;

        // Calculate change of concentration due to phase-transformation due to curvature
        for (auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
        {
            size_t PhaseIdx = Phase.FieldsStatistics[alpha->index].Phase;
            double locPhaseConcentration = -omega(PhaseIdx).dChemicalPotential(i,j,k,comp);
            double locPhaseDot1 = Phase.Dot1(i,j,k, dt).get_value(alpha->index);
            ChemicalPotentialDot(i,j,k)({comp}) -= locPhaseConcentration*locPhaseDot1*locInverseSusceptibility;
        }
    }
}
void Solver::CalculateLocalIncrements2(long i, long j, long k, const PhaseField& Phase, const Density& omega, double dt)
{
    for (size_t comp = 0; comp < Ncomp; comp++)
    {
        double locInverseSusceptibility = 1/Susceptibility(i,j,k,comp,Phase,omega);
        // Calculate change of concentration due to phase-transformation due to pressure difference
        for (auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
        {
            size_t PhaseIdx = Phase.FieldsStatistics[alpha->index].Phase;
            double locPhaseConcentration = -omega(PhaseIdx).dChemicalPotential(i,j,k,comp);
            double locPhaseDot2 = Phase.Dot2(i,j,k,dt).get_value(alpha->index);
            ChemicalPotentialDot2(i,j,k)({comp}) -= locPhaseConcentration*locPhaseDot2*locInverseSusceptibility;
        }
    }
}
void Solver::CalculateLocalPhaseFieldIncrements (long i, long j, long k, PhaseField& Phase, const Density& omega, const InterfaceProperties& IP)
{
    for (auto it = Phase.FieldsDot(i,j,k).begin(); it != Phase.FieldsDot(i,j,k).end(); ++it)
    {
        Phase.FieldsDot(i,j,k).set_sym2(it->indexA,it->indexB,0.0);
    }

    for(auto alpha = Phase.Fields(i, j, k).cbegin(); alpha != Phase.Fields(i, j, k).cend(); ++alpha)
    {
        size_t AlphaIdx = Phase.FieldsStatistics[alpha->index].Phase;
        double locOmegaAlpha = omega(AlphaIdx)(i,j,k);
        for(auto  beta = alpha + 1; beta != Phase.Fields(i, j, k).cend();  ++beta)
        {
            size_t BetaIdx = Phase.FieldsStatistics[beta->index].Phase;
            if (AlphaIdx != BetaIdx)
            {
                double locOmegaBeta  = omega(BetaIdx)(i,j,k);
                double locOmegaDelta = locOmegaBeta - locOmegaAlpha ;
                double Prefactor     = Pi*Pi/(8.0*Phase.Eta*Phase.LocalNumberOfPhaseFieldsSR(i,j,k));
                double loc_dPhi_dt   = locOmegaDelta*IP.get_mobility(i,j,k, alpha->index, beta->index)*Prefactor;
                Phase.FieldsDot(i,j,k).add_asym2(alpha->index, beta->index, loc_dPhi_dt);
            }
        }
    }
}
void Solver::MergeLocalIncrements(long i, long j, long k, double dt)
{
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        ChemicalPotential   (i,j,k)({comp}) += dt*ChemicalPotentialDot(i,j,k)({comp});
        ChemicalPotentialDot(i,j,k)({comp}) = 0.0;
    }
}
void Solver::MergeLocalIncrements1(long i, long j, long k, double dt)
{
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        ChemicalPotentialOld(i,j,k)({comp})  = ChemicalPotential(i,j,k)({comp});
        ChemicalPotential   (i,j,k)({comp}) += dt*ChemicalPotentialDot(i,j,k)({comp});
        ChemicalPotentialDot(i,j,k)({comp})  = 0.0;
    }
}
void Solver::MergeLocalIncrements2(long i, long j, long k, double dt)
{
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        ChemicalPotentialDot (i,j,k)({comp})  = dt*ChemicalPotentialDot2(i,j,k)({comp});
        ChemicalPotential    (i,j,k)({comp}) += ChemicalPotentialDot(i,j,k)({comp});
        ChemicalPotentialDot2(i,j,k)({comp})  = 0.0;
    }
}
void Solver::MergeLocalIncrements2Implicit(long i, long j, long k, PhaseField& Phase, const Density& omega, const InterfaceProperties& IP, const double dt)
{
    if (Ncomp == 1)
    {
        // This method uses the secant method (Newton's method with numeric calculation of the gradient) to solve the semi implicit time-step
        double ChemicalPotential_Curvature = ChemicalPotential(i,j,k)({0});
        auto residual = [this,i,j,k,&Phase,&omega,&IP,dt,ChemicalPotential_Curvature](double mu)
        {
            double mu0 = ChemicalPotential(i,j,k)({0});
            ChemicalPotential(i,j,k)({0}) = mu;
            ChemicalPotentialDot2(i,j,k)({0}) = 0.0;

            CalculateLocalConcentrations       (i,j,k,Phase,omega);
            CalculateLocalPhaseFieldIncrements (i,j,k,Phase,omega,IP);
            CalculateLocalIncrements2          (i,j,k,Phase,omega,dt);

            double res = ChemicalPotential(i,j,k)({0}) - ChemicalPotential_Curvature - ChemicalPotentialDot2(i,j,k)({0})*dt;
            ChemicalPotential(i,j,k)({0}) = mu0;
            return res;
        };

        // Use chemical potential of the previous time step as starting point for the root fining algorithm
        ChemicalPotential(i,j,k)({0}) = ChemicalPotentialOld(i,j,k)({0});

        try
        {
            FindRoot::Secant(residual, ChemicalPotential(i,j,k)({0}), ChemicalPotentialAccuracy, MaxIterations);
            //int error = FindRoot::Broyden(residual, ChemicalPotential(i,j,k)({0}), ChemicalPotentialAccuracy, MaxIterations);
        }
        catch (std::runtime_error& ecep)
        {
            Info::WriteWarning(ecep.what(),thisclassname,"MergeLocalIncrements2Implicit");
        }
    }
    else
    {
        // This method uses the Broyden's method to solve the semi implicit
        // time-step (Newton's method for multiple dimensions where the
        // inverse Jacobian is calculated numerically only for the first iteration).
        typedef std::vector<double> vector;
        vector ChemicalPotential_Curvature(Ncomp);
        for (std::size_t comp = 0; comp < Ncomp; comp++)
        {
            ChemicalPotential_Curvature[comp] = ChemicalPotential(i,j,k)({comp});
        }
        auto resiudal = [this,i,j,k,&Phase,&omega,&IP,dt,&ChemicalPotential_Curvature] (vector& x0)
        {
            std::vector<double> mu (Ncomp);
            std::vector<double> f0 (Ncomp);
            for (std::size_t comp = 0; comp < Ncomp; comp++)
            {
                mu[comp] = ChemicalPotential (i,j,k)({comp});
                ChemicalPotential     (i,j,k)({comp}) = x0[comp];
                ChemicalPotentialDot2 (i,j,k)({comp}) = 0;
            }

            CalculateLocalConcentrations       (i,j,k,Phase,omega);
            CalculateLocalPhaseFieldIncrements (i,j,k,Phase,omega,IP);
            CalculateLocalIncrements2          (i,j,k,Phase,omega,dt);

            for (std::size_t comp = 0; comp < Ncomp; comp++)
            {
                f0[comp] = ChemicalPotential(i,j,k)({comp}) - ChemicalPotential_Curvature[comp] - ChemicalPotentialDot2(i,j,k)({comp})*dt;
                ChemicalPotential (i,j,k)({comp}) = mu[comp];
            }
            return f0;
        };

        // Use chemical potential of the previous time step as starting point for the root fining algorithm
        std::vector<double> mu (Ncomp);
        for (std::size_t comp = 0; comp < Ncomp; comp++)
        {
            mu[comp] = ChemicalPotentialOld(i,j,k)({comp});
        }

        try
        {
            FindRoot::Broyden(resiudal, mu, ChemicalPotentialAccuracy, MaxIterations);
        }
        catch (std::runtime_error& ecep)
        {
            Info::WriteWarning(ecep.what(),thisclassname,"MergeLocalIncrements2Implicit");
        }

        for (std::size_t comp = 0; comp < Ncomp; comp++)
        {
            ChemicalPotential(i,j,k)({comp}) = mu[comp] ;
        }
    }
}
void Solver::CalculateLocalMobilities(long i, long j, long k, const PhaseField& Phase)
{
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        double locMobility = 0;

        for (auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
        {
            size_t PhaseIdx = Phase.FieldsStatistics[alpha->index].Phase;
            locMobility += alpha->value*PhaseMobilities({PhaseIdx,comp});
        }

        if (Use_InterfaceMobilities)
        for (auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
        for (auto beta = alpha+1; beta  != Phase.Fields(i,j,k).cend();  beta++)
        {
            size_t PhaseIdx = Phase.FieldsStatistics[alpha->index].Phase;
            {
                size_t BetaIdx = Phase.FieldsStatistics[beta->index].Phase;
                locMobility += 4*alpha->value*beta->value*InterfaceMobilities({PhaseIdx,BetaIdx,comp});
            }
        }

        if (Use_dPhaseMobility_dConcentration)
        for(size_t compB = 0; compB < Ncomp; compB++)
        for (auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
        {
            size_t PhaseIdx = Phase.FieldsStatistics[alpha->index].Phase;
            locMobility += alpha->value*dPhaseMobilities_dConcentration({PhaseIdx,comp,compB});
        }

        Mobilities(i,j,k)({comp}) = locMobility;
    }
}
void Solver::SolveExplicit(PhaseField& Phase, const Density& omega, const BoundaryConditions&  BC, const InterfaceProperties& IP, const double dt, [[maybe_unused]] const double Temp)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,ChemicalPotential.Bcells()-1,)
    {
        CalculateLocalConcentrations      (i,j,k,Phase,omega);
        CalculateLocalPhaseFieldIncrements(i,j,k,Phase,omega,IP);
        CalculateLocalMobilities          (i,j,k,Phase);
        CaclulateLocalDiffusionFlux       (i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,ChemicalPotential.Bcells()-2,)
    {
        CalculateLocalIncrements(i,j,k,Phase,omega,dt);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,ChemicalPotential.Bcells()-2,)
    {
        MergeLocalIncrements(i,j,k,dt);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    if (ConserveTOC) EnforceConservationOfTOC(Phase,omega);

    if (dNx) BC.SetX(ChemicalPotential);
    if (dNy) BC.SetY(ChemicalPotential);
    if (dNz) BC.SetZ(ChemicalPotential);
}
void Solver::SolveImplicit(PhaseField& Phase, const Density& omega, const BoundaryConditions& BC, const InterfaceProperties& IP, const double dt, [[maybe_unused]] const double Temp)
{
    // Apply first part of rhs explicitly
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,ChemicalPotential.Bcells()-1,)
    {
        CalculateLocalConcentrations      (i,j,k,Phase,omega);
        CalculateLocalPhaseFieldIncrements(i,j,k,Phase,omega,IP);
        CalculateLocalMobilities          (i,j,k,Phase);
        CaclulateLocalDiffusionFlux       (i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,ChemicalPotential.Bcells()-2,)
    {
        CalculateLocalIncrements1(i,j,k,Phase,omega,dt);
        //CalculateLocalIncrements2(i,j,k,Phase,omega,dt);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ChemicalPotential,ChemicalPotential.Bcells()-2,)
    {
        MergeLocalIncrements1(i,j,k,dt);
        MergeLocalIncrements2Implicit(i,j,k,Phase,omega,IP,dt);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    if (ConserveTOC) EnforceConservationOfTOC(Phase,omega);

    if (dNx) BC.SetX(ChemicalPotential);
    if (dNy) BC.SetY(ChemicalPotential);
    if (dNz) BC.SetZ(ChemicalPotential);
}
void Solver::Solve(PhaseField& Phase, const Density& omega, const BoundaryConditions& BC, const InterfaceProperties& IP, double dt, double Temp)
{
    if (Phase.Fields.Bcells() < 2)
    {
        Info::WriteExit("Too few phase-field boundary cell! Two are required!", thisclassname, "Initialize");
        std::abort();
    }

    if (UseImplicitSolver)
    {
        SolveImplicit(Phase, omega, BC, IP, dt, Temp);
    }
    else
    {
        SolveExplicit(Phase, omega, BC, IP, dt, Temp);
    }
}
double Solver::MolarVolume(long i, long j, long k, const PhaseField& Phase, const Density& omega) const
{
    double locTotalConcentration = 0.0;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        locTotalConcentration += Concentration(i,j,k,comp,Phase,omega);
    }
    if (locTotalConcentration < 0 ) locTotalConcentration = 0;
    double locMolarVolume = (locTotalConcentration > 0.0) ? 1.0/locTotalConcentration : std::numeric_limits<float>::max();
    return locMolarVolume;
}
double Solver::MoleFraction(long i, long j, long k, size_t comp, const PhaseField& Phase, const Density& omega) const
{
     double value = MolarVolume(i,j,k,Phase,omega)*Concentration(i,j,k,comp,Phase,omega);
     if (value > 1) value = 1;
     if (value < 0) value = 0;
     return value;
}
void Solver::WriteVTK(long tStep, const Settings& locSettings, const PhaseField& Phase, const Density& omega, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
         ListOfFields.push_back((VTK::Field_t){"Chemical Potential "+ElementNames[comp]+" [J]", [this, comp                ] (long i, long j, long k) {return ChemicalPotential (i,j,k)({comp});}});
         ListOfFields.push_back((VTK::Field_t){ElementNames[comp]+ " [mol/m^3]"               , [this, comp, &Phase, &omega] (long i, long j, long k) {return Concentration     (i,j,k,comp,Phase,omega);}});
         ListOfFields.push_back((VTK::Field_t){ElementNames[comp]+ " [mol/mol]"               , [this, comp, &Phase, &omega] (long i, long j, long k) {return MoleFraction      (i,j,k,comp,Phase,omega);}});
    }
    ListOfFields.push_back((VTK::Field_t){"Molar Volume [m^3/mol]",    [this, &Phase, &omega] (long i, long j, long k) {return MolarVolume(i,j,k,Phase,omega);}});
    ListOfFields.push_back((VTK::Field_t){"Grand Potential Density [J/m^3]", [&Phase, &omega] (long i, long j, long k) {return omega(i,j,k,Phase);}});
    std::string Filename = UserInterface::MakeFileName(VTKDir, thisclassname + "_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
void Solver::WriteVTKChemicalPotential(long tStep, const Settings& locSettings, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
         ListOfFields.push_back((VTK::Field_t){"Chemical Potential "+ElementNames[comp]+" [J]", [this, comp] (long i, long j, long k) {return ChemicalPotential (i,j,k)({comp});}});
    }
    std::string Filename = UserInterface::MakeFileName(VTKDir, "ChemicalPotential_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
void Solver::WriteVTKConcentration(long tStep, const Settings& locSettings, const PhaseField& Phase, const Density& omega, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
         ListOfFields.push_back((VTK::Field_t){ElementNames[comp]+" [1/m^3]" , [this, comp, &Phase, &omega] (long i, long j, long k) {return Concentration(i,j,k,comp,Phase,omega);}});
    }
    std::string Filename = UserInterface::MakeFileName(VTKDir, "Concentration_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void Solver::WriteVTKMoleFraction(long tStep, const Settings& locSettings, const PhaseField& Phase, const Density& omega, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
         ListOfFields.push_back((VTK::Field_t){ElementNames[comp]+ " [mol/mol]", [this, comp, &Phase, &omega] (long i, long j, long k) {return MoleFraction(i,j,k,comp,Phase,omega);}});
    }
    std::string Filename = UserInterface::MakeFileName(VTKDir, "MoleFraction_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
void Solver::WriteVTKMolarVolume(long tStep, const Settings& locSettings, const PhaseField& Phase, const Density& omega, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t){"Molar Volume [mole/m^3]" , [this, &Phase, &omega] (long i, long j, long k) {return MolarVolume(i,j,k,Phase,omega);}});
    std::string Filename = UserInterface::MakeFileName(VTKDir, "MolarVolume_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
void Solver::WriteVTKMobility(long tStep, const Settings& locSettings, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        ListOfFields.push_back((VTK::Field_t){"Mobility_"+ElementNames[comp]+" [(m^3 x)/kg]", [this,comp] (long i, long j, long k) {return Mobilities(i,j,k)({comp});}});
    }
    std::string Filename = UserInterface::MakeFileName(VTKDir, "Mobility_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
void Solver::WriteVTKDiffusionFlux(long tStep, const Settings& locSettings, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        ListOfFields.push_back((VTK::Field_t){"DiffusionFlux_"+ElementNames[comp]+" [1/(m^2 s)]", [this,comp] (long i, long j, long k) {return DiffusionFlux(i,j,k)({comp});}});
    }
    std::string Filename = UserInterface::MakeFileName(VTKDir, "DiffusionFlux_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
void Solver::WriteVTKVelocity(long tStep, const Settings& locSettings, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        ListOfFields.push_back((VTK::Field_t){"Velocity_"+ElementNames[comp]+" [1/(m^2 s)]", [this,comp] (long i, long j, long k) {return (Concentrations(i,j,k)({comp}) > 0.0) ? DiffusionFlux(i,j,k)({comp})/Concentrations(i,j,k)({comp}) : dVector3({0.0,0.0,0.0});}});
    }
    std::string Filename = UserInterface::MakeFileName(VTKDir, "Velocity_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
void Solver::WriteVTKSusceptibility(long tStep, const Settings& locSettings, const PhaseField& Phase, const Density& omega, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        ListOfFields.push_back((VTK::Field_t){"Susceptibility_"+ElementNames[comp]+" [(m^3 x)/kg]", [this,comp,&Phase,&omega] (long i, long j, long k) {return Susceptibility(i,j,k,comp,Phase,omega);}});
    }
    std::string Filename = UserInterface::MakeFileName(VTKDir, "Susceptibility_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
void Solver::WriteVTKDiffusivity(long tStep, const Settings& locSettings, const PhaseField& Phase, const Density& omega, long precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        ListOfFields.push_back((VTK::Field_t){"Diffusivity_"+ElementNames[comp]+" [(m^3 x)/kg]", [this,comp,&Phase,&omega] (long i, long j, long k) {return Mobilities(i,j,k)({comp})/Susceptibility(i,j,k,comp,Phase,omega);}});
    }
    std::string Filename = UserInterface::MakeFileName(VTKDir, "Diffusivity_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}
}// namespace openphase::GrandPotential
