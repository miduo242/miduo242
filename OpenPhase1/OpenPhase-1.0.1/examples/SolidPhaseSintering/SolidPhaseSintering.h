/*
 *   This file is part of the OpenPhase (R) software library.
 *
 *   Copyright (c) 2009-2020 Ruhr-Universitaet Bochum,
 *                 Universitaetsstrasse 150, D-44801 Bochum, Germany
 *             AND 2018-2020 OpenPhase Solutions GmbH,
 *                 Wasserstrasse 494, D-44795 Bochum, Germany.
 *
 *    All rights reserved.
 *
 *
 *    DEVELOPMENT VERSION, DO NOT PUBLISH OR DISTRIBUTE.
 *
 *
 *   OpenPhase (R) is a joint development of Interdisciplinary Centre for
 *   Advanced Materials Simulation (ICAMS), Ruhr University Bochum
 *   and OpenPhase Solutions GmbH.
 *
 *   File created :   2021
 *   Main contributors :   Raphael Schiedung
 *
 */

#ifndef SOLIDPHASESINTERING_H
#define SOLIDPHASESINTERING_H

#include "PhaseField.h"
#include "GrandPotential/Solver.h"
#include "GrandPotential/Density.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "RunTimeControl.h"

namespace op = openphase;

void ReadInitializationInuptParameters(op::Settings& locSettings,
        std::string InputFileName, double& MeanRadius,
        double& StdRadius, int& X0, int& XN, int& Y0, int& YN, int& Z0, int& ZN)
{
    // Read initialization input parameters
    op::Info::WriteBlankLine();
    op::Info::WriteLineInsert("ProjectInput");
    op::Info::WriteStandard("Source", InputFileName);
    std::fstream inp(InputFileName);
    if (!inp)
    {
        std::cerr << "File \"" <<  InputFileName << "\" could not be opened\n";
        std::exit(EXIT_FAILURE);
    };
    std::stringstream inp_data;
    inp_data << inp.rdbuf();
    inp.close();
    int moduleLocation = op::UserInterface::FindModuleLocation(inp_data, "ProjectInput");
    MeanRadius = op::UserInterface::ReadParameterD(inp_data, moduleLocation, "AvR" )/locSettings.dx;
    StdRadius  = op::UserInterface::ReadParameterD(inp_data, moduleLocation, "StdR")/locSettings.dx;
    double SX = op::UserInterface::ReadParameterD(inp_data, moduleLocation, "SX")/locSettings.dx;
    double SY = op::UserInterface::ReadParameterD(inp_data, moduleLocation, "SY")/locSettings.dx;
    double SZ = op::UserInterface::ReadParameterD(inp_data, moduleLocation, "SZ")/locSettings.dx;
    X0 = std::max(int((locSettings.TotalNx - SX)/2), 0);
    XN = std::min(int((locSettings.TotalNx - SX)/2 + SX), locSettings.TotalNx);
    Y0 = std::max(int((locSettings.TotalNy - SY)/2), 0);
    YN = std::min(int((locSettings.TotalNy - SY)/2 + SY), locSettings.TotalNx);
    Z0 = std::max(int((locSettings.TotalNz - SZ)/2), 0);
    ZN = std::min(int((locSettings.TotalNz - SZ)/2 + SZ), locSettings.TotalNx);
}

namespace CalculateDensity
{
double BoundingBox(const op::PhaseField &Phase,
        const op::GrandPotential::Density& omega,
        const op::GrandPotential::Solver& GPS,
        long int i0, long int i1, long int j0, long int j1, long int k0, long int k1)
{
    long int i_min = i0 - Phase.OffsetX;
    long int i_max = i1 - Phase.OffsetX;
    long int j_min = j0 - Phase.OffsetY;
    long int j_max = j1 - Phase.OffsetY;
    long int k_min = k0 - Phase.OffsetZ;
    long int k_max = k1 - Phase.OffsetZ;

    if (i_min < 0) i_min = 0;
    if (j_min < 0) j_min = 0;
    if (k_min < 0) k_min = 0;
    if (i_max > Phase.Fields.sizeX()) i_max = Phase.Fields.sizeX();
    if (j_max > Phase.Fields.sizeY()) j_max = Phase.Fields.sizeY();
    if (k_max > Phase.Fields.sizeZ()) k_max = Phase.Fields.sizeZ();

    double LocMass = 0;
    double Mass    = 0;
    if (i_min < Phase.Fields.sizeX())
    if (j_min < Phase.Fields.sizeY())
    if (k_min < Phase.Fields.sizeZ())
    if (i_max > 0)
    if (j_max > 0)
    if (k_max > 0)
    {
        #pragma omp parallel for collapse(3) reduction(+:LocMass)
        for (long int i = i_min; i < i_max; ++i)
        for (long int j = j_min; j < j_max; ++j)
        for (long int k = k_min; k < k_max; ++k)
        {
            LocMass += GPS.MassDensity(i,j,k,Phase,omega);
        }
    }
    #ifdef MPI_PARALLEL
        MPI_Allreduce(&LocMass, &Mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #else
        Mass = LocMass;
    #endif
    double Volume = (i1-i0)*(j1-j0)*(k1-k0);
    if (Volume > 0.0) return Mass/Volume;
    else return 0.0;
}
double LineIntersectionX(const op::PhaseField &Phase,
        const op::GrandPotential::Density& omega,
        const op::GrandPotential::Solver& GPS)
{
    double LocMass   = 0;
    double Mass      = 0;
    double LocVolume = 0;
    double Volume    = 0;

    if (Phase.Fields.sizeX() > 1)
    {
        #ifdef MPI_PARALLEL
        //TODO
        #else
            #pragma omp parallel for collapse(2) reduction(+:LocVolume,LocMass)
            for (long int j = 0; j < Phase.Fields.sizeY(); ++j)
            for (long int k = 0; k < Phase.Fields.sizeZ(); ++k)
            {
                long int i_min = Phase.Fields.sizeX()-1;
                for (int n = 0; n < Phase.Fields.sizeX(); ++n)
                if (Phase.Fractions(n,j,k)({0}) <= 0.05)
                {
                    i_min = n;
                    break;
                }
                long int i_max = 0;
                for (int n = Phase.Fields.sizeX()-1; n > 0; --n)
                if (Phase.Fractions(n,j,k)({0}) <= 0.05)
                {
                    i_max = n;
                    break;
                }
                for (long int i = i_min; i < i_max; ++i)
                {
                    LocMass   += GPS.MassDensity(i,j,k,Phase,omega);
                    LocVolume += 1.0;
                }
            }
        #endif
    }
    else
    {
        #pragma omp parallel for collapse(2) reduction(+:LocVolume,LocMass)
        for (long int j = 0; j < Phase.Fields.sizeY(); ++j)
        for (long int k = 0; k < Phase.Fields.sizeZ(); ++k)
        {
            if (Phase.Fractions(0,j,k)({0}) <= 0.05)
            {
                LocMass   += GPS.MassDensity(0,j,k,Phase,omega);
                LocVolume += 1.0;
            }
        }
    }

    #ifdef MPI_PARALLEL
        MPI_Allreduce(&LocMass,   &Mass,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&LocVolume, &Volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #else
        Mass   = LocMass;
        Volume = LocVolume;
    #endif
    if (Volume > 0.0) return Mass/Volume;
    else return 0.0;
}
double LineIntersectionY(const op::PhaseField &Phase,
        const op::GrandPotential::Density& omega,
        const op::GrandPotential::Solver& GPS)
{
    double LocMass   = 0;
    double Mass      = 0;
    double LocVolume = 0;
    double Volume    = 0;

    if (Phase.Fields.sizeY() > 1)
    {
        #pragma omp parallel for collapse(2) reduction(+:LocVolume,LocMass)
        for (long int i = 0; i < Phase.Fields.sizeX(); ++i)
        for (long int k = 0; k < Phase.Fields.sizeZ(); ++k)
        {
            long int j_min = Phase.Fields.sizeY()-1;
            for (int n = 0; n < Phase.Fields.sizeY(); ++n)
            if (Phase.Fractions(i,n,k)({0}) <= 0.05)
            {
                j_min = n;
                break;
            }
            long int j_max = 0;
            for (int n =  Phase.Fields.sizeY()-1; n > 0; --n)
            if (Phase.Fractions(i,n,k)({0}) <= 0.05)
            {
                j_max = n;
                break;
            }
            for (long int j = j_min; j < j_max; ++j)
            {
                LocMass   += GPS.MassDensity(i,j,k,Phase,omega);
                LocVolume += 1.0;
            }
        }
    }
    else
    {
        #pragma omp parallel for collapse(2) reduction(+:LocVolume,LocMass)
        for (long int i = 0; i < Phase.Fields.sizeX(); ++i)
        for (long int k = 0; k < Phase.Fields.sizeZ(); ++k)
        {
            if (Phase.Fractions(i,0,k)({0}) <= 0.05)
            {
                LocMass   += GPS.MassDensity(i,0,k,Phase,omega);
                LocVolume += 1.0;
            }
        }
    }
    #ifdef MPI_PARALLEL
        MPI_Allreduce(&LocMass,   &Mass,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&LocVolume, &Volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #else
        Mass   = LocMass;
        Volume = LocVolume;
    #endif
    if (Volume > 0.0) return Mass/Volume;
    else return 0.0;
}
double LineIntersectionZ(const op::PhaseField &Phase,
        const op::GrandPotential::Density& omega,
        const op::GrandPotential::Solver& GPS)
{
    double LocMass   = 0;
    double Mass      = 0;
    double LocVolume = 0;
    double Volume    = 0;

    if (Phase.Fields.sizeZ() > 1)
    {
        #pragma omp parallel for collapse(2) reduction(+:LocVolume,LocMass)
        for (long int i = 0; i < Phase.Fields.sizeX(); ++i)
        for (long int j = 0; j < Phase.Fields.sizeY(); ++j)
        {
            long int k_min = Phase.Fields.sizeZ()-1;
            for (int n = 0; n < Phase.Fields.sizeZ(); ++n)
            if (Phase.Fractions(i,j,n)({0}) <= 0.05)
            {
                k_min = n;
                break;
            }
            long int k_max = 0;
            for (int n = Phase.Fields.sizeZ()-1; n > 0; --n)
            if (Phase.Fractions(i,j,n)({0}) <= 0.05)
            {
                k_max = n;
                break;
            }
            for (long int k = k_min; k <= k_max; ++k)
            {
                LocMass   += GPS.MassDensity(i,j,k,Phase,omega);
                LocVolume += 1.0;
            }
        }
    }
    else
    {
        #pragma omp parallel for collapse(2) reduction(+:LocVolume,LocMass)
        for (long int i = 0; i < Phase.Fields.sizeX(); ++i)
        for (long int j = 0; j < Phase.Fields.sizeY(); ++j)
        {
            if (Phase.Fractions(i,j,0)({0}) <= 0.05)
            {
                LocMass   += GPS.MassDensity(i,j,0,Phase,omega);
                LocVolume += 1.0;
            }

        }
    }
    #ifdef MPI_PARALLEL
        MPI_Allreduce(&LocMass,   &Mass,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&LocVolume, &Volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #else
        Mass   = LocMass;
        Volume = LocVolume;
    #endif
    if (Volume > 0.0) return Mass/Volume;
    else return 0.0;
}

}
void DoSinteringDiagnostics(
        op::Settings& locSettings,
        op::PhaseField& Phase,
        op::DoubleObstacle& DO,
        op::InterfaceProperties& IP,
        op::GrandPotential::Density& omega,
        op::GrandPotential::Solver& GPS,
        op::RunTimeControl& RTC)
{
    std::fstream log(locSettings.TextDir + "TimeLog.csv", std::ios::app|std::ios::out);
    //double Density   = 0.0;
    double DensityLI  = 0.0;
    double DensityLIX = 0.0;
    double DensityLIY = 0.0;
    double DensityLIZ = 0.0;

    double dim = 0;
    if (locSettings.Nx > 1) dim+=1.0;
    if (locSettings.Ny > 1) dim+=1.0;
    if (locSettings.Nz > 1) dim+=1.0;

    //Density       = CalculateDensity::BoundingBox(Phase,omega,GPS,ii0,ii1,jj0,jj1,kk0,kk1);
    DensityLIX    = CalculateDensity::LineIntersectionX(Phase,omega,GPS);
    DensityLIY    = CalculateDensity::LineIntersectionY(Phase,omega,GPS);
    DensityLIZ    = CalculateDensity::LineIntersectionZ(Phase,omega,GPS);
    DensityLI     = 0.0;
    if (locSettings.Nx > 1) DensityLI += DensityLIX;
    if (locSettings.Ny > 1) DensityLI += DensityLIY;
    if (locSettings.Nz > 1) DensityLI += DensityLIZ;
    DensityLI    /= dim;

    std::array<std::stringstream,2> line;
    line[1] << std::scientific << std::setprecision(16);
    op::Info::WriteWithLog(line , RTC.tStep , "Time [s]", RTC.tStep*RTC.dt);
    op::Info::Write("-");
    op::Info::WriteWithLog(line , RTC.tStep , "Interface energy [J]"      , DO.Energy(Phase, IP));
    //op::Info::WriteWithLog(line , RTC.tStep , "Surface energy [J]"        , DO.Energy(Phase, IP,0,1));
    //op::Info::WriteWithLog(line , RTC.tStep , "Grain boundary energy [J]" , DO.Energy(Phase, IP,1,1));
    op::Info::WriteWithLog(line , RTC.tStep , "Bulk grand potential [J]"  , GPS.GrandPotential(Phase, omega));
    op::Info::Write("-");
    for (size_t comp = 0; comp < locSettings.Ncomp; comp++)
    {
        op::Info::WriteWithLog(line , RTC.tStep , "Amount of "             +locSettings.ElementNames[comp]+" [mol]"    , GPS.TotalAmountOfComponent(comp));
        op::Info::WriteWithLog(line , RTC.tStep , "Vapor concentration of "+locSettings.ElementNames[comp]+" [mol/m^3]", GPS.Concentration(0,0,0,comp,Phase,omega));
    }
    op::Info::Write("-");
    //op::Info::WriteWithLog(line , RTC.tStep , "Density BB [kg/m^3]"      , Density);
    op::Info::WriteWithLog(line , RTC.tStep , "Density LI [kg/m^3]"      , DensityLI);
    op::Info::WriteWithLog(line , RTC.tStep , "Density LIX [kg/m^3]"     , DensityLIX);
    op::Info::WriteWithLog(line , RTC.tStep , "Density LIY [kg/m^3]"     , DensityLIY);
    op::Info::WriteWithLog(line , RTC.tStep , "Density LIZ [kg/m^3]"     , DensityLIZ);
    op::Info::Write("-");
    for (size_t n = 1; n < locSettings.Nphases; n++)
    {
        size_t NGrains      = Phase.FieldsStatistics.NumberOfGrains();
        double dx           = locSettings.dx;
        double MeanVolume   = Phase.FieldsStatistics.MeanGrainVolume(n)*dx*dx*dx;
        double StddevVolume = Phase.FieldsStatistics.StandardDeviationOfGrainVolumes(n)*dx*dx*dx;
        double MeanRadius   = Phase.FieldsStatistics.MeanGrainRadius(n,locSettings.dNz != 0)*dx;
        double StddevRadius = Phase.FieldsStatistics.StandardDeviationOfGrainRadii(n,locSettings.dNz != 0)*dx;
        op::Info::WriteWithLog(line , RTC.tStep , "Number of grains of "       +locSettings.PhaseNames[n]         , NGrains     );
        op::Info::WriteWithLog(line , RTC.tStep , "Mean grain volume of "      +locSettings.PhaseNames[n]+" [m^3]", MeanVolume  );
        op::Info::WriteWithLog(line , RTC.tStep , "Stddev of grain volumes of "+locSettings.PhaseNames[n]+" [m^3]", StddevVolume);
        op::Info::WriteWithLog(line , RTC.tStep , "Mean grain radius of "      +locSettings.PhaseNames[n]+" [m]"  , MeanRadius  );
        op::Info::WriteWithLog(line , RTC.tStep , "Stddev of grain radii of "  +locSettings.PhaseNames[n]+" [m]"  , StddevRadius);
    }
    op::Info::WriteLineToLogfile(log, line, RTC.tStep);
}
#endif
