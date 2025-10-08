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
 *   File created :   2011
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */

#include "Info.h"
#include "DrivingForce.h"
#include "Settings.h"
#include "Base/UserInterface.h"
#include "VTK.h"
#include "GrainInfo.h"
#include "PhaseField.h"
#include "InterfaceProperties.h"
#include "BoundaryConditions.h"
#include "Base/CommonFunctions.h"

namespace openphase
{

using namespace std;
DrivingForce::DrivingForce(Settings& locSettings, const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void DrivingForce::Initialize(Settings& locSettings)
{
    thisclassname = "DrivingForce";

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

    Resolution = locSettings.Resolution;

    Nphases = locSettings.Nphases;
    Noise.Allocate(Nphases,Nphases);

    // Setting default values to be used if ReadInput() is not called
    dGcut = true;
    CutOff.Allocate(Nphases,Nphases);
    Averaging = false;

    if(locSettings.iWidth < 5.0)
    {
        PhiThreshold = locSettings.iWidth/15.0;
    }
    else
    {
        PhiThreshold = 1.0/3.0;
    }

    Range = max(locSettings.Bcells, size_t(locSettings.iWidth));

    if(Resolution == Resolutions::Double)
    {
        Range = locSettings.iWidth/2 + 1;
    }

    // End setting default values

    Raw.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, Range);

    if(Resolution == Resolutions::Double)
    {
        RawDR.Allocate((1+dNx)*Nx, (1+dNy)*Ny, (1+dNz)*Nz, dNx, dNy, dNz, Range*2);
    }

    OvershootCounter = 0;
    MAXOvershootPOS.Allocate(Nphases,Nphases);
    MAXOvershootNEG.Allocate(Nphases,Nphases);
    MAXDrivingForcePOS.Allocate(Nphases,Nphases);
    MAXDrivingForceNEG.Allocate(Nphases,Nphases);

    VTKDir = locSettings.VTKDir;
    RawDataDir = locSettings.RawDataDir;

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void DrivingForce::ReadInput(const string InputFileName)
{
    Info::WriteLineInsert("DrivingForce input");
    Info::WriteStandard("Source", InputFileName.c_str());

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };


    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);

    inp.close();
}

void DrivingForce::ReadInput(stringstream& inp)
{
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    Averaging    = UserInterface::ReadParameterB(inp, moduleLocation, string("Average"), true, true);
    Unifying     = UserInterface::ReadParameterB(inp, moduleLocation, string("bUnify"), false, false);
    Range        = UserInterface::ReadParameterI(inp, moduleLocation, string("Range"), false, Range);
    PhiThreshold = UserInterface::ReadParameterD(inp, moduleLocation, string("Threshold"), false, PhiThreshold);
    dGcut        = UserInterface::ReadParameterB(inp, moduleLocation, string("dGcut"), false, true);

    string tmp = UserInterface::ReadParameterK(inp, moduleLocation, string("NoiseMode"), false, "OFF");

    if(tmp == string("OFF")) NoiseMode = NoiseModes::Zero;
    if(tmp == string("RELATIVE")) NoiseMode = NoiseModes::Relative;
    if(tmp == string("ABSOLUTE")) NoiseMode = NoiseModes::Absolute;

 	tmp = UserInterface::ReadParameterK(inp, moduleLocation, string("dGClearingMode"), false, "AUTOMATIC");

    if(tmp == string("AUTOMATIC")) dGClearingMode = ClearingModes::Automatic;
    if(tmp == string("MANUAL")) dGClearingMode = ClearingModes::Manual;

	
    for(size_t m = 0; m < Nphases; ++m)
    for(size_t n = m; n < Nphases; ++n)
    {
        stringstream converter;
        converter << m << "_" << n;
        string counter = converter.str();
        Noise(m,n) = UserInterface::ReadParameterD(inp, moduleLocation, string("Noise_") + counter, false, 0.0);
        Noise(n,m) = Noise(m,n);

        CutOff(m,n) = UserInterface::ReadParameterD(inp, moduleLocation, string("CutOff_") + counter, false, CutOffDefault);
        if(CutOff(m,n) == 0.0)
        {
            stringstream message;
            message << "CutOff(" << m << ", " << n << ") = 0.0, please correct the driving force input!";
            Info::WriteExit(message.str(), thisclassname, "ReadInput()");
            exit(13);
        }
        CutOff(n,m) = CutOff(m,n);
    }

    Info::WriteLine();
    Info::WriteBlankLine();
}

void DrivingForce::Clear()
{
    ClearSR();
    if (Resolution == Resolutions::Double)
    {
        ClearDR();
    }
}

void DrivingForce::ClearSR()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,Raw.Bcells(),)
    {
        Raw(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DrivingForce::ClearDR()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,RawDR,RawDR.Bcells(),)
    {
        RawDR(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
/*
void DrivingForce::Refine(PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,RawDR,0,)
    {
        if(Phase.InterfaceDR(i,j,k))
        {
            RawDR(i,j,k) = Raw(i/2,j/2,k/2);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}*/

void DrivingForce::Refine(PhaseField& Phase)
{
    long int fx = 1 + dNx;
    long int fy = 1 + dNy;
    long int fz = 1 + dNz;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Fields(i,j,k).flag)
        for(int di = -dNx; di <= dNx; di+=2)
        for(int dj = -dNy; dj <= dNy; dj+=2)
        for(int dk = -dNz; dk <= dNz; dk+=2)
        {
            //RawDR(fx*i+(di+1)/2,fy*j+(dj+1)/2,fz*k+(dk+1)/2) = Raw.at(i+di*0.25,j+dj*0.25,k+dk*0.25);
            RawDR(fx*i+(di+1)/2,fy*j+(dj+1)/2,fz*k+(dk+1)/2) = Raw(i,j,k);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DrivingForce::SetBoundaryConditions(const BoundaryConditions& BC)
{
    if(dNx) BC.SetX(Raw);
    if(dNy) BC.SetY(Raw);
    if(dNz) BC.SetZ(Raw);
}

void DrivingForce::Average(const PhaseField& Phase, const BoundaryConditions& BC)
{
    if(Averaging)
    {
        SetBoundaryConditions(BC);
        CollectAverage(Phase);
        SetBoundaryConditions(BC);
        DistributeAverage(Phase);
    }
}

NodeAB DrivingForce::CalcUnified(const PhaseField& Phase, const BoundaryConditions& BC)
{
    NodeAB Unified;
    NodeAB Counter;
    //OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0, reduction(+:UnifiedDrivingForce),reduction(+:Counter))
    STORAGE_LOOP_BEGIN(i,j,k,Raw,0)
    {
        if (Phase.Interface(i,j,k))
        for(const auto& it : Raw(i,j,k))
        {
            const double PhiAlpha = Phase.Fields(i,j,k).get_value(it.indexA);
            const double PhiBeta  = Phase.Fields(i,j,k).get_value(it.indexB);
            const double weight   = std::sqrt(PhiAlpha*PhiBeta);

            Unified.add_asym1(it.indexA, it.indexB, it.value1*weight);
            Unified.add_asym2(it.indexA, it.indexB, it.value2*weight);
            Counter.add_sym1(it.indexA, it.indexB, weight);
        }
    }
    STORAGE_LOOP_END

    for(auto& it : Unified)
    {
        const double norm = Counter.get_sym1(it.indexA,it.indexB);
        if (norm > DBL_EPSILON)
        {
            it.value1 /= norm;
            it.value2 /= norm;
        }
    }
    return Unified;
}

void DrivingForce::Unify(const PhaseField& Phase, const BoundaryConditions& BC)
{
    if (Unifying)
    {
        NodeAB Unified = CalcUnified(Phase, BC);

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
        {
            if (Phase.Interface(i,j,k))
            {
                Raw(i,j,k) = Unified;
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        SetBoundaryConditions(BC);
    }
}

void DrivingForce::CollectAverage(const PhaseField& Phase)
{
    const int Xrange = min(Range, Nx-1)*dNx;
    const int Yrange = min(Range, Ny-1)*dNy;
    const int Zrange = min(Range, Nz-1)*dNz;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        if (Phase.Interface(i,j,k))
        for(auto it = Raw(i,j,k).begin();
                 it != Raw(i,j,k).end(); ++it)
        {
            double PhiAlpha = Phase.Fields(i,j,k).get_value(it->indexA);
            double PhiBeta  = Phase.Fields(i,j,k).get_value(it->indexB);
            const double PhiAlpha_PhiBeta = (PhiBeta > DBL_EPSILON) ? PhiAlpha/PhiBeta : DBL_MAX;

            if (PhiAlpha_PhiBeta > PhiThreshold/(1.0 - PhiThreshold) and
                PhiAlpha_PhiBeta < (1.0 - PhiThreshold)/PhiThreshold)
            {
                double value = 0.0;

                double SumWeights = 0.0;

                for(int ii = -Xrange; ii <= Xrange; ii++)
                for(int jj = -Yrange; jj <= Yrange; jj++)
                for(int kk = -Zrange; kk <= Zrange; kk++)
                if(Phase.Interface(i+ii,j+jj,k+kk))
                {
                    double dist = sqrt(ii*ii + jj*jj + kk*kk);

                    double locPhiAlphaValue = Phase.Fields(i+ii, j+jj, k+kk)[it->indexA];
                    double locPhiBetaValue = Phase.Fields(i+ii, j+jj, k+kk)[it->indexB];

                    double weight = sqrt(locPhiAlphaValue * locPhiBetaValue)*(Range - dist);
                    if (weight > DBL_EPSILON)
                    {
                        SumWeights += weight;
                        value += weight*Raw(i+ii, j+jj, k+kk).get_asym1(it->indexA, it->indexB);
                    }
                }
                if(SumWeights > 0.0) it->value2 = value/SumWeights;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DrivingForce::DistributeAverage(const PhaseField& Phase)
{
    const int Xrange = min(Range, Nx-1)*dNx;
    const int Yrange = min(Range, Ny-1)*dNy;
    const int Zrange = min(Range, Nz-1)*dNz;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        if (Phase.Interface(i,j,k))
        for(auto it = Raw(i,j,k).begin();
                 it != Raw(i,j,k).end(); ++it)
        {
            double PhiAlpha = Phase.Fields(i,j,k).get_value(it->indexA);
            double PhiBeta  = Phase.Fields(i,j,k).get_value(it->indexB);

            if (PhiAlpha*PhiBeta != 0.0)
            {
                double counter = 0.0;
                double value = 0.0;
                for(int ii = -Xrange; ii <= Xrange; ii++)
                for(int jj = -Yrange; jj <= Yrange; jj++)
                for(int kk = -Zrange; kk <= Zrange; kk++)
                if(Phase.Interface(i+ii,j+jj,k+kk))
                {
                    double dist = sqrt(ii*ii + jj*jj + kk*kk);
                    double locPhiAlpha = Phase.Fields(i+ii,j+jj,k+kk).get_value(it->indexA);
                    double locPhiBeta  = Phase.Fields(i+ii,j+jj,k+kk).get_value(it->indexB);
                    const double locPhiAlpha_locPhiBeta = (locPhiBeta > DBL_EPSILON) ? locPhiAlpha/locPhiBeta : DBL_MAX;

                    if (locPhiAlpha_locPhiBeta > PhiThreshold/(1.0 - PhiThreshold) and
                        locPhiAlpha_locPhiBeta < (1.0 - PhiThreshold)/PhiThreshold and (Range - dist) > 0.0)
                    {
                        counter ++;
                        value += Raw(i+ii, j+jj, k+kk).get_asym2(it->indexA, it->indexB);
                    }
                }
                if(counter > 0.0)
                {
                    it->value1 = value/counter;
                }
            }
            else if((!Phase.FieldsStatistics[it->indexA].Stage and
                     !(Phase.FieldsStatistics[it->indexA].MAXVolume == 0.0)) and
                    (!Phase.FieldsStatistics[it->indexB].Stage and
                     !(Phase.FieldsStatistics[it->indexB].MAXVolume == 0.0)))
            {
                it->value1 = 0.0;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DrivingForce::MergePhaseFieldIncrements(PhaseField& Phase, InterfaceProperties& IP)
{
    switch(Resolution)
    {
        case Resolutions::Single:
        {
            MergePhaseFieldIncrementsSR(Phase, IP);
            break;
        }
        case Resolutions::Double:
        {
            Refine(Phase);
            if (dGClearingMode == ClearingModes::Automatic) ClearSR();
            MergePhaseFieldIncrementsDR(Phase, IP);
            break;
        }
    }
}

void DrivingForce::MergePhaseFieldIncrementsSR(PhaseField& Phase,
                                               InterfaceProperties& IP)
{
    /** This function calculates a time independent phase-field increment by
    converting the local driving-forces. Additional noise-term can be applied
    and stability mechanism to strengthen the numerical phase-field profile. */

    double locMaxPsi = 0.0;
    long int locOvershootCounter = 0;
    Matrix<double> locMAXOvershootPOS = MAXOvershootPOS;
    Matrix<double> locMAXOvershootNEG = MAXOvershootNEG;
    Matrix<double> locMAXDrivingForcePOS = MAXDrivingForcePOS;
    Matrix<double> locMAXDrivingForceNEG = MAXDrivingForceNEG;

    /* If noise is applied to the driving force, noise distribution has to be
    initialized. */

    std::mt19937_64 NoiseGenerator(98310);
    std::uniform_real_distribution <double> NoiseDistribution(-1.0, 1.0);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Phase.Fields, 0, \
            reduction(+:locOvershootCounter) \
            reduction(MatrixDMAX:locMAXOvershootPOS, locMAXDrivingForcePOS) \
            reduction(MatrixDMIN:locMAXOvershootNEG, locMAXDrivingForceNEG) \
            reduction(max: locMaxPsi))
    {
        if (Phase.Interface(i,j,k))
        {
            double Prefactor = CalculatePrefactorSR(Phase,i,j,k);

            for(auto it = Raw(i,j,k).begin();
                     it != Raw(i,j,k).end(); ++it)
            {
                /* For each pair of grains the driving force is treated
                individually. */

                double locdG = it->value1;                                       // Local driving force of the grain pair.

                size_t indexA = it->indexA;
                size_t indexB = it->indexB;
                size_t pIndexA = Phase.FieldsStatistics[indexA].Phase;
                size_t pIndexB = Phase.FieldsStatistics[indexB].Phase;

                double alphaVal = Phase.Fields(i,j,k).get_value(indexA);
                double betaVal  = Phase.Fields(i,j,k).get_value(indexB);

                double norm = sqrt(alphaVal * betaVal);                         // Square root of the phase-fields of the local grain pair.

                if(Phase.FieldsStatistics[indexA].Stage == 2 or
                   Phase.FieldsStatistics[indexB].Stage == 2)
                {
                    /* In the case of a newly nucleated grain, a small value
                    is assigned to the normalization coefficient to allow the
                    initial growth of the nucleus. */

                    norm += 1.0e-6;
                }

                switch(NoiseMode)
                {
                    /* If noise is enabled, depending on the noise distribution
                    a random increment is added or subtracted from the local
                    driving force. Depending on the chosen mode, the noise
                    amplitude can be absolute or relative. */

                    case NoiseModes::Relative:
                    {
                        locdG *= (1.0 + Noise(pIndexA,pIndexB)
                                 *NoiseDistribution(NoiseGenerator));
                        break;
                    }
                    case NoiseModes::Absolute:
                    {
                        locdG += Noise(pIndexA,pIndexB)
                                *NoiseDistribution(NoiseGenerator);
                        break;
                    }
                    case NoiseModes::Zero:
                    default:
                    {
                        break;
                    }
                }

                if(dGcut)
                {
                    /* If the local driving force exceeds the interface stabilizing
                    force, the interface profile can distort leading to simulation
                    artifacts. Therefore the driving force is limited to a safe
                    value controlled by the user specified CutOff parameter.*/

                    double absDG = fabs(locdG);
                    double allowedDG = CutOff(pIndexA,pIndexB)*Prefactor
                                      *IP.InterfaceEnergy(pIndexA,pIndexB).MaxEnergy;

                    /* Collecting statistics for later output */
                    if(absDG > 0.4*allowedDG)
                    {
                        locOvershootCounter++;
                    }

                    double tmpMAXOvershoot  = locdG/allowedDG;

                    if(tmpMAXOvershoot > locMAXOvershootPOS(pIndexA,pIndexB))
                    {
                        locMAXOvershootPOS(pIndexA,pIndexB) = tmpMAXOvershoot;
                    }
                    if(tmpMAXOvershoot < locMAXOvershootNEG(pIndexA,pIndexB))
                    {
                        locMAXOvershootNEG(pIndexA,pIndexB) = tmpMAXOvershoot;
                    }
                    locMAXDrivingForcePOS(pIndexA,pIndexB) = max(locdG, locMAXDrivingForcePOS(pIndexA,pIndexB));
                    locMAXDrivingForceNEG(pIndexA,pIndexB) = min(locdG, locMAXDrivingForceNEG(pIndexA,pIndexB));
                    /* End collecting statistics for later output */
                    locdG = allowedDG*tanh(locdG/allowedDG);
                }

                double dPsi_dt = locdG*IP.get_mobility(i,j,k,indexA,indexB)*norm*Prefactor;

                locMaxPsi = max(fabs(dPsi_dt), locMaxPsi);

                Phase.FieldsDot(i,j,k).add_asym1(indexA,indexB,dPsi_dt);
            }
            if (dGClearingMode == ClearingModes::Automatic) Raw(i,j,k).clear();
        }
        else if (Phase.Fields(i,j,k).flag)
        {
            if (dGClearingMode == ClearingModes::Automatic) Raw(i,j,k).clear();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    MPI_Allreduce(&locMaxPsi, &maxPsi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    long int rlocOvershootCounter = locOvershootCounter;
    MPI_Allreduce(&rlocOvershootCounter, &locOvershootCounter, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    OvershootCounter  += locOvershootCounter;
    MPI_Allreduce(locMAXOvershootPOS.data(), MAXOvershootPOS.data(), Nphases*Nphases, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(locMAXOvershootNEG.data(), MAXOvershootNEG.data(), Nphases*Nphases, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(locMAXDrivingForcePOS.data(), MAXDrivingForcePOS.data(), Nphases*Nphases, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(locMAXDrivingForceNEG.data(), MAXDrivingForceNEG.data(), Nphases*Nphases, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

#else
    maxPsi = locMaxPsi;

    OvershootCounter   += locOvershootCounter;
    MAXOvershootPOS     = locMAXOvershootPOS;
    MAXOvershootNEG     = locMAXOvershootNEG;
    MAXDrivingForcePOS  = locMAXDrivingForcePOS;
    MAXDrivingForceNEG  = locMAXDrivingForceNEG;
#endif
}

void DrivingForce::MergePhaseFieldIncrementsDR(PhaseField& Phase,
                                               InterfaceProperties& IP)
{
    /** This function calculates a time independent phase-field increment by
    converting the local driving-forces. Additional noise-term can be applied
    and stability mechanism to strengthen the numerical phase-field profile. */

    double locMaxPsi = 0.0;
    long int locOvershootCounter = 0;
    Matrix<double> locMAXOvershootPOS = MAXOvershootPOS;
    Matrix<double> locMAXOvershootNEG = MAXOvershootNEG;
    Matrix<double> locMAXDrivingForcePOS = MAXDrivingForcePOS;
    Matrix<double> locMAXDrivingForceNEG = MAXDrivingForceNEG;
    /* If noise is applied to the driving force, noise distribution has to be
    initialized. */

    std::mt19937_64 NoiseGenerator(98310);
    std::uniform_real_distribution <double> NoiseDistribution(-1.0, 1.0);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Phase.FieldsDR, 0, \
            reduction(+:locOvershootCounter) \
            reduction(MatrixDMAX:locMAXOvershootPOS, locMAXDrivingForcePOS) \
            reduction(MatrixDMIN:locMAXOvershootNEG, locMAXDrivingForceNEG) \
            reduction(max: locMaxPsi))
    {
        if (Phase.InterfaceDR(i,j,k))
        {
            double Prefactor = CalculatePrefactorDR(Phase,i,j,k);

            for(auto it = RawDR(i,j,k).begin();
                     it != RawDR(i,j,k).end(); ++it)
            {
                /* For each pair of grains the driving force is treated
                individually. */

                double locdG = it->value1;                                      // Local driving force of the grain pair.

                size_t indexA = it->indexA;
                size_t indexB = it->indexB;
                size_t pIndexA = Phase.FieldsStatistics[indexA].Phase;
                size_t pIndexB = Phase.FieldsStatistics[indexB].Phase;

                double alphaVal = Phase.FieldsDR(i,j,k).get_value(indexA);
                double betaVal  = Phase.FieldsDR(i,j,k).get_value(indexB);

                double norm = sqrt(alphaVal * betaVal);                         // Square root of the phase-fields of the local grain pair.

                if(Phase.FieldsStatistics[indexA].Stage == 2 or
                   Phase.FieldsStatistics[indexB].Stage == 2)
                {
                    /* In the case of a newly nucleated grain, a small value
                    is assigned to the normalization coefficient to allow the
                    initial growth of the nucleus. */

                    norm += 1.0e-6;
                }

                switch(NoiseMode)
                {
                    /* If noise is enabled, depending on the noise distribution
                    a random increment is added or subtracted from the local
                    driving force. Depending on the chosen mode, the noise
                    amplitude can be absolute or relative. */

                    case NoiseModes::Relative:
                    {
                        locdG *= (1.0 + Noise(pIndexA,pIndexB)
                                 *NoiseDistribution(NoiseGenerator));
                        break;
                    }
                    case NoiseModes::Absolute:
                    {
                        locdG += Noise(pIndexA,pIndexB)
                                *NoiseDistribution(NoiseGenerator);
                        break;
                    }
                    case NoiseModes::Zero:
                    default:
                    {
                        break;
                    }
                }

                if(dGcut)
                {
                    /* If the local driving force exceeds the interface stabilizing
                    force, the interface profile can distort leading to simulation
                    artifacts. Therefore the driving force is limited to a safe
                    value controlled by the user specified CutOff parameter.*/

                    double absDG = fabs(locdG);
                    double allowedDG = CutOff(pIndexA,pIndexB)*Prefactor
                                      *IP.InterfaceEnergy(pIndexA,pIndexB).MaxEnergy;

                    /* Collecting statistics for later output */

                    if(absDG > 0.4*allowedDG)
                    {
                        locOvershootCounter++;
                    }

                    double tmpMAXOvershoot  = locdG/allowedDG;

                    if(tmpMAXOvershoot > locMAXOvershootPOS(pIndexA,pIndexB))
                    {
                        locMAXOvershootPOS(pIndexA,pIndexB) = tmpMAXOvershoot;
                    }
                    if(tmpMAXOvershoot < locMAXOvershootNEG(pIndexA,pIndexB))
                    {
                        locMAXOvershootNEG(pIndexA,pIndexB) = tmpMAXOvershoot;
                    }

                    locMAXDrivingForcePOS(pIndexA,pIndexB) = max(locdG, locMAXDrivingForcePOS(pIndexA,pIndexB));
                    locMAXDrivingForceNEG(pIndexA,pIndexB) = min(locdG, locMAXDrivingForceNEG(pIndexA,pIndexB));
                    /* End collecting statistics for later output */

                    locdG = allowedDG*tanh(locdG/allowedDG);
                }

                double dPsi_dt = locdG*IP.get_mobility_DR(i,j,k,indexA,indexB)*norm*Prefactor;

                locMaxPsi = max(fabs(dPsi_dt), locMaxPsi);

                Phase.FieldsDotDR(i,j,k).add_asym1(indexA,indexB,dPsi_dt);
            }
            RawDR(i,j,k).clear();
        }
        else if (Phase.FieldsDR(i,j,k).flag)
        {
            RawDR(i,j,k).clear();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    MPI_Allreduce(&locMaxPsi, &maxPsi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    long int rlocOvershootCounter = locOvershootCounter;
    MPI_Allreduce(&rlocOvershootCounter, &locOvershootCounter, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    OvershootCounter  += locOvershootCounter;
    MPI_Allreduce(locMAXOvershootPOS.data(), MAXOvershootPOS.data(), Nphases*Nphases, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(locMAXOvershootNEG.data(), MAXOvershootNEG.data(), Nphases*Nphases, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(locMAXDrivingForcePOS.data(), MAXDrivingForcePOS.data(), Nphases*Nphases, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(locMAXDrivingForceNEG.data(), MAXDrivingForceNEG.data(), Nphases*Nphases, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

#else
    maxPsi = locMaxPsi;

    OvershootCounter  += locOvershootCounter;
    MAXOvershootPOS    = locMAXOvershootPOS;
    MAXOvershootNEG    = locMAXOvershootNEG;
    MAXDrivingForcePOS = locMAXDrivingForcePOS;
    MAXDrivingForceNEG = locMAXDrivingForceNEG;
#endif
}

double DrivingForce::MaxTimeStep(PhaseField& Phase, InterfaceProperties& IP,
                     Settings& OPSettings, double TheorLimit, double NumerLimit)
{
    /** This function will calculate the maximum time step the phase field
    should be solved with.
    This function uses two different stability criteria. A theoretical stability
    function as 0.5*dx^2/(mu*Sigma) (if TheorLimit is set to 0.5), and a
    numerical criteria, where the maximum phase field increment in the whole
    simulation domain has to be lower than 1E-3 (if NumerLimit is set to 1E-3).
    If phase field gets unstable, decrease TheorLimit, if diffusion field
    becomes unstable due to too high phase field increments, decrease
    NumerLimit.*/

    double maxmu = IP.maxMu;
    double maxsg = IP.maxSigma;
    double dx2 = OPSettings.dx*OPSettings.dx;
    double maxTheorTimeStep = 0.0;
    double maxNumerTimeStep = 0.0;

    // Calculate theoretical maximum allowed time step
    if((maxmu > DBL_EPSILON) and (maxsg > DBL_EPSILON))
    {
        maxTheorTimeStep = TheorLimit*dx2/(maxmu*maxsg);
    }

    // Calculate maximum numerical allowed time step
    if(maxPsi > DBL_EPSILON)
    {
        maxNumerTimeStep = NumerLimit/maxPsi;
    }

    /* the equation would look like: maxNumerTimeStep=dt*NumerLimit/maxPsi
    if maxPsi would be the phase field increment. But as Psi it has to be
    multiplied with dt, it cancels out in this equation. So no dt needed!*/

    // Calculate total maximum allowed time step
    double maxTimeStep = min(maxTheorTimeStep, maxNumerTimeStep);

    return maxTimeStep;
}

void DrivingForce::PrintDiagnostics()
{
    if (OvershootCounter)
    {
        std::string message = "The driving force has been limited " +
                               std::to_string(OvershootCounter) + " times!\n";

        for(size_t n = 0; n < Nphases; n++)
        for(size_t m = n; m < Nphases; m++)
        if(MAXOvershootPOS(n, m) > DBL_EPSILON or MAXOvershootNEG(n, m) < -DBL_EPSILON or
           MAXOvershootPOS(m, n) > DBL_EPSILON or MAXOvershootNEG(m, n) < -DBL_EPSILON)
        {
            message += "   Phase pair (" + std::to_string(n) + ", " + std::to_string(m) + "): \n";
            if(n != m)
            {
                message += "       Max (positive) driving force value (overshoot ratio): "
                        + std::to_string(max(MAXDrivingForcePOS(n, m), -MAXDrivingForceNEG(m, n)))
                        + " ("
                        + std::to_string(max(MAXOvershootPOS(n, m), -MAXOvershootNEG(m, n)))
                        + ")\n";
                message += "       Min (negative) driving force value (overshoot ratio): "
                        + std::to_string(min(MAXDrivingForceNEG(n, m), -MAXDrivingForcePOS(m, n)))
                        + " ("
                        + std::to_string(min(MAXOvershootNEG(n, m), -MAXOvershootPOS(m, n)))
                        + ")\n";
            }
            else
            {
                message += "      Max (positive) driving force value (overshoot ratio): "
                        + std::to_string(MAXDrivingForcePOS(n, m))
                        + " ("
                        + std::to_string(MAXOvershootPOS(n, m))
                        + ")\n";
                message += "      Min (negative) driving force value (overshoot ratio): "
                        + std::to_string(MAXDrivingForceNEG(n, m))
                        + " ("
                        + std::to_string(MAXOvershootNEG(n, m))
                        + ")\n";
            }
        }
        Info::WriteStandard(thisclassname + "::PrintDiagnostics()", message);
        OvershootCounter = 0;
        for(size_t n = 0; n < Nphases; n++)
        for(size_t m = 0; m < Nphases; m++)
        {
            MAXOvershootPOS(n,m) = 0.0;
            MAXOvershootNEG(n,m) = 0.0;
            MAXDrivingForcePOS(n,m) = 0.0;
            MAXDrivingForceNEG(n,m) = 0.0;
        }
    }
}

void DrivingForce::PrintPointStatistics(const int x, const int y, const int z) const
{
    std::stringstream pointstat;
    pointstat << "DrivingForce Indices:\t";

    for (auto alpha = Raw(x,y,z).cbegin();
              alpha != Raw(x,y,z).cend(); ++alpha)
    {
        pointstat << alpha->indexA << ", " << alpha->indexB << "\t\t";
    }
    pointstat << endl;
    pointstat << "DrivingForce Values:\t";

    for (auto alpha = Raw(x,y,z).cbegin();
              alpha != Raw(x,y,z).cend(); ++alpha)
    {
        pointstat << alpha->value1 << "\t\t";
    }
    pointstat << endl;

    Info::WriteSimple(pointstat.str());
}

void DrivingForce::WriteVTK(const int tStep, const Settings& locSettings,
                            const size_t indexA, const size_t indexB,
                            const int precision) const
{
    stringstream converter;
    converter << indexA << "," << indexB;
    string phases = converter.str();

    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"dGavg(" + phases + ")", [indexA,indexB,this](int i,int j,int k){return Raw(i,j,k).get_asym1(indexA, indexB);}});
    ListOfFields.push_back((VTK::Field_t) {"dGraw(" + phases + ")", [indexA,indexB,this](int i,int j,int k){return Raw(i,j,k).get_asym2(indexA, indexB);}});
    std::string Filename = UserInterface::MakeFileName(locSettings.VTKDir, "DrivingForce_", tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void DrivingForce::WriteVTKforPhases(const int tStep, const Settings& locSettings,
                                     const PhaseField& Phase,
                                     const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;

    for(size_t alpha = 0; alpha < Phase.Nphases; alpha++)
    for(size_t beta = alpha; beta < Phase.Nphases; beta++)
    {
        stringstream converter;
        converter << alpha << "," << beta;
        string phases = converter.str();

        ListOfFields.push_back((VTK::Field_t) {"dGavg(" + phases + ")", [alpha,beta,&Phase,this](int i,int j,int k){
            double tempdG = 0.0;
            int counter = 0;
            for(auto it1 = Phase.Fields(i,j,k).cbegin();
                     it1 != Phase.Fields(i,j,k).cend(); ++it1)
            for(auto it2 = Phase.Fields(i,j,k).cbegin();
                     it2 != Phase.Fields(i,j,k).cend(); ++it2)
            if(it1 != it2 and
               Phase.FieldsStatistics[it1->index].Phase == alpha and
               Phase.FieldsStatistics[it2->index].Phase == beta)
            {
                tempdG += Raw(i,j,k).get_asym1(it1->index, it2->index);
                counter ++;
            }
            if(counter > 1.0)
            {
                tempdG /= counter;
            }
            return tempdG;
            }});
    }
    std::string Filename = UserInterface::MakeFileName(locSettings.VTKDir, "DrivingForcePhases_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void DrivingForce::Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)
{
    Raw.Reallocate(newNx, newNy, newNz);
    if(Resolution == Resolutions::Double)
    {
        RawDR.Reallocate((1+dNx)*newNx, (1+dNy)*newNy, (1+dNz)*newNz);
    }

    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    Info::WriteStandard(thisclassname, "Remeshed");
}

DrivingForce& DrivingForce::operator= (const DrivingForce& rhs)
{
    // protect against invalid self-assignment and copy of uninitialized object
    if (this != &rhs and rhs.thisclassname == "DrivingForce")
    {
        thisclassname = rhs.thisclassname;

        TotalNx = rhs.TotalNx;
        OffsetX = rhs.OffsetX;
        TotalNy = rhs.TotalNy;
        OffsetY = rhs.OffsetY;
        TotalNz = rhs.TotalNz;
        OffsetZ = rhs.OffsetZ;

        Nx = rhs.Nx;
        Ny = rhs.Ny;
        Nz = rhs.Nz;

        dNx = rhs.dNx;
        dNy = rhs.dNy;
        dNz = rhs.dNz;

        Nphases = rhs.Nphases;
        Range = rhs.Range;
        PhiThreshold = rhs.PhiThreshold;

        OvershootCounter = rhs.OvershootCounter;
        MAXOvershootPOS = rhs.MAXOvershootPOS;
        MAXOvershootNEG = rhs.MAXOvershootNEG;
        MAXDrivingForcePOS = rhs.MAXDrivingForcePOS;
        MAXDrivingForceNEG = rhs.MAXDrivingForceNEG;

        CutOff = rhs.CutOff;
        Averaging = rhs.Averaging;
        Resolution = rhs.Resolution;

        if (Raw.IsNotAllocated())
        {
            Raw.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, rhs.Raw.Bcells());
        }
        else if (not Raw.IsSize(rhs.Nx, rhs.Ny, rhs.Nz))
        {
            Raw.Reallocate(Nx, Ny, Nz);
        }

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,Raw.Bcells(),)
        {
            Raw(i,j,k) = rhs.Raw(i,j,k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    return *this;
}

NodeAB DrivingForce::Raw_at(const double x, const double y, const double z) const
{
#ifdef DEBUG
    if(x > Nx + Raw.Bcells()*dNx - 1 or
       y > Ny + Raw.Bcells()*dNy - 1 or
       z > Nz + Raw.Bcells()*dNz - 1 or
       x < -Raw.Bcells()*dNx or
       y < -Raw.Bcells()*dNy or
       z < -Raw.Bcells()*dNz)
    {
        std::stringstream message;
        message << "ERROR: DrivingForce::Raw_at()\n"
                << "Access beyond storage range -> ("
                << x << "," << y << "," << z << ")" << " is outside of storage bounds ["
                << -Raw.Bcells()*dNx << ", " << -Raw.Bcells()*dNy << ", " << -Raw.Bcells()*dNz << "] and ("
                << Nx + Raw.Bcells()*dNx << ", " << Ny + Raw.Bcells()*dNy << ", " << Nz + Raw.Bcells()*dNz << ")\n"
                << "Terminating!!!\n";
        throw std::logic_error(message.str());
    }
#endif

    long int x0 = floor(x)*dNx;
    long int y0 = floor(y)*dNy;
    long int z0 = floor(z)*dNz;
    double dx = fabs(x - x0)*dNx;
    double dy = fabs(y - y0)*dNy;
    double dz = fabs(z - z0)*dNz;

    NodeAB loc_dG;
    double sum_of_weights = 0.0;
    if(Raw(x0    ,y0    ,z0    ).size())
    {
        loc_dG.add_asym1(Raw(x0    ,y0    ,z0    )*((1.0 - dx)*(1.0 - dy)*(1.0 - dz)));
        sum_of_weights += ((1.0 - dx)*(1.0 - dy)*(1.0 - dz));
    }
    if(Raw(x0+dNx,y0    ,z0    ).size())
    {
        loc_dG.add_asym1(Raw(x0+dNx,y0    ,z0    )*(dx*(1.0 - dy)*(1.0 - dz)));
        sum_of_weights += (dx*(1.0 - dy)*(1.0 - dz));
    }
    if(Raw(x0    ,y0+dNy,z0    ).size())
    {
        loc_dG.add_asym1(Raw(x0    ,y0+dNy,z0    )*((1.0 - dx)*dy*(1.0 - dz)));
        sum_of_weights += ((1.0 - dx)*dy*(1.0 - dz));
    }
    if(Raw(x0    ,y0    ,z0+dNz).size())
    {
        loc_dG.add_asym1(Raw(x0    ,y0    ,z0+dNz)*((1.0 - dx)*(1.0 - dy)*dz));
        sum_of_weights += ((1.0 - dx)*(1.0 - dy)*dz);
    }
    if(Raw(x0+dNx,y0+dNy,z0    ).size())
    {
        loc_dG.add_asym1(Raw(x0+dNx,y0+dNy,z0    )*(dx*dy*(1.0 - dz)));
        sum_of_weights += (dx*dy*(1.0 - dz));
    }
    if(Raw(x0+dNx,y0    ,z0+dNz).size())
    {
        loc_dG.add_asym1(Raw(x0+dNx,y0    ,z0+dNz)*(dx*(1.0 - dy)*dz));
        sum_of_weights += (dx*(1.0 - dy)*dz);
    }
    if(Raw(x0    ,y0+dNy,z0+dNz).size())
    {
        loc_dG.add_asym1(Raw(x0    ,y0+dNy,z0+dNz)*((1.0 - dx)*dy*dz));
        sum_of_weights += ((1.0 - dx)*dy*dz);
    }
    if(Raw(x0+dNx,y0+dNy,z0+dNz).size())
    {
        loc_dG.add_asym1(Raw(x0+dNx,y0+dNy,z0+dNz)*(dx*dy*dz));
        sum_of_weights += (dx*dy*dz);
    }
    if(sum_of_weights >= DBL_EPSILON)
    {
        loc_dG *= 1.0/sum_of_weights;
    }
    return loc_dG;
}

}// namespace openphase

