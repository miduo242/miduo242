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
 *   File created :   2014
 *   Main contributors :   Oleg Shchyglo; Raphael Schiedung
 *
 */

#include "Base/CommonFunctions.h"
#include "Base/EulerAngles.h"
#include "FluidDynamics/InteractionSolidFluid.h"
#include "GrainInfo.h"
#include "Info.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Velocities.h"

namespace openphase
{
void InteractionSolidFluid::CollectGrainsStatistics(PhaseField& Phase, const BoundaryConditions& BC, const Settings& locSettings)
{
    if (BC.BC0X == BoundaryConditionTypes::Periodic and
        BC.BCNX == BoundaryConditionTypes::Periodic and
        BC.BC0Y == BoundaryConditionTypes::Periodic and
        BC.BCNY == BoundaryConditionTypes::Periodic and
        BC.BC0Z == BoundaryConditionTypes::Periodic and
        BC.BCNX == BoundaryConditionTypes::Periodic)
    {
       CalculateCenterOfMassWithPeriodicBoundaryConditions(Phase);
    }
    //TODO
    //else if ((BC.BC0X == BoundaryConditionTypes::Periodic) or (BC.BCNX == BoundaryConditionTypes::Periodic) or
    //         (BC.BC0Y == BoundaryConditionTypes::Periodic) or (BC.BCNY == BoundaryConditionTypes::Periodic) or
    //         (BC.BC0Z == BoundaryConditionTypes::Periodic) or (BC.BCNX == BoundaryConditionTypes::Periodic) or
    //         (BC.BC0X == BoundaryConditionTypes::NoFlux  ) or (BC.BCNX == BoundaryConditionTypes::NoFlux  ) or
    //         (BC.BC0Y == BoundaryConditionTypes::NoFlux  ) or (BC.BCNY == BoundaryConditionTypes::NoFlux  ) or
    //         (BC.BC0Z == BoundaryConditionTypes::NoFlux  ) or (BC.BCNX == BoundaryConditionTypes::NoFlux  ))
    //{
    //    //Info::WriteWarning("Chosen boundary conditions are not yet entirety supported!", "Interaction Solid Fluid", "CollectGrainsStatistics");
    //    CalculateCenterOfMass(Phase);
    //}
    else
    {
        CalculateCenterOfMass(Phase);
    }

    CollectGrainsStatisticsStepTwo(Phase, BC, locSettings);
}

void InteractionSolidFluid::CalculateCenterOfMass(PhaseField& Phase)
{

    for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    {
        Phase.FieldsStatistics[idx].Rcm.set_to_zero();
    }

    size_t Nthreads = 1;
    #ifdef _OPENMP
    Nthreads = omp_get_max_threads();
    #endif

    const size_t size = Phase.FieldsStatistics.size();
    std::vector<GrainInfo> locFieldsStatistics(Nthreads);
    for (size_t t = 0; t < Nthreads; t++)
    {
        locFieldsStatistics[t].Allocate(size);
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        size_t thread = 0;
        #ifdef _OPENMP
        thread = omp_get_thread_num();
        #endif

        for(auto it = Phase.Fields(i,j,k).cbegin();
                 it != Phase.Fields(i,j,k).cend(); ++it)
        if(it->value != 0.0)
        {
            locFieldsStatistics[thread][it->index].Rcm[0] += i*it->value;
            locFieldsStatistics[thread][it->index].Rcm[1] += j*it->value;
            locFieldsStatistics[thread][it->index].Rcm[2] += k*it->value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    {
        if (Phase.FieldsStatistics[idx].Volume != 0.0)
        {
            for (auto it : locFieldsStatistics) Phase.FieldsStatistics[idx].Rcm += it[idx].Rcm;

            #ifdef MPI_PARALLEL
            auto local_Rcm = Phase.FieldsStatistics[idx].Rcm;
            MPI_Allreduce(&local_Rcm[0], &(Phase.FieldsStatistics[idx].Rcm[0]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&local_Rcm[1], &(Phase.FieldsStatistics[idx].Rcm[1]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&local_Rcm[2], &(Phase.FieldsStatistics[idx].Rcm[2]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            #endif

            Phase.FieldsStatistics[idx].Rcm *= 1.0/double(Phase.FieldsStatistics[idx].Volume);
        }
        else
        {
            Phase.FieldsStatistics[idx].Rcm.set_to_zero();
        }
    }
}

void InteractionSolidFluid::CalculateCenterOfMassWithPeriodicBoundaryConditions(PhaseField& Phase)
{
    int Nx = Phase.Nx;
    int Ny = Phase.Ny;
    int Nz = Phase.Nz;
    GrainInfo gloFieldsStatistics1;
    GrainInfo gloFieldsStatistics2;
    gloFieldsStatistics1.Allocate(Phase.FieldsStatistics.size());
    gloFieldsStatistics2.Allocate(Phase.FieldsStatistics.size());

    for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    {
        gloFieldsStatistics1[idx].Rcm.set_to_zero();
        gloFieldsStatistics2[idx].Rcm.set_to_zero();
    }

    size_t Nthreads = 1;
    #ifdef _OPENMP
    Nthreads = omp_get_max_threads();
    #endif

    const size_t size = Phase.FieldsStatistics.size();
    std::vector<GrainInfo> locFieldsStatistics1(Nthreads);
    std::vector<GrainInfo> locFieldsStatistics2(Nthreads);
    for (auto& it : locFieldsStatistics1) it.Allocate(size);
    for (auto& it : locFieldsStatistics2) it.Allocate(size);

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
        {
            size_t thread = 0;
            #ifdef _OPENMP
            thread = omp_get_thread_num();
            #endif

            for(auto it = Phase.Fields(i,j,k).cbegin();
                     it != Phase.Fields(i,j,k).cend(); ++it)
            if(it->value != 0.0)
            {
                double theta = (double)i/(double)Nx*2.0*Pi;
                locFieldsStatistics1[thread][it->index].Rcm[0] += Nx/(2.0*Pi)*cos(theta)*it->value;
                locFieldsStatistics2[thread][it->index].Rcm[0] += Nx/(2.0*Pi)*sin(theta)*it->value;
                theta = (double)j/(double)Ny*2.0*Pi;
                locFieldsStatistics1[thread][it->index].Rcm[1] += Ny/(2.0*Pi)*cos(theta)*it->value;
                locFieldsStatistics2[thread][it->index].Rcm[1] += Ny/(2.0*Pi)*sin(theta)*it->value;
                theta = (double)k/(double)Nz*2.0*Pi;
                locFieldsStatistics1[thread][it->index].Rcm[2] += Nz/(2.0*Pi)*cos(theta)*it->value;
                locFieldsStatistics2[thread][it->index].Rcm[2] += Nz/(2.0*Pi)*sin(theta)*it->value;
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
        {
            for (auto it : locFieldsStatistics1) gloFieldsStatistics1[idx].Rcm += it[idx].Rcm;
            for (auto it : locFieldsStatistics2) gloFieldsStatistics2[idx].Rcm += it[idx].Rcm;

            if(Phase.FieldsStatistics[idx].Volume > 0.0)
            {
                gloFieldsStatistics1[idx].Rcm *= 1.0/double(Phase.FieldsStatistics[idx].Volume);
                gloFieldsStatistics2[idx].Rcm *= 1.0/double(Phase.FieldsStatistics[idx].Volume);
            }

            Phase.FieldsStatistics[idx].Rcm[0] = atan2(-gloFieldsStatistics2[idx].Rcm[0],-gloFieldsStatistics1[idx].Rcm[0]) + Pi;
            Phase.FieldsStatistics[idx].Rcm[0] = Nx*Phase.FieldsStatistics[idx].Rcm[0]/(2.0*Pi);
            Phase.FieldsStatistics[idx].Rcm[1] = atan2(-gloFieldsStatistics2[idx].Rcm[1],-gloFieldsStatistics1[idx].Rcm[1]) + Pi;
            Phase.FieldsStatistics[idx].Rcm[1] = Nx*Phase.FieldsStatistics[idx].Rcm[1]/(2.0*Pi);
            Phase.FieldsStatistics[idx].Rcm[2] = atan2(-gloFieldsStatistics2[idx].Rcm[2],-gloFieldsStatistics1[idx].Rcm[2]) + Pi;
            Phase.FieldsStatistics[idx].Rcm[2] = Nx*Phase.FieldsStatistics[idx].Rcm[2]/(2.0*Pi);

            #ifdef MPI_PARALLEL
            auto local_Rcm = Phase.FieldsStatistics[idx].Rcm;
            MPI_Allreduce(&local_Rcm[0], &(Phase.FieldsStatistics[idx].Rcm[0]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&local_Rcm[1], &(Phase.FieldsStatistics[idx].Rcm[1]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&local_Rcm[2], &(Phase.FieldsStatistics[idx].Rcm[2]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            #endif
        }
}

void InteractionSolidFluid::CollectGrainsStatisticsStepTwo(PhaseField& Phase, const BoundaryConditions& BC, const Settings& locSettings)
{
    for(unsigned int idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    {
        Phase.FieldsStatistics[idx].InertiaM.set_to_zero();
    }

    size_t Nthreads = 1;
    #ifdef _OPENMP
    Nthreads = omp_get_max_threads();
    #endif

    std::vector<GrainInfo> locFieldsStatistics(Nthreads);
    for (size_t t = 0; t < Nthreads; t++) locFieldsStatistics[t].Allocate(Phase.FieldsStatistics.size());

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        size_t thread = 0;
        #ifdef _OPENMP
        thread = omp_get_thread_num();
        #endif

        for(auto it = Phase.Fields(i,j,k).cbegin();
                 it != Phase.Fields(i,j,k).cend(); ++it)
        if(it->value != 0.0)
        {
            dVector3 locRad = CommonFunctions::Distance(dVector3({double(i),double(j),double(k)}),Phase.FieldsStatistics[it->index].Rcm, BC, locSettings);
            const double locRadX = locRad[0]*Phase.dx;
            const double locRadY = locRad[1]*Phase.dx;
            const double locRadZ = locRad[2]*Phase.dx;

            locFieldsStatistics[thread][it->index].InertiaM(0,0) += it->value*(locRadY*locRadY + locRadZ*locRadZ);
            locFieldsStatistics[thread][it->index].InertiaM(1,1) += it->value*(locRadX*locRadX + locRadZ*locRadZ);
            locFieldsStatistics[thread][it->index].InertiaM(2,2) += it->value*(locRadX*locRadX + locRadY*locRadY);

            locFieldsStatistics[thread][it->index].InertiaM(0,1) -= it->value*(locRadX*locRadY);
            locFieldsStatistics[thread][it->index].InertiaM(1,0) -= it->value*(locRadY*locRadX);
            locFieldsStatistics[thread][it->index].InertiaM(0,2) -= it->value*(locRadX*locRadZ);
            locFieldsStatistics[thread][it->index].InertiaM(2,0) -= it->value*(locRadZ*locRadX);
            locFieldsStatistics[thread][it->index].InertiaM(1,2) -= it->value*(locRadY*locRadZ);
            locFieldsStatistics[thread][it->index].InertiaM(2,1) -= it->value*(locRadZ*locRadY);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    if(Phase.FieldsStatistics[idx].Exist)
    {
        for(auto it : locFieldsStatistics) Phase.FieldsStatistics[idx].InertiaM += it [idx].InertiaM;

        #ifdef MPI_PARALLEL
        auto locInertiaM = Phase.FieldsStatistics[idx].InertiaM;
        MPI_Allreduce(&locInertiaM(0,0), &(Phase.FieldsStatistics[idx].InertiaM(0,0)), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locInertiaM(1,1), &(Phase.FieldsStatistics[idx].InertiaM(1,1)), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locInertiaM(2,2), &(Phase.FieldsStatistics[idx].InertiaM(2,2)), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        MPI_Allreduce(&locInertiaM(0,1), &(Phase.FieldsStatistics[idx].InertiaM(0,1)), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locInertiaM(1,0), &(Phase.FieldsStatistics[idx].InertiaM(1,0)), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locInertiaM(0,2), &(Phase.FieldsStatistics[idx].InertiaM(0,2)), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locInertiaM(2,0), &(Phase.FieldsStatistics[idx].InertiaM(2,0)), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locInertiaM(1,2), &(Phase.FieldsStatistics[idx].InertiaM(1,2)), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locInertiaM(2,1), &(Phase.FieldsStatistics[idx].InertiaM(2,1)), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        #endif
    }

    double dx = Phase.dx;
    double dV = dx*dx*dx;
    for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    {
        Phase.FieldsStatistics[idx].InertiaM *= dV*Phase.FieldsStatistics[idx].Density;
    }
}

size_t InteractionSolidFluid::CalculateSolidVelocities(PhaseField& Phase,
       Velocities& Vel, const BoundaryConditions& BC,
       const Settings& locSettings, const double dt, double* VLimit,
       bool EnforceZeroTotalMomentum,
       bool EnforceZeroTotalAngularMomentum)
{
    const double dx = Phase.dx;
    const double dV = dx*dx*dx;

    double stdVLimit = dx/dt;
    if (VLimit == nullptr) VLimit = &stdVLimit;
    size_t VelocityLimitApplied = 0;

    CollectGrainsStatistics(Phase, BC, locSettings);
    for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    {
        Grain& grain = Phase.FieldsStatistics[idx];
        if(grain.Exist and grain.State == AggregateStates::Solid and
                grain.Mobile and grain.Volume > DBL_EPSILON)
        {
            double SolidMass_1 = 1.0/(grain.Volume*dV*grain.Density);
            grain.Acm  += grain.Force * SolidMass_1;
            grain.Vcm  += grain.Acm * dt;

            if (grain.Vcm[0] > *VLimit or
                grain.Vcm[1] > *VLimit or
                grain.Vcm[2] > *VLimit)
            {
                VelocityLimitApplied++;
                grain.Vcm /= grain.Vcm.abs();
                grain.Vcm *= *VLimit;
            }

            if (grain.Volume > 10)
            {
                grain.aAcc += grain.InertiaM.inverted() * grain.Torque;
                grain.aVel += grain.aAcc * dt;
                EulerAngles locAngles({grain.aVel[0] * dt,
                                       grain.aVel[1] * dt,
                                       grain.aVel[2] * dt}, XYZ);
                grain.Orientation += locAngles.getQuaternion();
            }
        }
        else
        {
            grain.Acm  = {0,0,0};
            grain.Vcm  = {0,0,0};
            grain.aAcc = {0,0,0};
            grain.aVel = {0,0,0};
        }
    }

    if (EnforceZeroTotalMomentum)
    {
        dVector3 Momentum = {0.0,0.0,0.0};
        double   Mass     = 0.0;

        for (size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
        if (Phase.FieldsStatistics[idx].State == AggregateStates::Solid and
            Phase.FieldsStatistics[idx].Mobile)
        {
            const double GrainMass = Phase.FieldsStatistics[idx].Volume*dV*
                                     Phase.FieldsStatistics[idx].Density;
            Mass     += GrainMass;
            Momentum += Phase.FieldsStatistics[idx].Vcm * GrainMass ;
        }

        if (Mass > DBL_EPSILON)
        {
            const dVector3 VcmFix = Momentum/Mass;
            for (size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
            if (Phase.FieldsStatistics[idx].State == AggregateStates::Solid and
                Phase.FieldsStatistics[idx].Mobile)
            {
                Phase.FieldsStatistics[idx].Vcm -= VcmFix;
            }
        }

        if (EnforceZeroTotalAngularMomentum)
        {
            // TODO consider Avel
            dVector3 AMomentum = {0.0,0.0,0.0};
            double   Mass      = 0.0;
            double   Radius    = 0.0;

            for (size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
            if (Phase.FieldsStatistics[idx].State == AggregateStates::Solid and
                Phase.FieldsStatistics[idx].Mobile)
            {
                const double   GrainMass = Phase.FieldsStatistics[idx].Volume*dV*
                                           Phase.FieldsStatistics[idx].Density;
                const dVector3 Momentum  = Phase.FieldsStatistics[idx].Vcm * GrainMass;
                Mass      += GrainMass;
                AMomentum += Phase.FieldsStatistics[idx].Rcm.cross(Momentum);
                Radius    += Phase.FieldsStatistics[idx].Rcm.abs();
            }

            if ((Radius*Mass > DBL_EPSILON) and (AMomentum.abs() > DBL_EPSILON))
            {
                const double   Fix  = AMomentum.abs()/(Radius*Mass);
                const dVector3 e_AM = AMomentum/AMomentum.abs();

                for (size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
                if (Phase.FieldsStatistics[idx].State == AggregateStates::Solid and
                    Phase.FieldsStatistics[idx].Mobile)
                if (Phase.FieldsStatistics[idx].Rcm.abs() > DBL_EPSILON)
                {
                    const dVector3 e_R    = Phase.FieldsStatistics[idx].Rcm/
                                            Phase.FieldsStatistics[idx].Rcm.abs();
                    const dVector3 e_Fix  = e_AM.cross(e_R);

                    Phase.FieldsStatistics[idx].Vcm -= e_Fix * Fix;
                }
            }
        }
    }

    if (VelocityLimitApplied)
    {
        std::stringstream message;
        message << "Solid velocity limit as applied!";
        Info::WriteWarning(message.str(), "InteractionSolidFluid", "CalculateSolidVelocities");
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    for(auto alpha = Phase.Fields(i,j,k).cbegin();
             alpha != Phase.Fields(i,j,k).cend(); ++alpha)
    {
        Grain& grain = Phase.FieldsStatistics[alpha->index];
        if (grain.State == AggregateStates::Solid and grain.Mobile)
        {
            dMatrix3x3 W;
            W(0,0) = 0.0;
            W(1,1) = 0.0;
            W(2,2) = 0.0;
            W(0,1) = -grain.aVel[2];
            W(0,2) =  grain.aVel[1];
            W(1,2) = -grain.aVel[0];
            W(1,0) =  grain.aVel[2];
            W(2,0) = -grain.aVel[1];
            W(2,1) =  grain.aVel[0];

            const dVector3 pos = {double(i), double(j), double(k)};
            dVector3 distanceCM;
            CommonFunctions::CalculateDistancePeriodic(pos, grain.Rcm, distanceCM, Phase.Nx, Phase.Ny, Phase.Nz);
            const dVector3 locR = distanceCM * dx;

            Vel.Phase(i,j,k)({grain.Phase}) = grain.Vcm + W*locR;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Vel.CalculateAverage(Phase);

    return VelocityLimitApplied;
}

} //namespace openphase
