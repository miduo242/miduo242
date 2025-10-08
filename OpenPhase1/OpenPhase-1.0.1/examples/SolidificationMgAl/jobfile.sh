#!/bin/sh

#SBATCH --partition=compute
#SBATCH --output=job%j.out
#SBATCH --error=job%j.e
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --cpu_bind=cores
#SBATCH --threads-per-core=1
#SBATCH --hint=nomultithread
#SBATCH --time=8-00:00:00
#SBATCH --job-name=Mart_C03
#SBATCH --mail-type=ALL
#SBATCH --mail-user=oleg.shchyglo@rub.de

./SolidificationMgAl
