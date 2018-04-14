#!/bin/csh
#SBATCH --time=00:10:00
#SBATCH -p main
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16

srun ./heat


