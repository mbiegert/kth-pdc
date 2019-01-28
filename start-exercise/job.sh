#!/bin/bash -l

# job name
#SBATCH -J myjob
# account
#SBATCH -A edu19.SF2568
# email notification
#SBATCH --mail-type=BEGIN,END
# 10 minutes wall-clock time will be given to this job
#SBATCH -t 00:10:00
# Number of nodes
#SBATCH --nodes=2
# set tasks per node to 24 in order to disable hyperthreading
#SBATCH --ntasks-per-node=24

module add i-compilers intelmpi

mpirun -np 5 ./hello

