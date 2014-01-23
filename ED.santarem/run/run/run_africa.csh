#!/bin/bash
#PBS -q medvigy
#PBS -N dd4.lat5lon9
#PBS -l nodes=1:ppn=1,walltime=160:00:00
#
module load openmpi
cd $PBS_O_WORKDIR
mpirun -np 1 -hostfile $PBS_NODEFILE ./ed2 -f ED2IN.dd4 -lat 05 -lon 09
