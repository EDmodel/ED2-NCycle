#!/bin/bash
#PBS -q medvigy
#PBS -N ed2_postp
#PBS -l nodes=1:ppn=1,walltime=160:00:00
#
module load openmpi
cd $PBS_O_WORKDIR
mpirun -np 1 -hostfile $PBS_NODEFILE ./timeseries-M
