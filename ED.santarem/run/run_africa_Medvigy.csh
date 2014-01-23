
#!/bin/bash
#PBS -N Medvigy.ED2INsappd04cycFix
#PBS -l nodes=1:ppn=1,walltime=24:00:00
#
module load openmpi

cd $PBS_O_WORKDIR
mpirun -np 1 -hostfile $PBS_NODEFILE ./ed2 -f  ED2IN.medvigy2




