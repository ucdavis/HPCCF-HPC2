#!/bin/sh

#PBS -N WOPWOP3
#PBS -e errorfile
#PBS -o outputfile
#PBS -l nodes=10:ppn=2
#PBS -l walltime=8:00:00

# This job's working directory
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo This jobs runs on the following processors:
echo `cat $PBS_NODEFILE`
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS processors.

/usr/local/mvapich_intel_10-1.0.1/bin/mpirun -np $NPROCS -machinefile $PBS_NODEFILE ./wopwop3 > OUTPUT.txt
