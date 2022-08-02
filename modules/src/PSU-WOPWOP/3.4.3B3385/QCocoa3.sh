#!/bin/sh

#PBS -N WOPWOP3.1
#PBS -e r_wopwop.err
#PBS -o r_wopwop.out
#PBS -l nodes=15:ppn=2

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
echo This job has allocated $NPROCS nodes

/usr/local/mpich-1.2.5/bin/mpirun -np $NPROCS -machinefile $PBS_NODEFILE ~/WOPWOP/wopwop3 > OUTPUT.txt
#PBS -q debug