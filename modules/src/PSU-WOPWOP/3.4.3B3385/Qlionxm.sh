#!/bin/sh

#PBS -N PSU_WOPWOP
#PBS -e wopwop.err
#PBS -o wopwop.out
#PBS -l nodes=1:myrinet:ppn=2
#PBS -q lionxm-aero
#PBS -l walltime=3:00:00
# This job's working directory
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo This job runs on the following processors:
echo `cat $PBS_NODEFILE`
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

/usr/global/bin/mmpirun ./wopwop3 > OUTPUT.txt

echo Job Complete

