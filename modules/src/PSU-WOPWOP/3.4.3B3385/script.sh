#!/bin/sh

#PBS -N PSU_WOPWOP 
#PBS -e r_wopwop.err
#PBS -o r_wopwop.out
#PBS -q lionxm-aero
#PBS -l nodes=5:myrinet:ppn=2
#PBS -l walltime=48:00:00

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

/usr/global/bin/mmpirun ~/wopwop3/wopwop3 > OUTPUT.txt
#PBS -q debug
