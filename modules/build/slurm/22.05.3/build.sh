#!/bin/bash

source ./module.sh

if [ ! -f "slurm-${VERSION}.tar.bz2" ]; then
		wget https://download.schedmd.com/slurm/slurm-22.05.3.tar.bz2
		tar -xvf slurm-$VERSION.tar.bz2
fi

cd $SRC_DIR

./configure --prefix=/share/software/slurm/22.05.3/ucdhpc-20.04/ \
            --enable-salloc-kill-cmd \
            --with-pmix=/share/software/pmix/3.2.3/ucdhpc-20.04/ \
            --with-hwloc=/share/software/hwloc/2.4.1/ucdhpc-20.04/ \
            --with-hdf5=yes \
            --with-ucx=/share/software/ucx/1.9.0/ucdhpc-20.04/ \
            --with-bpf=yes \
            --with-munge=yes
make -j 32
