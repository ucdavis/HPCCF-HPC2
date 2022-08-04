#!/bin/bash

source ./module.sh

if [ ! -f "${SRC_DIR}.tar.gz" ]; then
		wget https://ftp.gromacs.org/gromacs/gromacs-$VERSION.tar.gz
		tar -xvzf gromacs-$VERSION.tar.gz
fi

cd $SRC_DIR
mkdir -p build
cd build

module load openmpi/4.1.0

cmake .. -DGMX_MPI=on -DGMX_HWLOC=ON -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}
make -j 16

