#!/bin/bash
set -e

source ./module.sh

git clone --recursive git@github.com:LLNL/MemSurfer.git
cp setup.py $SRC_DIR/

cd $SRC_DIR

module load conda3/4.13.0
source activate MemSurfer-$VERSION

CC=`which gcc` CXX=`which g++` LDCXXSHARED="`which g++` -bundle -undefined dynamic_lookup" python setup.py bdist_wheel
