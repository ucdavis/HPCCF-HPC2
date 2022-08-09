#!/bin/bash
set -e

source ./module.sh

mkdir -p $INSTALL_PREFIX
mkdir -p `dirname $MODULEFILE_INSTALL_DEST`
cp $MODULEFILE_SOURCE $MODULEFILE_INSTALL_DEST

echo "# NOTE: Now you must switch to the conda3 user and create"
echo "  a conda environment named MemSurfer-$VERSION using the "
echo "  environment.yml file in this directory: conda env create -f environment.yml"
echo "  Then, run install_conda.sh as the conda3 user."


