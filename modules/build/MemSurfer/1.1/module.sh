#!/bin/bash

source ../../../variables.sh

MODULE_NAME=MemSurfer
VERSION="1.1"

SRC_DIR=MemSurfer
INSTALL_PREFIX=$MODULES_ROOT/${MODULE_NAME}/${VERSION}/$SYSTEM_NAME
MODULEFILE_SOURCE=../../../modulefiles/${MODULE_NAME}/${VERSION}
MODULEFILE_INSTALL_DEST=$MODULEFILES_PATH/$MODULE_NAME/${VERSION}

module load SWIG/4.0.2
module load CGAL/4.13

export BOOST_ROOT=/usr/
export CGAL_ROOT=$CGAL_HOME
export EIGEN_ROOT=/usr/

echo "Using INSTALL_PREFIX=$INSTALL_PREFIX"
echo "Modulefile source: $MODULEFILE_SOURCE"
echo "Modulefile destination: $MODULEFILE_INSTALL_DEST"
