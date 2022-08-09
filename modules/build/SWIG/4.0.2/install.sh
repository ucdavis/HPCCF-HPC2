#!/bin/bash

source ./module.sh

mkdir -p $INSTALL_PREFIX
mkdir -p `dirname $MODULEFILE_INSTALL_DEST`
cp $MODULEFILE_SOURCE $MODULEFILE_INSTALL_DEST

cd $SRC_DIR
make install

