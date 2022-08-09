#!/bin/bash

source ./module.sh

if [ ! -f "swig-${VERSION}.tar.gz" ]; then
		wget http://prdownloads.sourceforge.net/swig/swig-4.0.2.tar.gz
		tar -xvzf swig-$VERSION.tar.gz
fi

cd $SRC_DIR

./configure --prefix=$INSTALL_PREFIX
make -j 16
