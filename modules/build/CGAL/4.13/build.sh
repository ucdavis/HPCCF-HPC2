#!/bin/bash

source ./module.sh

if [ ! -f "CGAL-${VERSION}.tar.gz" ]; then
		wget "https://github.com/CGAL/cgal/archive/releases/CGAL-${VERSION}.tar.gz"
		tar -xvzf CGAL-$VERSION.tar.gz
fi

cd $SRC_DIR
mkdir -p build
pushd build

cmake	-DCMAKE_INSTALL_PREFIX:STRING=$INSTALL_PREFIX \
			-DCMAKE_BUILD_TYPE:STRING=Release \
			-DCMAKE_CXX_FLAGS:STRING="-Wno-dev -Wno-unknown-warning-option" \
			-DBUILD_SHARED_LIBS:BOOL=ON \
			-DWITH_CGAL_Qt5:BOOL=OFF \
			-DWITH_GMP:BOOL=OFF \
			-DWITH_MPFR:BOOL=OFF \
			..
make -j 16
