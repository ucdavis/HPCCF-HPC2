#!/bin/bash

if [ $# -ne 1 ]; then
  echo "ERROR: compile_fftw.sh requires an argument."
  exit
fi

if [ ! -r "./fftw_installed/lib/libfftw3_$1.a" ]; then
  pushd $PWD  > /dev/null 2>&1
  cd fftw/
  echo "Compiling FFTW for $1 (this will take a few minutes)..."
  ./configure --prefix=$PWD/../fftw_installed/ --enable-fortran --enable-double > /dev/null 2>&1
  make > make.log
  make install >> make.log
  mv "../fftw_installed/lib/libfftw3.a" "../fftw_installed/lib/libfftw3_$1.a"
  mv "../fftw_installed/lib/libfftw3.la" "../fftw_installed/lib/libfftw3_$1.la"
  popd > /dev/null 2>&1
fi
