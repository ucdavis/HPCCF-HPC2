#!/bin/bash

if [ -r .compiler ]; then
  target=`cat .compiler`
  make $target
else
  echo "No .compiler file: using intel"
  make clean
  make intel
fi
