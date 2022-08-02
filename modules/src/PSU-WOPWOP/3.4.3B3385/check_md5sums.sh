#!/bin/bash

md5sum *.f90 *.f > .new_md5sums
difference=`diff .new_md5sums .md5sums`
if [ -n "$difference" ]; then
  echo "Local modifications have been made."
fi
