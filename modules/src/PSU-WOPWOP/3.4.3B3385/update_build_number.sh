#!/bin/bash

# The purpose of this script is twofold - first, it updates the "build number" 
# compiled into the code by creating a file called build_number.f90 that 
# returns a string with the current build (i.e. the current Subversion 
# revision number). It then checks to see if any of the files under subversion 
# control have been locally modified. In that case, it appends the string 
# "locally modified" to the build number.
#
# The build_number.f90 file is created on the fly so that at no point is a
# file with the revision number committed to Subversion - if it were then
# that file would by modified all the time, adding needless updates to the
# repository.

# Function: WriteBuildFile
# Arguments:
#   $1 - The new build number
#   $2 - The file to write to
#   $3 - The number of locally modified files
function WriteBuildFile
{
  echo "! DO NOT EDIT -- Automatically generated file. See update_build_number.sh" > $2
  echo "! " >> $2
  echo "! " >> $2
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!! " >> $2
  echo "! DO NOT ADD TO SUBVERSION ! " >> $2
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!! " >> $2
  echo "! " >> $2
  echo "! " >> $2
  echo "! " >> $2
  echo "function GetBuildNumber () result (build)" >> $2
  echo "  character(len=128) :: build" >> $2
  if [ "$3" -gt "0" ]; then
    echo "  build=\"$1 [locally modified]\"" >> $2
  else
    echo "  build=\"$1\"" >> $2
  fi
  echo "end function GetBuildNumber" >> $2
}

findsvn=`which svn 2>&1 | grep -i "which: no"`
if [ -n "$findsvn" ]; then
  echo "Setting build number to (no build number: no svn client)."
  WriteBuildFile '-1' 'build_number.f90' 0
else

  # Get the previous modification status:

  # Determine if anything has been modified locally: recall that the regular
  # expression "^M" searches for lines that begin with the letter "M", i.e.
  # files that subversion says have been locally modified.
  modified=`svn status | grep "^M" | wc -l` # Returns the number of modified files

  # Get the SVN rev:
  new_rev=`svn info ./ | grep 'Revision:' | awk '{ print $2 }'`

  # Get the current rev:
  old_rev=0
  old_modified=0
  if [ -e 'build_number.f90' ]; then
    old_rev=`cat build_number.f90 | grep 'build="' | sed 's/\(^  build="\)\([0-9]\+\).*/\2/'`
    check=`cat build_number.f90 | grep modified`
    if [ -n "$check" ]; then
      # The build number already indicates that the file was modified - no need to
      # change it.
      old_modified=$modified
    fi
  fi

  # Check the revs:
  if [ "$new_rev" -eq "$old_rev" -a "$old_modified" -eq "$modified" ]; then
    echo "Rebuilding $new_rev"
  else 
    echo "Setting build number to $new_rev."
    WriteBuildFile $new_rev 'build_number.f90' $modified
  fi

 # found=`which md5sum 2>&1 | grep "which: no"`
 # if [ -z "$found" ]; then
    #md5sum *.f90 *.f > .md5sums
 # fi

fi
