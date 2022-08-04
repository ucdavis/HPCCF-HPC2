#!/bin/bash 
echo $1
# A script to check out a clean version of PSU-WOPWOP set the build number, 
# clean it up, and create a tarball. It takes one non-optional argument: the
# build number to create. If no build number is given, the latest version
# is used.

# Check to see if a build number was given
if [ "$#" -eq  "2" ]; then
  # Try to check out that build number
  echo "Checking out build $1 from $2"
  svn -q co -r $1 $2 PSU-WOPWOP_Build_$1
  # Check to see if the return status was OK:
  if [ "$?" -ne "0" ]; then
    echo "Checkout failed: probably a bad build number."
    exit
  fi
  #Go into the directory and create the build number
  cd  PSU-WOPWOP_Build_$1
  if [ -e "update_build_number.sh" ]; then
    ./update_build_number.sh
  else
    echo "Couldn't update the build number: this build probably pre-dates that code."
  fi
  #Remove the version control files:
  /bin/rm -rf `find -name .svn`
  /bin/rm -f update_build_number.sh
  echo -e "#!/bin/bash\n\n./check_md5sums.sh\n">update_build_number.sh
  chmod +x update_build_number.sh
  /bin/rm -f make_build.sh
  # Tar up the build
  cd ..
  echo "Creating tarball of the build..."
  tar czf PSU-WOPWOP_Build_$1.tar.gz PSU-WOPWOP_Build_$1/
  rm -rf PSU-WOPWOP_Build_$1/
  echo "Created build tarball PSU-WOPWOP_Build_$1.tar.gz"
else
  echo "A build number bust be specified."
fi

