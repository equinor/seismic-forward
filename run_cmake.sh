#!/bin/bash

# Use this script as follows:
# Create the directory for your wanted architecture, for instance
#
# $ mkdir hostname-release    (hostname = name of host)
# $ cd hostname-release
# $ ../seismic-forward/run_cmake.sh
#
# Use -d or --debug for debug compilation
# Use -m or --matlab for MATLAB compilation

matlab=-DBUILD_MATLAB_MEX=OFF
if [ "$1" = "-m" ] || [ "$1" = "--matlab" ] || [ "$2" = "-m" ] || [ "$2" = "--matlab" ]
then
    echo "Build Matlab executables"
    matlab=-DBUILD_MATLAB_MEX=ON
fi

if [ "$1" = "-d" ] || [ "$1" = "--debug" ] || [ "$2" = "-d" ] || [ "$2" = "--debug" ]
then
    echo "Build debug executables"
    debug=-DCMAKE_BUILD_TYPE=DEBUG
fi

echo "Running cmake ../seismic-forward-4.1 ${matlab} ${debug}"
cmake -Wno-dev ../seismic-forward-4.1 ${matlab} ${debug}
