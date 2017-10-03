#!/bin/sh

ROOT_DIR=../Geo2Seis
CUR_DIR=`basename "$PWD"`
EXE_DIR=compare_segy_traces
EXE=compare_traces

if [ ${CUR_DIR} = "Geo2Seis" ]; then
    echo "The test program cannot be compiled in the source directory. Please,"
    echo "copy the script to the top directory where object files are stores"
    echo "and execute from there."
    exit;
fi

if [ -e ${EXE_DIR} ]; then
    echo "Directory compare_segy_traces exists"
else
    echo "Making directory compare_segy_traces"
    mkdir ${EXE_DIR}
fi

echo "Removing old executable and object file"
rm ${EXE_DIR}/compare_traces.o
rm ${EXE}

echo "Compiling compare_segy_traces/compare.cpp"
g++ -c -O2 \
    -I${ROOT_DIR}  \
    ${ROOT_DIR}/${EXE_DIR}/compare_traces.cpp -o compare_segy_traces/compare_traces.o

echo "Linking compare.exe"
g++ -I${ROOT_DIR} \
    -I${ROOT_DIR}/nr \
    -I${ROOT_DIR}/geo2seis \
    -I/opt/software/intel/compilers_and_libraries_2016.2.181/linux/tbb/include \
    ${EXE_DIR}/compare_traces.o                                             \
    CMakeFiles/geo2seis_lib.dir/nr/nrlib/iotools/fileio.cpp.o               \
    CMakeFiles/geo2seis_lib.dir/nr/nrlib/segy/segy.cpp.o                    \
    CMakeFiles/geo2seis_lib.dir/nr/nrlib/segy/segytrace.cpp.o               \
    CMakeFiles/geo2seis_lib.dir/nr/nrlib/segy/segygeometry.cpp.o            \
    CMakeFiles/geo2seis_lib.dir/nr/nrlib/segy/traceheader.cpp.o             \
    CMakeFiles/geo2seis_lib.dir/nr/nrlib/segy/commonheaders.cpp.o           \
    CMakeFiles/geo2seis_lib.dir/nr/nrlib/math/constants.cpp.o               \
    CMakeFiles/geo2seis_lib.dir/nr/nrlib/iotools/stringtools.cpp.o          \
    CMakeFiles/geo2seis_lib.dir/nr/nrlib/iotools/logkit.cpp.o               \
    CMakeFiles/geo2seis_lib.dir/nr/nrlib/volume/volume.cpp.o                \
    CMakeFiles/geo2seis_lib.dir/nr/nrlib/surface/surfaceio.cpp.o            \
    CMakeFiles/geo2seis_lib.dir/3rd-party/boost/system/error_code.cpp.o     \
    CMakeFiles/geo2seis_lib.dir/3rd-party/boost/filesystem/operations.cpp.o \
    CMakeFiles/geo2seis_lib.dir/3rd-party/boost/filesystem/path.cpp.o       \
    -o ${EXE}
