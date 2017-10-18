#!/bin/sh

ROOT_DIR=../Geo2Seis
CUR_DIR=`basename "$PWD"`
OBJ_DIR=CMakeFiles/geo2seis_lib.dir
EXE=compare_traces

Boost_INCLUDE_DIR=/usr/include
Boost_LIBRARY_DIR=/usr/lib/x86_64-linux-gnu

if [ ${CUR_DIR} = "Geo2Seis" ]; then
    echo "The test program cannot be compiled in the source directory. Please,"
    echo "copy the script to the top directory where object files are stores"
    echo "and execute from there."
    exit;
fi

if  [ -e ${OBJ_DIR}/geo2seis/compare_traces.o ] ; then
  echo "Removing old object file compare_traces.o"
  rm ${OBJ_DIR}/geo2seis/compare_traces.o
fi
if [ -e ${EXE} ] ; then
  echo "Removing old executable"
  rm ${EXE}
fi

echo "Compiling compare_traces.cpp"
g++ -c -O2 \
    -I${ROOT_DIR} \
    -I${ROOT_DIR}/nr \
    ${ROOT_DIR}/geo2seis/compare_traces.cpp -o ${OBJ_DIR}/geo2seis/compare_traces.o

if  [ -e ${OBJ_DIR}/geo2seis/compare_traces.o ] ; then
  echo "Linking compare_traces"
  g++ ${OBJ_DIR}/geo2seis/compare_traces.o              \
      ${OBJ_DIR}/nr/nrlib/iotools/fileio.cpp.o          \
      ${OBJ_DIR}/nr/nrlib/iotools/logstream.cpp.o       \
      ${OBJ_DIR}/nr/nrlib/iotools/bigfile.cpp.o         \
      ${OBJ_DIR}/nr/nrlib/iotools/stringtools.cpp.o     \
      ${OBJ_DIR}/nr/nrlib/iotools/logkit.cpp.o          \
      ${OBJ_DIR}/nr/nrlib/segy/segy.cpp.o               \
      ${OBJ_DIR}/nr/nrlib/segy/segytrace.cpp.o          \
      ${OBJ_DIR}/nr/nrlib/segy/segygeometry.cpp.o       \
      ${OBJ_DIR}/nr/nrlib/segy/traceheader.cpp.o        \
      ${OBJ_DIR}/nr/nrlib/segy/commonheaders.cpp.o      \
      ${OBJ_DIR}/nr/nrlib/stormgrid/stormcontgrid.cpp.o \
      ${OBJ_DIR}/nr/nrlib/volume/volume.cpp.o           \
      ${OBJ_DIR}/nr/nrlib/surface/surfaceio.cpp.o       \
      ${OBJ_DIR}/nr/nrlib/math/constants.cpp.o          \
      -lboost_filesystem                                \
      -lboost_system                                    \
      -o ${EXE}
fi
