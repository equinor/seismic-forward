cmake_minimum_required(VERSION 2.8)

set(C_FILES
    boost/filesystem/config.hpp
    boost/filesystem/convenience.hpp
    boost/filesystem/fstream.hpp
    boost/filesystem/operations.cpp
    boost/filesystem/path.cpp
    boost/filesystem/portability.cpp
    boost/system/config.hpp
    boost/system/cygwin_error.hpp
    boost/system/error_code.cpp
    boost/system/linux_error.hpp
    boost/system/system_error.hpp
    boost/system/windows_error.hpp
    boost/cerrno.hpp
    boost/filesystem.hpp

    flens/aux_complex.cc
    flens/aux_malloc.cc
    flens/blas.cc
    flens/id.cc
    flens/lapack.cc
    flens/range.cc
    flens/scalarclosures.cc
    flens/snapshot.cc
    flens/sparse_blas.cc
    flens/underscore.cc
)

set(H_FILES
    flens/aux_complex.h
    flens/aux_malloc.h
    flens/blas.h
    flens/id.h
    flens/lapack.h
    flens/range.h
    flens/scalarclosures.h
    flens/snapshot.h
    flens/sparse_blas.h
    flens/underscore.h
)

append_source_files_to_property(${CMAKE_CURRENT_SOURCE_DIR} "${C_FILES}")
add_include_root(${CMAKE_CURRENT_SOURCE_DIR})