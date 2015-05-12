cmake_minimum_required(VERSION 2.6)

set(C_FILES
    boost/filesystem/operations.cpp
    boost/filesystem/path.cpp
    boost/filesystem/portability.cpp
    boost/system/error_code.cpp

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
    boost/filesystem/config.hpp
    boost/filesystem/convenience.hpp
    boost/filesystem/fstream.hpp
    boost/filesystem/operations.hpp
    boost/filesystem/path.hpp
    boost/system/config.hpp
    boost/system/cygwin_error.hpp
    boost/system/error_code.hpp
    boost/system/linux_error.hpp
    boost/system/system_error.hpp
    boost/system/windows_error.hpp
    boost/cerrno.hpp
    boost/filesystem.hpp

    flens/array.h
    flens/aux_cmath.h
    flens/aux_complex.h
    flens/aux_malloc.h
    flens/bandstorage.h
    flens/blas.h
    flens/blas_flens.h
    flens/blas_generic.h
    flens/cg.h
    flens/crs.h
    flens/densevector.h
    flens/evalclosure.h
    flens/fixedsizearray.h
    flens/flens.h
    flens/fullstorage.h
    flens/generalmatrix.h
    flens/hacksforgmpxx.h
    flens/hermitianmatrix.h
    flens/id.h
    flens/lapack.h
    flens/lapack_flens.h
    flens/lapack_flens_aux.h
    flens/lapack_generic.h
    flens/listinitializer.h
    flens/macros.h
    flens/matvec.h
    flens/matvecclosures.h
    flens/matvecio.h
    flens/matvecoperations.h
    flens/multigrid.h
    flens/operationtypes.h
    flens/packedstorage.h
    flens/polynom.h
    flens/range.h
    flens/result.h
    flens/scalarclosures.h
    flens/scalaroperations.h
    flens/snapshot.h
    flens/sparsematrix.h
    flens/sparse_blas.h
    flens/sparse_blas_flens.h
    flens/storage.h
    flens/symmetricmatrix.h
    flens/tinymatrix.h
    flens/tinyvector.h
    flens/tiny_blas.h
    flens/traits.h
    flens/triangularmatrix.h
    flens/underscore.h
    flens/uplo.h
)


append_source_files_to_property(${CMAKE_CURRENT_SOURCE_DIR} "${C_FILES}")
append_source_files_to_property(${CMAKE_CURRENT_SOURCE_DIR} "${H_FILES}")
add_include_root(${CMAKE_CURRENT_SOURCE_DIR})

