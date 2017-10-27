cmake_minimum_required(VERSION 2.6)

set(C_FILES
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
    boost_1_65_1/boost/filesystem/config.hpp
    boost_1_65_1/boost/filesystem/convenience.hpp
    boost_1_65_1/boost/filesystem/fstream.hpp
    boost_1_65_1/boost/filesystem/operations.hpp
    boost_1_65_1/boost/filesystem/path.hpp
    boost_1_65_1/boost/system/config.hpp
    boost_1_65_1/boost/system/cygwin_error.hpp
    boost_1_65_1/boost/system/error_code.hpp
    boost_1_65_1/boost/system/linux_error.hpp
    boost_1_65_1/boost/system/system_error.hpp
    boost_1_65_1/boost/system/windows_error.hpp
    boost_1_65_1/boost/cerrno.hpp
    boost_1_65_1/boost/filesystem.hpp

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
