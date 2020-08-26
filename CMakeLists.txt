cmake_minimum_required(VERSION 2.6)

project(G2S)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

option( BUILD_MATLAB_MEX  "Build .mex files for use with Matlab." ON)

#
# Deactivated by PAAL 2017.10.25 because of problems to get this to work with new boost
#
#option( BUILD_UNIT_TESTS  "Build unit tests (Requires boost)." ON)
#enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

if(APPLE)
    cmake_policy(SET CMP0042 OLD) #OSX RPATH
endif(APPLE)

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RELEASE)
endif()

include(cmake/AddSource.cmake)

#
# Boost things added by PAAL
#
set(BOOST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/3rd-party/boost_1_65_1")
#find_package(Boost 1.65.1 REQUIRED COMPONENTS system filesystem)

find_package(Boost 1.65.1 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
#link_directories(${Boost_LIBRARY_DIRS})

option( USE_INTEL_STUDIO "Use libraries from Intel Studio" ON )
find_package(IntelStudio)
if(NOT INTEL_STUDIO_FOUND)
 option(USE_INTEL_STUDIO OFF)
endif(NOT INTEL_STUDIO_FOUND)

option( USE_OPENMP "Use OpenMP for parallelisation" ON )
if (USE_OPENMP)
	find_package(OpenMP)
	if (OPENMP_FOUND)
		add_definitions(-DWITH_OMP)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif()
endif()


if(USE_INTEL_STUDIO)
    include_directories(${MKL_INCLUDE_DIRS} ${MKL_FFTW_INCLUDE_DIRS} ${TBB_INCLUDE_DIRS})
    link_directories(${TBB_LIB_DIR})
    set(FFT_LIBRARIES ${MKL_LIBRARIES})
    set(MATH_LIBRARIES ${MKL_LIBRARIES} ${TBB_LIBRARIES})
    set(FFT_MATLAB_LIBRARIES -L${MKL_LIB_PATH} -lmkl_rt)
    add_definitions(-DMKL)
else()
    # Out-of-the-box on Unix Statoil systems this will not work and a path is required.
    # For example /prog/sdpsoft/gsl/lib/libgslcblas.a
    set(CBLAS_LIBRARY "blas" CACHE FILEPATH "Path to a static or dynamic version of CBLAS")

    find_package(FFTW REQUIRED)
    if(FFTW_FOUND)
        include_directories(${FFTW_INCLUDES})
    else(FFTW_FOUND)
        message(FATAL_ERROR "Unable to find FFTW library.")
    endif(FFTW_FOUND)

    set(FFT_LIBRARIES ${FFTW_STATIC_LIBRARIES} ${FFTWF_STATIC_LIBRARIES})
    set(FFT_MATLAB_LIBRARIES ${FFT_LIBRARIES})

    find_package(TBB REQUIRED)
    if(TBB_FOUND)
        include_directories(${TBB_INCLUDE_DIRS})
    else(TBB_FOUND)
        message(FATAL_ERROR "Unable to find Intel TBB library.")
    endif(TBB_FOUND)

    set(MATH_LIBRARIES ${CBLAS_LIBRARY} lapack ${FFT_LIBRARIES} ${TBB_LIBRARIES})
endif(USE_INTEL_STUDIO)

set_property(GLOBAL PROPERTY SOURCE_FILES)
set_property(GLOBAL PROPERTY INCLUDE_ROOTS)

add_definitions(-DFLENS_FIRST_INDEX=0)
if(MSVC)
    add_definitions(-wd4996)
else()
    add_definitions(-fPIC)
endif(MSVC)

add_subdirectory(3rd-party)
add_subdirectory(nr)
add_subdirectory(geo2seis)

get_property(SOURCE_FILES GLOBAL PROPERTY SOURCE_FILES)
get_property(INCLUDE_ROOTS GLOBAL PROPERTY INCLUDE_ROOTS)
foreach(INCLUDE_ROOT ${INCLUDE_ROOTS})
  include_directories(BEFORE "${INCLUDE_ROOT}")
endforeach()

add_library(geo2seis_lib STATIC ${SOURCE_FILES})

foreach(FILE ${SOURCE_FILES})
  get_filename_component(PARENT_DIR "${FILE}" PATH)

  # skip src or include and changes /'s to \\'s
  string(REGEX REPLACE "${CMAKE_CURRENT_SOURCE_DIR}" "" GROUP "${PARENT_DIR}")
  string(REPLACE "/" "\\" GROUP "${GROUP}")

  # group into "Source Files" and "Header Files"
  if ("${FILE}" MATCHES ".*\\.cpp")
    set(GROUP "Source Files\\${GROUP}")
  elseif ("${FILE}" MATCHES ".*\\.cc")
    set(GROUP "Source Files\\${GROUP}")
  elseif("${FILE}" MATCHES ".*\\.h")
    set(GROUP "Header Files\\${GROUP}")
  endif()

  source_group("${GROUP}" FILES "${FILE}")
endforeach()

add_executable(seismic_forward geo2seis/main.cpp)

#message(STATUS "Boost libraries to be linked with: ${Boost_LIBRARIES}")

#target_link_libraries(seismic_forward geo2seis_lib ${MATH_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(seismic_forward geo2seis_lib ${MATH_LIBRARIES})

if(BUILD_MATLAB_MEX)
  add_subdirectory(mex)
endif(BUILD_MATLAB_MEX)

if(BUILD_UNIT_TESTS)
add_subdirectory(tests)
endif(BUILD_UNIT_TESTS)
