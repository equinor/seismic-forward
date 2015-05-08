cmake_minimum_required(VERSION 2.8)

project(G2S)

option( BUILD_MATLAB_MEX  "Build .mex files for use with Matlab." ON)
option( BUILD_UNIT_TESTS  "Build unit tests (Requires boost)." ON)

enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

if(APPLE)
    cmake_policy(SET CMP0042 OLD) #OSX RPATH
endif(APPLE)

set(CMAKE_BUILD_TYPE RELEASE)

include(cmake/AddSource.cmake)

find_package(MKL)
if(MKL_FOUND)
    include_directories(${MKL_INCLUDE_DIRS} ${MKL_FFTW_INCLUDE_DIRS})
    set(MATH_LIBRARIES ${MKL_LIBRARIES})
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

    set(MATH_LIBRARIES ${CBLAS_LIBRARY} lapack ${FFTW_STATIC_LIBRARIES} ${FFTWF_STATIC_LIBRARIES})
endif(MKL_FOUND)

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
include_directories(BEFORE "${INCLUDE_ROOTS}")

add_library(geo2seis_lib ${SOURCE_FILES})
#add_library(geo2seis_static STATIC $<TARGET_OBJECTS:geo2seis_object>)

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

target_link_libraries(seismic_forward ${MATH_LIBRARIES} geo2seis_lib)

if(BUILD_MATLAB_MEX)
  add_subdirectory(mex)
endif(BUILD_MATLAB_MEX)

if(BUILD_UNIT_TESTS)
add_subdirectory(tests)
endif(BUILD_UNIT_TESTS)

