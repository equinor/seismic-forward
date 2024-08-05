# - Find Intel Studio libraries (MKL and TBB)
#
# Currently only supports static linking, sequential running and compilation
# for the Intel 64 architecture.
#
#  INTEL_STUDIO_ROOT     - Intel Studio XE root library
#  INTEL_STUDIO_FOUND    - True if Intel Studio is found.
#  MKL_ROOT              - Intel MKL root directory
#  MKL_INCLUDE_DIRS      - where to find MKL headers
#  MKL_FFTW_INCLUDE_DIRS - where to find MKL implementation of FFTW headers
#  MKL_LIBRARIES         - MKL link options
#  TBB_ROOT              - Intel TBB root director
#  TBB_INCLUDE_DIRS      - where to find TBB headers
#  TBB_LIBRARY           - TBB link options

if(UNIX)
    find_path(INTEL_STUDIO_ROOT
              NAMES
	           mkl/latest/include/mkl.h
	           mkl/include/mkl.h
              PATHS
		   ENV ONEAPI_ROOT
                   /prog/Intel/studioxe2016
                   /prog/Intel/studioxe2015
                   /prog/Intel/studioxe2013  # Statoil Linux
                   /nr/prog/intel/Compiler   # NR Linux
                   /opt/software/intel       # Paals laptop fra 2018
              DOC "Root directory for Intel MKL."
              )
elseif(WIN32)
    set(PROG_FILES_ENV "PROGRAMFILES(X86)") # Workaround to handle parenthesis
    find_path(INTEL_STUDIO_ROOT
              NAMES mkl/include/mkl.h
              PATHS
                   "$ENV{${PROG_FILES_ENV}}/IntelSWTools/compilers_and_libraries/windows"
                   "%ProgramFiles%/IntelSWTools/compilers_and_libraries/windows"
                   "$ENV{${PROG_FILES_ENV}}/Intel/Composer XE"
                   "%ProgramFiles%/Intel/Composer XE"
              DOC "Root directory for Intel MKL."
              )
#else()
#    find_path(MKL_ROOT
#              NAMES include/mkl.h
#              PATHS
#                   ENV MKLROOT
#              DOC "Root directory for Intel MKL."
#              )
endif(UNIX)

if(INTEL_STUDIO_ROOT)
   # Intel MKL
   find_path(MKL_ROOT
             NAMES include/mkl.h
             PATHS
		${INTEL_STUDIO_ROOT}/mkl/latest
		${INTEL_STUDIO_ROOT}/mkl
             DOC "Root directory for Intel MKL."
             )

   find_path(MKL_INCLUDE_DIRS
             NAMES mkl.h
             PATHS ${MKL_ROOT}/include
             DOC "Intel MKL include directory."
             )
   find_path(MKL_FFTW_INCLUDE_DIRS
             NAMES fftw.h
             PATHS ${MKL_ROOT}/include/fftw
             DOC "Include directory for FFTW interface to Intel MKL."
             )

   if(UNIX)
      if(APPLE)
         find_library(MKL_INTERFACE_LIBRARY
                      NAME libmkl_intel_lp64.a
                      PATHS ${MKL_ROOT}/lib)
         find_library(MKL_THREADING_LIBRARY
                      NAME libmkl_sequential.a
                      PATHS ${MKL_ROOT}/lib)
         find_library(MKL_COMPUTATIONAL_LIBRARY
                      NAME libmkl_core.a
                      PATHS ${MKL_ROOT}/lib)
         set(MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_COMPUTATIONAL_LIBRARY} pthread m)
      else() # Linux
         set(MKL_LIB_PATH ${MKL_ROOT}/lib/intel64)
         find_library(MKL_INTERFACE_LIBRARY
                      NAME libmkl_intel_lp64.a
                      PATHS ${MKL_LIB_PATH})
         find_library(MKL_THREADING_LIBRARY
                      NAME libmkl_sequential.a
                      PATHS ${MKL_LIB_PATH})
         find_library(MKL_COMPUTATIONAL_LIBRARY
                      NAME libmkl_core.a
                      PATHS ${MKL_LIB_PATH})
         set(MKL_LINK_GROUP "-Wl,--start-group ${MKL_INTERFACE_LIBRARY} ${MKL_COMPUTATIONAL_LIBRARY} ${MKL_THREADING_LIBRARY} -Wl,--end-group")
         set(MKL_LIBRARIES ${MKL_LINK_GROUP} pthread m dl)
      endif(APPLE)
   elseif(WIN32)
      set(MKL_LIB_PATH ${MKL_ROOT}/lib/intel64)
      find_library(MKL_INTERFACE_LIBRARY
                   NAME mkl_intel_lp64.lib
                   PATHS ${MKL_LIB_PATH})
      find_library(MKL_THREADING_LIBRARY
                   NAME mkl_sequential.lib
                   PATHS ${MKL_LIB_PATH})
      find_library(MKL_COMPUTATIONAL_LIBRARY
                   NAME mkl_core.lib
                   PATHS ${MKL_LIB_PATH})
      set(MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_COMPUTATIONAL_LIBRARY})
   endif(UNIX)

   # Intel TBB
   find_path(TBB_ROOT
             NAMES include/tbb/compat/thread include/tbb/tbb.h
             PATHS
		${INTEL_STUDIO_ROOT}/tbb/latest
		${INTEL_STUDIO_ROOT}/tbb
             DOC "Root directory for Intel TBB."
             )

   find_path(TBB_INCLUDE_DIRS
             NAMES tbb/compat/thread tbb/tbb.h
             PATHS ${TBB_ROOT}/include
             DOC "Intel TBB include directory."
             )

   if(UNIX)
     execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                     OUTPUT_VARIABLE GCC_VERSION)

     if(GCC_VERSION VERSION_LESS 4.4)
       set(TBB_GCC_VER gcc4.1)
     else()
       set(TBB_GCC_VER gcc4.4)
     endif()

     find_path(TBB_LIB_DIR
               NAMES libtbb.so
               PATHS ${TBB_ROOT}/lib/intel64/${TBB_GCC_VER}
               DOC "Library directory for Intel TBB"
               )

   elseif(WIN32)
     if(MSVC_VERSION EQUAL 1800)
       set(TBB_MSVC_VER vc12)
     elseif(MSVC_VERSION EQUAL 1700)
       set(TBB_MSVC_VER vc11)
     elseif(MSVC_VERSION EQUAL 1600)
       set(TBB_MSVC_VER vc10)
     else()
       set(TBB_MSVC_VER vc_mt)
     endif()

     find_path(TBB_LIB_DIR
               NAMES tbb.lib
               PATHS ${TBB_ROOT}/lib/intel64/${TBB_MSVC_VER}
               DOC "Library directory for Intel TBB"
               )
   endif(UNIX)

   find_library(TBB_LIBRARIES_RELEASE
                NAMES tbb libtbb
                HINTS ${TBB_LIB_DIR})

   #find_library(TBB_LIBRARIES_DEBUG
   #             NAMES tbb_debug libtbb_debug
   #             HINTS ${TBB_LIB_DIR})

   set(TBB_LIBRARIES
   #   debug ${TBB_LIBRARIES_DEBUG}
      optimized ${TBB_LIBRARIES_RELEASE})


endif(INTEL_STUDIO_ROOT)


# handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(INTEL_STUDIO DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIRS TBB_INCLUDE_DIRS TBB_LIBRARIES)

mark_as_advanced(MKL_LIBRARIES MKL_INCLUDE_DIRS MKL_FFTW_INCLUDE_DIRS)
mark_as_advanced(MKL_INTERFACE_LIBRARY MKL_THREADING_LIBRARY MKL_COMPUTATIONAL_LIBRARY)
mark_as_advanced(TBB_LIBRARIES TBB_INCLUDE_DIRS TBB_LIB_DIR)
mark_as_advanced(TBB_LIBRARIES_RELEASE)
#mark_as_advanced(TBB_LIBRARIES_RELEASE TBB_LIBRARIES_DEBUG)

# message(STATUS "TBB_ROOT = ${TBB_ROOT}")
# message(STATUS "TBB_INCLUDE_DIRS = ${TBB_INCLUDE_DIRS}")
# message(STATUS "TBB_LIBRARIES = ${TBB_LIBRARIES}")
# message(STATUS "TBB_LIB_DIR = ${TBB_LIB_DIR}")
