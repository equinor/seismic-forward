# - Find Intel MKL
# Find the native FFTW includes and library
#
# Currently only supports static linking, sequential running and compilation
# for the Intel 64 architecture.
#
#  MKLROOT               - Intel MKL root directory
#  MKL_INCLUDE_DIRS      - where to find mkl headers
#  MKL_FFTW_INCLUDE_DIRS - where to find MKL implementation of FFTW headers
#  MKL_LIBRARIES         - MKL link options
#  MKL_FOUND             - True if FFTW found.

if(UNIX)
    find_path(MKLROOT
              NAMES include/mkl.h
              PATHS 
                   ENV MKLROOT
                   /nr/prog/intel/Compiler/mkl
              DOC "Root directory for Intel MKL."
              )
elseif(WIN32)
    set(PROG_FILES_ENV "PROGRAMFILES(X86)") # Needed to handle parenthesis
    find_path(MKLROOT
              NAMES include/mkl.h
              PATHS 
                   ENV MKLROOT
                   "$ENV{${PROG_FILES_ENV}}/Intel/Composer XE/mkl"
                   "%ProgramFiles%/Intel/Composer XE/mkl"
              DOC "Root directory for Intel MKL."
              )
else()
    find_path(MKLROOT
              NAMES include/mkl.h
              PATHS 
                   ENV MKLROOT
              DOC "Root directory for Intel MKL."
              )
endif(UNIX)

if(MKLROOT)
   find_path(MKL_INCLUDE_DIRS
             NAMES mkl.h
             PATHS ${MKLROOT}/include
             DOC "Intel MKL include directory."
             )
   find_path(MKL_FFTW_INCLUDE_DIRS
             NAMES fftw.h
             PATHS ${MKLROOT}/include/fftw
             DOC "Include directory for FFTW interface to Intel MKL."
             )

   if(UNIX)
      if(APPLE)
         find_library(MKL_INTERFACE_LIBRARY
                      NAME libmkl_intel_lp64.a
                      PATHS ${MKLROOT}/lib)
         find_library(MKL_THREADING_LIBRARY
                      NAME libmkl_sequential.a
                      PATHS ${MKLROOT}/lib)
         find_library(MKL_COMPUTATIONAL_LIBRARY
                      NAME libmkl_core.a
                      PATHS ${MKLROOT}/lib)
         set(MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_COMPUTATIONAL_LIBRARY} pthread m)
      else() # Linux
         set(MKL_LIB_PATH ${MKLROOT}/lib/intel64) 
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
      set(MKL_LIB_PATH ${MKLROOT}/lib/intel64) 
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
endif(MKLROOT)


# handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIRS)

mark_as_advanced(MKL_LIBRARIES MKL_INCLUDE_DIRS MKL_FFTW_INCLUDE_DIRS) 
mark_as_advanced(MKL_INTERFACE_LIBRARY MKL_THREADING_LIBRARY MKL_COMPUTATIONAL_LIBRARY)