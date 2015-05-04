# - Find Intel MKL
# Find the native FFTW includes and library
#
#  MKLROOT               - Intel MKL root directory
#  MKL_INCLUDE_DIRS      - where to find mkl headers
#  MKL_FFTW_INCLUDE_DIRS - where to find MKL implementation of FFTW headers
#  MKL_LIBRARIES         - MKL link options
#  MKL_FOUND             - True if FFTW found.

if(NOT MKLROOT)
   if(DEFINED ENV{MKLROOT})
      set(MKLROOT $ENV{MKLROOT})
   endif()
endif()

if(MKLROOT)
   set(MKL_INCLUDE_DIRS ${MKLROOT}/include)
   set(MKL_FFTW_INCLUDE_DIRS ${MKLROOT}/include/fftw)

   if(UNIX)
      if(APPLE)
         set(MKL_LIBRARIES "${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_core.a ${MKLROOT}/lib/libmkl_sequential.a -lpthread -lm")
      else()
         set(MKL_LIBRARIES "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group")
      endif(APPLE)
   endif(UNIX)
   if(WIN32)
      set(MKL_LIBRARIES mkl_intel_lp64.lib mkl_core.lib mkl_sequential.lib)
   endif(WIN32)
endif(MKLROOT)


# handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIRS)

mark_as_advanced (MKL_LIBRARIES MKL_INCLUDE_DIRS MKL_FFTW_INCLUDE_DIRS)