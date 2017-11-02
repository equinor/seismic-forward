SRC += $(NRLIB_BASE_DIR)fft/fft.cpp \
       $(NRLIB_BASE_DIR)fft/fftgrid2d.cpp \
	   $(NRLIB_BASE_DIR)fft/fftgrid3d.cpp

# CXXFLAGS += -I /nr/prog/intel/Compiler/mkl/include/fftw

# If not FLENS is used, MKL must be linked in:

# Static linking
#MKLLFLAGS += -L$(MKLPATH) -I$(MKLINCLUDE) -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a \
	$(MKLPATH)/libmkl_sequential.a $(MKLPATH)/libmkl_core.a -Wl,--end-group \
	-lpthread

# Dynamic linking
#LFLAGS += -L /nr/prog/mkl21/lib/em64t -lmkl_sequential -lmkl_core -lmkl_intel_lp64
