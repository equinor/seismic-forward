SRC += flens/aux_complex.cc \
       flens/aux_malloc.cc \
       flens/blas.cc \
       flens/id.cc \
       flens/lapack.cc \
       flens/range.cc \
       flens/scalarclosures.cc \
       flens/snapshot.cc \
       flens/sparse_blas.cc \
       flens/underscore.cc

#--------------------------------------------
# Using ATLAS:
# ATLASLFLAGS += -L/usr/lib64/atlas -L/usr/lib64/atlas/atlas -L/lib64 -llapack -lblas -lcblas -latlas -lgfortran -lpthread
# CXXFLENS    += -DFLENS_FIRST_INDEX=0
# Compile with: -DWITH_ATLAS


#--------------------------------------------
# Using Intel MKL
#MKLINCLUDE := /nr/prog/intel/Compiler/mkl/include
#MKLPATH :=    /nr/prog/intel/Compiler/mkl/lib/intel64

#CXXFLAGS += -DFLENS_FIRST_INDEX=0 -DMKL -I $(MKLINCLUDE)

## Dynamic linking
##LFLAGS += -L /nr/prog/mkl21/lib/em64t -lmkl_sequential -lmkl_core -lmkl_intel_lp64

# Static linking
# MKLLFLAGS += -L$(MKLPATH) -I$(MKLINCLUDE) -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a \
#	$(MKLPATH)/libmkl_sequential.a $(MKLPATH)/libmkl_core.a -Wl,--end-group \
#	-lpthread -ldl

##-static -Wl,--start-group $(MKLLIBPATH)/libmkl_intel_lp64.a $(MKLLIBPATH)/libmkl_sequential.a $(MKLLIBPATH)/libmkl_core.a -Wl,--end-group
#--------------------------------------------

