SRC += $(NRLIB_BASE_DIR)segy/commonheaders.cpp \
       $(NRLIB_BASE_DIR)segy/segy.cpp \
       $(NRLIB_BASE_DIR)segy/traceheader.cpp \
       $(NRLIB_BASE_DIR)segy/segytrace.cpp \
       $(NRLIB_BASE_DIR)segy/segygeometry.cpp

TEST_SRC += $(NRLIB_BASE_DIR)segy/unittests/segy_findformat_test.cpp \
            $(NRLIB_BASE_DIR)segy/unittests/segygeometry_test.cpp \
            $(NRLIB_BASE_DIR)segy/unittests/segygeometry_line_test.cpp \
            $(NRLIB_BASE_DIR)segy/unittests/segyio_test.cpp
