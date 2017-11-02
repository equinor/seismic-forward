SRC += $(NRLIB_BASE_DIR)iotools/bigfile.cpp \
       $(NRLIB_BASE_DIR)iotools/fileio.cpp \
       $(NRLIB_BASE_DIR)iotools/logkit.cpp \
       $(NRLIB_BASE_DIR)iotools/logstream.cpp \
       $(NRLIB_BASE_DIR)iotools/stringtools.cpp \
       $(NRLIB_BASE_DIR)iotools/tabularfile.cpp

TEST_SRC += $(NRLIB_BASE_DIR)iotools/unittests/bigfile_test.cpp \
            $(NRLIB_BASE_DIR)iotools/unittests/fileio_binaryreadwrite_test.cpp \
            $(NRLIB_BASE_DIR)iotools/unittests/fileio_parsefrombuffer_test.cpp \
            $(NRLIB_BASE_DIR)iotools/unittests/fileio_writetobuffer_test.cpp \
            $(NRLIB_BASE_DIR)iotools/unittests/stringtools_test.cpp
