cmake_minimum_required(VERSION 2.8)

set(SOURCE_FILES
    nrlib/eclipsegrid/eclipsecarfin.cpp
    nrlib/eclipsegrid/eclipsecarfingroup.cpp
    nrlib/eclipsegrid/eclipseeditnnc.cpp
    nrlib/eclipsegrid/eclipsefaceindexer.cpp
    nrlib/eclipsegrid/eclipsefault.cpp
    nrlib/eclipsegrid/eclipsegeometry.cpp
    nrlib/eclipsegrid/eclipsegrid.cpp
    nrlib/eclipsegrid/eclipsetools.cpp
    nrlib/eclipsegrid/eclipsetransmult.cpp

    nrlib/fft/fft.cpp
    nrlib/fft/fftgrid2d.cpp
    nrlib/fft/fftgrid3d.cpp

    nrlib/flens/nrlib_flens.cpp

    nrlib/geometry/bilinearsurface.cpp
    nrlib/geometry/box.cpp
    nrlib/geometry/line.cpp
    nrlib/geometry/plane.cpp
    nrlib/geometry/point.cpp
    nrlib/geometry/polygon.cpp
    nrlib/geometry/transformation.cpp
    nrlib/geometry/triangle.cpp

    nrlib/iotools/fileio.cpp
    nrlib/iotools/logkit.cpp
    nrlib/iotools/stringtools.cpp
    nrlib/iotools/tabularfile.cpp

    nrlib/math/constants.cpp

    nrlib/pointset/pointset.cpp

    nrlib/random/beta.cpp
    nrlib/random/betawithendmass.cpp
    nrlib/random/binomial.cpp
    nrlib/random/chisquared.cpp
    nrlib/random/delta.cpp
    nrlib/random/dSFMT.cpp
    nrlib/random/fractal.cpp
    nrlib/random/functions.cpp
    nrlib/random/gamma.cpp
    nrlib/random/generaldiscrete.cpp
    nrlib/random/lognormal.cpp
    nrlib/random/normal.cpp
    nrlib/random/poisson.cpp
    nrlib/random/random.cpp
    nrlib/random/randomgenerator.cpp
    nrlib/random/splituniform.cpp
    nrlib/random/triangular.cpp
    nrlib/random/uniform.cpp

    nrlib/segy/commonheaders.cpp
    nrlib/segy/segy.cpp
    nrlib/segy/segygeometry.cpp
    nrlib/segy/segytrace.cpp
    nrlib/segy/traceheader.cpp

    nrlib/stormgrid/stormcontgrid.cpp
    nrlib/stormgrid/stormfaciesgrid.cpp

    nrlib/surface/surfaceio.cpp

    nrlib/tinyxml/tinyxml.cpp
    nrlib/tinyxml/tinyxmlerror.cpp
    nrlib/tinyxml/tinyxmlparser.cpp

    nrlib/volume/volume.cpp
)

set(H_FILES
    nrlib/eclipsegrid/eclipsecarfin.hpp
    nrlib/eclipsegrid/eclipsecarfingroup.hpp
    nrlib/eclipsegrid/eclipseeditnnc.hpp
    nrlib/eclipsegrid/eclipsefaceindexer.hpp
    nrlib/eclipsegrid/eclipsefault.hpp
    nrlib/eclipsegrid/eclipsegeometry.hpp
    nrlib/eclipsegrid/eclipsegrid.hpp
    nrlib/eclipsegrid/eclipsetools.hpp
    nrlib/eclipsegrid/eclipsetransmult.hpp

    nrlib/exception/exception.hpp

    nrlib/fft/fft.hpp
    nrlib/fft/fftgrid2d.hpp
    nrlib/fft/fftgrid3d.hpp

    nrlib/flens/nrlib_flens.hpp

    nrlib/geometry/bilinearsurface.hpp
    nrlib/geometry/box.hpp
    nrlib/geometry/line.hpp
    nrlib/geometry/plane.hpp
    nrlib/geometry/point.hpp
    nrlib/geometry/polygon.hpp
    nrlib/geometry/transformation.hpp
    nrlib/geometry/triangle.hpp

    nrlib/grid/grid.hpp
    nrlib/grid/grid2d.hpp
    nrlib/grid/grid4d.hpp
    nrlib/grid/grid5d.hpp

    nrlib/iotools/fileio.hpp
    nrlib/iotools/logkit.hpp
    nrlib/iotools/stringtools.hpp
    nrlib/iotools/tabularfile.hpp

    nrlib/math/constants.hpp
    nrlib/math/mathutility.hpp

    nrlib/pointset/pointset.hpp

    nrlib/random/beta.hpp
    nrlib/random/betawithendmass.hpp
    nrlib/random/binomial.hpp
    nrlib/random/chisquared.hpp
    nrlib/random/delta.hpp
    nrlib/random/distribution.hpp
    nrlib/random/dSFMT.h
    nrlib/random/dSFMT-common.h
    nrlib/random/dSFMT-params.h
    nrlib/random/dSFMT-params19937.h
    nrlib/random/fractal.hpp
    nrlib/random/functions.hpp
    nrlib/random/gamma.hpp
    nrlib/random/generaldiscrete.hpp
    nrlib/random/lognormal.hpp
    nrlib/random/normal.hpp
    nrlib/random/poisson.hpp
    nrlib/random/random.hpp
    nrlib/random/randomgenerator.hpp
    nrlib/random/splituniform.hpp
    nrlib/random/triangular.hpp
    nrlib/random/uniform.hpp

    nrlib/segy/commonheaders.hpp
    nrlib/segy/segy.hpp
    nrlib/segy/segygeometry.hpp
    nrlib/segy/segytrace.hpp
    nrlib/segy/traceheader.hpp

    nrlib/stormgrid/stormcontgrid.hpp
    nrlib/stormgrid/stormfaciesgrid.hpp

    nrlib/surface/regularsurface.hpp
    nrlib/surface/regularsurfacerotated.hpp
    nrlib/surface/surface.hpp
    nrlib/surface/surfaceio.hpp

    nrlib/tinyxml/tinyxml.h

    nrlib/volume/volume.hpp
)


set(TEST_SOURCES
    nrlib/iotools/unittests/fileio_test.cpp
)

append_source_files_to_property(${CMAKE_CURRENT_SOURCE_DIR} "${SOURCE_FILES}")
append_source_files_to_property(${CMAKE_CURRENT_SOURCE_DIR} "${H_FILES}")
add_include_root(${CMAKE_CURRENT_SOURCE_DIR})
