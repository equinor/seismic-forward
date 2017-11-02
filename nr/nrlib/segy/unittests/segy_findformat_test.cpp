// $Id: segy_findformat_test.cpp 1690 2017-09-13 09:49:32Z perroe $

/// \file Unit tests for finding SEGY format from input SEGY file.

#include "nrlib/segy/segy.hpp"
#include "unittests/util.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

using namespace boost::filesystem;
using namespace NRLib;


BOOST_AUTO_TEST_CASE(SegyFindTraceHeaderFormat) {
  path filename = GetTestDir() / "lnea32_mr_jura_IL_3461_4410_XL_1617_2378.segy";

  TraceHeaderFormat format;
  BOOST_CHECK_NO_THROW(format = SegY::FindTraceHeaderFormat(filename.string()));

  BOOST_CHECK(format.GetStandardType());
  BOOST_CHECK_EQUAL(format.GetFormatName(), "SeisWorks");
}


BOOST_AUTO_TEST_CASE(SegyFindTraceHeaderFormatLine) {
  path filename = GetTestDir() / "seismic" / "single_line.sgy";

  TraceHeaderFormat format;
  BOOST_CHECK_NO_THROW(format = SegY::FindTraceHeaderFormat(filename.string()));

  BOOST_CHECK(format.GetStandardType());
  BOOST_CHECK_EQUAL(format.GetFormatName(), "SeisWorks");
}
