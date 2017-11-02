// $Id: rmswell_test.cpp 1642 2017-07-14 10:54:31Z perroe $

/// \file Unit tests testing reading LAS well.

#include <nrlib/well/rmswell.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

#include <unittests/util.hpp>

using namespace boost::filesystem;
using namespace NRLib;

BOOST_AUTO_TEST_CASE(RMSWellReadTest)
{
  path filename = GetTestDir() / "wells" / "Well_A.rmswell";

  RMSWell well(filename.string());

  BOOST_CHECK_EQUAL(well.GetNContLog(), 8U);
  BOOST_CHECK_EQUAL(well.GetNlog(), 14U);
  BOOST_CHECK_EQUAL(well.GetWellName(), "Well_A");
  BOOST_CHECK_CLOSE(well.GetXPos0(), 463174.625, 0.001);
  BOOST_CHECK_CLOSE(well.GetYPos0(), 5933349.000, 0.001);
}
