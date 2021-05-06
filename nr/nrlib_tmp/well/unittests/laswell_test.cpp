// $Id: laswell_test.cpp 1642 2017-07-14 10:54:31Z perroe $

/// \file Unit tests testing reading LAS well.

#include <nrlib/well/laswell.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

#include <unittests/util.hpp>

using namespace boost::filesystem;
using namespace NRLib;

BOOST_AUTO_TEST_CASE(LASWellReadTest)
{
  path filename = GetTestDir() / "wells" / "example_3_0.las";

  LasWell well(filename.string());

  BOOST_CHECK_EQUAL(well.GetNContLog(), 19U);
  BOOST_CHECK_EQUAL(well.GetWellName(), "\"Example LAS well\"");
  BOOST_CHECK_CLOSE(well.GetXPos0(), 672591.73, 0.001);
  BOOST_CHECK_CLOSE(well.GetYPos0(), 9825042.82, 0.001);
}
