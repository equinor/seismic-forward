// $Id: norsarwell_test.cpp 1642 2017-07-14 10:54:31Z perroe $

/// \file Unit tests testing reading Norsar well format.

#include <nrlib/well/norsarwell.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

#include <unittests/util.hpp>

using namespace boost::filesystem;
using namespace NRLib;

BOOST_AUTO_TEST_CASE(NorsarWellReadTest)
{
  path filename = GetTestDir() / "wells" / "norsar_well_w_twt.nwh";

  NorsarWell well(filename.string());

  BOOST_CHECK_EQUAL(well.GetNContLog(), 13U);
  // BOOST_CHECK_EQUAL(well.GetWellName(), "NORSAR_WELL");
  BOOST_CHECK_CLOSE(well.GetXPos0(), 320.59173, 0.001);
  BOOST_CHECK_CLOSE(well.GetYPos0(), 7345.04282, 0.001);
}
