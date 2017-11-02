
/// \file Unit tests for Box class in the NRLib Geometry library.

#include <boost/test/unit_test.hpp>

#include "nrlib/geometry/box.hpp"
#include "nrlib/geometry/line.hpp"

BOOST_AUTO_TEST_CASE(CornerIntersectionsTest)
{
  NRLib::Box unit_box(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0);
  NRLib::Line line = NRLib::Line::InfiniteLine(NRLib::Point(0.0, 0.0, 0.0), NRLib::Point(1.0, 1.0, 1.0));

  // This does not work properly:
  // std::vector<NRLib::Point> intersections = unit_box.FindIntersections(line);
  // The results should be corners (-1, -1, -1) and (1, 1, 1)
}
