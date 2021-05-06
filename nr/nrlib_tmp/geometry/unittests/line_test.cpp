/// \file Unit tests for the line function in geometry

#include <nrlib/geometry/line.hpp>

#include <boost/test/unit_test.hpp>
//#include <math.h>
//#include <algorithm>

using namespace NRLib;

BOOST_AUTO_TEST_CASE ( OneLineNormalToXYPlane )
{
  Point pt1(2, 2, -1);
  Point pt2(2, 2, 1);
  Point pt3(-2, 1, 0);
  Point pt4(2, 1, 0);
  Line def_line(pt1, pt2);
  Line line_in(pt3, pt4);
  Point intersect_pt;
  bool is_intersect = def_line.IntersectXY(line_in, intersect_pt);
  BOOST_CHECK_EQUAL(is_intersect, false);
  BOOST_CHECK_EQUAL(intersect_pt, Point(0, 0, 0));
}

BOOST_AUTO_TEST_CASE( LinesAreCrossing )
{
  Point pt1(0, 0, 0);
  Point pt2(0, 2, 0);
  Point pt3(-2, 1, 0);
  Point pt4(2, 1, 0);
  Line def_line(pt1, pt2);
  Line line_in(pt3, pt4);
  Point intersect_pt;
  bool is_intersect = def_line.IntersectXY(line_in, intersect_pt);
  BOOST_CHECK_EQUAL(is_intersect, true);
  BOOST_CHECK_EQUAL(intersect_pt, Point(0, 1, 0));
}

BOOST_AUTO_TEST_CASE( LinesAreNotCrossing )
{
  Point pt1(0, 0, 0);
  Point pt2(0, 2, 0);
  Point pt3(-2, 3, 0);
  Point pt4(2, 3, 0);
  Line def_line(pt1, pt2);
  Line line_in(pt3, pt4);
  Point intersect_pt;
  bool is_intersect = def_line.IntersectXY(line_in, intersect_pt);
  BOOST_CHECK_EQUAL(is_intersect, false);
  BOOST_CHECK_EQUAL(intersect_pt, Point(0, 0, 0));
}
