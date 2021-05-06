// $Id: segygeometry_test.cpp 1433 2017-02-20 09:28:00Z perroe $

/// \file Unit tests for the SegyGeometry functions in the NRLib SEGY library.

#include <nrlib/segy/segygeometry.hpp>
#include <nrlib/segy/segytrace.hpp>
#include <boost/test/unit_test.hpp>

using namespace NRLib;

struct PosInput
{
  int il;
  int xl;
  double x;
  double y;
};

SegyGeometry TestGeometry1()
{
  PosInput input[] = { { 5, 7, 0.0, 0.0 },
                       { 5, 8, 0.0, 25.0 },
                       { 5, 9, 0.0, 50.0 },
                       { 6, 7, 25.0, 0.0 },
                       { 6, 8, 25.0, 25.0 },
                       { 6, 9, 25.0, 50.0 } };

  std::vector<SegYTrace *> traces(6);

  for (int i = 0; i < 6; ++i) {
    TraceHeader header;

    header.SetInline(input[i].il);
    header.SetCrossline(input[i].xl);
    header.SetUtmx(input[i].x);
    header.SetUtmy(input[i].y);

    traces[i] = new SegYTrace(header);
  }

  SegyGeometry geometry(traces);

  for (size_t i = 0; i < traces.size(); ++i)
    delete traces[i];

  return geometry;
}

SegyGeometry TestGeometry2()
{
  PosInput input[] = { { 10, 7, 0.0, 0.0 },
                       { 10, 8, 0.0, 25.0 },
                       { 10, 9, 0.0, 50.0 },
                       { 8, 7, 25.0, 0.0 },
                       { 8, 8, 25.0, 25.0 },
                       { 8, 9, 25.0, 50.0 } };

  std::vector<SegYTrace *> traces(6);

  for (int i = 0; i < 6; ++i) {
    TraceHeader header;

    header.SetInline(input[i].il);
    header.SetCrossline(input[i].xl);
    header.SetUtmx(input[i].x);
    header.SetUtmy(input[i].y);

    traces[i] = new SegYTrace(header);
  }

  SegyGeometry geometry(traces);

  for (size_t i = 0; i < traces.size(); ++i)
    delete traces[i];

  return geometry;
}


SegyGeometry TestGeometry3()
{
  PosInput input[] = { { 2, 1, 0.0, 0.0 },
                       { 3, 1, 0.0, 25.0 },
                       { 4, 1, 0.0, 50.0 },
                       { 2, 3, 25.0, 0.0 },
                       { 3, 3, 25.0, 25.0 },
                       { 4, 3, 25.0, 50.0 } };

  std::vector<SegYTrace *> traces(6);

  for (int i = 0; i < 6; ++i) {
    TraceHeader header;

    header.SetInline(input[i].il);
    header.SetCrossline(input[i].xl);
    header.SetUtmx(input[i].x);
    header.SetUtmy(input[i].y);

    traces[i] = new SegYTrace(header);
  }

  SegyGeometry geometry(traces);

  for (size_t i = 0; i < traces.size(); ++i)
    delete traces[i];

  return geometry;
}


SegyGeometry TestGeometry4()
{
  PosInput input[] = { { 20, 1, 0.0, 0.0 },
                       { 10, 1, 0.0, 25.0 },
                       {  0, 1, 0.0, 50.0 },
                       { 20, 3, 25.0, 0.0 },
                       { 10, 3, 25.0, 25.0 },
                       {  0, 3, 25.0, 50.0 } };

  std::vector<SegYTrace *> traces(6);

  for (int i = 0; i < 6; ++i) {
    TraceHeader header;

    header.SetInline(input[i].il);
    header.SetCrossline(input[i].xl);
    header.SetUtmx(input[i].x);
    header.SetUtmy(input[i].y);

    traces[i] = new SegYTrace(header);
  }

  SegyGeometry geometry(traces);

  for (size_t i = 0; i < traces.size(); ++i)
    delete traces[i];

  return geometry;
}


SegyGeometry TestGeometry5()
{
  PosInput input[] = { { 7, 1,  0.0, 0.0 },
                       { 6, 3, 25.0, 25.0 },
                       { 5, 9, 100.0, 50.0 },
                       { 5, 5, 50.0, 50.0 },
                       { 7, 9, 100.0, 0.0 },
                       { 3, 9, 100.0, 100.0 } };

  std::vector<SegYTrace *> traces(6);

  for (int i = 0; i < 6; ++i) {
    TraceHeader header;

    header.SetInline(input[i].il);
    header.SetCrossline(input[i].xl);
    header.SetUtmx(input[i].x);
    header.SetUtmy(input[i].y);

    traces[i] = new SegYTrace(header);
  }

  SegyGeometry geometry(traces);

  for (size_t i = 0; i < traces.size(); ++i)
    delete traces[i];

  return geometry;
}

SegyGeometry TestGeometry6()          //Single line at atan(2) = 63.4 deg relative to UTM x
{
  PosInput input[] = { { 1, 1, 0.0, 0.0 },
                       { 2, 1, 25.0, 50.0 },
                       { 3, 1, 50.0, 100.0 },
                       { 4, 1, 75.0, 150.0 },
                       { 5, 1, 100.0, 200.0 },
                       { 6, 1, 125.0, 250.0 } };

  std::vector<SegYTrace *> traces(6);

  for (int i = 0; i < 6; ++i) {
    TraceHeader header;

    header.SetInline(input[i].il);
    header.SetCrossline(input[i].xl);
    header.SetUtmx(input[i].x);
    header.SetUtmy(input[i].y);

    traces[i] = new SegYTrace(header);
  }

  SegyGeometry geometry(traces);

  for (size_t i = 0; i < traces.size(); ++i)
    delete traces[i];

  return geometry;
}

SegyGeometry TestGeometry7()          //Single line at atan(0.5) = 26.6 deg relative to UTM x
{
  PosInput input[] = { { 1, 1, 0.0, 0.0 },
                       { 2, 1, 50.0, 25.0 },
                       { 3, 1, 100.0, 50.0 },
                       { 4, 1, 150.0, 75.0 },
                       { 5, 1, 200.0, 100.0 },
                       { 6, 1, 250.0, 125.0 } };

  std::vector<SegYTrace *> traces(6);

  for (int i = 0; i < 6; ++i) {
    TraceHeader header;

    header.SetInline(input[i].il);
    header.SetCrossline(input[i].xl);
    header.SetUtmx(input[i].x);
    header.SetUtmy(input[i].y);

    traces[i] = new SegYTrace(header);
  }

  SegyGeometry geometry(traces);

  for (size_t i = 0; i < traces.size(); ++i)
    delete traces[i];

  return geometry;
}

SegyGeometry TestGeometry8()          //Single trace
{
  PosInput input[] = { { 1, 1, 0.0, 0.0 } };

  std::vector<SegYTrace *> traces(1);

  for (int i = 0; i < 1; ++i) {
    TraceHeader header;

    header.SetInline(input[i].il);
    header.SetCrossline(input[i].xl);
    header.SetUtmx(input[i].x);
    header.SetUtmy(input[i].y);

    traces[i] = new SegYTrace(header);
  }

  SegyGeometry geometry(traces);

  for (size_t i = 0; i < traces.size(); ++i)
    delete traces[i];

  return geometry;
}

SegyGeometry TestGeometry9()          //Single line at atan(1) = 45 deg relative to UTM x
{
  PosInput input[] = { { 1, 1, 0.0, 0.0 },
                       { 2, 1, 25.0, 25.0 },
                       { 3, 1, 50.0, 50.0 },
                       { 4, 1, 75.0, 75.0 },
                       { 5, 1, 100.0, 100.0 },
                       { 6, 1, 125.0, 125.0 } };

  std::vector<SegYTrace *> traces(6);

  for (int i = 0; i < 6; ++i) {
    TraceHeader header;

    header.SetInline(input[i].il);
    header.SetCrossline(input[i].xl);
    header.SetUtmx(input[i].x);
    header.SetUtmy(input[i].y);

    traces[i] = new SegYTrace(header);
  }

  SegyGeometry geometry(traces);

  for (size_t i = 0; i < traces.size(); ++i)
    delete traces[i];

  return geometry;
}

BOOST_AUTO_TEST_CASE( FindILXLFromIJTest1 )
{
  SegyGeometry geometry = TestGeometry1();

  int IL, XL;
  geometry.FindILXL(size_t(0), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 5);
  BOOST_CHECK_EQUAL(XL, 7);

  geometry.FindILXL(size_t(0), size_t(1), IL, XL);
  BOOST_CHECK_EQUAL(IL, 5);
  BOOST_CHECK_EQUAL(XL, 8);

  geometry.FindILXL(size_t(0), size_t(2), IL, XL);
  BOOST_CHECK_EQUAL(IL, 5);
  BOOST_CHECK_EQUAL(XL, 9);

  geometry.FindILXL(size_t(1), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 6);
  BOOST_CHECK_EQUAL(XL, 7);

  geometry.FindILXL(size_t(1), size_t(1), IL, XL);
  BOOST_CHECK_EQUAL(IL, 6);
  BOOST_CHECK_EQUAL(XL, 8);

  geometry.FindILXL(size_t(1), size_t(2), IL, XL);
  BOOST_CHECK_EQUAL(IL, 6);
  BOOST_CHECK_EQUAL(XL, 9);
}


BOOST_AUTO_TEST_CASE(FindILXLFromIJTest2)
{
  SegyGeometry geometry = TestGeometry2();

  int IL, XL;
  geometry.FindILXL(size_t(0), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 10);
  BOOST_CHECK_EQUAL(XL, 7);

  geometry.FindILXL(size_t(1), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 8);
  BOOST_CHECK_EQUAL(XL, 7);

  geometry.FindILXL(size_t(0), size_t(1), IL, XL);
  BOOST_CHECK_EQUAL(IL, 10);
  BOOST_CHECK_EQUAL(XL, 8);

  geometry.FindILXL(size_t(1), size_t(1), IL, XL);
  BOOST_CHECK_EQUAL(IL, 8);
  BOOST_CHECK_EQUAL(XL, 8);

  geometry.FindILXL(size_t(0), size_t(2), IL, XL);
  BOOST_CHECK_EQUAL(IL, 10);
  BOOST_CHECK_EQUAL(XL, 9);

  geometry.FindILXL(size_t(1), size_t(2), IL, XL);
  BOOST_CHECK_EQUAL(IL, 8);
  BOOST_CHECK_EQUAL(XL, 9);
}


BOOST_AUTO_TEST_CASE(FindILXLFromIJTest3)
{
  SegyGeometry geometry = TestGeometry3();

  int IL, XL;
  geometry.FindILXL(size_t(0), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 2);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(0), size_t(1), IL, XL);
  BOOST_CHECK_EQUAL(IL, 3);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(0), size_t(2), IL, XL);
  BOOST_CHECK_EQUAL(IL, 4);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(1), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 2);
  BOOST_CHECK_EQUAL(XL, 3);

  geometry.FindILXL(size_t(1), size_t(1), IL, XL);
  BOOST_CHECK_EQUAL(IL, 3);
  BOOST_CHECK_EQUAL(XL, 3);

  geometry.FindILXL(size_t(1), size_t(2), IL, XL);
  BOOST_CHECK_EQUAL(IL, 4);
  BOOST_CHECK_EQUAL(XL, 3);
}


BOOST_AUTO_TEST_CASE(FindILXLFromIJTest4)
{
  SegyGeometry geometry = TestGeometry4();

  int IL, XL;
  geometry.FindILXL(size_t(0), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 20);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(0), size_t(1), IL, XL);
  BOOST_CHECK_EQUAL(IL, 10);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(0), size_t(2), IL, XL);
  BOOST_CHECK_EQUAL(IL, 0);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(1), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 20);
  BOOST_CHECK_EQUAL(XL, 3);

  geometry.FindILXL(size_t(1), size_t(1), IL, XL);
  BOOST_CHECK_EQUAL(IL, 10);
  BOOST_CHECK_EQUAL(XL, 3);

  geometry.FindILXL(size_t(1), size_t(2), IL, XL);
  BOOST_CHECK_EQUAL(IL, 0);
  BOOST_CHECK_EQUAL(XL, 3);
}


BOOST_AUTO_TEST_CASE(FindILXLFromIJTest5)
{
  SegyGeometry geometry = TestGeometry5();

  int IL, XL;
  geometry.FindILXL(size_t(0), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 7);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(1), size_t(1), IL, XL);
  BOOST_CHECK_EQUAL(IL, 6);
  BOOST_CHECK_EQUAL(XL, 3);

  geometry.FindILXL(size_t(2), size_t(2), IL, XL);
  BOOST_CHECK_EQUAL(IL, 5);
  BOOST_CHECK_EQUAL(XL, 5);

  geometry.FindILXL(size_t(4), size_t(2), IL, XL);
  BOOST_CHECK_EQUAL(IL, 5);
  BOOST_CHECK_EQUAL(XL, 9);

  geometry.FindILXL(size_t(4), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 7);
  BOOST_CHECK_EQUAL(XL, 9);

  geometry.FindILXL(size_t(4), size_t(4), IL, XL);
  BOOST_CHECK_EQUAL(IL, 3);
  BOOST_CHECK_EQUAL(XL, 9);
}


BOOST_AUTO_TEST_CASE(FindILXLFromIJTest6)
{
  SegyGeometry geometry = TestGeometry6();

  int IL, XL;
  geometry.FindILXL(size_t(0), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 1);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(0), size_t(1), IL, XL);
  BOOST_CHECK_EQUAL(IL, 2);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(0), size_t(2), IL, XL);
  BOOST_CHECK_EQUAL(IL, 3);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(0), size_t(3), IL, XL);
  BOOST_CHECK_EQUAL(IL, 4);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(0), size_t(4), IL, XL);
  BOOST_CHECK_EQUAL(IL, 5);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(0), size_t(5), IL, XL);
  BOOST_CHECK_EQUAL(IL, 6);
  BOOST_CHECK_EQUAL(XL, 1);
}

BOOST_AUTO_TEST_CASE(FindILXLFromIJTest7)
{
  SegyGeometry geometry = TestGeometry7();

  int IL, XL;
  geometry.FindILXL(size_t(0), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 1);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(1), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 2);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(2), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 3);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(3), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 4);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(4), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 5);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(5), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 6);
  BOOST_CHECK_EQUAL(XL, 1);
}

BOOST_AUTO_TEST_CASE(FindILXLFromIJTest8)
{
  SegyGeometry geometry = TestGeometry8();

  int IL, XL;
  geometry.FindILXL(size_t(0), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 1);
  BOOST_CHECK_EQUAL(XL, 1);
}

BOOST_AUTO_TEST_CASE(FindILXLFromIJTest9)
{
  SegyGeometry geometry = TestGeometry9();

  int IL, XL;
  geometry.FindILXL(size_t(0), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 1);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(1), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 2);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(2), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 3);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(3), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 4);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(4), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 5);
  BOOST_CHECK_EQUAL(XL, 1);

  geometry.FindILXL(size_t(5), size_t(0), IL, XL);
  BOOST_CHECK_EQUAL(IL, 6);
  BOOST_CHECK_EQUAL(XL, 1);
}

BOOST_AUTO_TEST_CASE(FindContIndexFromContILXLTest1)
{
  SegyGeometry geometry = TestGeometry1();

  double i, j;
  const double tolerance = 0.00001;

  geometry.FindContIndexFromContILXL(5.0, 7.0, i, j);
  BOOST_CHECK_CLOSE(i, 0.5, tolerance);
  BOOST_CHECK_CLOSE(j, 0.5, tolerance);

  geometry.FindContIndexFromContILXL(6.0, 9.0, i, j);
  BOOST_CHECK_CLOSE(i, 1.5, tolerance);
  BOOST_CHECK_CLOSE(j, 2.5, tolerance);

  geometry.FindContIndexFromContILXL(7.0, 10.0, i, j);
  BOOST_CHECK_CLOSE(i, 2.5, tolerance);
  BOOST_CHECK_CLOSE(j, 3.5, tolerance);
}


BOOST_AUTO_TEST_CASE(FindIndexXYToIJLargeGrid)
{
  SegyGeometry geometry(406059.73088598251, 6566867.8688883781, 12.504052816501131, 12.499931945166741,
                        3611, 7107, 15953.499590131056, 11050.500071772763,
                        0.032895441665889015, 0.072924340216817499, -0.072903081774711609,
                        0.032878458499908447, -0.42375998458917319);

  size_t n = 8;
  int i_indexes[] = { 0, 1, 42, 360, 1805, 3000, 3609, 3610 };
  int j_indexes[] = { 0, 1, 99, 1500, 3553, 5555, 7105, 7106 };

  BOOST_CHECK_EQUAL(geometry.GetILStep(), 1);
  BOOST_CHECK_EQUAL(geometry.GetXLStep(), 1);
  BOOST_CHECK_EQUAL(geometry.GetMinIL(), 15954);
  BOOST_CHECK_EQUAL(geometry.GetMaxIL(), 23060);
  BOOST_CHECK_EQUAL(geometry.GetMinXL(), 7440);
  BOOST_CHECK_EQUAL(geometry.GetMaxXL(), 11050);

  for (size_t ii = 0; ii < n; ++ii) {
    for (size_t jj = 0; jj < n; ++jj) {
      double x, y;
      int il, xl;
      size_t i_in = i_indexes[ii];
      size_t j_in = j_indexes[jj];

      geometry.FindXYFromIJ(i_in, j_in, x, y);

      geometry.FindILXL(x, y, il, xl);
      int il2, xl2;
      geometry.FindILXL(i_in, j_in, il2, xl2);
      BOOST_CHECK_EQUAL(il, il2);
      BOOST_CHECK_EQUAL(xl, xl2);

      size_t i, j;
      BOOST_CHECK_NO_THROW(geometry.FindIndex(x, y, i, j));
      BOOST_CHECK_EQUAL(i, i_in);
      BOOST_CHECK_EQUAL(j, j_in);
    }
  }

  size_t i, j;
  BOOST_CHECK_THROW(geometry.FindIndex(406056.0, 6568865.0, i, j), NRLib::Exception);
  BOOST_CHECK_THROW(geometry.FindIndex(442600.0, 6647845.0, i, j), NRLib::Exception);
}

