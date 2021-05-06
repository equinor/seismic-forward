// $Id: segygeometry_line_test.cpp 1446 2017-02-22 11:45:24Z perroe $

/// \file Unit tests for handling a single line in the SegyGeometry functions in the NRLib SEGY library.

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


SegyGeometry SetupGeometry(const std::vector<PosInput>& input)
{
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


/// \param il_dir If true, line goes in IL direction, else XL direction
/// \param dir    Direction:  0    1     2     3    4      5     ...
///               Degrees:    0   26.6   45   63.4  90   116.6   ...
void GenerateGeometry(bool il_dir, int dir,
                      SegyGeometry& geometry,
                      std::vector<size_t>& i_indexes, std::vector<size_t>& j_indexes)
{
  double utmx0 = 122;
  double utmy0 = 1654;
  int    il0 = 1;
  int    xl0 = 1;
  int    n = 6;
  int    i_start = 0;
  int    j_start = 0;
  int    i_inc = 1;
  int    j_inc = 1;
  double utmx_inc = 0.0;
  double utmy_inc = 0.0;

  std::vector<PosInput> inputs;
  i_indexes.resize(0);
  j_indexes.resize(0);

  int il_inc = 0;
  int xl_inc = 0;
  if (il_dir)
    il_inc = 1;
  else
    xl_inc = 1;

  if (dir >= 7 && dir <= 10)
    i_start = n - 1;

  if (dir >= 11 && dir <= 14)
    j_start = n - 1;

  switch (dir % 8) {
  case 0:
    i_inc = 1;
    j_inc = 0;
    utmx_inc = 50.0;
    utmy_inc = 0.0;
    break;
  case 1:
    i_inc = 1;
    j_inc = 0;
    utmx_inc = 50.0;
    utmy_inc = 25.0;
    break;
  case 2:
    i_inc = 1;
    j_inc = 0;
    utmx_inc = 50.0;
    utmy_inc = 50.0;
    break;
  case 3:
    i_inc = 0;
    j_inc = 1;
    utmx_inc = 25.0;
    utmy_inc = 50.0;
    break;
  case 4:
    i_inc = 0;
    j_inc = 1;
    utmx_inc = 0.0;
    utmy_inc = 50.0;
    break;
  case 5:
    i_inc = 0;
    j_inc = 1;
    utmx_inc = -25.0;
    utmy_inc = 50.0;
    break;
  case 6:
    i_inc = 0;
    j_inc = 1;
    utmx_inc = -50.0;
    utmy_inc = 50.0;
    break;
  case 7:
    i_inc = -1;
    j_inc = 0;
    utmx_inc = -50.0;
    utmy_inc = 25.0;
    break;
  }

  if (dir / 8 == 1) {
    i_inc = -i_inc;
    j_inc = -j_inc;
    utmx_inc = -utmx_inc;
    utmy_inc = -utmy_inc;
  }

  for (int i = 0; i < n; ++i)
  {
    PosInput input;
    input.il = il0 + il_inc * i;
    input.xl = xl0 + xl_inc * i;
    input.x = utmx0 + i * utmx_inc;
    input.y = utmy0 + i * utmy_inc;

    inputs.push_back(input);

    i_indexes.push_back(static_cast<size_t>(i_start + i*i_inc));
    j_indexes.push_back(static_cast<size_t>(j_start + i*j_inc));
  }

  geometry = SetupGeometry(inputs);
}


template<bool il_dir, int dir>
void test_function()
{
  SegyGeometry geometry;
  std::vector<size_t> i_indexes;
  std::vector<size_t> j_indexes;

  GenerateGeometry(il_dir, dir, geometry, i_indexes, j_indexes);

  int IL, XL;
  if (il_dir) {
    BOOST_CHECK_EQUAL(geometry.GetMinXL(), geometry.GetMaxXL());

    for (size_t i = 0; i < i_indexes.size(); ++i) {
      geometry.FindILXL(i_indexes[i], j_indexes[i], IL, XL);
      BOOST_CHECK_EQUAL(IL, static_cast<int>(i) + 1);
      BOOST_CHECK_EQUAL(XL, 1);
    }
  }
  else {
    BOOST_CHECK_EQUAL(geometry.GetMinIL(), geometry.GetMaxIL());

    for (size_t i = 0; i < i_indexes.size(); ++i) {
      geometry.FindILXL(i_indexes[i], j_indexes[i], IL, XL);
      BOOST_CHECK_EQUAL(IL, 1);
      BOOST_CHECK_EQUAL(XL, static_cast<int>(i) + 1);
    }
  }
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_0DegIL) {
  test_function<true, 0>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_0DegXL) {
  test_function<false, 0>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_23DegIL) {
  test_function<true, 1>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_23DegXL) {
  test_function<false, 1>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_45DegIL) {
  test_function<true, 2>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_45DegXL) {
  test_function<false, 2>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_67DegIL) {
  test_function<true, 3>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_67DegXL) {
  test_function<false, 3>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_90DegIL) {
  test_function<true, 4>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_90DegXL) {
  test_function<false, 4>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_113DegIL) {
  test_function<true, 5>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_113DegXL) {
  test_function<false, 5>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_135DegIL) {
  test_function<true, 6>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_135DegXL) {
  test_function<false, 6>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_157DegIL) {
  test_function<true, 7>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_157DegXL) {
  test_function<false, 7>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_180DegIL) {
  test_function<true, 8>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_180DegXL) {
  test_function<false, 8>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_203DegIL) {
  test_function<true, 9>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_203DegXL) {
  test_function<false, 9>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_225DegIL) {
  test_function<true, 10>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_225DegXL) {
  test_function<false, 10>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_247DegIL) {
  test_function<true, 11>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_247DegXL) {
  test_function<false, 11>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_270DegIL) {
  test_function<true, 12>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_270DegXL) {
  test_function<false, 12>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_293DegIL) {
  test_function<true, 13>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_293DegXL) {
  test_function<false, 13>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_315DegIL) {
  test_function<true, 14>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_315DegXL) {
  test_function<false, 14>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_337DegIL) {
  test_function<true, 15>();
}

BOOST_AUTO_TEST_CASE(SegygeometryLine_337DegXL) {
  test_function<false, 15>();
}
