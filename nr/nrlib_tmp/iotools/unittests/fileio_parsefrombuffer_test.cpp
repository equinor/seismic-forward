// $Id: fileio_parsefrombuffer_test.cpp 1584 2017-07-03 12:28:42Z perroe $

/// \file Unit tests for the FileIO functions in the NRLib IOTools library.

#include <nrlib/iotools/fileio.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <cmath>
#include <limits>

using namespace NRLib;
using namespace NRLib::NRLibPrivate;

BOOST_AUTO_TEST_CASE( ParseUInt16BETest )
{
  unsigned short us;

  // 0
  unsigned char buffer[2] = {0x00, 0x00};
  char* buf = reinterpret_cast<char*>(&buffer);
  ParseUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(us, 0U);

  // 1
  buffer[0] = 0x00; buffer[1] = 0x01;
  ParseUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(us, 1U);

  // 256
  buffer[0] = 0x01; buffer[1] = 0x00;
  ParseUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(us, 256U);

  // 65535
  buffer[0] = 0xff; buffer[1] = 0xff;
  ParseUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(us, 65535U);

  // 32767
  buffer[0] = 0x7f; buffer[1] = 0xff;
  ParseUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(us, 32767U);

  // 31415
  buffer[0] = 0x7a; buffer[1] = 0xb7;
  ParseUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(us, 31415U);
}


BOOST_AUTO_TEST_CASE( ParseUInt16LETest )
{
  unsigned short us;

  // 0
  unsigned char buffer[2] = {0x00, 0x00};
  char* buf = reinterpret_cast<char*>(&buffer);
  ParseUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(us, 0U);

  // 1
  buffer[0] = 0x01; buffer[1] = 0x00;
  ParseUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(us, 1U);

  // 256
  buffer[0] = 0x00; buffer[1] = 0x01;
  ParseUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(us, 256U);

  // 65535
  buffer[0] = 0xff; buffer[1] = 0xff;
  ParseUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(us, 65535U);

  // 32767
  buffer[0] = 0xff; buffer[1] = 0x7f;
  ParseUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(us, 32767U);

  // 31415
  buffer[0] = 0xb7; buffer[1] = 0x7a;
  ParseUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(us, 31415U);
}


BOOST_AUTO_TEST_CASE(ParseUInt32BETest)
{
  unsigned int ui;

  // 0
  unsigned char buffer[4] = { 0x00, 0x00, 0x00, 0x00 };
  char* buf = reinterpret_cast<char*>(&buffer);
  ParseUInt32BE(buf, ui);
  BOOST_CHECK_EQUAL(ui, 0U);

  // 1
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x01;
  ParseUInt32BE(buf, ui);
  BOOST_CHECK_EQUAL(ui, 1U);

  // 256
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x01; buffer[3] = 0x00;
  ParseUInt32BE(buf, ui);
  BOOST_CHECK_EQUAL(ui, 256U);

  // 65536
  buffer[0] = 0x00; buffer[1] = 0x01; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseUInt32BE(buf, ui);
  BOOST_CHECK_EQUAL(ui, 65536U);

  // 16777216
  buffer[0] = 0x01; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseUInt32BE(buf, ui);
  BOOST_CHECK_EQUAL(ui, 16777216U);

  // 4294967295
  buffer[0] = 0xff; buffer[1] = 0xff; buffer[2] = 0xff; buffer[3] = 0xff;
  ParseUInt32BE(buf, ui);
  BOOST_CHECK_EQUAL(ui, 4294967295U);

  // 2147483647
  buffer[0] = 0x7f; buffer[1] = 0xff; buffer[2] = 0xff; buffer[3] = 0xff;
  ParseUInt32BE(buf, ui);
  BOOST_CHECK_EQUAL(ui, 2147483647U);
}


BOOST_AUTO_TEST_CASE(ParseUInt32LETest)
{
  unsigned int ui;

  // 0
  unsigned char buffer[4] = { 0x00, 0x00, 0x00, 0x00 };
  char* buf = reinterpret_cast<char*>(&buffer);
  ParseUInt32LE(buf, ui);
  BOOST_CHECK_EQUAL(ui, 0U);

  // 1
  buffer[0] = 0x01; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseUInt32LE(buf, ui);
  BOOST_CHECK_EQUAL(ui, 1U);

  // 256
  buffer[0] = 0x00; buffer[1] = 0x01; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseUInt32LE(buf, ui);
  BOOST_CHECK_EQUAL(ui, 256U);

  // 65536
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x01; buffer[3] = 0x00;
  ParseUInt32LE(buf, ui);
  BOOST_CHECK_EQUAL(ui, 65536U);

  // 16777216
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x01;
  ParseUInt32LE(buf, ui);
  BOOST_CHECK_EQUAL(ui, 16777216U);

  // 4294967295
  buffer[0] = 0xff; buffer[1] = 0xff; buffer[2] = 0xff; buffer[3] = 0xff;
  ParseUInt32LE(buf, ui);
  BOOST_CHECK_EQUAL(ui, 4294967295U);

  // 2147483647
  buffer[0] = 0xff; buffer[1] = 0xff; buffer[2] = 0xff; buffer[3] = 0x7f;
  ParseUInt32LE(buf, ui);
  BOOST_CHECK_EQUAL(ui, 2147483647U);
}


BOOST_AUTO_TEST_CASE(ParseInt16BETest)
{
  short in;

  // 0
  unsigned char buffer[2] = { 0x00, 0x00 };
  char* buf = reinterpret_cast<char*>(&buffer);
  ParseInt16BE(buf, in);
  BOOST_CHECK_EQUAL(in, 0);

  // -1
  buffer[0] = 0xff; buffer[1] = 0xff;
  ParseInt16BE(buf, in);
  BOOST_CHECK_EQUAL(in, -1);

  // 1
  buffer[0] = 0x00; buffer[1] = 0x01;
  ParseInt16BE(buf, in);
  BOOST_CHECK_EQUAL(in, 1);

  // 256
  buffer[0] = 0x01; buffer[1] = 0x00;
  ParseInt16BE(buf, in);
  BOOST_CHECK_EQUAL(in, 256);

  // 32767
  buffer[0] = 0x7f; buffer[1] = 0xff;
  ParseInt16BE(buf, in);
  BOOST_CHECK_EQUAL(in, 32767);

  // -32768
  buffer[0] = 0x80; buffer[1] = 0x00;
  ParseInt16BE(buf, in);
  BOOST_CHECK_EQUAL(in, -32768);

  // 16777
  buffer[0] = 0x41; buffer[1] = 0x89;
  ParseInt16BE(buf, in);
  BOOST_CHECK_EQUAL(in, 16777);

  // 31415
  buffer[0] = 0x7a; buffer[1] = 0xb7;
  ParseInt16BE(buf, in);
  BOOST_CHECK_EQUAL(in, 31415);
}


BOOST_AUTO_TEST_CASE(ParseInt32BETest)
{
  int in;

  // 0
  unsigned char buffer[4] = { 0x00, 0x00, 0x00, 0x00 };
  char* buf = reinterpret_cast<char*>(&buffer);
  ParseInt32BE(buf, in);
  BOOST_CHECK_EQUAL(in, 0);

  // -1
  buffer[0] = 0xff; buffer[1] = 0xff; buffer[2] = 0xff; buffer[3] = 0xff;
  ParseInt32BE(buf, in);
  BOOST_CHECK_EQUAL(in, -1);

  // 1
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x01;
  ParseInt32BE(buf, in);
  BOOST_CHECK_EQUAL(in, 1);

  // 256
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x01; buffer[3] = 0x00;
  ParseInt32BE(buf, in);
  BOOST_CHECK_EQUAL(in, 256);

  // 2147483647
  buffer[0] = 0x7f; buffer[1] = 0xff; buffer[2] = 0xff; buffer[3] = 0xff;
  ParseInt32BE(buf, in);
  BOOST_CHECK_EQUAL(in, std::numeric_limits<int>::max());

  // -2147483648
  buffer[0] = 0x80; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseInt32BE(buf, in);
  BOOST_CHECK_EQUAL(in, std::numeric_limits<int>::min());

  // 16777216
  buffer[0] = 0x01; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseInt32BE(buf, in);
  BOOST_CHECK_EQUAL(in, 16777216);

  // 314159265
  buffer[0] = 0x12; buffer[1] = 0xb9; buffer[2] = 0xb0; buffer[3] = 0xa1;
  ParseInt32BE(buf, in);
  BOOST_CHECK_EQUAL(in, 314159265);
}


BOOST_AUTO_TEST_CASE(ParseIEEEFloatBETest)
{
  float f;
  // 0
  unsigned char buffer[4] = { 0x00, 0x00, 0x00, 0x00 };
  char* buf = reinterpret_cast<char*>(&buffer);
  ParseIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(f, 0.0F);

  // -0
  buffer[0] = 0x80; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(f, -0.0F);

  // 1
  buffer[0] = 0x3f; buffer[1] = 0x80; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(f, 1.0F);

  // -1
  buffer[0] = 0xbf; buffer[1] = 0x80; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(f, -1.0F);

  // Smallest denormalized number
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x01;
  ParseIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(f, std::numeric_limits<float>::denorm_min());

  // Middle denormalized number
  buffer[0] = 0x80; buffer[1] = 0x40; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseIEEEFloatBE(buf, f);
  BOOST_CHECK_CLOSE_FRACTION(f, -5.9e-39F, 0.1);

  // Largest denormalized number
  buffer[0] = 0x00; buffer[1] = 0x7f; buffer[2] = 0xff; buffer[3] = 0xff;
  ParseIEEEFloatBE(buf, f);
  BOOST_CHECK_CLOSE_FRACTION(f, 1.175494e-38F, 0.00001);

  // Smallest normalized number
  buffer[0] = 0x80; buffer[1] = 0x80; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseIEEEFloatBE(buf, f);
  BOOST_CHECK_CLOSE_FRACTION(f, -1.175494e-38, 0.00001);

  // Largest normalized number
  buffer[0] = 0x7f; buffer[1] = 0x7f; buffer[2] = 0xff; buffer[3] = 0xff;
  ParseIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(f, std::numeric_limits<float>::max());

  // +infinity
  buffer[0] = 0x7f; buffer[1] = 0x80; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(f, std::numeric_limits<float>::infinity());

  // -infinity
  buffer[0] = 0xff; buffer[1] = 0x80; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(f, -std::numeric_limits<float>::infinity());

  // Nan
  buffer[0] = 0x7f; buffer[1] = 0xff; buffer[2] = 0x01; buffer[3] = 0x00;
  ParseIEEEFloatBE(buf, f);
  BOOST_CHECK(boost::math::isnan(f));

  // Nan
  buffer[0] = 0xff; buffer[1] = 0xff; buffer[2] = 0xff; buffer[3] = 0xff;
  ParseIEEEFloatBE(buf, f);
  BOOST_CHECK(boost::math::isnan(f));

  // Nan
  buffer[0] = 0xff; buffer[1] = 0xc0; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseIEEEFloatBE(buf, f);
  BOOST_CHECK(boost::math::isnan(f));
}


BOOST_AUTO_TEST_CASE(ParseIEEEFloatLETest)
{
  float f;
  // 0
  unsigned char buffer[4] = { 0x00, 0x00, 0x00, 0x00 };
  char* buf = reinterpret_cast<char*>(&buffer);
  ParseIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(f, 0.0F);

  // -0
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x80;
  ParseIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(f, -0.0F);

  // 1
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x80; buffer[3] = 0x3f;
  ParseIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(f, 1.0F);

  // -1
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x80; buffer[3] = 0xbf;
  ParseIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(f, -1.0F);

  // Smallest denormalized number
  buffer[0] = 0x01; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(f, std::numeric_limits<float>::denorm_min());

  // Middle denormalized number
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x40; buffer[3] = 0x80;
  ParseIEEEFloatLE(buf, f);
  BOOST_CHECK_CLOSE_FRACTION(f, -5.9e-39F, 0.1);

  // Largest denormalized number
  buffer[0] = 0xff; buffer[1] = 0xff; buffer[2] = 0x7f; buffer[3] = 0x00;
  ParseIEEEFloatLE(buf, f);
  BOOST_CHECK_CLOSE_FRACTION(f, 1.17549421e-38F, 0.00001);

  // Smallest normalized number
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x80; buffer[3] = 0x80;
  ParseIEEEFloatLE(buf, f);
  BOOST_CHECK_CLOSE_FRACTION(f, -1.17549421e-38, 0.00001);

  // Largest normalized number
  buffer[0] = 0xff; buffer[1] = 0xff; buffer[2] = 0x7f; buffer[3] = 0x7f;
  ParseIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(f, std::numeric_limits<float>::max());

  // +infinity
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x80; buffer[3] = 0x7f;
  ParseIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(f, std::numeric_limits<float>::infinity());

  // -infinity
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x80; buffer[3] = 0xff;
  ParseIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(f, -std::numeric_limits<float>::infinity());

  // Nan
  buffer[0] = 0x00; buffer[1] = 0x01; buffer[2] = 0xff; buffer[3] = 0x7f;
  ParseIEEEFloatLE(buf, f);
  BOOST_CHECK(boost::math::isnan(f));

  // Nan
  buffer[0] = 0xff; buffer[1] = 0xff; buffer[2] = 0xff; buffer[3] = 0xff;
  ParseIEEEFloatLE(buf, f);
  BOOST_CHECK(boost::math::isnan(f));

  // Nan
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0xc0; buffer[3] = 0xff;
  ParseIEEEFloatLE(buf, f);
  BOOST_CHECK(boost::math::isnan(f));
}


BOOST_AUTO_TEST_CASE(ParseIEEEDoubleBETest)
{
  double d;
  // 0
  unsigned char buffer[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
  char* buf = reinterpret_cast<char*>(&buffer);
  ParseIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(d, 0.0);

  // -0
  buffer[0] = 0x80; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0x00; buffer[7] = 0x00;
  ParseIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(d, -0.0);

  // 1
  buffer[0] = 0x3f; buffer[1] = 0xf0; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0x00; buffer[7] = 0x00;
  ParseIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(d, 1.0);

  // -1
  buffer[0] = 0xbf; buffer[1] = 0xf0; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0x00; buffer[7] = 0x00;
  ParseIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(d, -1.0);

  // Smallest denormalized number
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0x00; buffer[7] = 0x01;
  ParseIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(d, std::numeric_limits<double>::denorm_min());
  //std::cout << "d = " << d << ", should be ~5e-324.\n";

  // Middle denormalized number
  buffer[0] = 0x80; buffer[1] = 0x08; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0x00; buffer[7] = 0x00;
  ParseIEEEDoubleBE(buf, d);
  BOOST_CHECK_CLOSE_FRACTION(d, -1.1e-308, 0.1);

  // Largest denormalized number
  buffer[0] = 0x00; buffer[1] = 0x0f; buffer[2] = 0xff; buffer[3] = 0xff;
  buffer[4] = 0xff; buffer[5] = 0xff; buffer[6] = 0xff; buffer[7] = 0xff;
  ParseIEEEDoubleBE(buf, d);
  BOOST_CHECK_CLOSE_FRACTION(d, 2.2e-308, 0.1);

  // Smallest normalized number
  buffer[0] = 0x80; buffer[1] = 0x10; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0x00; buffer[7] = 0x00;
  ParseIEEEDoubleBE(buf, d);
  BOOST_CHECK_CLOSE_FRACTION(d, -2.226e-308, 0.1);

  // Largest normalized number
  buffer[0] = 0x7f; buffer[1] = 0xef; buffer[2] = 0xff; buffer[3] = 0xff;
  buffer[4] = 0xff; buffer[5] = 0xff; buffer[6] = 0xff; buffer[7] = 0xff;
  ParseIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(d, std::numeric_limits<double>::max());

  // +infinity
  buffer[0] = 0x7f; buffer[1] = 0xf0; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0x00; buffer[7] = 0x00;
  ParseIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(d, std::numeric_limits<double>::infinity());

  // -infinity
  buffer[0] = 0xff; buffer[1] = 0xf0; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0x00; buffer[7] = 0x00;
  ParseIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(d, -std::numeric_limits<double>::infinity());

  // Nan
  buffer[0] = 0x7f; buffer[1] = 0xff; buffer[2] = 0x01; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0x00; buffer[7] = 0x00;
  ParseIEEEDoubleBE(buf, d);
  BOOST_CHECK(boost::math::isnan(d));

  // Nan
  buffer[0] = 0xff; buffer[1] = 0xff; buffer[2] = 0xff; buffer[3] = 0xff;
  buffer[4] = 0xff; buffer[5] = 0xff; buffer[6] = 0xff; buffer[7] = 0xff;
  ParseIEEEDoubleBE(buf, d);
  BOOST_CHECK(boost::math::isnan(d));

  // "Random" number
  buffer[0] = 0x3d; buffer[1] = 0x8b; buffer[2] = 0xa2; buffer[3] = 0x3c;
  buffer[4] = 0xdd; buffer[5] = 0x51; buffer[6] = 0x22; buffer[7] = 0xb5;
  ParseIEEEDoubleBE(buf, d);
  BOOST_CHECK_CLOSE_FRACTION(d, 3.141592653589793e-12, 1e-20);
}


BOOST_AUTO_TEST_CASE(ParseIEEEDoubleLETest)
{
  double d;
  // 0
  unsigned char buffer[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
  char* buf = reinterpret_cast<char*>(&buffer);
  ParseIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(d, 0.0);

  // -0
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0x00; buffer[7] = 0x80;
  ParseIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(d, -0.0);

  // 1
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0xf0; buffer[7] = 0x3f;
  ParseIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(d, 1.0);

  // -1
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0xf0; buffer[7] = 0xbf;
  ParseIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(d, -1.0);

  // Smallest denormalized number
  buffer[0] = 0x01; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0x00; buffer[7] = 0x00;
  ParseIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(d, std::numeric_limits<double>::denorm_min());

  // Middle denormalized number
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0x08; buffer[7] = 0x80;
  ParseIEEEDoubleLE(buf, d);
  BOOST_CHECK_CLOSE_FRACTION(d, -1.1e-308, 0.1);

  // Largest denormalized number
  buffer[0] = 0xff; buffer[1] = 0xff; buffer[2] = 0xff; buffer[3] = 0xff;
  buffer[4] = 0xff; buffer[5] = 0xff; buffer[6] = 0x0f; buffer[7] = 0x00;
  ParseIEEEDoubleLE(buf, d);
  BOOST_CHECK_CLOSE_FRACTION(d, 2.2e-308, 0.1);

  // Smallest normalized number
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0x10; buffer[7] = 0x80;
  ParseIEEEDoubleLE(buf, d);
  BOOST_CHECK_CLOSE_FRACTION(d, -2.226e-308, 0.1);

  // Largest normalized number
  buffer[0] = 0xff; buffer[1] = 0xff; buffer[2] = 0xff; buffer[3] = 0xff;
  buffer[4] = 0xff; buffer[5] = 0xff; buffer[6] = 0xef; buffer[7] = 0x7f;
  ParseIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(d, std::numeric_limits<double>::max());

  // +infinity
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0xf0; buffer[7] = 0x7f;
  ParseIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(d, std::numeric_limits<double>::infinity());

  // -infinity
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x00; buffer[6] = 0xf0; buffer[7] = 0xff;
  ParseIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(d, -std::numeric_limits<double>::infinity());

  // Nan
  buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  buffer[4] = 0x00; buffer[5] = 0x01; buffer[6] = 0xff; buffer[7] = 0x7f;
  ParseIEEEDoubleLE(buf, d);
  BOOST_CHECK(boost::math::isnan(d));

  // Nan
  buffer[0] = 0xff; buffer[1] = 0xff; buffer[2] = 0xff; buffer[3] = 0xff;
  buffer[4] = 0xff; buffer[5] = 0xff; buffer[6] = 0xff; buffer[7] = 0xff;
  ParseIEEEDoubleLE(buf, d);
  BOOST_CHECK(boost::math::isnan(d));

  // "Random" number
  buffer[0] = 0xb5; buffer[1] = 0x22; buffer[2] = 0x51; buffer[3] = 0xdd;
  buffer[4] = 0x3c; buffer[5] = 0xa2; buffer[6] = 0x8b; buffer[7] = 0x3d;
  ParseIEEEDoubleLE(buf, d);
  BOOST_CHECK_CLOSE_FRACTION(d, 3.141592653589793e-12, 1e-20);
}


BOOST_AUTO_TEST_CASE(ParseIBMFloatBETest)
{
  float f;
  // 0
  unsigned char buffer[4] = { 0x00, 0x00, 0x00, 0x00 };
  char* buf = reinterpret_cast<char*>(&buffer);
  ParseIBMFloatBE(buf, f);
  BOOST_CHECK_EQUAL(f, 0.0F);

  // -0
  //buffer[0] = 0x80; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x00;
  //ParseIBMFloatBE(buf, f);
  //BOOST_CHECK_EQUAL(f, -0.0F);

  // 1
  buffer[0] = 0x41; buffer[1] = 0x10; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseIBMFloatBE(buf, f);
  BOOST_CHECK_EQUAL(f, 1.0F);

  // -1
  buffer[0] = 0xc1; buffer[1] = 0x10; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseIBMFloatBE(buf, f);
  BOOST_CHECK_EQUAL(f, -1.0F);

  // Smallest denormalized number
  //buffer[0] = 0x00; buffer[1] = 0x00; buffer[2] = 0x00; buffer[3] = 0x01;
  //ParseIBMFloatBE(buf, f);
  //BOOST_CHECK_CLOSE_FRACTION(f, 1e-14, 0.00001);

  // Middle denormalized number
  // buffer[0] = 0x80; buffer[1] = 0x40; buffer[2] = 0x00; buffer[3] = 0x00;
  //ParseIBMFloatBE(buf, f);
  //BOOST_CHECK_CLOSE_FRACTION(f, -5.9e-39, 0.1);

  // Largest denormalized number
  //buffer[0] = 0x00; buffer[1] = 0x7f; buffer[2] = 0xff; buffer[3] = 0xff;
  //ParseIBMFloatBE(buf, f);
  //BOOST_CHECK_CLOSE_FRACTION(f, 1.17549421e-38, 0.00001);

  // Smallest normalized number
  //buffer[0] = 0x80; buffer[1] = 0x80; buffer[2] = 0x00; buffer[3] = 0x00;
  //ParseIBMFloatBE(buf, f);
  //BOOST_CHECK_CLOSE_FRACTION(f, -1.17549421e-38, 0.00001);

  // Largest normalized number
  //buffer[0] = 0xef; buffer[1] = 0xff; buffer[2] = 0xff; buffer[3] = 0xff;
  //ParseIBMFloatBE(buf, f);
  //BOOST_CHECK_CLOSE_FRACTION(f, 1e38, 0.00001);

  // -118.625
  buffer[0] = 0xc2; buffer[1] = 0x76; buffer[2] = 0xa0; buffer[3] = 0x00;
  ParseIBMFloatBE(buf, f);
  BOOST_CHECK_CLOSE_FRACTION(f, -118.625F, 0.00001);

  // +infinity (As in current implementation. Maybe not correct per standard. )
  buffer[0] = 0x61; buffer[1] = 0x10; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseIBMFloatBE(buf, f);
  BOOST_CHECK_EQUAL(f, std::numeric_limits<float>::infinity());

  // -infinity (As in current implementation. Maybe not correct per standard. )
  buffer[0] = 0xe1; buffer[1] = 0x10; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseIBMFloatBE(buf, f);
  BOOST_CHECK_EQUAL(f, -std::numeric_limits<float>::infinity());

  // Nan
  buffer[0] = 0x7f; buffer[1] = 0xff; buffer[2] = 0x01; buffer[3] = 0x00;
  ParseIBMFloatBE(buf, f);
  BOOST_CHECK(boost::math::isnan(f));

  // Nan
  buffer[0] = 0xff; buffer[1] = 0xff; buffer[2] = 0xff; buffer[3] = 0xff;
  ParseIBMFloatBE(buf, f);
  BOOST_CHECK(boost::math::isnan(f));

  // Nan
  buffer[0] = 0xff; buffer[1] = 0xc0; buffer[2] = 0x00; buffer[3] = 0x00;
  ParseIBMFloatBE(buf, f);
  BOOST_CHECK(boost::math::isnan(f));

  // 3.14159265e-12
  buffer[0] = 0x37; buffer[1] = 0x37; buffer[2] = 0x44; buffer[3] = 0x79;
  ParseIBMFloatBE(buf, f);
  BOOST_CHECK_CLOSE_FRACTION(f, 3.14159265e-12, 1e-6);
}


