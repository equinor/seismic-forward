// $Id: fileio_writetobuffer_test.cpp 1584 2017-07-03 12:28:42Z perroe $

/// \file Unit tests for the writing binary data to buffers in NRLib fileio.

#include <nrlib/iotools/fileio.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <limits>

using namespace NRLib;
using namespace NRLib::NRLibPrivate;

BOOST_AUTO_TEST_CASE(WriteUInt16BETest)
{
  unsigned short us;
  unsigned char buffer[2] = { 0xff, 0xff };
  char* buf = reinterpret_cast<char*>(&buffer);

  // 0
  us = 0U;
  WriteUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);

  us = 1U;
  WriteUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x01);

  us = 256U;
  WriteUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(buffer[0], 0x01);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);

  us = 65535U;
  WriteUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(buffer[0], 0xff);
  BOOST_CHECK_EQUAL(buffer[1], 0xff);

  us = 32767U;
  WriteUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(buffer[0], 0x7f);
  BOOST_CHECK_EQUAL(buffer[1], 0xff);

  us = 31415U;
  WriteUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(buffer[0], 0x7a);
  BOOST_CHECK_EQUAL(buffer[1], 0xb7);
}


BOOST_AUTO_TEST_CASE(WriteUInt16LETest)
{
  unsigned short us;
  unsigned char buffer[2] = { 0xff, 0xff };
  char* buf = reinterpret_cast<char*>(&buffer);

  // 0
  us = 0U;
  WriteUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);

  us = 1U;
  WriteUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(buffer[0], 0x01);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);

  us = 256U;
  WriteUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x01);

  us = 65535U;
  WriteUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(buffer[0], 0xff);
  BOOST_CHECK_EQUAL(buffer[1], 0xff);

  us = 32767U;
  WriteUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(buffer[0], 0xff);
  BOOST_CHECK_EQUAL(buffer[1], 0x7f);

  us = 31415U;
  WriteUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(buffer[0], 0xb7);
  BOOST_CHECK_EQUAL(buffer[1], 0x7a);
}


BOOST_AUTO_TEST_CASE(WriteUInt32BETest)
{
  unsigned int ui;
  unsigned char buffer[4] = { 0xff, 0xff, 0xff, 0xff };
  char* buf = reinterpret_cast<char*>(&buffer);

  // 0
  ui = 0U;
  WriteUInt32BE(buf, ui);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  ui = 1U;
  WriteUInt32BE(buf, ui);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x01);

  ui = 256U;
  WriteUInt32BE(buf, ui);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x01);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  ui = 65536U;
  WriteUInt32BE(buf, ui);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x01);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  ui = 16777216U;
  WriteUInt32BE(buf, ui);
  BOOST_CHECK_EQUAL(buffer[0], 0x01);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  ui = 4294967295U;
  WriteUInt32BE(buf, ui);
  BOOST_CHECK_EQUAL(buffer[0], 0xff);
  BOOST_CHECK_EQUAL(buffer[1], 0xff);
  BOOST_CHECK_EQUAL(buffer[2], 0xff);
  BOOST_CHECK_EQUAL(buffer[3], 0xff);

  ui = 2147483647U;
  WriteUInt32BE(buf, ui);
  BOOST_CHECK_EQUAL(buffer[0], 0x7f);
  BOOST_CHECK_EQUAL(buffer[1], 0xff);
  BOOST_CHECK_EQUAL(buffer[2], 0xff);
  BOOST_CHECK_EQUAL(buffer[3], 0xff);

  ui = 3141592653U;
  WriteUInt32BE(buf, ui);
  BOOST_CHECK_EQUAL(buffer[0], 0xbb);
  BOOST_CHECK_EQUAL(buffer[1], 0x40);
  BOOST_CHECK_EQUAL(buffer[2], 0xe6);
  BOOST_CHECK_EQUAL(buffer[3], 0x4d);
}


BOOST_AUTO_TEST_CASE(WriteUInt32LETest)
{
  unsigned int ui;
  unsigned char buffer[4] = { 0xff, 0xff, 0xff, 0xff };
  char* buf = reinterpret_cast<char*>(&buffer);

  // 0
  ui = 0U;
  WriteUInt32LE(buf, ui);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  ui = 1U;
  WriteUInt32LE(buf, ui);
  BOOST_CHECK_EQUAL(buffer[0], 0x01);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  ui = 256U;
  WriteUInt32LE(buf, ui);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x01);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  ui = 65536U;
  WriteUInt32LE(buf, ui);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x01);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  ui = 16777216U;
  WriteUInt32LE(buf, ui);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x01);

  ui = 4294967295U;
  WriteUInt32LE(buf, ui);
  BOOST_CHECK_EQUAL(buffer[0], 0xff);
  BOOST_CHECK_EQUAL(buffer[1], 0xff);
  BOOST_CHECK_EQUAL(buffer[2], 0xff);
  BOOST_CHECK_EQUAL(buffer[3], 0xff);

  ui = 2147483647U;
  WriteUInt32LE(buf, ui);
  BOOST_CHECK_EQUAL(buffer[0], 0xff);
  BOOST_CHECK_EQUAL(buffer[1], 0xff);
  BOOST_CHECK_EQUAL(buffer[2], 0xff);
  BOOST_CHECK_EQUAL(buffer[3], 0x7f);

  ui = 3141592653U;
  WriteUInt32LE(buf, ui);
  BOOST_CHECK_EQUAL(buffer[0], 0x4d);
  BOOST_CHECK_EQUAL(buffer[1], 0xe6);
  BOOST_CHECK_EQUAL(buffer[2], 0x40);
  BOOST_CHECK_EQUAL(buffer[3], 0xbb);
}


BOOST_AUTO_TEST_CASE(WriteInt16BETest)
{
  short n;
  unsigned char buffer[2] = { 0xff, 0xff };
  char* buf = reinterpret_cast<char*>(&buffer);

  // 0
  n = 0;
  WriteInt16BE(buf, n);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);

  n = -1;
  WriteInt16BE(buf, n);
  BOOST_CHECK_EQUAL(buffer[0], 0xff);
  BOOST_CHECK_EQUAL(buffer[1], 0xff);

  n = 1;
  WriteInt16BE(buf, n);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x01);

  n = 256;
  WriteInt16BE(buf, n);
  BOOST_CHECK_EQUAL(buffer[0], 0x01);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);

  n = 32767;
  WriteInt16BE(buf, n);
  BOOST_CHECK_EQUAL(buffer[0], 0x7f);
  BOOST_CHECK_EQUAL(buffer[1], 0xff);

  n = -32768;
  WriteInt16BE(buf, n);
  BOOST_CHECK_EQUAL(buffer[0], 0x80);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);

  n = 16777;
  WriteInt16BE(buf, n);
  BOOST_CHECK_EQUAL(buffer[0], 0x41);
  BOOST_CHECK_EQUAL(buffer[1], 0x89);

  n = 31415;
  WriteInt16BE(buf, n);
  BOOST_CHECK_EQUAL(buffer[0], 0x7a);
  BOOST_CHECK_EQUAL(buffer[1], 0xb7);
}


BOOST_AUTO_TEST_CASE(WriteInt32BETest)
{
  int n;
  unsigned char buffer[4] = { 0xff, 0xff, 0xff, 0xff };
  char* buf = reinterpret_cast<char*>(&buffer);

  // 0
  n = 0;
  WriteInt32BE(buf, n);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  n = -1;
  WriteInt32BE(buf, n);
  BOOST_CHECK_EQUAL(buffer[0], 0xff);
  BOOST_CHECK_EQUAL(buffer[1], 0xff);
  BOOST_CHECK_EQUAL(buffer[2], 0xff);
  BOOST_CHECK_EQUAL(buffer[3], 0xff);

  n = 1;
  WriteInt32BE(buf, n);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x01);

  n = 256;
  WriteInt32BE(buf, n);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x01);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  n = std::numeric_limits<int>::max();
  WriteInt32BE(buf, n);
  BOOST_CHECK_EQUAL(buffer[0], 0x7f);
  BOOST_CHECK_EQUAL(buffer[1], 0xff);
  BOOST_CHECK_EQUAL(buffer[2], 0xff);
  BOOST_CHECK_EQUAL(buffer[3], 0xff);

  n = std::numeric_limits<int>::min();
  WriteInt32BE(buf, n);
  BOOST_CHECK_EQUAL(buffer[0], 0x80);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  n = 16777216;
  WriteInt32BE(buf, n);
  BOOST_CHECK_EQUAL(buffer[0], 0x01);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  n = 314159265;
  WriteInt32BE(buf, n);
  BOOST_CHECK_EQUAL(buffer[0], 0x12);
  BOOST_CHECK_EQUAL(buffer[1], 0xb9);
  BOOST_CHECK_EQUAL(buffer[2], 0xb0);
  BOOST_CHECK_EQUAL(buffer[3], 0xa1);
}


BOOST_AUTO_TEST_CASE(WriteIEEEFloatBETest)
{
  float f;
  unsigned char buffer[4] = { 0xff, 0xff, 0xff, 0xff };
  char* buf = reinterpret_cast<char*>(&buffer);

  // 0
  f = 0.0F;
  WriteIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // -0
  f = -0.0F;
  WriteIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0x80);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // 1
  f = 1.0F;
  WriteIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0x3f);
  BOOST_CHECK_EQUAL(buffer[1], 0x80);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // -1
  f = -1.0F;
  WriteIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0xbf);
  BOOST_CHECK_EQUAL(buffer[1], 0x80);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // Smallest denormalized number
  f = std::numeric_limits<float>::denorm_min();
  WriteIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x01);

  // Largest normalized number
  f = std::numeric_limits<float>::max();
  WriteIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0x7f);
  BOOST_CHECK_EQUAL(buffer[1], 0x7f);
  BOOST_CHECK_EQUAL(buffer[2], 0xff);
  BOOST_CHECK_EQUAL(buffer[3], 0xff);

  // Inf
  f = std::numeric_limits<float>::infinity();
  WriteIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0x7f);
  BOOST_CHECK_EQUAL(buffer[1], 0x80);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // -Inf
  f = -std::numeric_limits<float>::infinity();
  WriteIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0xff);
  BOOST_CHECK_EQUAL(buffer[1], 0x80);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // Signalling NaN
  f = std::numeric_limits<float>::signaling_NaN();
  WriteIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0x7f);
  unsigned char flag = buffer[1] & 0x80;
  BOOST_CHECK_EQUAL(flag, 0x80);
  flag = (buffer[1] & 0x7f) | buffer[2] | buffer[3];
  BOOST_CHECK_NE(flag, 0x00);

  // -Signalling NaN
  f = -std::numeric_limits<float>::signaling_NaN();
  WriteIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0xff);
  flag = buffer[1] & 0x80;
  BOOST_CHECK_EQUAL(flag, 0x80);
  flag = (buffer[1] & 0x7f) | buffer[2] | buffer[3];
  BOOST_CHECK_NE(flag, 0x00);

  // Quiet NaN
  f = std::numeric_limits<float>::quiet_NaN();
  WriteIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0x7f);
  flag = buffer[1] & 0x80;
  BOOST_CHECK_EQUAL(flag, 0x80);
  flag = (buffer[1] & 0x7f) | buffer[2] | buffer[3];
  BOOST_CHECK_NE(flag, 0x00);

  // -Quiet NaN
  f = -std::numeric_limits<float>::quiet_NaN();
  WriteIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0xff);
  flag = buffer[1] & 0x80;
  BOOST_CHECK_EQUAL(flag, 0x80);
  flag = (buffer[1] & 0x7f) | buffer[2] | buffer[3];
  BOOST_CHECK_NE(flag, 0x00);

  // "Random" number
  // Maybe machine dependent...
  f = 3.14e-12F;
  WriteIEEEFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0x2c);
  BOOST_CHECK_EQUAL(buffer[1], 0x5c);
  BOOST_CHECK_EQUAL(buffer[2], 0xf5);
  BOOST_CHECK_EQUAL(buffer[3], 0x36);
}


BOOST_AUTO_TEST_CASE(WriteIEEEFloatLETest)
{
  float f;
  unsigned char buffer[4] = { 0xff, 0xff, 0xff, 0xff };
  char* buf = reinterpret_cast<char*>(&buffer);

  // 0
  f = 0.0F;
  WriteIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);

  // -0
  f = -0.0F;
  WriteIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(buffer[3], 0x80);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);

  // 1
  f = 1.0F;
  WriteIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(buffer[3], 0x3f);
  BOOST_CHECK_EQUAL(buffer[2], 0x80);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);

  // -1
  f = -1.0F;
  WriteIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(buffer[3], 0xbf);
  BOOST_CHECK_EQUAL(buffer[2], 0x80);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);

  // Smallest denormalized number
  f = std::numeric_limits<float>::denorm_min();
  WriteIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[0], 0x01);

  // Largest normalized number
  f = std::numeric_limits<float>::max();
  WriteIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(buffer[3], 0x7f);
  BOOST_CHECK_EQUAL(buffer[2], 0x7f);
  BOOST_CHECK_EQUAL(buffer[1], 0xff);
  BOOST_CHECK_EQUAL(buffer[0], 0xff);

  // Inf
  f = std::numeric_limits<float>::infinity();
  WriteIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(buffer[3], 0x7f);
  BOOST_CHECK_EQUAL(buffer[2], 0x80);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);

  // -Inf
  f = -std::numeric_limits<float>::infinity();
  WriteIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(buffer[3], 0xff);
  BOOST_CHECK_EQUAL(buffer[2], 0x80);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);

  // Signalling NaN
  f = std::numeric_limits<float>::signaling_NaN();
  WriteIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(buffer[3], 0x7f);
  unsigned char flag = buffer[2] & 0x80;
  BOOST_CHECK_EQUAL(flag, 0x80);
  flag = (buffer[2] & 0x7f) | buffer[1] | buffer[0];
  BOOST_CHECK_NE(flag, 0x00);

  // -Signalling NaN
  f = -std::numeric_limits<float>::signaling_NaN();
  WriteIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(buffer[3], 0xff);
  flag = buffer[2] & 0x80;
  BOOST_CHECK_EQUAL(flag, 0x80);
  flag = (buffer[2] & 0x7f) | buffer[1] | buffer[0];
  BOOST_CHECK_NE(flag, 0x00);

  // Quiet NaN
  f = std::numeric_limits<float>::quiet_NaN();
  WriteIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(buffer[3], 0x7f);
  flag = buffer[2] & 0x80;
  BOOST_CHECK_EQUAL(flag, 0x80);
  flag = (buffer[2] & 0x7f) | buffer[1] | buffer[0];
  BOOST_CHECK_NE(flag, 0x00);

  // -Quiet NaN
  f = -std::numeric_limits<float>::quiet_NaN();
  WriteIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(buffer[3], 0xff);
  flag = buffer[2] & 0x80;
  BOOST_CHECK_EQUAL(flag, 0x80);
  flag = (buffer[2] & 0x7f) | buffer[1] | buffer[0];
  BOOST_CHECK_NE(flag, 0x00);

  // "Random" number
  // Maybe machine dependent...
  f = 3.14e-12F;
  WriteIEEEFloatLE(buf, f);
  BOOST_CHECK_EQUAL(buffer[3], 0x2c);
  BOOST_CHECK_EQUAL(buffer[2], 0x5c);
  BOOST_CHECK_EQUAL(buffer[1], 0xf5);
  BOOST_CHECK_EQUAL(buffer[0], 0x36);
}


BOOST_AUTO_TEST_CASE(WriteIEEEDoubleBETest)
{
  double d;
  unsigned char buffer[8] = { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff };
  char* buf = reinterpret_cast<char*>(&buffer);

  // 0
  d = 0.0;
  WriteIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[6], 0x00);
  BOOST_CHECK_EQUAL(buffer[7], 0x00);

  // -0
  d = -0.0;
  WriteIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(buffer[0], 0x80);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[6], 0x00);
  BOOST_CHECK_EQUAL(buffer[7], 0x00);

  // 1
  d = 1.0;
  WriteIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(buffer[0], 0x3f);
  BOOST_CHECK_EQUAL(buffer[1], 0xf0);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[6], 0x00);
  BOOST_CHECK_EQUAL(buffer[7], 0x00);

  // -1
  d = -1.0;
  WriteIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(buffer[0], 0xbf);
  BOOST_CHECK_EQUAL(buffer[1], 0xf0);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[6], 0x00);
  BOOST_CHECK_EQUAL(buffer[7], 0x00);

  // Smallest denormalized number
  d = std::numeric_limits<double>::denorm_min();
  WriteIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[6], 0x00);
  BOOST_CHECK_EQUAL(buffer[7], 0x01);

  // Largest normalized number
  d = std::numeric_limits<double>::max();
  WriteIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(buffer[0], 0x7f);
  BOOST_CHECK_EQUAL(buffer[1], 0xef);
  BOOST_CHECK_EQUAL(buffer[2], 0xff);
  BOOST_CHECK_EQUAL(buffer[3], 0xff);
  BOOST_CHECK_EQUAL(buffer[4], 0xff);
  BOOST_CHECK_EQUAL(buffer[5], 0xff);
  BOOST_CHECK_EQUAL(buffer[6], 0xff);
  BOOST_CHECK_EQUAL(buffer[7], 0xff);

  // Inf
  d = std::numeric_limits<double>::infinity();
  WriteIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(buffer[0], 0x7f);
  BOOST_CHECK_EQUAL(buffer[1], 0xf0);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[6], 0x00);
  BOOST_CHECK_EQUAL(buffer[7], 0x00);

  // -Inf
  d = -std::numeric_limits<double>::infinity();
  WriteIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(buffer[0], 0xff);
  BOOST_CHECK_EQUAL(buffer[1], 0xf0);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[6], 0x00);
  BOOST_CHECK_EQUAL(buffer[7], 0x00);

  // Signalling NaN
  d = std::numeric_limits<double>::signaling_NaN();
  WriteIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(buffer[0], 0x7f);
  unsigned char flag = buffer[1] & 0xf0;
  BOOST_CHECK_EQUAL(flag, 0xf0);
  flag = (buffer[1] & 0x0f) | buffer[2] | buffer[3] | buffer[4] | buffer[5]
          | buffer[6] | buffer[7];
  BOOST_CHECK_NE(flag, 0x00);

  // -Signalling NaN
  d = -std::numeric_limits<double>::signaling_NaN();
  WriteIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(buffer[0], 0xff);
  flag = buffer[1] & 0xf0;
  BOOST_CHECK_EQUAL(flag, 0xf0);
  flag = (buffer[1] & 0x0f) | buffer[2] | buffer[3] | buffer[4] | buffer[5]
    | buffer[6] | buffer[7];
  BOOST_CHECK_NE(flag, 0x00);

  // "Random" number
  // Maybe machine dependent...
  d = 3.1415926535897932384e-12;
  WriteIEEEDoubleBE(buf, d);
  BOOST_CHECK_EQUAL(buffer[0], 0x3d);
  BOOST_CHECK_EQUAL(buffer[1], 0x8b);
  BOOST_CHECK_EQUAL(buffer[2], 0xa2);
  BOOST_CHECK_EQUAL(buffer[3], 0x3c);
  BOOST_CHECK_EQUAL(buffer[4], 0xdd);
  BOOST_CHECK_EQUAL(buffer[5], 0x51);
  BOOST_CHECK_EQUAL(buffer[6], 0x22);
  BOOST_CHECK_EQUAL(buffer[7], 0xb5);
}


BOOST_AUTO_TEST_CASE(WriteIEEEDoubleLETest)
{
  double d;
  unsigned char buffer[8] = { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff };
  char* buf = reinterpret_cast<char*>(&buffer);

  // 0
  d = 0.0;
  WriteIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[6], 0x00);
  BOOST_CHECK_EQUAL(buffer[7], 0x00);

  // -0
  d = -0.0;
  WriteIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(buffer[7], 0x80);
  BOOST_CHECK_EQUAL(buffer[6], 0x00);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);

  // 1
  d = 1.0;
  WriteIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(buffer[7], 0x3f);
  BOOST_CHECK_EQUAL(buffer[6], 0xf0);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);

  // -1
  d = -1.0;
  WriteIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(buffer[7], 0xbf);
  BOOST_CHECK_EQUAL(buffer[6], 0xf0);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);

  // Smallest denormalized number
  d = std::numeric_limits<double>::denorm_min();
  WriteIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(buffer[7], 0x00);
  BOOST_CHECK_EQUAL(buffer[6], 0x00);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[0], 0x01);

  // Largest normalized number
  d = std::numeric_limits<double>::max();
  WriteIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(buffer[7], 0x7f);
  BOOST_CHECK_EQUAL(buffer[6], 0xef);
  BOOST_CHECK_EQUAL(buffer[5], 0xff);
  BOOST_CHECK_EQUAL(buffer[4], 0xff);
  BOOST_CHECK_EQUAL(buffer[3], 0xff);
  BOOST_CHECK_EQUAL(buffer[2], 0xff);
  BOOST_CHECK_EQUAL(buffer[1], 0xff);
  BOOST_CHECK_EQUAL(buffer[0], 0xff);

  // Inf
  d = std::numeric_limits<double>::infinity();
  WriteIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(buffer[7], 0x7f);
  BOOST_CHECK_EQUAL(buffer[6], 0xf0);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);

  // -Inf
  d = -std::numeric_limits<double>::infinity();
  WriteIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(buffer[7], 0xff);
  BOOST_CHECK_EQUAL(buffer[6], 0xf0);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);

  // Signalling NaN
  d = std::numeric_limits<double>::signaling_NaN();
  WriteIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(buffer[7], 0x7f);
  unsigned char flag = buffer[6] & 0xf0;
  BOOST_CHECK_EQUAL(flag, 0xf0);
  flag = (buffer[6] & 0x0f) | buffer[5] | buffer[4] | buffer[3] | buffer[2]
    | buffer[1] | buffer[0];
  BOOST_CHECK_NE(flag, 0x00);

  // -Signalling NaN
  d = -std::numeric_limits<double>::signaling_NaN();
  WriteIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(buffer[7], 0xff);
  flag = buffer[6] & 0xf0;
  BOOST_CHECK_EQUAL(flag, 0xf0);
  flag = (buffer[6] & 0x0f) | buffer[5] | buffer[4] | buffer[3] | buffer[2]
    | buffer[1] | buffer[0];
  BOOST_CHECK_NE(flag, 0x00);

  // "Random" number
  // Maybe machine dependent...
  d = 3.1415926535897932384e-12;
  WriteIEEEDoubleLE(buf, d);
  BOOST_CHECK_EQUAL(buffer[7], 0x3d);
  BOOST_CHECK_EQUAL(buffer[6], 0x8b);
  BOOST_CHECK_EQUAL(buffer[5], 0xa2);
  BOOST_CHECK_EQUAL(buffer[4], 0x3c);
  BOOST_CHECK_EQUAL(buffer[3], 0xdd);
  BOOST_CHECK_EQUAL(buffer[2], 0x51);
  BOOST_CHECK_EQUAL(buffer[1], 0x22);
  BOOST_CHECK_EQUAL(buffer[0], 0xb5);
}


BOOST_AUTO_TEST_CASE(WriteIBMFloatBETest)
{
  float f;
  unsigned char buffer[4] = { 0xff, 0xff, 0xff, 0xff };
  char* buf = reinterpret_cast<char*>(&buffer);

  // 0
  f = 0.0F;
  WriteIBMFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // -0    (Not supported by NRLib)
  //f = -0.0F;
  //WriteIBMFloatBE(buf, f);
  // BOOST_CHECK_EQUAL(buffer[0], 0x80);
  //BOOST_CHECK_EQUAL(buffer[1], 0x00);
  //BOOST_CHECK_EQUAL(buffer[2], 0x00);
  //BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // 1
  f = 1.0F;
  WriteIBMFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0x41);
  BOOST_CHECK_EQUAL(buffer[1], 0x10);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // -1
  f = -1.0F;
  WriteIBMFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0xc1);
  BOOST_CHECK_EQUAL(buffer[1], 0x10);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // Inf
  f = std::numeric_limits<float>::infinity();
  WriteIBMFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0x61);
  BOOST_CHECK_EQUAL(buffer[1], 0x10);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // -Inf
  f = -std::numeric_limits<float>::infinity();
  WriteIBMFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0xe1);
  BOOST_CHECK_EQUAL(buffer[1], 0x10);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // nan
  f = std::numeric_limits<float>::signaling_NaN();
  WriteIBMFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0x61);
  unsigned char flag = buffer[1] & 0x10;
  BOOST_CHECK_EQUAL(flag, 0x10);
  // flag = (buffer[1] & 0xef) | buffer[2] | buffer[3];
  // BOOST_TEST(flag != 0x00);

  // Wikipedia example number.
  f = -118.625F;
  WriteIBMFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0xc2);
  BOOST_CHECK_EQUAL(buffer[1], 0x76);
  BOOST_CHECK_EQUAL(buffer[2], 0xa0);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // "Random" number
  f = 3.14159265e-12F;
  WriteIBMFloatBE(buf, f);
  BOOST_CHECK_EQUAL(buffer[0], 0x37);
  BOOST_CHECK_EQUAL(buffer[1], 0x37);
  BOOST_CHECK_EQUAL(buffer[2], 0x44);
  BOOST_CHECK_EQUAL(buffer[3], 0x79);
}

