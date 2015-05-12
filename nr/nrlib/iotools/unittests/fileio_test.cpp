// $Id: fileio_test.cpp 1255 2014-02-28 15:07:01Z vigsnes $

/// \file Unit tests for the FileIO functions in the NRLib IOTools library.

#include <nrlib/iotools/fileio.hpp>
#include <boost/test/unit_test.hpp>

using namespace NRLib;
using namespace NRLib::NRLibPrivate;

BOOST_AUTO_TEST_CASE( ParseUInt16BETest )
{
  unsigned short us;

  // 0
  unsigned char buffer[2] = {0x00, 0x00};
  char* buf = reinterpret_cast<char*>(&buffer);
  ParseUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(us, 0);

  // 1
  buffer[0] = 0x00; buffer[1] = 0x01;
  ParseUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(us, 1);

  // 256
  buffer[0] = 0x01; buffer[1] = 0x00;
  ParseUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(us, 256);

  // 65535
  buffer[0] = 0xff; buffer[1] = 0xff;
  ParseUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(us, 65535);

  // 32767
  buffer[0] = 0x7f; buffer[1] = 0xff;
  ParseUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(us, 32767);

  // 31415
  buffer[0] = 0x7a; buffer[1] = 0xb7;
  ParseUInt16BE(buf, us);
  BOOST_CHECK_EQUAL(us, 31415);
}


BOOST_AUTO_TEST_CASE( ParseUInt16LETest )
{
  unsigned short us;

  // 0
  unsigned char buffer[2] = {0x00, 0x00};
  char* buf = reinterpret_cast<char*>(&buffer);
  ParseUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(us, 0);

  // 1
  buffer[0] = 0x01; buffer[1] = 0x00;
  ParseUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(us, 1);

  // 256
  buffer[0] = 0x00; buffer[1] = 0x01;
  ParseUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(us, 256);

  // 65535
  buffer[0] = 0xff; buffer[1] = 0xff;
  ParseUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(us, 65535);

  // 32767
  buffer[0] = 0xff; buffer[1] = 0x7f;
  ParseUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(us, 32767);

  // 31415
  buffer[0] = 0xb7; buffer[1] = 0x7a;
  ParseUInt16LE(buf, us);
  BOOST_CHECK_EQUAL(us, 31415);
}
