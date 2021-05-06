// $Id: bigfile_test.cpp 1606 2017-07-09 20:46:14Z perroe $

/// \file Unit tests for BigFile class in the NRLib IOTools library.

#include <nrlib/iotools/bigfile.hpp>

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

using namespace NRLib;
using namespace boost::filesystem;

BOOST_AUTO_TEST_CASE(BigFileReadWriteTest)
{
  path test_file_name = temp_directory_path() / "test_bigfile.bin";

  BigFile file_out;
  file_out.OpenWrite(test_file_name);

  char buffer[256];
  for (int i = 0; i < 256; ++i) {
    buffer[i] = static_cast<char>(i);
  }

  file_out.write(buffer, 256);
  file_out.close();

  BigFile file_in;
  file_in.OpenRead(test_file_name);

  char buffer2[256];
  size_t nread = file_in.read(buffer2, 256);
  BOOST_CHECK_EQUAL(nread, 256U);

  for (int i = 0; i < 256; ++i) {
    BOOST_CHECK_EQUAL(buffer[i], static_cast<char>(i));
  }

  char dummy;

  nread = file_in.read(&dummy, 1);
  BOOST_CHECK_EQUAL(nread, 0U);
  BOOST_CHECK(file_in.eof());
  BOOST_CHECK(!file_in.good());

  file_in.clear();
  BOOST_CHECK(!file_in.eof());
  BOOST_CHECK(file_in.good());

  file_in.seek(0, SEEK_SET);
  file_in.read(&dummy, 1);
  BOOST_CHECK_EQUAL(dummy, '\x00');
  BOOST_CHECK_EQUAL(file_in.tell(), 1);

  file_in.seek(-1, SEEK_END);
  file_in.read(&dummy, 1);
  BOOST_CHECK_EQUAL(dummy, '\xff');
  BOOST_CHECK_EQUAL(file_in.tell(), 256);

  file_in.seek(-128, SEEK_CUR);
  file_in.read(&dummy, 1);
  BOOST_CHECK_EQUAL(dummy, '\x80');
  BOOST_CHECK_EQUAL(file_in.tell(), 129);

  BOOST_CHECK_EQUAL(file_in.GetFileName(), test_file_name);
  BOOST_CHECK_EQUAL(file_in.FileSize(), 256);

  BOOST_CHECK(file_in.good());

  file_in.close();

  remove(test_file_name);
}
