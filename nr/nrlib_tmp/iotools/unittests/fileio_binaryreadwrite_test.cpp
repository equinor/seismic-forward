// $Id: fileio_binaryreadwrite_test.cpp 1537 2017-06-23 09:36:41Z perroe $

/// \file Unit tests for reading and writing binary files using
///       the FileIO functions in the NRLib IOTools library.

#include <nrlib/iotools/fileio.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/test/unit_test.hpp>

#include <fstream>

using namespace NRLib;
using namespace boost::filesystem;

BOOST_AUTO_TEST_CASE(ReadWriteBinaryShortTest)
{
  path filepath = temp_directory_path() / "test_shorts.bin";

  std::vector<short> ints(8);
  ints[0] = 0;
  ints[1] = -1;
  ints[2] = 1;
  ints[3] = 256;
  ints[4] = 32767;
  ints[5] = -32768;
  ints[6] = 16777;
  ints[7] = 31415;

  std::ofstream file(filepath.c_str(), std::ios::binary | std::ios::out);

  WriteBinaryShort(file, ints[0]);
  WriteBinaryShort(file, ints[1]);
  WriteBinaryShortArray(file, ints.begin() + 2, ints.end() - 1);
  WriteBinaryShort(file, ints[7]);

  file.close();

  std::ifstream ifile(filepath.c_str(), std::ios::binary | std::ios::in);
  unsigned char buffer[16];
  ifile.read(reinterpret_cast<char*>(buffer), 16);
  ifile.close();

  // 0
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);

  // -1
  BOOST_CHECK_EQUAL(buffer[2], 0xff);
  BOOST_CHECK_EQUAL(buffer[3], 0xff);

  // 1
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[5], 0x01);

  // 256
  BOOST_CHECK_EQUAL(buffer[6], 0x01);
  BOOST_CHECK_EQUAL(buffer[7], 0x00);

  // 32767
  BOOST_CHECK_EQUAL(buffer[8], 0x7f);
  BOOST_CHECK_EQUAL(buffer[9], 0xff);

  // -32768
  BOOST_CHECK_EQUAL(buffer[10], 0x80);
  BOOST_CHECK_EQUAL(buffer[11], 0x00);

  // 16777
  BOOST_CHECK_EQUAL(buffer[12], 0x41);
  BOOST_CHECK_EQUAL(buffer[13], 0x89);

  // 31415
  BOOST_CHECK_EQUAL(buffer[14], 0x7a);
  BOOST_CHECK_EQUAL(buffer[15], 0xb7);

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  for (int i = 0; i < 8; ++i) {
    BOOST_CHECK_EQUAL(ReadBinaryShort(ifile), ints[i]);
  }

  ifile.close();

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), false);

  std::vector<short> in_i(8);
  ReadBinaryShortArray(ifile, in_i.begin(), 8);

  for (int i = 0; i < 8; ++i)
    BOOST_CHECK_EQUAL(in_i[i], ints[i]);

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), true);

  ifile.close();

  remove(filepath);
}


BOOST_AUTO_TEST_CASE(ReadWriteBinaryShortLETest)
{
  path filepath = temp_directory_path() / "test_shorts_le.bin";

  std::vector<short> ints(8);
  ints[0] = 0;
  ints[1] = -1;
  ints[2] = 1;
  ints[3] = 256;
  ints[4] = 32767;
  ints[5] = -32768;
  ints[6] = 16777;
  ints[7] = 31415;

  std::ofstream file(filepath.c_str(), std::ios::binary | std::ios::out);

  WriteBinaryShort(file, ints[0], END_LITTLE_ENDIAN);
  WriteBinaryShort(file, ints[1], END_LITTLE_ENDIAN);
  WriteBinaryShortArray(file, ints.begin() + 2, ints.end() - 1, END_LITTLE_ENDIAN);
  WriteBinaryShort(file, ints[7], END_LITTLE_ENDIAN);

  file.close();

  std::ifstream ifile(filepath.c_str(), std::ios::binary | std::ios::in);
  unsigned char buffer[16];
  ifile.read(reinterpret_cast<char*>(buffer), 16);
  ifile.close();

  // 0
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);

  // -1
  BOOST_CHECK_EQUAL(buffer[2], 0xff);
  BOOST_CHECK_EQUAL(buffer[3], 0xff);

  // 1
  BOOST_CHECK_EQUAL(buffer[4], 0x01);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);

  // 256
  BOOST_CHECK_EQUAL(buffer[6], 0x00);
  BOOST_CHECK_EQUAL(buffer[7], 0x01);

  // 32767
  BOOST_CHECK_EQUAL(buffer[8], 0xff);
  BOOST_CHECK_EQUAL(buffer[9], 0x7f);

  // -32768
  BOOST_CHECK_EQUAL(buffer[10], 0x00);
  BOOST_CHECK_EQUAL(buffer[11], 0x80);

  // 16777
  BOOST_CHECK_EQUAL(buffer[12], 0x89);
  BOOST_CHECK_EQUAL(buffer[13], 0x41);

  // 31415
  BOOST_CHECK_EQUAL(buffer[14], 0xb7);
  BOOST_CHECK_EQUAL(buffer[15], 0x7a);


  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  for (int i = 0; i < 8; ++i) {
    BOOST_CHECK_EQUAL(ReadBinaryShort(ifile, END_LITTLE_ENDIAN), ints[i]);
  }

  ifile.close();

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), false);

  std::vector<short> in_i(8);
  ReadBinaryShortArray(ifile, in_i.begin(), 8, END_LITTLE_ENDIAN);

  for (int i = 0; i < 8; ++i)
    BOOST_CHECK_EQUAL(in_i[i], ints[i]);

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), true);

  ifile.close();

  remove(filepath);
}


BOOST_AUTO_TEST_CASE(ReadWriteBinaryIntTest)
{
  path filepath = temp_directory_path() / "test_ints.bin";

  std::vector<int> ints(8);
  ints[0] = 0;
  ints[1] = -1;
  ints[2] = 1;
  ints[3] = 256;
  ints[4] = std::numeric_limits<int>::max();
  ints[5] = std::numeric_limits<int>::min();
  ints[6] = 16777216;
  ints[7] = 314159265;

  std::ofstream file(filepath.c_str(), std::ios::binary | std::ios::out);

  WriteBinaryInt(file, ints[0]);
  WriteBinaryInt(file, ints[1]);
  WriteBinaryIntArray(file, ints.begin() + 2, ints.end() - 1);
  WriteBinaryInt(file, ints[7]);

  file.close();

  std::ifstream ifile(filepath.c_str(), std::ios::binary | std::ios::in);
  unsigned char buffer[32];
  ifile.read(reinterpret_cast<char*>(buffer), 32);
  ifile.close();

  // 0
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // -1
  BOOST_CHECK_EQUAL(buffer[4], 0xff);
  BOOST_CHECK_EQUAL(buffer[5], 0xff);
  BOOST_CHECK_EQUAL(buffer[6], 0xff);
  BOOST_CHECK_EQUAL(buffer[7], 0xff);

  // 1
  BOOST_CHECK_EQUAL(buffer[8], 0x00);
  BOOST_CHECK_EQUAL(buffer[9], 0x00);
  BOOST_CHECK_EQUAL(buffer[10], 0x00);
  BOOST_CHECK_EQUAL(buffer[11], 0x01);

  // 256
  BOOST_CHECK_EQUAL(buffer[12], 0x00);
  BOOST_CHECK_EQUAL(buffer[13], 0x00);
  BOOST_CHECK_EQUAL(buffer[14], 0x01);
  BOOST_CHECK_EQUAL(buffer[15], 0x00);

  // 2147483647
  BOOST_CHECK_EQUAL(buffer[16], 0x7f);
  BOOST_CHECK_EQUAL(buffer[17], 0xff);
  BOOST_CHECK_EQUAL(buffer[18], 0xff);
  BOOST_CHECK_EQUAL(buffer[19], 0xff);

  // -2147483648
  BOOST_CHECK_EQUAL(buffer[20], 0x80);
  BOOST_CHECK_EQUAL(buffer[21], 0x00);
  BOOST_CHECK_EQUAL(buffer[22], 0x00);
  BOOST_CHECK_EQUAL(buffer[23], 0x00);

  // 16777216
  BOOST_CHECK_EQUAL(buffer[24], 0x01);
  BOOST_CHECK_EQUAL(buffer[25], 0x00);
  BOOST_CHECK_EQUAL(buffer[26], 0x00);
  BOOST_CHECK_EQUAL(buffer[27], 0x00);

  // 314159265
  BOOST_CHECK_EQUAL(buffer[28], 0x12);
  BOOST_CHECK_EQUAL(buffer[29], 0xb9);
  BOOST_CHECK_EQUAL(buffer[30], 0xb0);
  BOOST_CHECK_EQUAL(buffer[31], 0xa1);


  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  for (int i = 0; i < 8; ++i) {
    BOOST_CHECK_EQUAL(ReadBinaryInt(ifile), ints[i]);
  }

  ifile.close();

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), false);

  std::vector<int> in_i(8);
  ReadBinaryIntArray(ifile, in_i.begin(), 8);

  for (int i = 0; i < 8; ++i)
    BOOST_CHECK_EQUAL(in_i[i], ints[i]);

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), true);

  ifile.close();

  remove(filepath);
}


BOOST_AUTO_TEST_CASE(ReadWriteBinaryIntLETest)
{
  path filepath = temp_directory_path() / "test_ints_le.bin";

  std::vector<int> ints(8);
  ints[0] = 0;
  ints[1] = -1;
  ints[2] = 1;
  ints[3] = 256;
  ints[4] = std::numeric_limits<int>::max();
  ints[5] = std::numeric_limits<int>::min();
  ints[6] = 16777216;
  ints[7] = 314159265;

  std::ofstream file(filepath.c_str(), std::ios::binary | std::ios::out);

  WriteBinaryInt(file, ints[0], END_LITTLE_ENDIAN);
  WriteBinaryInt(file, ints[1], END_LITTLE_ENDIAN);
  WriteBinaryIntArray(file, ints.begin() + 2, ints.end() - 1, END_LITTLE_ENDIAN);
  WriteBinaryInt(file, ints[7], END_LITTLE_ENDIAN);

  file.close();

  std::ifstream ifile(filepath.c_str(), std::ios::binary | std::ios::in);
  unsigned char buffer[32];
  ifile.read(reinterpret_cast<char*>(buffer), 32);
  ifile.close();

  // 0
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // -1
  BOOST_CHECK_EQUAL(buffer[4], 0xff);
  BOOST_CHECK_EQUAL(buffer[5], 0xff);
  BOOST_CHECK_EQUAL(buffer[6], 0xff);
  BOOST_CHECK_EQUAL(buffer[7], 0xff);

  // 1
  BOOST_CHECK_EQUAL(buffer[8], 0x01);
  BOOST_CHECK_EQUAL(buffer[9], 0x00);
  BOOST_CHECK_EQUAL(buffer[10], 0x00);
  BOOST_CHECK_EQUAL(buffer[11], 0x00);

  // 256
  BOOST_CHECK_EQUAL(buffer[12], 0x00);
  BOOST_CHECK_EQUAL(buffer[13], 0x01);
  BOOST_CHECK_EQUAL(buffer[14], 0x00);
  BOOST_CHECK_EQUAL(buffer[15], 0x00);

  // 2147483647
  BOOST_CHECK_EQUAL(buffer[16], 0xff);
  BOOST_CHECK_EQUAL(buffer[17], 0xff);
  BOOST_CHECK_EQUAL(buffer[18], 0xff);
  BOOST_CHECK_EQUAL(buffer[19], 0x7f);

  // -2147483648
  BOOST_CHECK_EQUAL(buffer[20], 0x00);
  BOOST_CHECK_EQUAL(buffer[21], 0x00);
  BOOST_CHECK_EQUAL(buffer[22], 0x00);
  BOOST_CHECK_EQUAL(buffer[23], 0x80);

  // 16777216
  BOOST_CHECK_EQUAL(buffer[24], 0x00);
  BOOST_CHECK_EQUAL(buffer[25], 0x00);
  BOOST_CHECK_EQUAL(buffer[26], 0x00);
  BOOST_CHECK_EQUAL(buffer[27], 0x01);

  // 314159265
  BOOST_CHECK_EQUAL(buffer[28], 0xa1);
  BOOST_CHECK_EQUAL(buffer[29], 0xb0);
  BOOST_CHECK_EQUAL(buffer[30], 0xb9);
  BOOST_CHECK_EQUAL(buffer[31], 0x12);


  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  for (int i = 0; i < 8; ++i) {
    BOOST_CHECK_EQUAL(ReadBinaryInt(ifile, END_LITTLE_ENDIAN), ints[i]);
  }

  ifile.close();

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), false);

  std::vector<int> in_i(8);
  ReadBinaryIntArray(ifile, in_i.begin(), 8, END_LITTLE_ENDIAN);

  for (int i = 0; i < 8; ++i)
    BOOST_CHECK_EQUAL(in_i[i], ints[i]);

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), true);

  ifile.close();

  remove(filepath);
}


BOOST_AUTO_TEST_CASE(ReadWriteBinaryFloatTest)
{
  path filepath = temp_directory_path() / "test_floats.bin";

  std::vector<float> f(8);
  f[0] = 0.0F;
  f[1] = -1.0F;
  f[2] = 1.0F;
  f[3] = std::numeric_limits<float>::infinity();
  f[4] = -std::numeric_limits<float>::infinity();
  f[5] = 1.5622F;
  f[6] = std::numeric_limits<float>::signaling_NaN();
  f[7] = 3.141592653589793e-12F;

  std::ofstream file(filepath.c_str(), std::ios::binary | std::ios::out);

  WriteBinaryFloat(file, f[0]);
  WriteBinaryFloat(file, f[1]);
  WriteBinaryFloatArray(file, f.begin() + 2, f.end() - 1);
  WriteBinaryFloat(file, f[7]);

  file.close();

  std::ifstream ifile(filepath.c_str(), std::ios::binary | std::ios::in);
  unsigned char buffer[32];
  ifile.read(reinterpret_cast<char*>(buffer), 32);
  ifile.close();

  // 0.0F
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // -1.0F
  BOOST_CHECK_EQUAL(buffer[4], 0xbf);
  BOOST_CHECK_EQUAL(buffer[5], 0x80);
  BOOST_CHECK_EQUAL(buffer[6], 0x00);
  BOOST_CHECK_EQUAL(buffer[7], 0x00);

  // 1.0F
  BOOST_CHECK_EQUAL(buffer[8], 0x3f);
  BOOST_CHECK_EQUAL(buffer[9], 0x80);
  BOOST_CHECK_EQUAL(buffer[10], 0x00);
  BOOST_CHECK_EQUAL(buffer[11], 0x00);

  // +Infinity
  BOOST_CHECK_EQUAL(buffer[12], 0x7f);
  BOOST_CHECK_EQUAL(buffer[13], 0x80);
  BOOST_CHECK_EQUAL(buffer[14], 0x00);
  BOOST_CHECK_EQUAL(buffer[15], 0x00);

  // -Infinity
  BOOST_CHECK_EQUAL(buffer[16], 0xff);
  BOOST_CHECK_EQUAL(buffer[17], 0x80);
  BOOST_CHECK_EQUAL(buffer[18], 0x00);
  BOOST_CHECK_EQUAL(buffer[19], 0x00);

  // 1.5622F
  BOOST_CHECK_EQUAL(buffer[20], 0x3f);
  BOOST_CHECK_EQUAL(buffer[21], 0xc7);
  BOOST_CHECK_EQUAL(buffer[22], 0xf6);
  BOOST_CHECK_EQUAL(buffer[23], 0x2b);

  // NaN
  BOOST_CHECK_EQUAL(buffer[24], 0x7f);
  unsigned char flag = buffer[25] & 0x80;
  BOOST_CHECK_EQUAL(flag, 0x80);
  flag = (buffer[25] & 0x7f) | buffer[26] | buffer[27];
  BOOST_CHECK_NE(flag, 0x00);

  // 3.141592653589793e-12F
  BOOST_CHECK_EQUAL(buffer[28], 0x2c);
  BOOST_CHECK_EQUAL(buffer[29], 0x5d);
  BOOST_CHECK_EQUAL(buffer[30], 0x11);
  BOOST_CHECK_EQUAL(buffer[31], 0xe7);

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  for (int i = 0; i < 8; ++i) {
    if (i == 6)
      BOOST_CHECK(boost::math::isnan(ReadBinaryFloat(ifile)));
    else
      BOOST_CHECK_EQUAL(ReadBinaryFloat(ifile), f[i]);
  }
  ifile.close();

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), false);

  std::vector<float> in_f(8);
  ReadBinaryFloatArray(ifile, in_f.begin(), 8);

  BOOST_CHECK(boost::math::isnan(in_f[6]));

  for (int i = 0; i < 8; ++i) {
    if (i != 6)
      BOOST_CHECK_EQUAL(in_f[i], f[i]);
  }

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), true);

  ifile.close();

  remove(filepath);
}


BOOST_AUTO_TEST_CASE(ReadWriteBinaryFloatLETest)
{
  path filepath = temp_directory_path() / "test_floats_le.bin";

  std::vector<float> f(8);
  f[0] = 0.0F;
  f[1] = -1.0F;
  f[2] = 1.0F;
  f[3] = std::numeric_limits<float>::infinity();
  f[4] = -std::numeric_limits<float>::infinity();
  f[5] = 1.5622F;
  f[6] = std::numeric_limits<float>::signaling_NaN();
  f[7] = 3.141592653589793e-12F;

  std::ofstream file(filepath.c_str(), std::ios::binary | std::ios::out);

  WriteBinaryFloat(file, f[0], END_LITTLE_ENDIAN);
  WriteBinaryFloat(file, f[1], END_LITTLE_ENDIAN);
  WriteBinaryFloatArray(file, f.begin() + 2, f.end() - 1, END_LITTLE_ENDIAN);
  WriteBinaryFloat(file, f[7], END_LITTLE_ENDIAN);

  file.close();

  std::ifstream ifile(filepath.c_str(), std::ios::binary | std::ios::in);
  unsigned char buffer[32];
  ifile.read(reinterpret_cast<char*>(buffer), 32);
  ifile.close();

  // 0.0F
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // -1.0F
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[6], 0x80);
  BOOST_CHECK_EQUAL(buffer[7], 0xbf);

  // 1.0F
  BOOST_CHECK_EQUAL(buffer[8], 0x00);
  BOOST_CHECK_EQUAL(buffer[9], 0x00);
  BOOST_CHECK_EQUAL(buffer[10], 0x80);
  BOOST_CHECK_EQUAL(buffer[11], 0x3f);

  // +Infinity
  BOOST_CHECK_EQUAL(buffer[12], 0x00);
  BOOST_CHECK_EQUAL(buffer[13], 0x00);
  BOOST_CHECK_EQUAL(buffer[14], 0x80);
  BOOST_CHECK_EQUAL(buffer[15], 0x7f);

  // -Infinity
  BOOST_CHECK_EQUAL(buffer[16], 0x00);
  BOOST_CHECK_EQUAL(buffer[17], 0x00);
  BOOST_CHECK_EQUAL(buffer[18], 0x80);
  BOOST_CHECK_EQUAL(buffer[19], 0xff);

  // 1.5622F
  BOOST_CHECK_EQUAL(buffer[20], 0x2b);
  BOOST_CHECK_EQUAL(buffer[21], 0xf6);
  BOOST_CHECK_EQUAL(buffer[22], 0xc7);
  BOOST_CHECK_EQUAL(buffer[23], 0x3f);

  // NaN
  BOOST_CHECK_EQUAL(buffer[27], 0x7f);
  unsigned char flag = buffer[26] & 0x80;
  BOOST_CHECK_EQUAL(flag, 0x80);
  flag = (buffer[26] & 0x7f) | buffer[25] | buffer[24];
  BOOST_CHECK_NE(flag, 0x00);

  // 3.141592653589793e-12F
  BOOST_CHECK_EQUAL(buffer[28], 0xe7);
  BOOST_CHECK_EQUAL(buffer[29], 0x11);
  BOOST_CHECK_EQUAL(buffer[30], 0x5d);
  BOOST_CHECK_EQUAL(buffer[31], 0x2c);

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  for (int i = 0; i < 8; ++i) {
    if (i == 6)
      BOOST_CHECK(boost::math::isnan(ReadBinaryFloat(ifile, END_LITTLE_ENDIAN)));
    else
      BOOST_CHECK_EQUAL(ReadBinaryFloat(ifile, END_LITTLE_ENDIAN), f[i]);
  }
  ifile.close();

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), false);

  std::vector<float> in_f(8);
  ReadBinaryFloatArray(ifile, in_f.begin(), 8, END_LITTLE_ENDIAN);

  BOOST_CHECK(boost::math::isnan(in_f[6]));

  for (int i = 0; i < 8; ++i) {
    if (i != 6)
      BOOST_CHECK_EQUAL(in_f[i], f[i]);
  }

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), true);

  ifile.close();

  remove(filepath);
}


BOOST_AUTO_TEST_CASE(ReadWriteBinaryDoubleTest)
{
  path filepath = temp_directory_path() / "test_doubles.bin";

  std::vector<double> f(8);
  f[0] = 0.0;
  f[1] = -1.0;
  f[2] = 1.0;
  f[3] = std::numeric_limits<double>::infinity();
  f[4] = -std::numeric_limits<double>::infinity();
  f[5] = 1.5622;
  f[6] = std::numeric_limits<double>::signaling_NaN();
  f[7] = 3.141592653589793e-12;

  std::ofstream file(filepath.c_str(), std::ios::binary | std::ios::out);

  WriteBinaryDouble(file, f[0]);
  WriteBinaryDouble(file, f[1]);
  WriteBinaryDoubleArray(file, f.begin() + 2, f.end() - 1);
  WriteBinaryDouble(file, f[7]);

  file.close();

  std::ifstream ifile(filepath.c_str(), std::ios::binary | std::ios::in);
  unsigned char buffer[64];
  ifile.read(reinterpret_cast<char*>(buffer), 64);
  ifile.close();

  // 0
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[6], 0x00);
  BOOST_CHECK_EQUAL(buffer[7], 0x00);

  // -1.0
  BOOST_CHECK_EQUAL(buffer[8], 0xbf);
  BOOST_CHECK_EQUAL(buffer[9], 0xf0);
  BOOST_CHECK_EQUAL(buffer[10], 0x00);
  BOOST_CHECK_EQUAL(buffer[11], 0x00);
  BOOST_CHECK_EQUAL(buffer[12], 0x00);
  BOOST_CHECK_EQUAL(buffer[13], 0x00);
  BOOST_CHECK_EQUAL(buffer[14], 0x00);
  BOOST_CHECK_EQUAL(buffer[15], 0x00);

  // 1.0
  BOOST_CHECK_EQUAL(buffer[16], 0x3f);
  BOOST_CHECK_EQUAL(buffer[17], 0xf0);
  BOOST_CHECK_EQUAL(buffer[18], 0x00);
  BOOST_CHECK_EQUAL(buffer[19], 0x00);
  BOOST_CHECK_EQUAL(buffer[20], 0x00);
  BOOST_CHECK_EQUAL(buffer[21], 0x00);
  BOOST_CHECK_EQUAL(buffer[22], 0x00);
  BOOST_CHECK_EQUAL(buffer[23], 0x00);

  // 3.141592653589793e-12
  BOOST_CHECK_EQUAL(buffer[56], 0x3d);
  BOOST_CHECK_EQUAL(buffer[57], 0x8b);
  BOOST_CHECK_EQUAL(buffer[58], 0xa2);
  BOOST_CHECK_EQUAL(buffer[59], 0x3c);
  BOOST_CHECK_EQUAL(buffer[60], 0xdd);
  BOOST_CHECK_EQUAL(buffer[61], 0x51);
  BOOST_CHECK_EQUAL(buffer[62], 0x22);
  BOOST_CHECK_EQUAL(buffer[63], 0xb5);

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  for (int i = 0; i < 8; ++i) {
    if (i == 6)
      BOOST_CHECK(boost::math::isnan(ReadBinaryDouble(ifile)));
    else
      BOOST_CHECK_EQUAL(ReadBinaryDouble(ifile), f[i]);
  }
  ifile.close();

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), false);

  std::vector<double> in_f(8);
  ReadBinaryDoubleArray(ifile, in_f.begin(), 8);

  BOOST_CHECK(boost::math::isnan(in_f[6]));

  for (int i = 0; i < 8; ++i) {
    if (i != 6)
      BOOST_CHECK_EQUAL(in_f[i], f[i]);
  }

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), true);

  ifile.close();

  remove(filepath);
}


BOOST_AUTO_TEST_CASE(ReadWriteBinaryDoubleLETest)
{
  path filepath = temp_directory_path() / "test_doubles.bin";

  std::vector<double> f(8);
  f[0] = 0.0;
  f[1] = -1.0;
  f[2] = 1.0;
  f[3] = std::numeric_limits<double>::infinity();
  f[4] = -std::numeric_limits<double>::infinity();
  f[5] = 1.5622;
  f[6] = std::numeric_limits<double>::signaling_NaN();
  f[7] = 3.141592653589793e-12;

  std::ofstream file(filepath.c_str(), std::ios::binary | std::ios::out);

  WriteBinaryDouble(file, f[0], END_LITTLE_ENDIAN);
  WriteBinaryDouble(file, f[1], END_LITTLE_ENDIAN);
  WriteBinaryDoubleArray(file, f.begin() + 2, f.end() - 1, END_LITTLE_ENDIAN);
  WriteBinaryDouble(file, f[7], END_LITTLE_ENDIAN);

  file.close();

  std::ifstream ifile(filepath.c_str(), std::ios::binary | std::ios::in);
  unsigned char buffer[64];
  ifile.read(reinterpret_cast<char*>(buffer), 64);
  ifile.close();

  // 0
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[6], 0x00);
  BOOST_CHECK_EQUAL(buffer[7], 0x00);

  // -1.0
  BOOST_CHECK_EQUAL(buffer[8], 0x00);
  BOOST_CHECK_EQUAL(buffer[9], 0x00);
  BOOST_CHECK_EQUAL(buffer[10], 0x00);
  BOOST_CHECK_EQUAL(buffer[11], 0x00);
  BOOST_CHECK_EQUAL(buffer[12], 0x00);
  BOOST_CHECK_EQUAL(buffer[13], 0x00);
  BOOST_CHECK_EQUAL(buffer[14], 0xf0);
  BOOST_CHECK_EQUAL(buffer[15], 0xbf);

  // 1.0
  BOOST_CHECK_EQUAL(buffer[16], 0x00);
  BOOST_CHECK_EQUAL(buffer[17], 0x00);
  BOOST_CHECK_EQUAL(buffer[18], 0x00);
  BOOST_CHECK_EQUAL(buffer[19], 0x00);
  BOOST_CHECK_EQUAL(buffer[20], 0x00);
  BOOST_CHECK_EQUAL(buffer[21], 0x00);
  BOOST_CHECK_EQUAL(buffer[22], 0xf0);
  BOOST_CHECK_EQUAL(buffer[23], 0x3f);

  // 3.141592653589793e-12
  BOOST_CHECK_EQUAL(buffer[56], 0xb5);
  BOOST_CHECK_EQUAL(buffer[57], 0x22);
  BOOST_CHECK_EQUAL(buffer[58], 0x51);
  BOOST_CHECK_EQUAL(buffer[59], 0xdd);
  BOOST_CHECK_EQUAL(buffer[60], 0x3c);
  BOOST_CHECK_EQUAL(buffer[61], 0xa2);
  BOOST_CHECK_EQUAL(buffer[62], 0x8b);
  BOOST_CHECK_EQUAL(buffer[63], 0x3d);

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  for (int i = 0; i < 8; ++i) {
    if (i == 6)
      BOOST_CHECK(boost::math::isnan(ReadBinaryDouble(ifile, END_LITTLE_ENDIAN)));
    else
      BOOST_CHECK_EQUAL(ReadBinaryDouble(ifile, END_LITTLE_ENDIAN), f[i]);
  }
  ifile.close();

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), false);

  std::vector<double> in_f(8);
  ReadBinaryDoubleArray(ifile, in_f.begin(), 8, END_LITTLE_ENDIAN);

  BOOST_CHECK(boost::math::isnan(in_f[6]));

  for (int i = 0; i < 8; ++i) {
    if (i != 6)
      BOOST_CHECK_EQUAL(in_f[i], f[i]);
  }

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), true);

  ifile.close();

  remove(filepath);
}


BOOST_AUTO_TEST_CASE(ReadWriteBinaryIbmFloatTest)
{
  path filepath = temp_directory_path() / "test_IBM_float.bin";

  std::vector<float> f(8);
  f[0] = 0.0F;
  f[1] = -1.0F;
  f[2] = 1.0F;
  f[3] = std::numeric_limits<float>::infinity();
  f[4] = -std::numeric_limits<float>::infinity();
  f[5] = 1.5622F;
  f[6] = std::numeric_limits<float>::signaling_NaN();
  f[7] = 3.141592653589793e-12F;

  std::ofstream file(filepath.c_str(), std::ios::binary | std::ios::out);

  WriteBinaryIbmFloat(file, f[0]);
  WriteBinaryIbmFloat(file, f[1]);
  WriteBinaryIbmFloatArray(file, f.begin() + 2, f.end() - 1);
  WriteBinaryIbmFloat(file, f[7]);

  file.close();

  std::ifstream ifile(filepath.c_str(), std::ios::binary | std::ios::in);
  unsigned char buffer[32];
  ifile.read(reinterpret_cast<char*>(buffer), 32);
  ifile.close();

  // 0.0F
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // -1.0F
  BOOST_CHECK_EQUAL(buffer[4], 0xc1);
  BOOST_CHECK_EQUAL(buffer[5], 0x10);
  BOOST_CHECK_EQUAL(buffer[6], 0x00);
  BOOST_CHECK_EQUAL(buffer[7], 0x00);

  // 1.0F
  BOOST_CHECK_EQUAL(buffer[8], 0x41);
  BOOST_CHECK_EQUAL(buffer[9], 0x10);
  BOOST_CHECK_EQUAL(buffer[10], 0x00);
  BOOST_CHECK_EQUAL(buffer[11], 0x00);


  // 1.5622F
  BOOST_CHECK_EQUAL(buffer[20], 0x41);
  BOOST_CHECK_EQUAL(buffer[21], 0x18);
  BOOST_CHECK_EQUAL(buffer[22], 0xfe);
  BOOST_CHECK_EQUAL(buffer[23], 0xc5);

  // 3.141592653589793e-12F
  BOOST_CHECK_EQUAL(buffer[28], 0x37);
  BOOST_CHECK_EQUAL(buffer[29], 0x37);
  BOOST_CHECK_EQUAL(buffer[30], 0x44);
  BOOST_CHECK_EQUAL(buffer[31], 0x79);

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  for (int i = 0; i < 8; ++i) {
    if (i == 6)
      // Not true on VS2012
      // BOOST_CHECK(boost::math::isnan(ReadBinaryIbmFloat(ifile)));
      ReadBinaryIbmFloat(ifile);
    else if (i == 5 || i == 7) // IBM float has lower precision than IEEE float:
      BOOST_CHECK_CLOSE_FRACTION(ReadBinaryIbmFloat(ifile), f[i], 3e-7);
    else
      BOOST_CHECK_EQUAL(ReadBinaryIbmFloat(ifile), f[i]);
  }
  ifile.close();

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), false);

  std::vector<float> in_f(8);
  ReadBinaryIbmFloatArray(ifile, in_f.begin(), 8);

  // Not true on VS2012
  // BOOST_CHECK(boost::math::isnan(in_f[6]));

  for (int i = 0; i < 8; ++i) {
    if (i == 5 || i == 7) // IBM float has lower precision than IEEE float:
      BOOST_CHECK_CLOSE_FRACTION(in_f[i], f[i], 3e-7);
    else if (i != 6)
      BOOST_CHECK_EQUAL(in_f[i], f[i]);
  }

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), true);

  ifile.close();

  remove(filepath);

}


BOOST_AUTO_TEST_CASE(ReadWriteBinaryIbmFloatLETest)
{
  path filepath = temp_directory_path() / "test_IBM_float_le.bin";

  std::vector<float> f(8);
  f[0] = 0.0F;
  f[1] = -1.0F;
  f[2] = 1.0F;
  f[3] = std::numeric_limits<float>::infinity();
  f[4] = -std::numeric_limits<float>::infinity();
  f[5] = 1.5622F;
  f[6] = std::numeric_limits<float>::signaling_NaN();
  f[7] = 3.141592653589793e-12F;

  std::ofstream file(filepath.c_str(), std::ios::binary | std::ios::out);

  WriteBinaryIbmFloat(file, f[0], END_LITTLE_ENDIAN);
  WriteBinaryIbmFloat(file, f[1], END_LITTLE_ENDIAN);
  WriteBinaryIbmFloatArray(file, f.begin() + 2, f.end() - 1, END_LITTLE_ENDIAN);
  WriteBinaryIbmFloat(file, f[7], END_LITTLE_ENDIAN);

  file.close();

  std::ifstream ifile(filepath.c_str(), std::ios::binary | std::ios::in);
  unsigned char buffer[32];
  ifile.read(reinterpret_cast<char*>(buffer), 32);
  ifile.close();

  // 0.0F
  BOOST_CHECK_EQUAL(buffer[0], 0x00);
  BOOST_CHECK_EQUAL(buffer[1], 0x00);
  BOOST_CHECK_EQUAL(buffer[2], 0x00);
  BOOST_CHECK_EQUAL(buffer[3], 0x00);

  // -1.0F
  BOOST_CHECK_EQUAL(buffer[4], 0x00);
  BOOST_CHECK_EQUAL(buffer[5], 0x00);
  BOOST_CHECK_EQUAL(buffer[6], 0x10);
  BOOST_CHECK_EQUAL(buffer[7], 0xc1);

  // 1.0F
  BOOST_CHECK_EQUAL(buffer[8], 0x00);
  BOOST_CHECK_EQUAL(buffer[9], 0x00);
  BOOST_CHECK_EQUAL(buffer[10], 0x10);
  BOOST_CHECK_EQUAL(buffer[11], 0x41);

  // 1.5622F
  BOOST_CHECK_EQUAL(buffer[20], 0xc5);
  BOOST_CHECK_EQUAL(buffer[21], 0xfe);
  BOOST_CHECK_EQUAL(buffer[22], 0x18);
  BOOST_CHECK_EQUAL(buffer[23], 0x41);

  // 3.141592653589793e-12F
  BOOST_CHECK_EQUAL(buffer[28], 0x79);
  BOOST_CHECK_EQUAL(buffer[29], 0x44);
  BOOST_CHECK_EQUAL(buffer[30], 0x37);
  BOOST_CHECK_EQUAL(buffer[31], 0x37);

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  for (int i = 0; i < 8; ++i) {
    if (i == 6)
      // Not true on VS2012
      // BOOST_CHECK(boost::math::isnan(ReadBinaryIbmFloat(ifile, END_LITTLE_ENDIAN)));
      ReadBinaryIbmFloat(ifile, END_LITTLE_ENDIAN);
    else if (i == 5 || i == 7) // IBM float has lower precision than IEEE float:
      BOOST_CHECK_CLOSE_FRACTION(ReadBinaryIbmFloat(ifile, END_LITTLE_ENDIAN), f[i], 3e-7);
    else
      BOOST_CHECK_EQUAL(ReadBinaryIbmFloat(ifile, END_LITTLE_ENDIAN), f[i]);
  }
  ifile.close();

  ifile.open(filepath.c_str(), std::ios::binary | std::ios::in);

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), false);

  std::vector<float> in_f(8);
  ReadBinaryIbmFloatArray(ifile, in_f.begin(), 8, END_LITTLE_ENDIAN);

  // Not true on VS2012
  // BOOST_CHECK(boost::math::isnan(in_f[6]));

  for (int i = 0; i < 8; ++i) {
    if (i == 5 || i == 7) // IBM float has lower precision than IEEE float:
      BOOST_CHECK_CLOSE_FRACTION(in_f[i], f[i], 3e-7);
    else if (i != 6)
      BOOST_CHECK_EQUAL(in_f[i], f[i]);
  }

  BOOST_CHECK_EQUAL(CheckEndOfFile(ifile), true);

  ifile.close();

  remove(filepath);
}
