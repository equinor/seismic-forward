// $Id: fileio.hpp 1594 2017-07-06 14:06:55Z perroe $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// •  Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// •  Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
// SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef NRLIB_FILEIO_HPP
#define NRLIB_FILEIO_HPP

// #include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <boost/endian/conversion.hpp>

#include "stringtools.hpp"
#include "../exception/exception.hpp"

namespace NRLib {
  enum Endianess {END_LITTLE_ENDIAN,
                  END_BIG_ENDIAN};
  enum GridFileFormat {UNKNOWN             = -1,
                       STORM_PETRO_BINARY  = 0,
                       STORM_PETRO_ASCII   = 1,
                       STORM_FACIES_BINARY = 2,
                       STORM_FACIES_ASCII  = 3,
                       SGRI                = 4,
                       SEGY                = 5,
                       PLAIN_ASCII         = 6,
                       BINARY_TREND        = 7};

  /// \brief Open file for reading.
  void OpenRead(std::ifstream&          stream,
                const std::string&      filename,
                std::ios_base::openmode mode = std::ios_base::in);
  void OpenRead(std::fstream&          stream,
                const std::string&      filename,
                std::ios_base::openmode mode = std::ios_base::in);

  /// \brief Open file for writing.
  /// \param create_dir If true the directory is created if it does not exist.
  void OpenWrite(std::ofstream&          stream,
                 const std::string&      filename,
                 std::ios_base::openmode mode = std::ios_base::out,
                 bool                    create_dir = true);
  void OpenWrite(std::fstream&           stream,
                 const std::string&      filename ,
                 std::ios_base::openmode mode= std::ios_base::out ,
                 bool                    create_dir = true);

  void CopyFile(const std::string & from_path,
                const std::string & to_path,
                bool                allow_overwrite);

  bool FileExists(const std::string& filename);

  void CreateDirIfNotExists(const std::string & filename);

  void RemoveFile(const std::string & filename);

  /// \brief Finds the size of a file. Throws IOError if file not found.
  /// \return Size of file in bytes.
  unsigned long long FindFileSize(const std::string & filename);

  /// \brief Find type of file, for 3D grid files.
  /// \todo Move to a suitable place.
  int FindGridFileType(const std::string& filename);

  // ---------------------------------
  // ASCII read and write
  // ---------------------------------

  /// \brief Reads and discard all characters until and including
  ///        next newline.
  /// \param throw_if_non_whitespace If true, throw exception if non-whitespace
  ///                                character is encountered.
  void DiscardRestOfLine(std::istream& stream,
                         int&          line_num,
                         bool          throw_if_non_whitespace);

  /// \brief Gets the next line that not only contains whitespace.
  /// Returns the stream. Stream will be bad if end-of-file is reached.
  std::istream& GetNextNonEmptyLine(std::istream & stream,
                                    int          & line_num,
                                    std::string  & line);

  /// \brief Returns last non-empty line in file.
  /// Stream will be at end of file.
  /// \note Does not work on Windows for files larger than 2 GB. (Uses seekg + tellg)
  std::string  FindLastNonEmptyLine(std::istream             & stream,
                                    const std::ios::pos_type & max_line_len = 1000);

  /// \brief Returns next line that not only contains whitespace, or starts with the
  ///        given comment token.
  /// \throws EndOfFile if end of file is reached before next non-empty line.
  void SkipComments(std::istream & stream,
                    char           comment_token,
                    int          & line_num);

  /// \brief Check if end of file is reached.
  /// Discards all whitespace before end of file or next token.
  /// \return True if end of file is reached, else false.
  bool CheckEndOfFile(std::istream& stream);

  /// \brief Read next white-space seperated token from file.
  ///        Updates line number when new line is read in.
  ///        If end of file is reached the state of stream is set to eof.
  ///        If ReadNextToken is called with a stream that is not good, the state
  ///        is set to fail.
  std::istream& ReadNextToken(std::istream& stream, std::string& s, int& line);

  /// \brief Gets next token from file, and parses it as type T
  ///        Might be relatively slow, due to typechecking, and counting line
  ///        numbers.
  /// \throws EndOfFile if end of file is reached.
  /// \throws Exception if the token could not be parsed as the given type.
  template <typename T>
  T ReadNext(std::istream& stream, int& line);

  /// \brief Reads the next, possibly quoted, string.
  void ReadNextQuoted(std::istream& stream, char quote, std::string& s, int& line);

  /// \brief Writes array
  template <typename I>
  void WriteAsciiArray(std::ostream& stream,
                       I             begin,
                       I             end,
                       int           n_per_line = 6);

  /// \brief Gets sequence with elements of type T from input stream.
  ///        Might be relatively slow, due to typechecking, and counting line
  ///        numbers.
  /// \note  The container must already be big enough to read all n
  ///        elements.
  template <typename I>
  I ReadAsciiArray(std::istream& stream, I begin, size_t n, int& line);

  /// \brief Gets sequence with elements of type T from input stream.
  ///        Does no type checking, and does not count line numbers.
  /// \note  The container must already be big enough to read all n
  ///        elements.
  template <typename I>
  I ReadAsciiArrayFast(std::istream& stream, I begin, size_t n);

  /// \brief Gets sequence with elements of type T from input stream.
  ///        Does no type checking, and does not count line numbers.
  ///        Reads the rest of the file, and does no type checking at all.
  ///        If a double is attempted read as an int, the first integer part
  ///        of the double is read, while the fraction part is discarded.
  ///        For int, float or double this function is 10x as fast as
  ///        GetAsciiArrayFast on Windows.
  /// \note  The container must already be big enough to read all n
  ///        elements.
  template <typename I>
  I ReadAsciiArrayFastRestOfFile(std::istream& stream, I begin, size_t n);

  // ---------------------------------
  // Binary read and write
  //
  // 2-byte integer
  // ---------------------------------

  /// \brief Write a 2-byte integer to a binary file.
  template <typename F>
  void WriteBinaryShort(F& stream,
                        short s,
                        Endianess number_representation = END_BIG_ENDIAN);

  /// \brief Read a 2-byte integer from a binary file.
  template <typename F>
  short ReadBinaryShort(F& stream,
                        Endianess number_representation = END_BIG_ENDIAN);

  /// \brief Write an array of binary 2-byte integers.
  template <typename F, typename I>
  void WriteBinaryShortArray(F& stream,
                             I begin,
                             I end,
                             Endianess number_representation = END_BIG_ENDIAN);

  /// \brief Read an array of binary 2-byte integers.
  /// \note  The container must already be big enough to read all n
  ///        elements.
  template<typename F, typename I>
  I ReadBinaryShortArray(F& stream,
                         I begin,
                         size_t n,
                         Endianess number_representation = END_BIG_ENDIAN);

  // ---------------------------------
  // 4-byte integer
  // ---------------------------------

  /// \brief Write a 4-byte integer to a binary file.
  template <typename F>
  void WriteBinaryInt(F& stream,
                      int i,
                      Endianess number_representation = END_BIG_ENDIAN);

  /// \brief Read a 4-byte integer from a binary file.
  template <typename F>
  int ReadBinaryInt(F& stream,
                    Endianess number_representation = END_BIG_ENDIAN);

  /// \brief Write an array of binary 4-byte integers.
  template<typename F, typename I>
  void WriteBinaryIntArray(F& stream,
                           I begin,
                           I end,
                           Endianess number_representation = END_BIG_ENDIAN);

  /// \brief Read an array of binary 4-byte integers.
  /// \note  The container must already be big enough to read all n
  ///        elements.
  template<typename F, typename I>
  I ReadBinaryIntArray(F& stream,
                       I begin,
                       size_t n,
                       Endianess number_representation = END_BIG_ENDIAN);

  // ---------------------------------
  // 4-byte IEEE floating point number
  // ---------------------------------

  /// \brief Write a 4-byte float on standard IEEE format.
  template <typename F>
  void WriteBinaryFloat(F& stream,
                        float f,
                        Endianess number_representation = END_BIG_ENDIAN);

  /// \brief Read a 4-byte float on standard IEEE format.
  template <typename F>
  float ReadBinaryFloat(F& stream,
                        Endianess number_representation = END_BIG_ENDIAN);

  /// \brief Write an array of 4-byte floats on standard IEEE format.
  template<typename F, typename I>
  void WriteBinaryFloatArray(F& stream,
                             I begin,
                             I end,
                             Endianess number_representation = END_BIG_ENDIAN);

  /// \brief Read an array of 4-byte floats on standard IEEE format.
  /// \note  The container must already be big enough to read all n
  ///        elements.
  template<typename F, typename I>
  I ReadBinaryFloatArray(F& stream,
                         I begin,
                         size_t n,
                         Endianess number_representation = END_BIG_ENDIAN);

  // ---------------------------------
  // 8-byte IEEE floating point number
  // ---------------------------------

  /// \brief Write a 8-byte float on standard IEEE format.
  template <typename F>
  void WriteBinaryDouble(F& stream,
                         double d,
                         Endianess number_representation = END_BIG_ENDIAN);

  /// \brief Read a 8-byte float on standard IEEE format.
  template <typename F>
  double ReadBinaryDouble(F& stream,
                          Endianess number_representation = END_BIG_ENDIAN);

  /// \brief Write an array of 8-byte floats on standard IEEE format.
  template<typename F, typename I>
  void WriteBinaryDoubleArray(F& stream,
                              I begin,
                              I end,
                              Endianess number_representation = END_BIG_ENDIAN);

  /// \brief Read an array of 8-byte floats on standard IEEE format.
  /// \note  The container must already be big enough to read all n
  ///        elements.
  template<typename F, typename I>
  I ReadBinaryDoubleArray(F& stream,
                          I begin,
                          size_t n,
                          Endianess number_representation = END_BIG_ENDIAN);

  // ---------------------------------
  // 4-byte IBM floating point number
  // ---------------------------------

  /// \brief Write a 4-byte float on standard IEEE format.
  template <typename F>
  void WriteBinaryIbmFloat(F& stream,
                           float f,
                           Endianess number_representation = END_BIG_ENDIAN);

  /// \brief Read a 4-byte float on standard IEEE format.
  template <typename F>
  float ReadBinaryIbmFloat(F& stream,
                           Endianess number_representation = END_BIG_ENDIAN);

  /// \brief Write an array of 4-byte floats on standard IEEE format.
  template<typename F, typename I>
  void WriteBinaryIbmFloatArray(F& stream,
                                I begin,
                                I end,
                                Endianess number_representation = END_BIG_ENDIAN);

  /// \brief Read an array of 4-byte floats on standard IEEE format.
  /// \note  The container must already be big enough to read all n
  ///        elements.
  template<typename F, typename I>
  I ReadBinaryIbmFloatArray(F& stream,
                            I begin,
                            size_t n,
                            Endianess number_representation = END_BIG_ENDIAN);

bool IgnoreComment(std::ifstream& file,
                   char chin);


  // ----------------------------------------------------------
  // Generic seek-function.
  // Requires FILE *, since streams fail after 2G on windows.
  // ----------------------------------------------------------

  /// Seek to given position in file, works for large files.
  int Seek(FILE * file, long long offset, int origin);

  /// Tell position in file, works for large files.
  long long Tell(FILE * file);


  // -------------------------------------------
  // Parsing of binary data from char buffer
  // -------------------------------------------

  /// Parse unsigned 32-bit integer from big-endian buffer.
  inline void ParseInt16BE(const char* buffer, /*int16_t*/ short& ui);

  /// Parse unsigned 32-bit integer from big-endian buffer.
  inline void ParseInt32BE(const char* buffer, /*int32_t*/ int& ui);

  /// Parse unsigned 32-bit integer from big-endian buffer.
  inline void ParseUInt16BE(const char* buffer, /*uint16_t*/ unsigned short& ui);

  /// Parse unsigned 32-bit integer from big-endian buffer.
  inline void ParseUInt32BE(const char* buffer, /*uint32_t*/ unsigned int& ui);

  /// Parse IEEE single-precision float from big-endian buffer.
  inline void ParseIEEEFloatBE(const char* buffer, float& f);

  /// Parse IEEE double-precision float from big-endian buffer.
  inline void ParseIBMFloatBE(const char* buffer, float& f);

  // -------------------------------------------
  // Writing binary data to char buffer
  // -------------------------------------------

  /// Parse unsigned 32-bit integer from big-endian buffer.
  inline void WriteInt16BE(char* buffer, /*int16_t*/ short ui);

  /// Parse unsigned 32-bit integer from big-endian buffer.
  inline void WriteInt32BE(char* buffer, /*int32_t*/ int ui);

  /// Write unsigned 32-bit integer to big-endian buffer.
  inline void WriteUInt16BE(char* buffer, /*uint16_t*/ unsigned short us);

  /// Write unsigned 32-bit integer to big-endian buffer.
  inline void WriteUInt32BE(char* buffer, /*uint32_t*/ unsigned int ui);


namespace NRLibPrivate {
  /// \todo Use stdint.h if available.
  // typedef unsigned int uint32_t;
  // typedef unsigned long long uint64_t;

  union FloatAsInt {
    /*uint32_t*/ unsigned int ui;
    float f;
  };

  union DoubleAsLongLong {
    /*uint64_t*/ unsigned long long ull;
    double d;
  };


  /// Parse unsigned 32-bit integer from little-endian buffer.
  inline void ParseUInt16LE(const char* buffer, /*uint16_t*/ unsigned short& ui);

  /// Write unsigned 32-bit integer to little-endian buffer.
  inline void WriteUInt16LE(char* buffer, /*uint16_t*/ unsigned short us);

  /// Parse unsigned 32-bit integer from little-endian buffer.
  inline void ParseUInt32LE(const char* buffer, /*uint32_t*/ unsigned int& ui);

  /// Write unsigned 32-bit integer to little-endian buffer.
  inline void WriteUInt32LE(char* buffer, /*uint32_t*/ unsigned int ui);

  /// Parse IEEE single-precision float from little-endian buffer.
  inline void ParseIEEEFloatLE(const char* buffer, float& f);

  /// Write IEEE single-precision float to big-endian buffer.
  inline void WriteIEEEFloatBE(char* buffer, float f);

  /// Write IEEE single-precision float to little-endian buffer.
  inline void WriteIEEEFloatLE(char* buffer, float f);

  /// Parse IEEE double-precision float from big-endian buffer.
  inline void ParseIEEEDoubleBE(const char* buffer, double& d);

  /// Parse IEEE double-precision float from little-endian buffer.
  inline void ParseIEEEDoubleLE(const char* buffer, double& d);

  /// Write IEEE double-precision float to big-endian buffer.
  inline void WriteIEEEDoubleBE(char* buffer, double d);

  /// Write IEEE double-precision float to little-endian buffer.
  inline void WriteIEEEDoubleLE(char* buffer, double d);

  /// Parse IEEE double-precision float from little-endian buffer.
  inline void ParseIBMFloatLE(const char* buffer, float& f);

  /// Write IEEE double-precision float to big-endian buffer.
  inline void WriteIBMFloatBE(char* buffer, float f);

  /// Write IEEE double-precision float to little-endian buffer.
  inline void WriteIBMFloatLE(char* buffer, float f);

  /// Inplace translation from IEEE to IBM 4-byte floating point format.
  inline void Ieee2Ibm(boost::uint32_t& in);

  /// Inplace translation from IEEE to IBM 4-byte floating point format.
  inline void Ibm2Ieee(boost::uint32_t& in);
} // namespace NRLibPrivate

} // namespace NRLib


// ========== INLINE AND TEMPLATE FUNCTION DEFINITIONS =========

namespace NRLib { // Needed to prevent gcc compilation error.

template <>
inline std::string ReadNext<std::string>(std::istream& stream, int& line)
{
  std::string s;
  ReadNextToken(stream, s, line);
  if (s == "")
    throw EndOfFile();
  return s;
}


template <typename T>
T ReadNext(std::istream& stream, int& line)
{
  std::string s;
  ReadNextToken(stream, s, line);
  if (s == "")
    throw EndOfFile();
  return ParseType<T>(s);
}


} // namespace NRLib


template <typename I>
void NRLib::WriteAsciiArray(std::ostream& stream,
                             I             begin,
                             I             end,
                             int           n_per_line)
{
  int count = 1;
  for (I it = begin; it != end; ++it, ++count) {
    stream << *it << ' ';
    if (count % n_per_line == 0) {
      stream << '\n';
    }
  }
}


template <typename I>
I NRLib::ReadAsciiArray(std::istream& stream, I begin, size_t n, int& line)
{
  typedef typename std::iterator_traits<I>::value_type T;
  for (size_t i = 0; i < n; ++i) {
    *begin = ReadNext<T>(stream, line);
    ++begin;
  }
  return begin;
}


template <typename I>
I NRLib::ReadAsciiArrayFast(std::istream& stream, I begin, size_t n)
{
  for (size_t i = 0; i < n; ++i) {
    stream >> *begin;
    ++begin;
    if ( stream.eof() ) {
      throw EndOfFile();
    }
    if ( stream.fail() ) {
      stream.clear();
      std::string nextToken;
      stream >> nextToken;
      throw Exception("Failure during reading element " + ToString(static_cast<unsigned int>(i)) + " of array. "
        + "Next token is " + nextToken + "\n");
    }
  }
  return begin;
}


template <typename I>
I NRLib::ReadAsciiArrayFastRestOfFile(std::istream& stream, I begin, size_t n)
{
  std::streampos pos = stream.tellg();
  stream.seekg(0, std::ios_base::end);
  std::streampos end = stream.tellg();
  stream.seekg(pos);
  size_t len = static_cast<size_t>(end - pos);
  std::string buffer(len, ' ');
  stream.read(&buffer[0], len);

  return ParseAsciiArrayFast(buffer, begin, n);
}


template <typename F>
void NRLib::WriteBinaryShort(F& stream,
                             short s,
                             Endianess file_format)
{
  using namespace NRLib::NRLibPrivate;

  char buffer[2];

  switch (file_format) {
  case END_BIG_ENDIAN:
    WriteUInt16BE(buffer, static_cast<unsigned short>(s));
    break;
  case END_LITTLE_ENDIAN:
    NRLibPrivate::WriteUInt16LE(buffer, static_cast<unsigned short>(s));
    break;
  default:
    assert(false); // Invalid endianess descriptor
  }

  if (!stream.write(buffer, 2)) {
    throw Exception("Error writing to stream.");
  }
}


template <typename F>
short NRLib::ReadBinaryShort(F& stream,
                             Endianess file_format)
{
  using namespace NRLib::NRLibPrivate;

  unsigned short us;
  char buffer[2];

  if (!stream.read(buffer, 2)) {
    if (stream.eof())
      throw EndOfFile();
    else
      throw Exception("Error reading from stream (a).");
  }

  switch (file_format) {
  case END_BIG_ENDIAN:
    ParseUInt16BE(buffer, us);
    break;
  case END_LITTLE_ENDIAN:
    ParseUInt16LE(buffer, us);
    break;
  default:
    throw Exception("Invalid file format.");
  }

  return static_cast<short>(us);
}


template <typename F, typename I>
void NRLib::WriteBinaryShortArray(F& stream,
                                  I begin,
                                  I end,
                                  NRLib::Endianess number_representation)
{
  typename I::difference_type n_in = std::distance(begin, end);

  std::vector<boost::uint16_t> buffer(n_in);
  memcpy(&buffer[0], &begin[0], 2 * n_in);

  switch (number_representation) {
  case END_BIG_ENDIAN:
    for (int i = 0; i < n_in; ++i)
      boost::endian::native_to_big_inplace(buffer[i]);
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; i < n_in; ++i)
      boost::endian::native_to_little_inplace(buffer[i]);
    break;
  default:
    assert(false); // Invalid endianess descriptor
  }

  if (!stream.write(reinterpret_cast<const char*>(&buffer[0]),
                    static_cast<std::streamsize>(2 * n_in))) {
    throw Exception("Error writing to stream.");
  }
}


template <typename F, typename I>
I NRLib::ReadBinaryShortArray(F& stream,
                              I begin,
                              size_t n,
                              NRLib::Endianess number_representation)
{
  std::vector<boost::uint16_t> buffer(n);

  if (!stream.read(reinterpret_cast<char *>(&buffer[0]),
                   static_cast<std::streamsize>(2 * n))) {
    throw Exception("Error reading from stream (f).");
  }

  switch (number_representation) {
  case END_BIG_ENDIAN:
    for (size_t i = 0; i < n; ++i)
      boost::endian::big_to_native_inplace(buffer[i]);
    break;
  case END_LITTLE_ENDIAN:
    for (size_t i = 0; i < n; ++i)
      boost::endian::little_to_native_inplace(buffer[i]);
    break;
  default:
    assert(false); // Invalid endianess descriptor
  }

  std::memcpy(&begin[0], &buffer[0], 2 * n);

  return begin;
}


template <typename F>
void NRLib::WriteBinaryInt(F& stream,
                           int i,
                           Endianess file_format)
{
  using namespace NRLib::NRLibPrivate;

  char buffer[4];

  switch (file_format) {
  case END_BIG_ENDIAN:
    WriteUInt32BE(buffer, static_cast<unsigned int>(i));
    break;
  case END_LITTLE_ENDIAN:
    WriteUInt32LE(buffer, static_cast<unsigned int>(i));
    break;
  default:
    throw Exception("Invalid file format.");
  }

  if (!stream.write(buffer, 4)) {
    throw Exception("Error writing to stream.");
  }
}


template <typename F>
int NRLib::ReadBinaryInt(F& stream,
                         Endianess file_format)
{
  using namespace NRLib::NRLibPrivate;

  unsigned int ui;
  char buffer[4];

  if (!stream.read(buffer, 4)) {
    if (stream.eof())
      throw EndOfFile();
    else
      throw Exception("Error reading from stream (b).");
  }

  switch (file_format) {
  case END_BIG_ENDIAN:
    ParseUInt32BE(buffer, ui);
    break;
  case END_LITTLE_ENDIAN:
    ParseUInt32LE(buffer, ui);
    break;
  default:
    throw Exception("Invalid file format.");
  }

  return static_cast<int>(ui);
}


template <typename F, typename I>
void NRLib::WriteBinaryIntArray(F& stream,
                                I begin,
                                I end,
                                NRLib::Endianess number_representation)
{
  typename I::difference_type n_in = std::distance(begin, end);

  std::vector<boost::uint32_t> buffer(n_in);
  memcpy(&buffer[0], &begin[0], 4 * n_in);

  switch (number_representation) {
  case END_BIG_ENDIAN:
    for (int i = 0; i < n_in; ++i)
      boost::endian::native_to_big_inplace(buffer[i]);
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; i < n_in; ++i)
      boost::endian::native_to_little_inplace(buffer[i]);
    break;
  default:
    assert(false); // Invalid endianess descriptor
  }

  if (!stream.write(reinterpret_cast<const char*>(&buffer[0]),
                    static_cast<std::streamsize>(4 * n_in))) {
    throw Exception("Error writing to stream.");
  }
}


template <typename F, typename I>
I NRLib::ReadBinaryIntArray(F& stream,
                            I begin,
                            size_t n,
                            NRLib::Endianess number_representation)
{
  std::vector<boost::uint32_t> buffer(n);

  if (!stream.read(reinterpret_cast<char *>(&buffer[0]),
                   static_cast<std::streamsize>(4 * n))) {
    throw Exception("Error reading from stream (g).");
  }

  switch (number_representation) {
  case END_BIG_ENDIAN:
    for (size_t i = 0; i < n; ++i)
      boost::endian::big_to_native_inplace(buffer[i]);
    break;
  case END_LITTLE_ENDIAN:
    for (size_t i = 0; i < n; ++i)
      boost::endian::little_to_native_inplace(buffer[i]);
    break;
  default:
    assert(false); // Invalid endianess descriptor
  }

  std::memcpy(&begin[0], &buffer[0], 4 * n);

  return begin;
}


template <typename F>
void NRLib::WriteBinaryFloat(F& stream,
                             float f,
                             Endianess file_format)
{
  using namespace NRLib::NRLibPrivate;

  char buffer[4];

  switch (file_format) {
  case END_BIG_ENDIAN:
    WriteIEEEFloatBE(buffer, f);
    break;
  case END_LITTLE_ENDIAN:
    WriteIEEEFloatLE(buffer, f);
    break;
  default:
    throw Exception("Invalid file format.");
  }

  if (!stream.write(buffer, 4)) {
    throw Exception("Error writing to stream.");
  }
}


template <typename F>
float NRLib::ReadBinaryFloat(F& stream,
                             Endianess file_format)
{
  using namespace NRLib::NRLibPrivate;

  float f;
  char buffer[4];

  if (!stream.read(buffer, 4)) {
    if (stream.eof())
      throw EndOfFile();
    else
      throw Exception("Error reading from stream (c).");
  }

  switch (file_format) {
  case END_BIG_ENDIAN:
    ParseIEEEFloatBE(buffer, f);
    break;
  case END_LITTLE_ENDIAN:
    ParseIEEEFloatLE(buffer, f);
    break;
  default:
    throw Exception("Invalid file format.");
  }

  return f;
}


template <typename F, typename I>
void NRLib::WriteBinaryFloatArray(F& stream,
                                  I begin,
                                  I end,
                                  NRLib::Endianess number_representation)
{
  typename I::difference_type n_in = std::distance(begin, end);

  std::vector<boost::uint32_t> buffer(n_in);

  memcpy(&buffer[0], &begin[0], 4 * n_in);

  switch (number_representation) {
  case END_BIG_ENDIAN:
    for (int i = 0; i < n_in; ++i)
      boost::endian::native_to_big_inplace(buffer[i]);
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; i < n_in; ++i)
      boost::endian::native_to_little_inplace(buffer[i]);
    break;
  default:
    assert(false); // Invalid endianess descriptor
  }

  if (!stream.write(reinterpret_cast<const char*>(&buffer[0]),
                    static_cast<std::streamsize>(4 * n_in))) {
    throw Exception("Error writing to stream.");
  }
}


template <typename F, typename I>
I NRLib::ReadBinaryFloatArray(F& stream,
                              I begin,
                              size_t n,
                              NRLib::Endianess number_representation)
{
  std::vector<boost::uint32_t> buffer(n);

  if (!stream.read(reinterpret_cast<char *>(&buffer[0]), static_cast<std::streamsize>(4*n))) {
    throw Exception("Error reading from stream (h).");
  }

  switch (number_representation) {
  case END_BIG_ENDIAN:
    for (size_t i = 0; i < n; ++i)
      boost::endian::big_to_native_inplace(buffer[i]);
    break;
  case END_LITTLE_ENDIAN:
    for (size_t i = 0; i < n; ++i)
      boost::endian::little_to_native_inplace(buffer[i]);
    break;
  default:
    assert(false); // Invalid endianess descriptor
  }

  std::memcpy(&begin[0], &buffer[0], 4 * n);

  return begin;
}


template <typename F>
void NRLib::WriteBinaryDouble(F& stream,
                              double d,
                              Endianess file_format)
{
  using namespace NRLib::NRLibPrivate;

  char buffer[8];

  switch (file_format) {
  case END_BIG_ENDIAN:
    WriteIEEEDoubleBE(buffer, d);
    break;
  case END_LITTLE_ENDIAN:
    WriteIEEEDoubleLE(buffer, d);
    break;
  default:
    throw Exception("Invalid file format.");
  }

  if (!stream.write(buffer, 8)) {
    throw Exception("Error writing to stream.");
  }
}


template <typename F>
double NRLib::ReadBinaryDouble(F& stream,
                               Endianess file_format)
{
  using namespace NRLib::NRLibPrivate;

  double d;
  char buffer[8];

  if (!stream.read(buffer, 8)) {
    if (stream.eof())
      throw EndOfFile();
    else
      throw Exception("Error reading from stream (d).");
  }

  switch (file_format) {
  case END_BIG_ENDIAN:
    ParseIEEEDoubleBE(buffer, d);
    break;
  case END_LITTLE_ENDIAN:
    ParseIEEEDoubleLE(buffer, d);
    break;
  default:
    throw Exception("Invalid file format.");
  }

  return d;
}


template <typename F, typename I>
void NRLib::WriteBinaryDoubleArray(F& stream,
                                   I begin,
                                   I end,
                                   NRLib::Endianess number_representation)
{
  typename I::difference_type n_in = std::distance(begin, end);
  std::vector<boost::uint64_t> buffer(n_in);

  memcpy(&buffer[0], &begin[0], 8 * n_in);

  switch (number_representation) {
  case END_BIG_ENDIAN:
    for (int i = 0; i < n_in; ++i)
      boost::endian::native_to_big_inplace(buffer[i]);
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; i < n_in; ++i)
      boost::endian::native_to_little_inplace(buffer[i]);
    break;
  default:
    assert(false); // Invalid endianess descriptor
  }

  if (!stream.write(reinterpret_cast<const char*>(&buffer[0]),
                    static_cast<std::streamsize>(8 * n_in))) {
    throw Exception("Error writing to stream.");
  }
}


template <typename F, typename I>
I NRLib::ReadBinaryDoubleArray(F& stream,
                               I begin,
                               size_t n,
                               NRLib::Endianess number_representation)
{
  std::vector<boost::uint64_t> buffer(n);

  if (!stream.read(reinterpret_cast<char *>(&buffer[0]),
                   static_cast<std::streamsize>(8 * n))) {
    throw Exception("Error reading from stream (i).");
  }

  switch (number_representation) {
  case END_BIG_ENDIAN:
    for (size_t i = 0; i < n; ++i)
      boost::endian::big_to_native_inplace(buffer[i]);
    break;
  case END_LITTLE_ENDIAN:
    for (size_t i = 0; i < n; ++i)
      boost::endian::little_to_native_inplace(buffer[i]);
    break;
  default:
    assert(false); // Invalid endianess descriptor
  }

  std::memcpy(&begin[0], &buffer[0], 8 * n);

  return begin;
}


template <typename F, typename I>
void NRLib::WriteBinaryIbmFloatArray(F& stream,
                                     I begin,
                                     I end,
                                     NRLib::Endianess number_representation)
{
  typename I::difference_type n_in = std::distance(begin, end);
  std::vector<boost::uint32_t> buffer(n_in);

  memcpy(&buffer[0], &begin[0], 4 * n_in);

  for (int i = 0; i < n_in; ++i)
    NRLibPrivate::Ieee2Ibm(buffer[i]);

  switch (number_representation) {
  case END_BIG_ENDIAN:
    for (int i = 0; i < n_in; ++i)
      boost::endian::native_to_big_inplace(buffer[i]);
    break;
  case END_LITTLE_ENDIAN:
    for (int i = 0; i < n_in; ++i)
      boost::endian::native_to_little_inplace(buffer[i]);
    break;
  default:
    assert(false); // Invalid endianess descriptor
  }

  if (!stream.write(reinterpret_cast<const char*>(&buffer[0]),
                    static_cast<std::streamsize>(4 * n_in))) {
    throw Exception("Error writing to stream.");
  }
}


template <typename F, typename I>
I NRLib::ReadBinaryIbmFloatArray(F& stream,
                                 I      begin,
                                 size_t n,
                                 NRLib::Endianess number_representation)
{
  std::vector<boost::uint32_t> buffer(n);

  std::string error;
  if (!stream.read(reinterpret_cast<char *>(&buffer[0]), static_cast<std::streamsize>(4*n))) {
    if (stream.eof())
      error = "Error reading binary IBM float array. Trying to read 4*" + NRLib::ToString(n) + " elements when end-of-file was reached.\n";
    else {
      error = "Error reading binary IBM float array. Hardware error? Full disk?\n";
    }
  }

  switch (number_representation) {
  case END_BIG_ENDIAN:
    for (size_t i = 0; i < n; ++i)
      boost::endian::big_to_native_inplace(buffer[i]);
    break;
  case END_LITTLE_ENDIAN:
    for (size_t i = 0; i < n; ++i)
      boost::endian::little_to_native_inplace(buffer[i]);
    break;
  default:
    assert(false); // Invalid endianess descriptor
  }

  for (size_t i = 0; i < n; ++i)
    NRLibPrivate::Ibm2Ieee(buffer[i]);

  memcpy(&begin[0], &buffer[0], 4 * n);

  if (error != "") {
    throw Exception(error);
  }

  return begin;
}


template <typename F>
void NRLib::WriteBinaryIbmFloat(F& stream,
                                float f,
                                Endianess file_format)
{
  char buffer[4];

  switch (file_format) {
  case END_BIG_ENDIAN:
    NRLibPrivate::WriteIBMFloatBE(buffer, f);
    break;
  case END_LITTLE_ENDIAN:
    NRLibPrivate::WriteIBMFloatLE(buffer, f);
    break;
  default:
    assert(false); // Invalid endianess descriptor
  }

  if (!stream.write(buffer, 4)) {
    throw Exception("Error writing to stream.");
  }
}


template <typename F>
float NRLib::ReadBinaryIbmFloat(F& stream,
                                Endianess file_format)
{
  float f = 0.0;
  char buffer[4];

  if (!stream.read(buffer, 4)) {
    if (stream.eof())
      throw EndOfFile();
    else
      throw Exception("Error reading from stream (e).");
  }

  switch (file_format) {
  case END_BIG_ENDIAN:
    ParseIBMFloatBE(buffer, f);
    break;
  case END_LITTLE_ENDIAN:
    NRLibPrivate::ParseIBMFloatLE(buffer, f);
    break;
  default:
    assert(false); // Invalid endianess descriptor
  }

  return f;
}




static const unsigned short smasks[2] = {0x00ff, 0xff00};

static const unsigned int masks[4] =
{0x000000ff, 0x0000ff00, 0x00ff0000, 0xff000000};

static const unsigned long long lmasks[8] =
{0x00000000000000ffULL, 0x000000000000ff00ULL,
 0x0000000000ff0000ULL, 0x00000000ff000000ULL,
 0x000000ff00000000ULL, 0x0000ff0000000000ULL,
 0x00ff000000000000ULL, 0xff00000000000000ULL};


// -----------------  16 bit integer ---------------------

void NRLib::ParseInt16BE(const char* buffer,
                         /*uint16_t*/ short& us)
{
  us = static_cast<unsigned char>(buffer[0]);
  us <<= 8;
  us |= static_cast<unsigned char>(buffer[1]);
}


// -----------------  16 bit integer ---------------------

void NRLib::ParseUInt16BE(const char* buffer,
                          /*uint16_t*/ unsigned short& us)
{
  us = static_cast<unsigned char>(buffer[0]);
  us <<= 8;
  us |= static_cast<unsigned char>(buffer[1]);
}


void NRLib::NRLibPrivate::ParseUInt16LE(const char* buffer,
                                    /*uint32_t*/ unsigned short& us)
{
  us = static_cast<unsigned char>(buffer[1]);
  us <<= 8;
  us |= static_cast<unsigned char>(buffer[0]);
}


void NRLib::WriteUInt16BE(char* buffer,
                          /*uint32_t*/ unsigned short us)
{
  buffer[0] = static_cast<char>( (us & smasks[1]) >> 8 );
  buffer[1] = static_cast<char>( us & smasks[0] );
}


void NRLib::NRLibPrivate::WriteUInt16LE(char* buffer,
                                    /*uint32_t*/ unsigned short us)
{
  buffer[0] = static_cast<char>(us & smasks[0]);
  buffer[1] = static_cast<char>( (us & smasks[1]) >> 8 );
}


// -----------------  32 bit integer ---------------------


void NRLib::ParseInt32BE(const char* buffer,
                         /*int32_t*/ int& ui)
{
  ui = static_cast<unsigned char>(buffer[0]);
  for (int i = 1; i < 4; ++i) {
    ui <<= 8;
    ui |= static_cast<unsigned char>(buffer[i]);
  }
}

// --------------  32 bit unsigned integer -----------------


void NRLib::ParseUInt32BE(const char* buffer,
                                    /*uint32_t*/ unsigned int& ui)
{
  ui = static_cast<unsigned char>(buffer[0]);
  for (int i = 1; i < 4; ++i) {
    ui <<= 8;
    ui |= static_cast<unsigned char>(buffer[i]);
  }
}


void NRLib::NRLibPrivate::ParseUInt32LE(const char* buffer,
                                    /*uint32_t*/ unsigned int& ui)
{
  ui = static_cast<unsigned char>(buffer[3]);
  for (int i = 1; i < 4; ++i) {
    ui <<= 8;
    ui |= static_cast<unsigned char>(buffer[3-i]);
  }
}

void NRLib::WriteInt32BE(char* buffer,
                          /*uint32_t*/ int ui)
{
  WriteUInt32BE(buffer, static_cast<unsigned int>(ui));
}

void NRLib::WriteInt16BE(char* buffer,
                         /*uint32_t*/ short ui)
{
  WriteUInt16BE(buffer, static_cast<unsigned short>(ui));
}


void NRLib::WriteUInt32BE(char* buffer,
                          /*uint32_t*/ unsigned int ui)
{
  buffer[0] = static_cast<char>( (ui & masks[3]) >> 24);
  buffer[1] = static_cast<char>( (ui & masks[2]) >> 16 );
  buffer[2] = static_cast<char>( (ui & masks[1]) >> 8 );
  buffer[3] = static_cast<char>( ui & masks[0] );
}


void NRLib::NRLibPrivate::WriteUInt32LE(char* buffer,
                                    /*uint32_t*/ unsigned int ui)
{
  buffer[0] = static_cast<char>(ui & masks[0]);
  buffer[1] = static_cast<char>( (ui & masks[1]) >> 8 );
  buffer[2] = static_cast<char>( (ui & masks[2]) >> 16 );
  buffer[3] = static_cast<char>( (ui & masks[3]) >> 24 );
}


// -----------------  32 bit IEEE floating point -------------------


// Big endian number representation.
void NRLib::ParseIEEEFloatBE(const char* buffer, float& f)
{
  NRLibPrivate::FloatAsInt tmp;
  tmp.ui = static_cast<unsigned char>(buffer[0]);
  for (int i = 1; i < 4; ++i) {
    tmp.ui <<= 8;
    tmp.ui |= static_cast<unsigned char>(buffer[i]);
  }
  f = tmp.f;
}

// Little endian number representation.
void NRLib::NRLibPrivate::ParseIEEEFloatLE(const char* buffer, float& f)
{
  FloatAsInt tmp;
  tmp.ui = static_cast<unsigned char>(buffer[3]);
  for (int i = 1; i < 4; ++i) {
    tmp.ui <<= 8;
    tmp.ui |= static_cast<unsigned char>(buffer[3-i]);
  }
  f = tmp.f;
}


// Big endian number representation.
void NRLib::NRLibPrivate::WriteIEEEFloatBE(char* buffer, float f)
{
  FloatAsInt tmp;
  tmp.f = f;
  buffer[0] = static_cast<char>( (tmp.ui & masks[3]) >> 24);
  buffer[1] = static_cast<char>( (tmp.ui & masks[2]) >> 16 );
  buffer[2] = static_cast<char>( (tmp.ui & masks[1]) >> 8 );
  buffer[3] = static_cast<char>( tmp.ui & masks[0] );
}

// Little endian number representation.
void NRLib::NRLibPrivate::WriteIEEEFloatLE(char* buffer, float f)
{
  FloatAsInt tmp;
  tmp.f = f;
  buffer[0] = static_cast<char>(tmp.ui & masks[0]);
  buffer[1] = static_cast<char>( (tmp.ui & masks[1]) >> 8 );
  buffer[2] = static_cast<char>( (tmp.ui & masks[2]) >> 16 );
  buffer[3] = static_cast<char>( (tmp.ui & masks[3]) >> 24 );
}


// -----------------  64 bit IEEE floating point --------------------


// Big endian number representation.
void NRLib::NRLibPrivate::ParseIEEEDoubleBE(const char* buffer, double& d)
{
  DoubleAsLongLong tmp;
  tmp.ull = static_cast<unsigned char>(buffer[0]);
  for (int i = 1; i < 8; ++i) {
    tmp.ull <<= 8;
    tmp.ull |= static_cast<unsigned char>(buffer[i]);
  }
  d = tmp.d;
}

// Little endian number representation.
void NRLib::NRLibPrivate::ParseIEEEDoubleLE(const char* buffer, double& d)
{
  DoubleAsLongLong tmp;
  tmp.ull = static_cast<unsigned char>(buffer[7]);
  for (int i = 1; i < 8; ++i) {
    tmp.ull <<= 8;
    tmp.ull |= static_cast<unsigned char>(buffer[7-i]);
  }
  d = tmp.d;
}


// Big endian number representation.
void NRLib::NRLibPrivate::WriteIEEEDoubleBE(char* buffer, double d)
{
  DoubleAsLongLong tmp;
  tmp.d = d;
  buffer[0] = static_cast<char>( (tmp.ull & lmasks[7]) >> 56);
  buffer[1] = static_cast<char>( (tmp.ull & lmasks[6]) >> 48);
  buffer[2] = static_cast<char>( (tmp.ull & lmasks[5]) >> 40);
  buffer[3] = static_cast<char>( (tmp.ull & lmasks[4]) >> 32);
  buffer[4] = static_cast<char>( (tmp.ull & lmasks[3]) >> 24);
  buffer[5] = static_cast<char>( (tmp.ull & lmasks[2]) >> 16);
  buffer[6] = static_cast<char>( (tmp.ull & lmasks[1]) >> 8 );
  buffer[7] = static_cast<char>(  tmp.ull & lmasks[0] );
}

// Little endian number representation.
void NRLib::NRLibPrivate::WriteIEEEDoubleLE(char* buffer, double d)
{
  DoubleAsLongLong tmp;
  tmp.d = d;
  buffer[0] = static_cast<char>(  tmp.ull & lmasks[0]);
  buffer[1] = static_cast<char>( (tmp.ull & lmasks[1]) >> 8 );
  buffer[2] = static_cast<char>( (tmp.ull & lmasks[2]) >> 16);
  buffer[3] = static_cast<char>( (tmp.ull & lmasks[3]) >> 24);
  buffer[4] = static_cast<char>( (tmp.ull & lmasks[4]) >> 32);
  buffer[5] = static_cast<char>( (tmp.ull & lmasks[5]) >> 40);
  buffer[6] = static_cast<char>( (tmp.ull & lmasks[6]) >> 48);
  buffer[7] = static_cast<char>( (tmp.ull & lmasks[7]) >> 56);
}


// -----------------  32 bit IBM floating point -------------------

static unsigned int IEEEMAX  = 0x7FFFFFFF;
static unsigned int IEMAXIB  = 0x611FFFFF;
static unsigned int IEMINIB  = 0x21200000;


/// Converts IBM float (represented as 32-bit int) to IEEE float.
void NRLib::NRLibPrivate::Ibm2Ieee(boost::uint32_t& in)
{
  static int it[8] = { 0x21800000, 0x21400000, 0x21000000, 0x21000000,
          0x20c00000, 0x20c00000, 0x20c00000, 0x20c00000 };
  static int mt[8] = { 8, 4, 2, 2, 1, 1, 1, 1 };
  /*uint32_t*/ unsigned int manthi, iexp, inabs;
  int ix;

  manthi = in & 0x00ffffff;
  ix     = manthi >> 21;
  iexp   = ( ( in & 0x7f000000 ) - it[ix] ) << 1;
  manthi = manthi * mt[ix] + iexp;
  inabs  = in & 0x7fffffff;
  if ( inabs > IEMAXIB ) manthi = IEEEMAX;
  manthi = manthi | ( in & 0x80000000 );
  in = ( inabs < IEMINIB ) ? 0 : manthi;
}


/// Converts IEEE float (represented as 32-bit int) to IBM float.
void NRLib::NRLibPrivate::Ieee2Ibm(boost::uint32_t& in)
{
  static int it[4] = { 0x21200000, 0x21400000, 0x21800000, 0x22100000 };
  static int mt[4] = { 2, 4, 8, 1 };
  /*uint32_t*/ unsigned int manthi, iexp, ix;

  ix     = ( in & 0x01800000 ) >> 23;
  iexp   = ( ( in & 0x7e000000 ) >> 1 ) + it[ix];
  manthi = ( mt[ix] * ( in & 0x007fffff) ) >> 3;
  manthi = (manthi + iexp) | ( in & 0x80000000 );
  in     = ( in & 0x7fffffff ) ? manthi : 0;
}


// Big endian number representation.
void NRLib::ParseIBMFloatBE(const char* buffer, float& f)
{
  NRLibPrivate::FloatAsInt tmp;
  tmp.ui = static_cast<unsigned char>(buffer[0]);
  for (int i = 1; i < 4; ++i) {
    tmp.ui <<= 8;
    tmp.ui |= static_cast<unsigned char>(buffer[i]);
  }

  NRLibPrivate::Ibm2Ieee(tmp.ui);
  f = tmp.f;
}


// Little endian number representation.
void NRLib::NRLibPrivate::ParseIBMFloatLE(const char* buffer, float& f)
{
  FloatAsInt tmp;
  tmp.ui = static_cast<unsigned char>(buffer[3]);
  for (int i = 1; i < 4; ++i) {
    tmp.ui <<= 8;
    tmp.ui |= static_cast<unsigned char>(buffer[3-i]);
  }

  NRLibPrivate::Ibm2Ieee(tmp.ui);
  f = tmp.f;
}


// Big endian number representation.
void NRLib::NRLibPrivate::WriteIBMFloatBE(char* buffer, float f)
{
  FloatAsInt tmp;
  tmp.f = f;
  Ieee2Ibm(tmp.ui);
  buffer[0] = static_cast<char>( (tmp.ui & masks[3]) >> 24);
  buffer[1] = static_cast<char>( (tmp.ui & masks[2]) >> 16 );
  buffer[2] = static_cast<char>( (tmp.ui & masks[1]) >> 8 );
  buffer[3] = static_cast<char>( tmp.ui & masks[0] );
}


// Little endian number representation.
void NRLib::NRLibPrivate::WriteIBMFloatLE(char* buffer, float f)
{
  FloatAsInt tmp;
  tmp.f = f;
  Ieee2Ibm(tmp.ui);
  buffer[0] = static_cast<char>(tmp.ui & masks[0]);
  buffer[1] = static_cast<char>( (tmp.ui & masks[1]) >> 8 );
  buffer[2] = static_cast<char>( (tmp.ui & masks[2]) >> 16 );
  buffer[3] = static_cast<char>( (tmp.ui & masks[3]) >> 24 );
}

#endif // NRLIB_FILEIO_HPP
