// $Id: bigfile.hpp 1601 2017-07-09 19:34:28Z perroe $

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

#ifndef NRLIB_IOTOOLS_BIGFILE_H
#define NRLIB_IOTOOLS_BIGFILE_H

#include "nrlib/exception/exception.hpp"
#include "nrlib/iotools/fileio.hpp"

#include <boost/filesystem.hpp>

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>

namespace NRLib {

/// \brief Class used for reading and writing of big files.
///        Should as far as possible have the same system as fstream, but will not support stream operators.
class BigFile {
public:
  BigFile();

  ~BigFile();

  void OpenRead(const boost::filesystem::path& filename, bool is_binary = true);

  void OpenWrite(const boost::filesystem::path& filename,
                 bool is_binary = true,
                 bool create_dir = true);

  /// Clear error codes.
  void clear();
  void close();

  bool good() const { return (file_ != NULL && !ferror(file_) && !feof(file_)); }
  bool eof() const  { return (feof(file_) != 0); }
  bool bad() const { return (ferror(file_) != 0); }

  inline size_t read(char* s, size_t n);
  inline size_t write(const char* s, size_t n);

  int seek(long long offset, int origin)  { return NRLib::Seek(file_, offset, origin); }
  long long tell() const                  { return NRLib::Tell(file_); }

  const boost::filesystem::path& GetFileName() const { return filename_; }

  long long FileSize() const { return boost::filesystem::file_size(filename_); }

  /// Error code. Only valid if bad().
  int GetErrorCode() const { return error_code_; }

  /// System error message, associated with error code.
  std::string GetErrorMessage() const;

private:
  boost::filesystem::path   filename_;
  FILE                    * file_;
  char                    * buffer_;
  int                       error_code_;
};


// ============= INLINE FUNCTION DEFINITIONS ===========================

#define NRLIB_BIGFILE_CHUNK_MAX 1048576

size_t BigFile::read(char* s, size_t n)
{
  size_t n_read = fread(s, 1, n, file_);

  if (n_read != n) {
    if (ferror(file_)) {
      error_code_ = errno;
      throw IOError("Error reading file " + filename_.string() + ": " + GetErrorMessage());
    }
  }

  return n_read;
}


size_t BigFile::write(const char* s, size_t n)
{
  size_t n_written = 0;
  if (n < NRLIB_BIGFILE_CHUNK_MAX)
    n_written = fwrite(s, 1, n, file_);
  else {
    // Some FILE* implementation are slow when writing very large chunks.
    size_t current_n = fwrite(s, 1, NRLIB_BIGFILE_CHUNK_MAX, file_);
    n_written = current_n;
    while (current_n > 0 && n_written < n) {
      size_t n_to_write = std::min(static_cast<size_t>(NRLIB_BIGFILE_CHUNK_MAX), n - n_written);
      current_n = fwrite(s + n_written, 1, n_to_write, file_);
      n_written += current_n;
    }
  }

  if (n_written != n) {
    if (ferror(file_)) {
      error_code_ = errno;
      throw IOError("Error writing to file " + filename_.string() + ": " + GetErrorMessage());
    }
  }

  return n_written;
}


}

#endif // NRLIB_IOTOOLS_BIGFILE_H
