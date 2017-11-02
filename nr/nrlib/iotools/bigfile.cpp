// $Id: bigfile.cpp 1607 2017-07-09 21:02:21Z perroe $

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

#include "bigfile.hpp"

#include "nrlib/exception/exception.hpp"
#include "nrlib/iotools/stringtools.hpp"

#include <cstdio>

#if defined(_MSC_VER)
#include <wchar.h>
#endif

#define NRLIB_BIGFILE_BUFFERSIZE 131072

using namespace NRLib;

namespace fs = boost::filesystem;


BigFile::BigFile()
  : filename_(""),
    file_(NULL),
    buffer_(0),
    error_code_(0)
{}


BigFile::~BigFile()
{
  close();
}


void BigFile::OpenRead(const fs::path& filename, bool is_binary)
{
  filename_ = filename;
#if defined(_MSC_VER)
  std::wstring mode = (is_binary) ? L"rb" : L"r";
#else
  std::string mode = (is_binary) ? "rb" : "r";
#endif

  if (file_)
    close();

  if (!fs::exists(filename)) {
    throw IOError("Failed to open " + filename.string() + " for reading: " +
                  "File does not exist.");
  }
  if (fs::is_directory(filename)) {
    throw IOError("Failed to open " + filename.string() + " for reading: " +
                  " It is a directory.");
  }
#if defined(_MSC_VER)
  file_ = _wfopen(filename.c_str(), mode.c_str());
#else
  file_ = fopen(filename.c_str(), mode.c_str());
#endif
  if (!file_) {
    error_code_ = errno;
    throw IOError("Failed to open " + filename.string() + " for reading: " +
                  GetErrorMessage());
  }

  buffer_ = new char[NRLIB_BIGFILE_BUFFERSIZE];
  if (!buffer_)
    throw Exception("Failed to open " + filename.string() +
                    ": Failed to allocate file buffer of size "
                    + NRLib::ToString(NRLIB_BIGFILE_BUFFERSIZE));

  if (setvbuf(file_, buffer_, _IOFBF, NRLIB_BIGFILE_BUFFERSIZE)) {
    error_code_ = errno;
    throw Exception("Failed to set buffer when opening " + filename.string()
                    + ": " + GetErrorMessage());
  }
}


void BigFile::OpenWrite(const boost::filesystem::path& filename, bool is_binary,
                        bool create_dir)
{
  filename_ = filename;
#if defined(_MSC_VER)
  std::wstring mode = (is_binary) ? L"wb" : L"w";
#else
  std::string mode = (is_binary) ? "wb" : "w";
#endif

  if (file_)
    close();

  try {
    if (fs::is_directory(filename)) {
      throw IOError("Failed to open " + filename.string() + " for writing: " +
                    " It is a directory.");
    }
    fs::path dir = filename.parent_path();
    if (!dir.empty()) {
      if (!fs::exists(dir)) {
        if (create_dir)
          fs::create_directories(dir);
        else {
          throw IOError("Failed to open " + filename.string() + " for writing: "
                        + "Parent directory does not exist.");
        }
      }
    }
  }
  catch (fs::filesystem_error& e) {
    throw IOError("Failed to open " + filename.string() + " for writing: " + e.what());
  }

#if defined(_MSC_VER)
  file_ = _wfopen(filename.c_str(), mode.c_str());
#else
  file_ = fopen(filename.string().c_str(), mode.c_str());
#endif

  if (!file_) {
    error_code_ = errno;
    throw IOError("Failed to open " + filename.string() + " for writing: " +
                  GetErrorMessage());
  }

  buffer_ = new char[NRLIB_BIGFILE_BUFFERSIZE];
  if (!buffer_)
    throw Exception("Failed to open " + filename.string() +
                    ": Failed to allocate file buffer of size "
                    + NRLib::ToString(NRLIB_BIGFILE_BUFFERSIZE));

  if (setvbuf(file_, buffer_, _IOFBF, NRLIB_BIGFILE_BUFFERSIZE)) {
    error_code_ = errno;
    throw Exception("Failed to set buffer when opening " + filename.string()
                    + ": " + GetErrorMessage());
  }
}


void BigFile::clear()
{
  if (file_ != NULL)
    clearerr(file_);
}


void BigFile::close()
{
  if (file_ != NULL)
    fclose(file_);
  file_ = NULL;

  delete[] buffer_;
  buffer_ = 0;
}


std::string BigFile::GetErrorMessage() const
{
#if defined(_MSC_VER)
  char msg[512];
  strerror_s(msg, 512, GetErrorCode());

  return std::string(msg);
#else
  return std::string(strerror(GetErrorCode()));
#endif
}
