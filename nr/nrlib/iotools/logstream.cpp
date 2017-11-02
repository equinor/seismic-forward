// $Id: logstream.cpp 1392 2016-11-29 08:40:07Z perroe $

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

#include <fstream>

#include "logstream.hpp"
#include "fileio.hpp"

using namespace NRLib;

// ============================= LogStream ==================================

LogStream::LogStream(int level) {
  fullLevel_ = level;
}


LogStream::LogStream(const std::vector<int> & levels, bool ignore_general) {
  unsigned int i;
  fullLevel_ = 0;
  for (i = 0; i < levels.size(); i++) {
    if (ignore_general == false)
      fullLevel_ = (fullLevel_ | levels[i]);
    levels_.push_back(levels[i]);
  }
}


LogStream::~LogStream() { }


bool
LogStream::ShouldLog(int level, int phase)
{
  if (phase < static_cast<int>(levels_.size())) {
    if ((level & levels_[phase]) > 0)
      return true;
  }
  else if ((level & fullLevel_) > 0)
    return true;

  return false;
}


// =========================== FileLogStream ================================


FileLogStream::FileLogStream(const std::string& filename, int level)
  : LogStream(level)
{
  NRLib::OpenWrite(file_, filename);
}


FileLogStream::FileLogStream(const std::string      & filename,
                             const std::vector<int> & levels,
                             bool                     ignore_general)
  : LogStream(levels, ignore_general)
{
  NRLib::OpenWrite(file_, filename);
}


FileLogStream::~FileLogStream()
{
  file_.close();
}


void
FileLogStream::LogMessage(int level, const std::string & message) {
  if (ShouldLog(level)) {
    file_ << message;
    file_.flush();
  }
}


void
FileLogStream::LogMessage(int level, int phase, const std::string & message) {
  if (ShouldLog(level, phase)) {
    file_ << message;
    file_.flush();
  }
}


// ========================== ScreenLogStream ===============================


ScreenLogStream::ScreenLogStream(int level)
  : LogStream(level),
    is_writing_progress_(false)
{ }


ScreenLogStream::ScreenLogStream(const std::vector<int> & levels,
                                 bool                     ignore_general)
  : LogStream(levels, ignore_general),
    is_writing_progress_(false)
{ }


ScreenLogStream::~ScreenLogStream()
{ }


void
ScreenLogStream::LogMessage(int level, const std::string & message) {
  if (ShouldLog(level)) {
    if (is_writing_progress_) {
      std::cout << "\n\n";
      is_writing_progress_ = false;
    }

    std::cout << message;
    std::cout.flush();
  }
}


void
ScreenLogStream::LogMessage(int level, int phase, const std::string & message) {
  if (ShouldLog(level, phase)) {
    if (is_writing_progress_) {
      std::cout << "\n\n";
      is_writing_progress_ = false;
    }

    std::cout << message;
    std::cout.flush();
  }
}


void
ScreenLogStream::UpdateProgress(double progress, const std::string & message) {
  if (!is_writing_progress_)
    WriteProgressHeader();

  int n_hats = 1 + static_cast<int>(progress * n_progress_hats_);

  std::cout << "\r  ";
  std::cout << std::string(n_hats, '^') << std::string(n_progress_hats_ - n_hats + 2, ' ');
  std::cout << message;
}


void
ScreenLogStream::WriteProgressHeader() {
  std::cout << "\n  0%       20%       40%       60%       80%      100%\n";
  std::cout << ("  |    |    |    |    |    |    |    |    |    |    |\n");
  is_writing_progress_ = true;
}
