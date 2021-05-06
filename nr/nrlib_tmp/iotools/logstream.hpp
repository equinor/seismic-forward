// $Id: logstream.hpp 1392 2016-11-29 08:40:07Z perroe $

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

#ifndef NRLIB_LOGSTREAM_HPP
#define NRLIB_LOGSTREAM_HPP

#include <fstream>
#include <string>
#include <vector>


namespace NRLib {

  /// Abstract base class for log streams.
  class LogStream {
  public:
    /// Creates log stream
    /// \param level LimitLevels flag that controls which messages are logged to this stream.
    LogStream(int level);

    /// Creates log stream
    /// \param levels         LimitLevels flags for each phase.
    /// \param ignore_general If true, general logging (from NRLib) is ignored.
    LogStream(const std::vector<int> & levels, bool ignore_general = false);

    virtual ~LogStream();

    virtual void LogMessage(int level, const std::string & message) = 0;
    virtual void LogMessage(int level, int phase, const std::string & message) = 0;
    virtual void UpdateProgress(double progress, const std::string & message) = 0;

  protected:
    /// Return true if message with given level should be logged by stream.
    bool ShouldLog(int level)             { return (fullLevel_ & level) != 0; }
    bool ShouldLog(int level, int phase);

  private:
    std::vector<int> levels_;
    int fullLevel_;
  };


  /// Log stream for logging to files.
  class FileLogStream : public LogStream {
  public:
    FileLogStream(const std::string& filename, int level);
    FileLogStream(const std::string& filename, const std::vector<int> & levels, bool ignore_general = false);
    ~FileLogStream();

    void LogMessage(int level, const std::string & message);
    void LogMessage(int level, int phase, const std::string & message);
    void UpdateProgress(double /*progress*/, const std::string & /*message*/) {}
  private:
    std::ofstream file_;
  };


  /// Log stream for logging to screen, using stdout.
  class ScreenLogStream : public LogStream {
  public:
    ScreenLogStream(int level);
    ScreenLogStream(const std::vector<int> & levels, bool ignore_general = false);
    ~ScreenLogStream();

    void LogMessage(int level, const std::string & message);
    void LogMessage(int level, int phase, const std::string & message);
    void UpdateProgress(double progress, const std::string & message);

  private:
    bool is_writing_progress_;
    static const int n_progress_hats_ = 50;

    void WriteProgressHeader();
  };
}

#endif // NRLIB_LOGSTREAM_HPP
