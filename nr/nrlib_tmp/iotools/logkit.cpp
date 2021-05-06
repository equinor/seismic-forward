// $Id: logkit.cpp 1429 2017-02-15 11:50:06Z perroe $

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

#include <algorithm>
#include <cstdarg>
#include <iostream>

#include "fileio.hpp"
#include "logkit.hpp"
#include "logstream.hpp"
#include "../exception/exception.hpp"

using namespace NRLib;

std::vector<LogStream*> LogKit::logstreams_(0);
int LogKit::screenLog_ = -1;
std::vector<BufferMessage *> * LogKit::buffer_ = NULL;

// Making a long table to allow direct access.
std::vector<int> LogKit::n_messages_(65, 0);
std::vector<std::string> LogKit::prefix_(65, "");


void
LogKit::SetFileLog(const std::string & fileName, int levels,
                   bool includeNRLibLogging)
{
  LogStream * curStream;
  if (includeNRLibLogging == true)
    curStream = new FileLogStream(fileName, levels);
  else {
    std::vector<int> phaseLevels;
    phaseLevels.push_back(0); //No logging in NRLib, phase 0.
    phaseLevels.push_back(levels); //This will be used for all other phases.
    curStream = new FileLogStream(fileName, phaseLevels);
  }
  logstreams_.push_back(curStream);
  DumpBuffer(curStream);
}


void
LogKit::SetFileLog(const std::string & fileName,
                   const std::vector<int> & levels,
                   bool ignore_general)
{
  LogStream * curStream = new FileLogStream(fileName, levels, ignore_general);
  logstreams_.push_back(curStream);
  DumpBuffer(curStream);
}


void
LogKit::SetFileLog(const std::string & fileName, int levels, int phase, bool ignore_general) {
  std::vector<int> phaseLevels(1000,0);
  phaseLevels[phase] = levels;
  LogStream * curStream = new FileLogStream(fileName, phaseLevels, ignore_general);
  logstreams_.push_back(curStream);
  DumpBuffer(curStream);
}


void
LogKit::SetScreenLog(int levels, bool includeNRLibLogging)
{
  LogStream * curStream;
  if (includeNRLibLogging == true)
    curStream = new ScreenLogStream(levels);
  else {
    std::vector<int> phaseLevels;
    phaseLevels.push_back(0); //No logging in NRLib, phase 0.
    phaseLevels.push_back(levels); //This will be used for all other phases.
    curStream = new ScreenLogStream(phaseLevels);
  }
  if (screenLog_ < 0) {
    screenLog_ = static_cast<int>(logstreams_.size());
    logstreams_.push_back(curStream);
  }
  else {
    delete logstreams_[screenLog_];
    logstreams_[screenLog_] = curStream;
  }
}


void
LogKit::SetScreenLog(const std::vector<int> & levels, bool ignore_general) {
  LogStream * curStream = new ScreenLogStream(levels, ignore_general);
  if (screenLog_ < 0) {
    screenLog_ = static_cast<int>(logstreams_.size());
    logstreams_.push_back(curStream);
  }
  else {
    delete logstreams_[screenLog_];
    logstreams_[screenLog_] = curStream;
  }
}


void
LogKit::AddLogStream(LogStream* stream) {
  logstreams_.push_back(stream);
  DumpBuffer(stream);
}


void
LogKit::RemoveLogStream(LogStream* stream) {
  std::vector<LogStream*>::iterator it
    = std::find(logstreams_.begin(), logstreams_.end(), stream);
  if (it < logstreams_.end()) {
    if (static_cast<int>(it - logstreams_.begin()) < screenLog_)
      screenLog_ = screenLog_ - 1;
    delete *it;
    logstreams_.erase(it);
  }
}


void
LogKit::LogMessage(int level, const std::string & message) {
  unsigned int i;
  n_messages_[level]++;
  std::string new_message = prefix_[level] + message;
  for (i=0;i<logstreams_.size();i++)
    logstreams_[i]->LogMessage(level, new_message);
  SendToBuffer(level,-1,new_message);
}


void
LogKit::LogMessage(int level, int phase, const std::string & message) {
  unsigned int i;
  n_messages_[level]++;
  std::string new_message = prefix_[level] + message;
  for (i=0;i<logstreams_.size();i++)
    logstreams_[i]->LogMessage(level, phase, new_message);
  SendToBuffer(level,phase,new_message);
}


void
LogKit::LogFormatted(int level, const char* format, ...) {
  va_list ap;
  char message[5000];
  va_start(ap, format);
  vsprintf(message, format, ap);
  va_end(ap);
  LogMessage(level, std::string(message));
}


void
LogKit::LogFormatted(int level, int phase, const char* format, ...) {
  va_list ap;
  char message[1000];
  va_start(ap, format);
  vsprintf(message, format, ap);
  va_end(ap);
  LogMessage(level, phase, std::string(message));
}


void
LogKit::UpdateProgress(double progress, const std::string& message) {
  for (size_t i = 0; i < logstreams_.size(); i++)
    logstreams_[i]->UpdateProgress(progress, message);
}


void
LogKit::EndLog() {
  unsigned int i;
  for (i=0;i<logstreams_.size();i++)
    delete logstreams_[i];

  logstreams_.resize(0);
  screenLog_ = -1;

  if (buffer_ != NULL)
    EndBuffering(); //Also deletes buffer.
}


void
LogKit::StartBuffering() {
  buffer_ = new std::vector<BufferMessage *>;
}


void
LogKit::EndBuffering() {
  if (buffer_ != NULL) {
    for (unsigned int i=0;i<buffer_->size();i++) {
      delete (*buffer_)[i];
      (*buffer_)[i] = NULL;
    }
    delete buffer_;
    buffer_ = NULL;
  }
}


void
LogKit::SendToBuffer(int level, int phase, const std::string & message) {
  if (buffer_ != NULL) {
    BufferMessage * bm = new BufferMessage;
    bm->level_ = level;
    bm->phase_ = phase;
    bm->text_  = message;
    buffer_->push_back(bm);
  }
}


void
LogKit::DumpBuffer(LogStream *logstream) {
  if (buffer_ != NULL) {
    for (unsigned int i=0;i<buffer_->size();i++) {
      if ((*buffer_)[i]->phase_ < 0)
        logstream->LogMessage((*buffer_)[i]->level_, (*buffer_)[i]->text_);
      else
        logstream->LogMessage((*buffer_)[i]->level_, (*buffer_)[i]->phase_, (*buffer_)[i]->text_);
    }
  }
}


void
LogKit::SetPrefix(const std::string & prefix, int level) {
  prefix_[level] = prefix;
}


void
LogKit::WriteHeader(const std::string & text,
                    MessageLevels       logLevel)
{
  int width = 100; // Total width of header
  std::string ruler(width, '*');
  std::string stars("*****");
  LogMessage(logLevel,"\n"+ruler+"\n");
  int starLength  = int(stars.length());
  int textLength  = int(text.length());
  int blankLength = std::max(0, width - textLength - 2*starLength);
  std::string blanks(blankLength/2, ' ');
  std::string center;
  if(blankLength % 2)
    center = stars + blanks + text + blanks + " " + stars;
  else
    center = stars + blanks + text + blanks +stars;
  LogMessage(logLevel, center+"\n");
  LogMessage(logLevel, ruler+"\n");
}
