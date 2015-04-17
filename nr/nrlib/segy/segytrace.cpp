// $Id: segytrace.cpp 1199 2013-10-02 08:24:02Z anner $

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

#include <cmath>
#include <ctime>
#include <cstdio>
#include <cassert>
#include <iostream>

#include "segy.hpp"
#include "commonheaders.hpp"
#include "traceheader.hpp"

#include "../iotools/logkit.hpp"
#include "../exception/exception.hpp"
#include "../iotools/fileio.hpp"
#include "../surface/surface.hpp"
#include "../iotools/stringtools.hpp"

const float segyRMISSING = -99999.000;
const int   segyIMISSING = -99999;

using namespace NRLib;


SegYTrace::SegYTrace(std::fstream & file, size_t jStart, size_t jEnd, int format,
                     size_t nz, const TraceHeader * trace_header)
{
  rmissing_      = segyRMISSING;
  imissing_      = segyIMISSING;
  j_start_       = jStart;
  j_end_         = jEnd;
  x_             = trace_header->GetUtmx();
  y_             = trace_header->GetUtmy();
  in_line_       = trace_header->GetInline();
  cross_line_    = trace_header->GetCrossline();
  coord1_        = trace_header->GetCoord1();
  coord2_        = trace_header->GetCoord2();
  trace_header_  = new TraceHeader(*trace_header);
  table_index_   = 0;
  file_position_ = 0;

  size_t nData = jEnd - jStart + 1;
  size_t i;
  std::vector<float> predata;
  predata.resize(nz);
  data_.resize(nData);

  try {
    if (format == 1)
    {
      //IBM
      ReadBinaryIbmFloatArray(file, predata.begin(), nz);
      for (i = 0; i < nData; i++)
        data_[i] = predata[jStart+i];
    }
    else if (format == 2)
    {
      std::vector<int> b(nz);
      ReadBinaryIntArray(file, b.begin(), nz);
      for (i = 0; i < nData; i++)
        data_[i] =  static_cast<float> (b[jStart+i]);
    }
    else if (format == 3)
    {
      std::vector<short> b(nz);
      ReadBinaryShortArray(file, b.begin(), nz);
      for (i = 0; i < nData; i++)
        data_[i] =  static_cast<float> (b[jStart+i]);
    }
    else if (format == 5)
    {
      ReadBinaryFloatArray(file, predata.begin(), nz);
      for (i = 0; i < nData; i++)
        data_[i] = predata[jStart + i];
    }
    else
      throw FileFormatError("Bad format");
  }
  catch (NRLib::Exception & e)
  {
    std::string text = e.what();
    text += "\nTrace header information for failing trace:";
    text += "\n  x         = " + ToString(x_,2);
    text += "\n  y         = " + ToString(y_,2);
    text += "\n  IL        = " + ToString(in_line_);
    text += "\n  XL        = " + ToString(cross_line_);
    text += "\n  coord1    = " + ToString(coord1_,2);
    text += "\n  coord2    = " + ToString(coord2_,2);
    text += "\n  nz        = " + ToString(nz) + "\n";
    text += "\nInformation about values read for last trace (float values may be unset):";
    text += "\n  buffer length                        = " + ToString(nz);
    text += "\n  inversion interval start index       = " + ToString(j_start_);
    text += "\n  inversion interval end index         = " + ToString(j_end_);
    text += "\n  number of data in inversion interval = " + ToString(nData);
    // This is an ugly way to get hold of the numbers, but the only(?) way
    // due to the exception throwing in read routines.
    if (e.what() != std::string("Bad format")) {
      text += "\n  trace values:\n";
      if (format == 1 || format == 5) {
        for (i = 0; i < nz; i++)
          text += std::string("    ") + ToString(i) + ":  " + ToString(predata[i]) + "\n";
      }
    }
    throw Exception(text);
  }
}

SegYTrace::SegYTrace(std::vector<float> indata, size_t jStart, size_t jEnd, double x, double y, int inLine, int crossLine)
{
  rmissing_   = segyRMISSING;
  imissing_   = segyIMISSING;
  j_start_    = jStart;
  j_end_      = jEnd;
  x_          = x;
  y_          = y;
  in_line_    = inLine;
  cross_line_ = crossLine;
  coord1_     = x;
  coord2_     = y;
  size_t nData = jEnd - jStart + 1;
  size_t i;
  data_.resize(nData);
  for (i = 0; i < nData; i++)
    data_[i] = indata[jStart + i];

  table_index_   = 0;
  file_position_ = 0;
  trace_header_  = NULL;
}

SegYTrace::SegYTrace(const TraceHeader& trace_header, bool keep_header)
{
  rmissing_      = segyRMISSING;
  imissing_      = segyIMISSING;
  j_start_       = 1;
  j_end_         = 0;
  x_             = trace_header.GetUtmx();
  y_             = trace_header.GetUtmy();
  in_line_       = trace_header.GetInline();
  cross_line_    = trace_header.GetCrossline();
  coord1_        = trace_header.GetCoord1();
  coord2_        = trace_header.GetCoord2();
  table_index_   = 0;
  file_position_ = 0;

  if(keep_header == true)
    trace_header_ = new TraceHeader(trace_header);
  else
    trace_header_ = NULL;

}

SegYTrace::~SegYTrace()
{
  delete trace_header_;
}

float
SegYTrace::GetValue(size_t j)const
{
  float value;
  if (j < j_start_ || j > j_end_)
    value = rmissing_;
  else
    value = data_[j - j_start_];
  return(value);
}

size_t
SegYTrace::GetLegalIndex(size_t index) const
{
  if (index < j_start_)
    return(j_start_);
  else if (index > j_end_)
    return(j_end_);
  else
    return(index);
}

const TraceHeader&
SegYTrace::GetTraceHeader() const
{
  if(trace_header_ != NULL)
    return (*trace_header_);
  else
    throw Exception("This trace has no header read from file.\n");
}
