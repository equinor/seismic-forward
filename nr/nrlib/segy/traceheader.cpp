// $Id: traceheader.cpp 1279 2014-06-03 14:47:05Z hauge $

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

#include "traceheader.hpp"

#include <cstdio>
#include <cstring>
#include <sstream>
#include <iostream>

#include "../iotools/fileio.hpp"
#include "../iotools/stringtools.hpp"
#include "../iotools/logkit.hpp"
const float RMISSING = -99999.000;
const int   IMISSING = -99999;

using namespace std;
using namespace NRLib;

enum standardLoc {
  SCALCO_LOC    = 71,
  SX_LOC        = 73,
  SY_LOC        = 77,
  NS_LOC        = 115,
  DT_LOC        = 117,
  INLINE_LOC    = 9, // 189,
  CROSSLINE_LOC = 21, // 193
  DELAY_REC_TIME = 109
};

TraceHeaderFormat::TraceHeaderFormat(int headerformat)
{
  Init(headerformat);
}

TraceHeaderFormat::TraceHeaderFormat(int headerformat,
                                     int bypassCoordScaling,
                                     int scalCoLoc,
                                     int utmxLoc,
                                     int utmyLoc,
                                     int inlineLoc,
                                     int crosslineLoc,
                                     int coordSys)
{
  Init(headerformat);
  //
  // Redefined format if parameters have been given values
  //
  if (scalCoLoc != IMISSING)
  {
    scal_co_loc_ = scalCoLoc;
    standard_type_ = false;
  }
  if (utmxLoc != IMISSING)
  {
    utmx_loc_ = utmxLoc;
    standard_type_  = false;
  }
  if (utmyLoc != IMISSING)
  {
    utmy_loc_ = utmyLoc;
  standard_type_ = false;
  }
  if (inlineLoc != IMISSING)
  {
    inline_loc_ = inlineLoc;
    standard_type_ = false;
  }
  if (crosslineLoc != IMISSING)
  {
    crossline_loc_ = crosslineLoc;
    standard_type_ = false;
  }
  if (coordSys != IMISSING)
  {
    coord_sys_ = static_cast<coordSys_t>(coordSys);
    standard_type_ = false;
  }
  if (bypassCoordScaling == 1)
  {
    scal_co_loc_ = -1;
  }
  CheckFormat();
}

TraceHeaderFormat::TraceHeaderFormat(int scalCoLoc,
                                     int utmxLoc,
                                     int utmyLoc,
                                     int inlineLoc,
                                     int crosslineLoc,
                                     coordSys_t coordSys)
  : format_name_("unnamed"),
    scal_co_loc_(scalCoLoc),
    utmx_loc_(utmxLoc),
    utmy_loc_(utmyLoc),
    inline_loc_(inlineLoc),
    crossline_loc_(crosslineLoc),
    coord_sys_(coordSys),
    delay_rec_time_loc_(),
    offset_loc_(),
    standard_type_(true)
{}

TraceHeaderFormat::TraceHeaderFormat()
{
  Init(0);
}

TraceHeaderFormat::TraceHeaderFormat(const TraceHeaderFormat & thf)
 : format_name_  (thf.GetFormatName()),
   scal_co_loc_   (thf.GetScalCoLoc()),
   utmx_loc_     (thf.GetUtmxLoc()),
   utmy_loc_     (thf.GetUtmyLoc()),
   inline_loc_   (thf.GetInlineLoc()),
   crossline_loc_(thf.GetCrosslineLoc()),
   coord_sys_    (thf.GetCoordSys()),
   delay_rec_time_loc_(thf.GetDelayRecTimeLoc()),
   offset_loc_   (thf.GetOffsetLoc()),
   standard_type_(thf.GetStandardType())
{
  CheckFormat();
}

void
TraceHeaderFormat::Init(int headerformat)
{
  standard_type_ = true;
  if (headerformat==SEISWORKS)
  {
    format_name_   = std::string("SeisWorks");
    scal_co_loc_   = SCALCO_LOC;
    utmx_loc_      = SX_LOC;
    utmy_loc_      = SY_LOC;
    inline_loc_    = INLINE_LOC;
    crossline_loc_ = CROSSLINE_LOC;
    coord_sys_     = UTM;
    delay_rec_time_loc_= IMISSING;
    offset_loc_    = IMISSING;
  }
  else if (headerformat==IESX)
  {
    format_name_   = std::string("IESX");
    scal_co_loc_   = SCALCO_LOC;
    utmx_loc_      = SX_LOC;
    utmy_loc_      = SY_LOC;
    inline_loc_    = 221;
    crossline_loc_ = CROSSLINE_LOC;
    coord_sys_     = UTM;
    delay_rec_time_loc_= IMISSING;
    offset_loc_    = IMISSING;
  }
  else if (headerformat==SIP)
  {
    format_name_   = std::string("SIP");
    scal_co_loc_   = SCALCO_LOC;
    utmx_loc_      = 181;
    utmy_loc_      = 185;
    inline_loc_    = 189;
    crossline_loc_ = 193;
    coord_sys_     = UTM;
    delay_rec_time_loc_= DELAY_REC_TIME;
    offset_loc_    = 37;
  }
  else if (headerformat == CHARISMA)
  {
    format_name_   = std::string("Charisma");
    scal_co_loc_   = SCALCO_LOC;
    utmx_loc_      = SX_LOC;
    utmy_loc_      = SY_LOC;
    inline_loc_    = 5;
    crossline_loc_ = CROSSLINE_LOC;
    coord_sys_     = UTM;
    delay_rec_time_loc_= IMISSING;
    offset_loc_    = IMISSING;
  }
  else if (headerformat == SIPX) // Sebn: SIP probably messed up when they made volumes with this header specification.
  {
    format_name_   = std::string("SIPX");
    scal_co_loc_   = SCALCO_LOC;
    utmx_loc_      = SX_LOC;
    utmy_loc_      = SY_LOC;
    inline_loc_    = 181;
    crossline_loc_ = 185;
    coord_sys_     = UTM;
    delay_rec_time_loc_= IMISSING;
    offset_loc_    = IMISSING;
  }
  else
  {
    std::string error;
    std::stringstream format;
    error += "\n\nERROR: Undefined trace header format encountered. The recognized";
    error += "\nformat names and their associated trace header locations are:\n\n";
    error += "Name             X     Y      IL    XL  CoorScal   CoorSys\n";
    error += "----------------------------------------------------------\n";
    format << "SeisWorks   "
           << std::right
           << std::setw(6)  << SX_LOC
           << std::setw(6)  << SY_LOC
           << std::setw(8)  << INLINE_LOC
           << std::setw(6)  << CROSSLINE_LOC
           << std::setw(10) << SCALCO_LOC
           << std::setw(10) << "UTM";
    error += format.str();
    error += "\n";
    format.str("");
    format << "Charisma    "
           << std::right
           << std::setw(6)  << SX_LOC
           << std::setw(6)  << SY_LOC
           << std::setw(8)  << 5
           << std::setw(6)  << CROSSLINE_LOC
           << std::setw(10) << SCALCO_LOC
           << std::setw(10) << "UTM";
    error += format.str();
    error += "\n";
    format.str("");
    format << "IESX        "
           << std::right
           << std::setw(6)  << SX_LOC
           << std::setw(6)  << SY_LOC
           << std::setw(8)  << 221
           << std::setw(6)  << CROSSLINE_LOC
           << std::setw(10) << SCALCO_LOC
           << std::setw(10) << "UTM";
    error += format.str();
    error += "\n";
    format.str("");
    format << "SIP         "
           << std::right
           << std::setw(6)  << 181
           << std::setw(6)  << 185
           << std::setw(8)  << 189
           << std::setw(6)  << 193
           << std::setw(10) << SCALCO_LOC
           << std::setw(10) << "UTM";
    error += format.str();
    error += "\n";
    format.str("");
    format << "SIPX        "
           << std::right
           << std::setw(6)  << SX_LOC
           << std::setw(6)  << SY_LOC
           << std::setw(8)  << 181
           << std::setw(6)  << 185
           << std::setw(10) << SCALCO_LOC
           << std::setw(10) << "UTM";
    error += format.str();
    error += "\n";
    format.str("");
    throw Exception(error);
  }
}


void
TraceHeaderFormat::CheckFormat()
{
  if(scal_co_loc_ > 0) {
    if(scal_co_loc_ == utmx_loc_) {
      std::string message = "Both Scaling Coefficient and UTMY in SegY format set to " + ToString(scal_co_loc_) + ".\n";
      throw(Exception(message));
    }
    if(scal_co_loc_ == utmy_loc_) {
      std::string message = "Both Scaling Coefficient and UTMY in SegY format set to " + ToString(scal_co_loc_) + ".\n";
      throw(Exception(message));
    }
    if(scal_co_loc_ == inline_loc_) {
      std::string message = "Both Scaling Coefficient and IL in SegY format set to " + ToString(scal_co_loc_) + ".\n";
      throw(Exception(message));
    }
    if(scal_co_loc_ == crossline_loc_) {
      std::string message = "Both Scaling Coefficient and XL in SegY format set to " + ToString(scal_co_loc_) + ".\n";
      throw(Exception(message));
    }
  }
  if(utmx_loc_ > 0) {
    if(utmx_loc_ == utmy_loc_) {
      std::string message = "Both UTMX and UTMY in SegY format set to " + ToString(utmx_loc_) + ".\n";
      throw(Exception(message));
    }
    if(utmx_loc_ == inline_loc_) {
      std::string message = "Both UTMX and IL in SegY format set to " + ToString(utmx_loc_) + ".\n";
      throw(Exception(message));
    }
    if(utmx_loc_ == crossline_loc_) {
      std::string message = "Both UTMX and XL in SegY format set to " + ToString(utmx_loc_) + ".\n";
      throw(Exception(message));
    }
  }
  if(utmy_loc_ > 0) {
    if(utmy_loc_ == inline_loc_) {
      std::string message = "Both UTMY and IL in SegY format set to " + ToString(utmy_loc_) + ".\n";
      throw(Exception(message));
    }
    if(utmy_loc_ == crossline_loc_) {
      std::string message = "Both UTMY and XL in SegY format set to " + ToString(utmy_loc_) + ".\n";
      throw(Exception(message));
    }
  }
  if(inline_loc_ > 0 && inline_loc_ == crossline_loc_) {
    std::string message = "Both IL and XL in SegY format set to " + ToString(inline_loc_) + ".\n";
    throw(Exception(message));
  }
}

std::string
TraceHeaderFormat::toString() const
{
  std::string output;
  std::string coordsys = (coord_sys_ == UTM ? "UTM" : "IL/XL");
  output = "Coord used: " + coordsys + "Location in trace header:" +
           " UTM-X: " + NRLib::ToString(utmx_loc_) +
           " UTM-Y: " + NRLib::ToString(utmy_loc_) +
           " IL: " + NRLib::ToString(inline_loc_) +
           " XL: " + NRLib::ToString(crossline_loc_);
  return output;
}

std::vector<TraceHeaderFormat*>
TraceHeaderFormat::GetListOfStandardHeaders()
{
  int n = numberOfFormats_*2;
  std::vector<TraceHeaderFormat*> thf(n);
  int j = 0;
  for (int i = 0 ; i < n ; i += 2)
  {
    // With coordinate scaling
    thf[i+0] = new TraceHeaderFormat(j);
    // Without coordinate scaling
    thf[i+1] = new TraceHeaderFormat(j, 1, IMISSING, IMISSING, IMISSING, IMISSING, IMISSING, IMISSING);
    j++;
  }
  return thf;
}


int
TraceHeaderFormat::IsDifferent(TraceHeaderFormat inFormat)
{
  int ok = 0;
  if ((scal_co_loc_ ==-1 && inFormat.GetScalCoLoc() != -1)
       || (scal_co_loc_ != -1 && inFormat.GetScalCoLoc() == -1))
    ok = 1;
  else if ( scal_co_loc_ != inFormat.GetScalCoLoc())
    ok = 2;
  if (utmx_loc_ != inFormat.GetUtmxLoc())
    ok = 2;
  if (utmy_loc_ != inFormat.GetUtmyLoc())
    ok = 2;
  if (inline_loc_ != inFormat.GetInlineLoc())
    ok = 2;
  if ( crossline_loc_ != inFormat.GetCrosslineLoc())
    ok = 2;

  return ok;

}

void TraceHeaderFormat::WriteValues() const
{
  LogKit::LogFormatted(LogKit::Medium,"This traceheader format has the following values:\n");
  LogKit::LogFormatted(LogKit::Medium," utmxLoc utmyLoc inlineLoc crosslineLoc scalcoLoc \n");
  LogKit::LogFormatted(LogKit::Medium,"--------------------------------------------------\n");
  LogKit::LogFormatted(LogKit::Medium,"%5d  %5d   %5d         %5d     %5d    \n",
                                       utmx_loc_, utmy_loc_, inline_loc_, crossline_loc_, scal_co_loc_);
}


TraceHeader::TraceHeader(const TraceHeaderFormat& format)
  : format_(format),
    status_(0),
    scal_co_(1),
    scalcoinitial_(1)
{
  memset(buffer_, 0, 240);
  rmissing_ = RMISSING;
  imissing_ = IMISSING;
}


void TraceHeader::Read(std::istream& inFile, int lineNo)
{
  if (!(inFile.read(buffer_,240))) {
    // end of file
    throw EndOfFile();
  }

  if (buffer_[0] == 'Ã' && buffer_[1] == '@' && buffer_[2] == 'ñ'
      && buffer_[80] == 'Ã' && buffer_[160] == 'Ã')
  {
    // This is not a trace header, but the start of an EDBDIC-header.
    // Set file pointer at end of EDBDIC header.
   // fseek(inFile, 2960, SEEK_CUR);
    inFile.seekg(2960);
    status_ = -1;
    return;
  }

  std::string buf_string(buffer_,240);
  std::istringstream header(buf_string, std::ios::in | std::ios::binary);

  int i = 0;
  while (i < 240)
  {
    if (i==(format_.GetScalCoLoc()-1))
    {
      scalcoinitial_ = ReadBinaryShort(header);
      i=i+2;
    }
    else if (i==(format_.GetUtmxLoc()-1))
    {
      utmx_ = float(ReadBinaryInt(header));
      i=i+4;
    }
    else if (i==(format_.GetUtmyLoc()-1))
    {
      utmy_ = float(ReadBinaryInt(header));
      i=i+4;
    }
    else if (i==(NS_LOC-1))
    {
      ns_ = ReadBinaryShort(header);
      i=i+2;
    }
    else if (i==(DT_LOC-1))
    {
      dt_ = ReadBinaryShort(header);
      i=i+2;
    }
    else if (i==(format_.GetInlineLoc()-1))
    {
      inline_ = ReadBinaryInt(header);
      i=i+4;
    }
    else if (i==(format_.GetCrosslineLoc()-1))
    {
      crossline_ = ReadBinaryInt(header);
      i=i+4;
    }
    else
    {
      ReadBinaryShort(header);
      i=i+2;
    }

  }

 // swapBuffer();

  if (format_.GetScalCoLoc() > 0) {
    switch (scalcoinitial_) {
      case -10000:
        scal_co_ = 0.0001f;
        break;
      case -1000:
        scal_co_ = 0.001f;
        break;
      case -100:
        scal_co_ = 0.01f;
        break;
      case -10:
        scal_co_ = 0.1f;
        break;
      case 10:
      case 100:
      case 1000:
      case 10000:
        scal_co_ = static_cast<float>(scalcoinitial_);
    }
  }
useBinaryInline = false;
  // if (format.useInlineFromBinHeader()) {
  if (lineNo > 0 && GetInline() == 0 && useBinaryInline==true) {
    SetInline(lineNo);
  }
}

int TraceHeader::Write(std::ostream& outFile)
{
  int errCode = 0;

 // write on correct locations. What to write between?
  int i = 0;

  while(i<240)
  {
    if (i==(format_.GetScalCoLoc()-1))
    {
      WriteBinaryShort(outFile, scalcoinitial_);
      i = i+2;
    }
    else if (i==(format_.GetUtmxLoc()-1))
    {
      WriteBinaryInt(outFile, static_cast<int>(utmx_));
      i=i+4;
    }
    else if (i==(format_.GetUtmyLoc()-1))
    {
      WriteBinaryInt(outFile, static_cast<int>(utmy_));
      i=i+4;
    }
    else if (i==(NS_LOC-1))
    {
      WriteBinaryShort(outFile, ns_);
      i=i+2;
    }
    else if (i==(DT_LOC-1))
    {
      WriteBinaryShort(outFile, dt_);
      i=i+2;
    }
    else if (i==(format_.GetInlineLoc()-1))
    {
      WriteBinaryInt(outFile, inline_);
      i=i+4;
    }
    else if (i==(format_.GetCrosslineLoc()-1))
    {
      WriteBinaryInt(outFile, crossline_);
      i=i+4;
    }
    else if (i==(format_.GetDelayRecTimeLoc()-1))
    {
      WriteBinaryShort(outFile, delay_rec_time_);
      i=i+2;
    }
     else if (i==(format_.GetOffsetLoc()-1))
    {
      WriteBinaryInt(outFile, static_cast<int>(offset_));
      i=i+4;
    }
   else
    {
      outFile.write(&(buffer_[i]), 2);
      i=i+2;
    }
  }

  return errCode;
}


void
TraceHeader::Dump(std::ostream& outFile, bool changeNs)
{
  if(changeNs == false)
    outFile.write(buffer_, 240);
  else {
    outFile.write(buffer_, NS_LOC-1);
    WriteBinaryShort(outFile, ns_);
    outFile.write(&(buffer_[NS_LOC+1]),239-NS_LOC);
  }
}


void TraceHeader::WriteValues()
{
  float dtms = static_cast<float>(dt_)/1000.0f;
  float lms  = static_cast<float>(ns_-1)*dtms;
  if (format_.GetScalCoLoc() == -1) {
    LogKit::LogFormatted(LogKit::High,"\n\nThe following header information was extracted from the first trace:\n\n");
    LogKit::LogFormatted(LogKit::High,"     UTMx         UTMy        IL    XL      Samples   dt(ms)  Length(ms)\n");
    LogKit::LogFormatted(LogKit::High,"-------------------------------------------------------------------------\n");
    LogKit::LogFormatted(LogKit::High,"%9.2f  %11.2f     %5d %5d       %6d     %4.2f     %7.2f\n",
                         utmx_, utmy_, inline_, crossline_, ns_, dtms, lms);
  }
  else {
    LogKit::LogFormatted(LogKit::High,"\n\nThe following header information was extracted from the first trace:\n\n");
    LogKit::LogFormatted(LogKit::High,"     UTMx         UTMy     CoScal        IL    XL      Samples   dt(ms)  Length(ms)\n");
    LogKit::LogFormatted(LogKit::High,"------------------------------------------------------------------------------------\n");
    LogKit::LogFormatted(LogKit::High,"%9.2f  %11.2f     %6.1f     %5d %5d       %6d     %4.2f     %7.2f\n",
                         utmx_, utmy_, scal_co_, inline_, crossline_, ns_, dtms, lms);
  }
}

double TraceHeader::GetUtmx() const
{
  int loc = format_.GetUtmxLoc();
  if (loc < 0) {
    return rmissing_;
  }
  return(scal_co_*utmx_);
}

void TraceHeader::SetUtmx(double utmx)
{
  int loc = format_.GetUtmxLoc();
  if (loc > 0) {
    utmx_ = utmx/scal_co_;
  }
}

double TraceHeader::GetUtmy() const
{
  int loc = format_.GetUtmyLoc();
  if (loc < 0) {
    return rmissing_;
  }
  return(scal_co_*utmy_);
}

void TraceHeader::SetUtmy(double utmy)
{
  int loc = format_.GetUtmyLoc();
  if (loc > 0) {
    utmy_ = utmy/scal_co_;
  }
}


int TraceHeader::GetInline() const
{
  int loc = format_.GetInlineLoc();
  if (loc < 0) {
    return imissing_;
  }
  return inline_;
}


void TraceHeader::SetInline(int inLine)
{
  int loc = format_.GetInlineLoc();
  if (loc > 0) {
    inline_ = inLine;
  }
}


int TraceHeader::GetCrossline() const
{
  int loc = format_.GetCrosslineLoc();
  if (loc < 0) {
    return imissing_;
  }
  return crossline_;
}


void TraceHeader::SetCrossline(int crossLine)
{
  int loc = format_.GetCrosslineLoc();
  if (loc > 0) {
    crossline_ = crossLine;
  }
}


double TraceHeader::GetCoord1() const
{
  if(format_.GetCoordSys() == TraceHeaderFormat::UTM)
    return(GetUtmx());
  else
    return(static_cast<double>(GetInline()));
}


double TraceHeader::GetCoord2() const
{
  if(format_.GetCoordSys() == TraceHeaderFormat::UTM)
    return(GetUtmy());
  else
    return(static_cast<double>(GetCrossline()));
}


void TraceHeader::SetNSamples(size_t ns)
{
  ns_ = static_cast<short>(ns);
}


void TraceHeader::SetDt(int dt)
{
  dt_ = static_cast<short>(dt);
}


short TraceHeader::GetScalCo() const
{
  int loc = format_.GetScalCoLoc();
  if (loc < 0) {
    return 0;
  }
  return scalcoinitial_;
}

void TraceHeader::SetScalCo(short scalcoinitial)
{
  scalcoinitial_ = scalcoinitial;
  switch (scalcoinitial) {
  case -10000:
    scal_co_ = 0.0001f;
    break;
  case -1000:
    scal_co_ = 0.001f;
    break;
  case -100:
    scal_co_ = 0.01f;
    break;
  case -10:
    scal_co_ = 0.1f;
    break;
  case 10:
  case 100:
  case 1000:
  case 10000:
    scal_co_ = static_cast<float>(scalcoinitial);
  }


}

void TraceHeader::SetDelayRecTime(short value)
{
  int loc = format_.GetDelayRecTimeLoc();
  if (loc > 0) {
    delay_rec_time_ = value;
  }
}

short TraceHeader::GetDelayRecTime() const
{
  int loc = format_.GetDelayRecTimeLoc();
  if (loc < 0) {
    return 0;
  }
  return delay_rec_time_;
}

