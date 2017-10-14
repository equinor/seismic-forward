// $Id: traceheader.hpp 1689 2017-09-13 09:37:00Z perroe $

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

#ifndef TRACEHEADER_HPP
#define TRACEHEADER_HPP

#include <string>
#include <vector>

namespace NRLib {
class BigFile;


///  The format of the trace header. Specifies the location
///  of the fields of interest. All locations are 1-based, i.e.
///  the location of the first field in the header is 1.
///  The locations are set to -1 if the fields are not given in the header.
class TraceHeaderFormat {
public:

  ///   Possible coordinate systems.
  enum coordSys_t {
    UTM  = 0,   ///< UTM coordinates.
    ILXL = 1    ///< Inline/Crossline coordinates.
  };
  enum headers {
    SEISWORKS = 0,
    IESX      = 1,
    SIP       = 2,
    CHARISMA  = 3,
    SIPX      = 4, // Sebn: SIP probably messed up when they made volumes with this header specification.
    HESS      = 5,
    numberOfFormats_
  };

  ///   Constructor, the standard format.
  TraceHeaderFormat(int headerformat);


  ///   Constructor, the standard format with modification.
  TraceHeaderFormat(int headerformat,
                    int byPassCoordScaling,
                    int scaleCoLoc,
                    int utmxLoc,
                    int utmyLoc,
                    int inlineLoc,
                    int crosslineLoc,
                    int startTimeLoc,
                    int coordSys);


  ///  Constructor. Locations should be set to -1 if not set.
  ///  \param[in] utmxLoc      Location scaling coefficient for UTM X and Y.
  ///  \param[in] utmxLoc      Location of the UTM-X field.
  ///  \param[in] utmyLoc      Location of the UTM-Y field.
  ///  \param[in] inlineLoc    Location of the inline coordinate field.
  ///  \param[in] crosslineLoc Location of crossline coordinate field.

  TraceHeaderFormat(int scaleCoLoc,
                    int utmxLoc,
                    int utmyLoc,
                    int inlineLoc,
                    int crosslineLoc,
                    int startTimeLoc,
                    coordSys_t coordSys);

  TraceHeaderFormat();

  TraceHeaderFormat(const TraceHeaderFormat& thf);

  /// Get format type
  std::string GetFormatName() const {return format_name_;}

  /// Get location of the UTM-X field. (-1 if non-existant)
  int GetUtmxLoc() const {return utmx_loc_;}

  /// Get location of the UTM-Y field. (-1 if non-existant)
  int GetUtmyLoc() const {return utmy_loc_;}

  /// Get location of the inline field. (-1 if non-existant)
  int GetInlineLoc() const {return inline_loc_;}

  /// Get location of the crossline field. (-1 if non-existant)
  int GetCrosslineLoc() const {return crossline_loc_;}

  /// Get location of the start time. (-1 if non-existant)
  int GetStartTimeLoc() const {return start_time_loc_;}

  /// Get location of the scaling cooefficient field. (-1 if non-existant)
  int GetScalCoLoc() const {return scal_co_loc_;}

  /// Get bypassCoodScaling status
  bool GetBypassCoordScaling() const {return scal_co_loc_ == -1;}

  /// Get offset location
  int GetOffsetLoc() const { return offset_loc_;}

  /// Get coordinate system.
  coordSys_t GetCoordSys() const {return coord_sys_;}

  /// Get is this a standard type
  bool GetStandardType() const {return standard_type_;}

  void SetScaleCoLoc(int loc)       {scal_co_loc_    = loc; standard_type_ = false;}
  void SetUtmxLoc(int loc)          {utmx_loc_       = loc; standard_type_ = false;}
  void SetUtmyLoc(int loc)          {utmy_loc_       = loc; standard_type_ = false;}
  void SetInlineLoc(int loc)        {inline_loc_     = loc; standard_type_ = false;}
  void SetCrosslineLoc(int loc)     {crossline_loc_  = loc; standard_type_ = false;}
  void SetStartTimeLoc(int loc)     {start_time_loc_ = loc; standard_type_ = false;}
  void SetCoordSys(coordSys_t type) {coord_sys_      = type;}

  /// Check that no two values point to the same byte. Throws if error.
  void CheckFormat();

  /// String representation.
  std::string toString() const;
  /// Return a list of available trace header formats.
  static std::vector<TraceHeaderFormat*> GetListOfStandardHeaders();

  /// Check if format equals input format. Returns 0 if equal, 1 if equal
  /// except bypasscoordiantescaling, 2 otherwise
  int IsDifferent(TraceHeaderFormat inFormat);

  void WriteValues() const;

private:
  /// Set default values for standard formats
  void Init(int headerformat);
  /// Format name
  std::string format_name_;
  /// Location of scaling coefficient. (-1 if non-existant)
  int scal_co_loc_;
  /// Location of the UTM-X field. (-1 if non-existant)
  int utmx_loc_;
  /// Location of the UTM-Y field. (-1 if non-existant)
  int utmy_loc_;
  /// Location of inline coordinate field. (-1 if non-existant)
  int inline_loc_;
  /// Location of crossline coordinate field. (-1 if non-existant)
  int crossline_loc_;
  /// Location of start time. (-1 if non-existant)
  int start_time_loc_;
  /// Coordinate system to use.
  coordSys_t coord_sys_;
  /// Location of offset
  int offset_loc_;
  /// Standard type
  bool standard_type_;

};


/// The trace header for a SEGY file.
class TraceHeader {
public:
  /// Constructor generating an empty header.
  /// \param[in] format  header format.
  TraceHeader(const TraceHeaderFormat& format = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS));

  /// Read in a new header.
  /// \param[in] in_file  input file.
  void Read(NRLib::BigFile& in_file);

  /// Write header to file.
  /// \param[in]  outFile output file.
  int Write(NRLib::BigFile& outFile);

  /// Dump what is stored in header buffer to file.
  /// May override the number of data. Intended for use when we copy headers
  /// from input to output.
  /// \param[in]  outFile output file.
  void Dump(NRLib::BigFile& outFile, bool changeNs = false);

  void WriteValues();

  /// Get UTM X coordinate.
  /// returns #RMISSING if UTM X location is not set in the format.
  double GetUtmx() const;

  /// Set UTM X coordinate.
  /// Does nothing if UTM X location is not set in the format.
  void SetUtmx(double utmx);

  /// Get UTM Y coordinate.
  /// returns #RMISSING if UTM Y location is not set in the format.
  double GetUtmy() const;

  /// Set UTM Y coordinate.
  /// Does nothing if UTM Y location is not set in the format.
  void SetUtmy(double utmy);

  /// Get inline coordinate.
  /// returns #IMISSING if inline location is not set in the format.
  int GetInline() const;

  /// Set inline coordinate.
  /// Does nothing if inline location is not set in the format.
  void SetInline(int inLine);

  /// Get crossline coordinate.
  /// returns #IMISSING if crossline location is not set in the format.
  int GetCrossline() const;

  /// Set crossline coordinate.
  /// Does nothing if crossline location is not set in the format.
  void SetCrossline(int crossLine);

  /// Get start time.
  /// returns 0.0 if start time is not set in the format.
  float GetStartTime() const;

  /// Set start time.
  /// Does nothing if start time is not set in the format.
  void SetStartTime(float start_time);

  /// Get current first coordinate (either UTM x or IL depending on coordsys).
  /// returns #RMISSING if value is not set in the format.
  double GetCoord1() const;

  /// Get current second coordinate (either UTM y or XL depending on coordsys).
  /// returns #RMISSING if value is not set in the format.
  double GetCoord2() const;


  /// Set number of samples.
  void SetNSamples(size_t ns);

  /// Set sample interval in mikroseconds.
  void SetDt(int dt);
  short GetDt() {return(dt_);}

  void SetScalCo(short scalcoinitial);

  void SetOffset(double offset) { offset_ = offset; };

  /// Get status code.
  ///  0 - everything went OK.
  /// -1 - not a trace header, but EDBDIC header.
  /// -2 - error reading from file.
  int GetStatus() const {return status_;}

  /// Returns false if header values are illegal.
  bool IsHeaderOK() const;

private:
  /// Header buffer in machine-specific byte order.
  char buffer_[240];

  /// Header format.
  TraceHeaderFormat format_;

  /// Status code.
  int status_;

  /// Scaling coefficient for UTM X and UTM Y.
  double scal_co_;

  short scalcoinitial_;
  double utmx_;
  double utmy_;
  int inline_;
  int crossline_;
  short start_time_;
  short ns_;
  short dt_;
  int imissing_;
  float rmissing_;
  double offset_;


  bool useBinaryInline;

  /// Get scaling coefficient for SX and SY from buffer.
  short GetScalCo() const;

};

} // namespace NRLib

#endif
