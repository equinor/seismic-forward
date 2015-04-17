// $Id: segy.hpp 1257 2014-03-06 13:55:54Z vigsnes $

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

#ifndef SEGY_HPP
#define SEGY_HPP

#include <fstream>
#include <string>
#include <vector>

#include "traceheader.hpp"
#include "commonheaders.hpp"
#include "../volume/volume.hpp"
#include "../segy/segygeometry.hpp"
#include "../segy/segytrace.hpp"

namespace NRLib {

const int segyIMISSING = -99999;
class SegYTrace;
class SegyGeometry;
class BinaryHeader;
class TextualHeader;


class SegY{
public:

  /// Constructor for reading
  /// Read only the headers on top of the file
  /// \param[in] fileName  Name of file to read data from
  /// \param[in] z0
  /// \throw IOError if the file can not be opened.
  SegY(const std::string       & fileName,
       float                     z0,
       const TraceHeaderFormat & traceHeaderFormat);


  /// Constructor for reading unknown format
  /// Read only the headers on top of the file
  /// \param[in] fileName  Name of file to read data from
  /// \param[in] z0
  /// \param[in] thf Vector of pointers to possible traceheaderformats. If NULL, default list is used.
  /// \throw IOError if the file can not be opened.
  /// \throw FileFormatError if the traceheaderformat can not be recognized.
  SegY(const std::string               & fileName,
       float                             z0,
       std::vector<TraceHeaderFormat *>  thf = std::vector<TraceHeaderFormat *>(0),
       bool                              searchStandardFormats = true);

  /// Constructor for writing
  /// \param[in] fileName  Name of file to write data to
  /// \throw IOError if the file can not be opened.
 SegY(const std::string       & fileName,
      float                     z0,
      size_t                    nz,
      float                     dz,
      const TextualHeader     & ebcdicHeader,
      const TraceHeaderFormat & traceHeaderFormat = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS));

  ~SegY();

  //>>>Begin read all traces mode
  void                      ReadAllTraces(const NRLib::Volume * volume,
                                          double                zPad,
                                          bool                  onlyVolume       = false,
                                          bool                  relative_padding = true); ///< Read all traces with header
  float                     GetValue(double x,
                                     double y,
                                     double z,
                                     int    outsideMode = segyIMISSING) const;

  size_t                    FindNumberOfSamplesInLongestTrace(void) const;

  void                      GetNearestTrace(std::vector<float> & trace_data,
                                            bool               & missing,
                                            float              & z0_data,
                                            float                x,
                                            float                y) const;

  std::vector<float>        GetAllValues();                           ///< Return vector with all values.

  void                      CreateRegularGrid();
  const SegyGeometry      * GetGeometry(void) const { return geometry_ ;} //Only makes sense after command above, or FindAndSetGeometry below.
  //<<<End read all trace mode

  //These functions are valid both after ReadAllTraces and FindAndSetGridGeometry
  const TraceHeader       & GetTraceHeader(int IL, int XL) const;
  const TraceHeader       & GetTraceHeader(float x, float y) const;

  //>>>Begin read single trace mode
  SegyGeometry            * FindGridGeometry(bool only_ilxl = false, bool keep_header = false);        //Note: This and function below also sets all traceheaders with file position.
  void                      FindAndSetGridGeometry(bool only_ilxl = false, bool keep_header = false);

  SegYTrace               * GetNextTrace(double                zPad = 0,
                                         const NRLib::Volume * volume = NULL,
                                         bool                  onlyVolume = false);
  void                      GetTraceData(int                  IL,
                                         int                  XL,
                                         std::vector<float> & result,
                                         const Volume       * volume = NULL);   // Only makes sense after FindAndSetGridGeometry.
  void                      GetTraceData(float                x,
                                         float                y,
                                         std::vector<float> & result,
                                         const Volume       * volume = NULL);   // Only makes sense after FindAndSetGridGeometry.
  void                      GetTraceData(std::streampos       pos,
                                         std::vector<float> & result,
                                         float                z_top = -1,       // Does not check pos, z_top or z_bot for validity.
                                         float                z_bot = -1);      // Top/bot = -1 means start from top/go to bottom
  std::streampos            GetFilePos(int IL, int XL) const;                   // Only makes sense after FindAndSetGridGeometry. Returns 0 for invalid trace.
  std::streampos            GetFilePos(float x, float y) const;                 // Only makes sense after FindAndSetGridGeometry. Returns 0 for invalid trace.
  //<<<End read single trace mode

  //>>>Begin write mode
  void                      SetGeometry(const SegyGeometry * geometry);

  void                      SetDelayRecTime(short value) { delay_rec_time_ = value;}

  void                      StoreTrace(double                     x,
                                       double                     y,
                                       const std::vector<float>   data,
                                       const NRLib::Volume      * volume,
                                       float                      topVal  = 0.0f,
                                       float                      baseVal = 0.0f);

  /// Write single trace to file
  void                      WriteTrace(const TraceHeader        & traceHeader,
                                       const std::vector<float> & data,
                                       const NRLib::Volume      * volume,
                                       float                      topVal  = 0.0f,
                                       float                      baseVal = 0.0f);
  /// Write single trace to internal memory

  void                      WriteTrace(double                   x,
                                       double                   y,
                                       const std::vector<float> data,
                                       const NRLib::Volume    * volume,
                                       float                    topVal  = 0.0f,
                                       float                    baseVal = 0.0f,
                                       short                    scalcoinitial = 1);
  void                      WriteAllTracesToFile(short scalcoinitial = 1); ///< Use only after writeTrace with x and y as input is used for the whole cube
  //<<<End write mode


  // int checkError(char * errText)
  //   {if(error_ > 0) strcpy(errText, errMsg_);return(error_);}
  /// Return (possibly upper limit for) number of traces

  size_t                    GetNTraces() const { return n_traces_ ;}
  size_t                    GetNz()      const { return nz_       ;}
  float                     GetDz()      const { return dz_       ;}

  enum                      OutsideModes{MISSING, ZERO, CLOSEST};

  size_t                    FindNumberOfTraces(void);
  static size_t             FindNumberOfTraces(const std::string       & fileName,
                                               const TraceHeaderFormat * traceHeaderFormat = NULL);


  static SegyGeometry     * FindGridGeometry(const std::string       & fileName,
                                             const TraceHeaderFormat * traceHeaderFormat = NULL);
  TraceHeaderFormat         GetTraceHeaderFormat(){return trace_header_format_;};
  static TraceHeaderFormat  FindTraceHeaderFormat(const std::string & fileName);

private:
  //void                      ebcdicHeader(std::string& outstring);               ///<
  bool                      ReadHeader(TraceHeader & header);                   ///< Trace header
  SegYTrace               * ReadTrace(const NRLib::Volume * volume,
                                      double                zPad,
                                      bool                & duplicateHeader,
                                      bool                  onlyVolume,
                                      bool                & outsideSurface,
                                      bool                  writevalues      = true,
                                      double              * outsideTopBot    = NULL,
                                      bool                  relative_padding = true);  ///< Read single trace from file
  //Note: If outsideTopBot == NULL, lack of data on top or bot will throw exception.
  //      Otherwise, outsideTopBot[0] will be top lack, [1] for bottom,
  //      [2] is x-coord, [3] is y-coord. Allocate outside.

  void                      WriteMainHeader(const TextualHeader& ebcdicHeader); ///< Quasi-dummy at the moment.
  void                      ReadDummyTrace(std::fstream & file, int format, size_t nz);
  /// Used to find correct trace header format.
  bool                      CompareTraces(TraceHeader *header1, TraceHeader *header2, int &delta, int &deltail, int &deltaxl);

  void                      SetBogusILXLUndefined(std::vector<NRLib::SegYTrace*> & traces);

  bool                      TraceHeaderOK(std::fstream &file, const TraceHeaderFormat *headerFormat);
  void                      FindDeltaILXL(TraceHeader *t1, TraceHeader *t2, TraceHeader *t3, double &dil, double &dxl, bool x);
  void                      CheckTopBotError(const double * tE, const double * bE); ///<Summarizes lack of data at top and bottom.

  TraceHeaderFormat         trace_header_format_;

  SegyGeometry            * geometry_;             ///< Parameters to find final index from i and j
  BinaryHeader            * binary_header_;         ///<

  bool                      single_trace_;          ///< Read one and one trace
  bool                      simbox_only_;           ///<
  bool                      check_simbox_;          ///<

  std::vector<SegYTrace*>   traces_;               ///< All traces
  size_t                    n_traces_;              ///< Holds the number of traces. May be an estimate if not all read.

  int                       datasize_;             ///< Bytes per datapoint in file.

  size_t                    nz_;                   ///< Number of time samples

  float                     z0_;                   ///< Top of segy cube
  float                     dz_;                   ///< Sampling density in time
  short                     delay_rec_time_;       ///< Delay recording time

  std::fstream              file_;
  std::string               file_name_;

  float                     rmissing_;

};




} // namespace NRLib

#endif
