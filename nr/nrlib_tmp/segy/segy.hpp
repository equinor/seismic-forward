// $Id: segy.hpp 1695 2017-09-29 12:02:24Z heidi $

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

#include <string>
#include <vector>

#include "traceheader.hpp"
#include "commonheaders.hpp"

#include "../iotools/bigfile.hpp"
#include "../segy/segygeometry.hpp"
#include "../segy/segytrace.hpp"
#include "../volume/volume.hpp"
#include "../stormgrid/stormcontgrid.hpp"

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
       double                    z0,
       const TraceHeaderFormat & traceHeaderFormat);


  /// Constructor for reading unknown format
  /// Read only the headers on top of the file
  /// \param[in] fileName  Name of file to read data from
  /// \param[in] z0
  /// \param[in] thf Vector of pointers to possible traceheaderformats. If NULL, default list is used.
  /// \throw IOError if the file can not be opened.
  /// \throw FileFormatError if the traceheaderformat can not be recognized.
  SegY(const std::string               & fileName,
       double                            z0,
       std::vector<TraceHeaderFormat *>  thf = std::vector<TraceHeaderFormat *>(0),
       bool                              searchStandardFormats = true);

  /// Constructor for writing
  /// \param[in] fileName  Name of file to write data to
  /// \throw IOError if the file can not be opened.
  SegY(const std::string       & fileName,
       double                    z0,
       size_t                    nz,
       double                    dz,
       const TextualHeader     & ebcdicHeader,
       const TraceHeaderFormat & traceHeaderFormat = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS));

  /// Constructor for writing
  /// Must be initialized
  SegY();

 /// Constructor for writing of segy from storm cube
 SegY(const StormContGrid     * storm_grid,
      const SegyGeometry      * geometry, //May be set to 0.
      double                    z0,
      double                    dz,
      int                       nz,
      const std::string       & file_name = "",
      bool                      write_to_file = true,
      const TraceHeaderFormat & trace_header_format = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS),
      bool                      is_seismic = false);

  ~SegY();

  /// Initializing empty constructor for writing
  void Initialize(const std::string       & fileName,
                  double                    z0,
                  size_t                    nz,
                  double                    dz,
                  const TextualHeader     & ebcdicHeader,
                  const TraceHeaderFormat & traceHeaderFormat = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS),
                  short                     n_traces_per_ensamble = 1);


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
  void                      ReportSizeOfVolume(void) const;

  void                      GetNearestTrace(std::vector<float> & trace_data,
                                            bool               & missing,
                                            double             & z0_data,
                                            double               x,
                                            double               y) const;

  std::vector<float>        GetAllValues();                           ///< Return vector with all values.

  void                      CreateRegularGrid(bool regularize_if_needed = true);
  const SegyGeometry      * GetGeometry(void) const { return geometry_ ;} //Only makes sense after command above, or FindAndSetGeometry below.
  //<<<End read all trace mode

  //These functions are valid both after ReadAllTraces and FindAndSetGridGeometry
  const TraceHeader       & GetTraceHeader(int IL, int XL) const;
  const TraceHeader       & GetTraceHeader(double x, double y) const;

  bool                      IsTraceDefined(int IL, int XL) const;

  //>>>Begin read single trace mode
  /// \param remove_bogus_traces   Remove traces that fall between the most-used steps in file.
  SegyGeometry            * FindGridGeometry(bool only_ilxl = false,
                                             bool keep_header = false,
                                             bool remove_bogus_traces = true);        //Note: This and function below also sets all traceheaders with file position.
  void                      FindAndSetGridGeometry(bool only_ilxl = false,
                                                   bool keep_header = false,
                                                   bool remove_bogus_traces = true);

  SegYTrace               * GetNextTrace(double                zPad = 0.0,
                                         const NRLib::Volume * volume = NULL,
                                         bool                  onlyVolume = false);
  void                      GetTraceData(int                  IL,
                                         int                  XL,
                                         std::vector<float> & result,
                                         const Volume       * volume = NULL);   // Only makes sense after FindAndSetGridGeometry.
  void                      GetTraceData(double               x,
                                         double               y,
                                         std::vector<float> & result,
                                         const Volume       * volume = NULL);   // Only makes sense after FindAndSetGridGeometry.
  void                      GetTraceData(std::streampos       pos,
                                         std::vector<float> & result,
                                         double               z_top = -1,       // Does not check pos, z_top or z_bot for validity.
                                         double               z_bot = -1);      // Top/bot = -1 means start from top/go to bottom
  std::streampos            GetFilePos(int IL, int XL) const;                   // Only makes sense after FindAndSetGridGeometry. Returns 0 for invalid trace.
  std::streampos            GetFilePos(double x, double y) const;                 // Only makes sense after FindAndSetGridGeometry. Returns 0 for invalid trace.
  //<<<End read single trace mode

  //>>>Begin write mode
  void                      SetGeometry(const SegyGeometry * geometry);

  void                      StoreTrace(double                     x,
                                       double                     y,
                                       const std::vector<float> & data,
                                       const NRLib::Volume      * volume,
                                       float                      topVal  = 0.0f,
                                       float                      baseVal = 0.0f);

  /// Write single trace to file
  /// Assumes traceheader read from input file.
  void                      WriteTrace(const TraceHeader        & traceHeader,
                                       const std::vector<float> & data,
                                       const NRLib::Volume      * volume,
                                       float                      topVal  = 0.0f,
                                       float                      baseVal = 0.0f);

  /// Write single trace to internal memory
  /// Assumes that segygeometry is set.
  void                      WriteTrace(double                     x,
                                       double                     y,
                                       const std::vector<float> & data,
                                       const NRLib::Volume      * volume,
                                       float                      topVal  = 0.0f,
                                       float                      baseVal = 0.0f,
                                       short                      scalcoinitial = 1,
                                       short                      offset = 0);

  /// Write single trace to file
  void                      WriteTrace(double                     x,
                                       double                     y,
                                       int                        IL,
                                       int                        XL,
                                       const std::vector<float> & data,
                                       short                      scalcoinitial = 1,
                                       short                      offset = 0);


  void                      WriteAllTracesToFile(short scalcoinitial = 1); ///< Use only after writeTrace with x and y as input is used for the whole cube
  //<<<End write mode


  // int checkError(char * errText)
  //   {if(error_ > 0) strcpy(errText, errMsg_);return(error_);}
  /// Return (possibly upper limit for) number of traces

  size_t                    GetNTraces() const { return n_traces_ ;}
  size_t                    GetNz()      const { return nz_       ;}
  double                    GetDz()      const { return dz_       ;}
  double                    GetTop()     const { return z0_       ;}

  bool                      GetSamplingInconsistency() const { return sampling_inconsistency_ ;}

  enum                      OutsideModes{MISSING, ZERO, CLOSEST};

  size_t                    FindNumberOfTraces(void);
  static size_t             FindNumberOfTraces(const std::string       & fileName,
                                               const TraceHeaderFormat * traceHeaderFormat = NULL);


  static SegyGeometry     * FindGridGeometry(const std::string       & fileName,
                                             const TraceHeaderFormat * traceHeaderFormat = NULL);
  TraceHeaderFormat         GetTraceHeaderFormat(){return trace_header_format_;};
  static TraceHeaderFormat  FindTraceHeaderFormat(const std::string & fileName);

  SegYTrace *              getTrace(int i) {return traces_[i];};

private:
  //void                      ebcdicHeader(std::string& outstring);               ///<
  bool                      ReadHeader(TraceHeader & header);                   ///< Trace header
  void                      ReadAndSetBinaryHeader(NRLib::BigFile& file);
  SegYTrace               * ReadTrace(const NRLib::Volume * volume,
                                      double                zPad,
                                      bool                & duplicateHeader,
                                      bool                  onlyVolume,
                                      bool                & outsideSurface,
                                      bool                  writevalues      = false,
                                      double              * outsideTopBot    = NULL,
                                      bool                  relative_padding = true);  ///< Read single trace from file
  //Note: If outsideTopBot == NULL, lack of data on top or bot will throw exception.
  //      Otherwise, outsideTopBot[0] will be top lack, [1] for bottom,
  //      [2] is x-coord, [3] is y-coord. Allocate outside.

  void                      WriteMainHeader(const TextualHeader& ebcdicHeader); ///< Quasi-dummy at the moment.
  static void               SkipTraceData(NRLib::BigFile & file, int datasize, size_t nz);
  /// Used to find correct trace header format.
  bool                      CompareTraces(const TraceHeader& header1, const TraceHeader& header2,
                                          int& delta, int& deltail, int& deltaxl);

  void                      SetBogusILXLUndefined(std::vector<NRLib::SegYTrace*> & traces);

  bool                      TraceHeaderOK(NRLib::BigFile & file, const TraceHeaderFormat & header_format);
  void                      FindDeltaILXL(const TraceHeader& t1, const TraceHeader& t2, const TraceHeader& t3,
                                          double &dil, double &dxl, bool x);
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

  double                    z0_;                   ///< Top of segy cube
  double                    dz_;                   ///< Sampling density in time

  NRLib::BigFile            file_;

  float                     rmissing_;
  size_t                    n_traces_per_ensamble_;

  bool                      sampling_inconsistency_;   ///< If sampling in trace header is inconsistent with sampling in binary header

};




} // namespace NRLib

#endif
