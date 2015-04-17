// $Id: segytrace.hpp 1199 2013-10-02 08:24:02Z anner $

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

#ifndef SEGYTRACE_HPP
#define SEGYTRACE_HPP

#include <fstream>
#include <string>
#include <vector>

#include "traceheader.hpp"
#include "commonheaders.hpp"

namespace NRLib
{
class SegYTrace
{
public:

 ///  Constructor
 /// \param[in] data           Data in trace
 /// \param[in] jStart         Start index in trace
 /// \param[in] jEnd           End index in trace
 /// \param[in] format
 /// \param[in] x              x coordinate of trace
 /// \param[in] y              y coordinate of trace
 /// \param[in] inLine         inline number
 /// \param[in] crossLine      crossline number


  SegYTrace(std::fstream      & file,
            size_t              jStart,
            size_t              jEnd,
            int                 format,
            size_t              nz,
            const TraceHeader * trace_header = NULL);                                     ///< Standard reading constructor.

  SegYTrace(std::vector<float> indata,
            size_t             jStart,
            size_t             jEnd,
            double              x,
            double              y,
            int                inLine,
            int                crossLine);                                                ///< Standard writing constructor

  SegYTrace(const TraceHeader & trace_header,
            bool                keep_header = true);                                      ///< Constructor for handling only headers.

  ~SegYTrace();

  void SetTableIndex(size_t index) {table_index_ = index;}                                ///< Set table index

  const std::vector<float> & GetTrace(void)              const { return data_       ;}
  float                      GetValue(size_t j)          const;                           ///< get trace value at index j
  size_t                     GetLegalIndex(size_t index) const;
  size_t                     GetStart()                  const { return j_start_    ;}    ///< Get start index
  size_t                     GetEnd()                    const { return j_end_      ;}    ///< Get end index
  double                     GetX()                      const { return x_          ;}    ///< Get UTM x coordinate
  double                     GetY()                      const { return y_          ;}    ///< Get UTM y coordinate
  int                        GetInline()                 const { return in_line_    ;}    ///< Get inline number
  int                        GetCrossline()              const { return cross_line_ ;}    ///< Get crossline number
  double                     GetCoord1()                 const { return coord1_     ;}    ///< Get first coordinate current system (either UTM or ILXL, defined in header)
  double                     GetCoord2()                 const { return(coord2_);}        ///< Get second coordinate current system (either UTM or ILXL, defined in header)
  size_t                     GetTableIndex()             const { return(table_index_);}   ///< Get tableIndex
  const TraceHeader        & GetTraceHeader()            const;                           ///< Traceheader
  std::streampos             GetFilePos()                const { return(file_position_);} ///< Get file position

  void SetFilePos(std::streampos pos) {file_position_ = pos;} /// Set file position


  void RemoveXY() { /// Void invalid x and y coordinates
    x_      = static_cast<float>(in_line_);
    y_      = static_cast<float>(cross_line_);
    coord1_ = x_;
    coord2_ = y_;}

  void SetILXLUndefined() {
    in_line_ = imissing_ ;
    cross_line_ = imissing_;
  }

private:
  ///(note that this class can live without trace_header, hence duplicates of information
  ///that may also be stored there.)

  std::vector<float> data_;         ///< Data in trace
  size_t             j_start_;      ///< Start index
  size_t             j_end_;        ///< End index
  double             x_;            ///< UTM x coord
  double             y_;            ///< UTM y coord
  int                in_line_;      ///< inline index
  int                cross_line_;   ///< crossline index
  double             coord1_;       ///< Coordinate 1 (either x or IL)
  double             coord2_;       ///< Coordinate 2 (either y or XL)
  float              rmissing_;
  int                imissing_;

  TraceHeader      * trace_header_; ///< The trace header as read from file.
  size_t             table_index_;  ///< Index in table, handy for some transitions, used for sorting in WriteAllTraces in SegY.
  std::streampos     file_position_;///< Position in file, used for GetRandomTrace. 0 means undefined (no trace starts at 0).
};

} // namespace NRLib

#endif
