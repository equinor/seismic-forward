// $Id: eclipsefault.hpp 1068 2012-09-18 11:21:53Z perroe $

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

#ifndef NRLIB_ECLIPSEGRID_ECLIPSEFAULT_HPP
#define NRLIB_ECLIPSEGRID_ECLIPSEFAULT_HPP

#include "eclipsetools.hpp"
#include <vector>
#include <string>

namespace NRLib {

class EclipseFault {

public:

  struct Segment{
    int i_from;
    int i_to;
    int j_from;
    int j_to;
    int k_from;
    int k_to;
    Face face;
  };

  typedef std::vector<Segment>::iterator       SegmentIterator;
  typedef std::vector<Segment>::const_iterator SegmentConstIterator;

  EclipseFault();
  explicit EclipseFault(const std::string& name);

  const std::string& GetName() const         { return fault_name_; }
  void SetName(const std::string& name)      { fault_name_ = name; }

  SegmentConstIterator SegmentsBegin() const { return segments_.begin(); }
  SegmentConstIterator SegmentsEnd()   const { return segments_.end(); }

  void AddSegment(Segment segment)           { segments_.push_back(segment); }

  bool IsStairStepped() const;

  void ReadSegments(const std::vector<std::string>& line_segment);
  void WriteFault(std::ofstream& out_file, const std::string& fault_name) const;

  // =============== Segment handeling methods ===================
  static void GetSegmentLowerCoord(const Segment& segment, int& i, int& j);
  static void GetSegmentUpperCoord(const Segment& segment, int& i, int& j);
  static void IncreaseSegmentLength(Segment& segment, int increment);
  static bool SegmentsCross(const Segment & segment1,
                            const Segment & segment2,
                            int           & i,
                            int           & j);

private:
  std::string fault_name_;

  std::vector<Segment> segments_;

};

}

#endif
