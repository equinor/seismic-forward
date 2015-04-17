// $Id: eclipsefault.cpp 1105 2012-11-01 07:33:35Z georgsen $

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

#include "eclipsefault.hpp"
#include "../iotools/stringtools.hpp"

#include <cassert>
#include <fstream>
#include <string>
#include <vector>

using namespace NRLib;


EclipseFault::EclipseFault()
{}


EclipseFault::EclipseFault(const std::string & name)
  : fault_name_(name)
{}


bool EclipseFault::IsStairStepped() const {
  for (SegmentConstIterator it = SegmentsBegin(); it != SegmentsEnd(); ++it) {
    if (it->face == PosZ || it->face == NegZ)
      return true;
  }
  return false;
}


void EclipseFault::ReadSegments(const std::vector<std::string>& line_segment) {
  Segment segment;

  segment.i_from = ParseType<int>(line_segment[1]) - 1;
  segment.i_to   = ParseType<int>(line_segment[2]) - 1;
  segment.j_from = ParseType<int>(line_segment[3]) - 1;
  segment.j_to   = ParseType<int>(line_segment[4]) - 1;
  segment.k_from = ParseType<int>(line_segment[5]) - 1;
  segment.k_to   = ParseType<int>(line_segment[6]) - 1;

  std::string face_str = line_segment[7];
  // Remove the 's
  face_str.erase(face_str.size() - 1, 1);
  face_str.erase(0, 1);
  if (face_str == "X" || face_str == "I") {
    segment.face = PosX;
  }
  else if (face_str == "Y" || face_str == "J") {
    segment.face = PosY;
  }
  else if (face_str == "Z" || face_str == "K") {
    segment.face = PosZ;
  }
  else if (face_str == "X-" || face_str == "I-") {
    segment.face = NegX;
  }
  else if (face_str == "Y-" || face_str == "J-") {
    segment.face = NegY;
  }
  else if (face_str == "Z-" || face_str == "K-") {
    segment.face = NegZ;
  }
  segments_.push_back(segment);
}


void EclipseFault::WriteFault(std::ofstream& out_file, const std::string& fault_name) const
{
  size_t i;
  std::string quoted_name = "'" + fault_name + "'";
  for(i = 0; i < segments_.size(); i++){
    out_file << quoted_name << setw(18)
             << segments_[i].i_from + 1 << std::setw(8) << segments_[i].i_to + 1 << std::setw(10)
             << segments_[i].j_from + 1 << std::setw(8) << segments_[i].j_to + 1 << std::setw(10)
             << segments_[i].k_from + 1 << std::setw(8) << segments_[i].k_to + 1 << std::setw(8);
    if (segments_[i].face == PosX){
      out_file  << "'X'" << std::setw(5) << "/\n";
    }
    else if (segments_[i].face == NegX){
      out_file  << "'X-'" << std::setw(5) << "/\n";
    }
    else if (segments_[i].face == PosY){
      out_file  << "'Y'" << std::setw(5) << "/\n";
    }
    else if (segments_[i].face == NegY){
      out_file  << "'Y-'" << std::setw(5) << "/\n";
    }
    else if (segments_[i].face == PosZ){
      out_file  << "'Z'" << std::setw(5) << "/\n";
    }
    else if (segments_[i].face == NegZ){
      out_file  << "'Z-'" << std::setw(5) << "/\n";
    }
  }
}


void EclipseFault::GetSegmentLowerCoord(const Segment& segment, int& i, int& j)
{
  switch(segment.face) {
  case PosX:
    i = segment.i_from + 1;
    j = segment.j_from;
    break;
  case PosY:
    i = segment.i_from;
    j = segment.j_from + 1;
    break;
  case NegX:
  case NegY:
  case PosZ:
  case NegZ:
    i = segment.i_from;
    j = segment.j_from;
    break;
  default:
    assert(0); // All faces are accounted for.
  }
}


void EclipseFault::GetSegmentUpperCoord(const Segment& segment, int& i, int& j)
{
  switch(segment.face) {
  case NegX:
    i = segment.i_to;
    j = segment.j_to + 1;
    break;
  case NegY:
    i = segment.i_to + 1;
    j = segment.j_to;
    break;
  case PosX:
  case PosY:
  case PosZ:
  case NegZ:
    i = segment.i_to + 1;
    j = segment.j_to + 1;
    break;
  default:
    assert(0); // All faces are accounted for.
  }
}


void EclipseFault::IncreaseSegmentLength(Segment& segment, int increment)
{
  switch (segment.face) {
  case NRLib::PosX:
  case NRLib::NegX:
    segment.j_from -= increment;
    segment.j_to   += increment;
    break;
  case NRLib::PosY:
  case NRLib::NegY:
    segment.i_from -= increment;
    segment.i_to   += increment;
    break;
  default:
    assert(0);    // Horizontal segments are not supported.
  }
}


bool EclipseFault::SegmentsCross(const Segment & segment1,
                                 const Segment & segment2,
                                 int           & i,
                                 int           & j)
{
  // Horizontal fault segments are not supported.
  assert(GetDirection(segment1.face) != DirectionZ &&
         GetDirection(segment2.face) != DirectionZ);

  // parallel segments don't cross.
  if ( GetDirection(segment1.face) == GetDirection(segment2.face) ) {
    return false;
  }

  int i_from_1, i_to_1, j_from_1, j_to_1, i_from_2, i_to_2, j_from_2, j_to_2;

  GetSegmentLowerCoord(segment1, i_from_1, j_from_1);
  GetSegmentUpperCoord(segment1, i_to_1, j_to_1);
  GetSegmentLowerCoord(segment2, i_from_2, j_from_2);
  GetSegmentUpperCoord(segment2, i_to_2, j_to_2);

  if (i_from_1 <= i_from_2 && i_to_1 >= i_to_2) {
    if (j_from_2 <= j_from_1 && j_to_2 >= j_to_1) {
      if ( (segment1.face == PosX || segment1.face == NegX) ) {
        i = i_from_1;
        j = j_from_2;
      } else {
        i = i_from_2;
        j = j_from_1;
      }
      return true;
    } else {
      return false;
    }
  }

  if ( j_from_1 <= j_from_2 && j_to_1 >= j_to_2) {
    if (i_from_2 <= i_from_1 && i_to_2 >= i_to_1) {
      if ( (segment1.face == PosX || segment1.face == NegX) ) {
        i = i_from_1;
        j = j_from_2;
      } else {
        i = i_from_2;
        j = j_from_1;
      }
      return true;
    } else {
      return false;
    }
  }
  return false;
}

