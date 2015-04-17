// $Id: eclipsefaceindexer.cpp 882 2011-09-23 13:10:16Z perroe $

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

#include "eclipsefaceindexer.hpp"
#include "eclipsetools.hpp"


using namespace std;

namespace NRLib {

FaceIndexer::FaceIndexer(const EclipseFault &fault)
{
  segment_iterator_ = fault.SegmentsBegin();
  end_iterator_ = fault.SegmentsEnd();
  i_ = segment_iterator_->i_from;
  j_ = segment_iterator_->j_from;
  k_ = segment_iterator_->k_from;
}

CellFace FaceIndexer::GetFace()
{
  CellFace cellface;
  cellface.i_ = i_;
  cellface.j_ = j_;
  cellface.k_ = k_;
  cellface.face_ = segment_iterator_->face;

  return cellface;
}

void FaceIndexer::operator ++()
{
  if ( this->IsFinished() == false ) {
    if (i_ == segment_iterator_->i_to) {
      if (j_ == segment_iterator_->j_to) {
        if (k_ == segment_iterator_->k_to) {
          segment_iterator_++;
          if (this->IsFinished() == false) {
            i_ = segment_iterator_->i_from;
            j_ = segment_iterator_->j_from;
            k_ = segment_iterator_->k_from;
          }
        }
        else{
          k_++;
          i_ = segment_iterator_->i_from;
          j_ = segment_iterator_->j_from;
        }
      }
      else{
        j_++;
        i_ = segment_iterator_->i_from;
      }
    }
    else{
      i_++;
    }
  }//end of IsFinished
}

bool FaceIndexer::IsFinished()
{
  if (segment_iterator_ == end_iterator_)
    return true;
  else
    return false;
}

} // namespace NRLib

