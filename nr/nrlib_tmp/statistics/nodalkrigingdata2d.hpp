
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


#ifndef NRLIB_STATISTICS_NODALKRIGINGDATA2D_HPP
#define NRLIB_STATISTICS_NODALKRIGINGDATA2D_HPP

#include <vector>
#include "nrlib/statistics/krigingdata2d.hpp"
#include "nrlib/geometry/point.hpp"
#include "nrlib/surface/regularsurface.hpp"

namespace NRLib {
class NodalKrigingData2D : public KrigingData2D
{
public:
  NodalKrigingData2D(int nData = 0);
  NodalKrigingData2D(const std::vector<Point* > & data_points, const RegularSurface<double> & trend_grid);
  ~NodalKrigingData2D(void);

  KrigingData2D* Clone()                           const { return new NodalKrigingData2D(*this); }
  KrigingData2D* NewInstance()                     const { return new NodalKrigingData2D(); }

  void                       AddData(int   i,
                                     int   j,
                                     float value);
  void                       FindMeanValues(void);

  const std::vector<float> & GetData(void)         const { return data_                          ;}
  const std::vector<int>   & GetIndexI(void)       const { return indexI_                        ;}
  const std::vector<int>   & GetIndexJ(void)       const { return indexJ_                        ;}

  inline int                 GetNumberOfData(void) const { return static_cast<int>(data_.size()) ;}
  inline float               GetData(int k)        const { return data_[k]                       ;}
  inline int                 GetIndexI(int k)      const { return indexI_[k]                     ;}
  inline int                 GetIndexJ(int k)      const { return indexJ_[k]                     ;}

  size_t                     CountDataInBlock(size_t imin, size_t imax,
                                              size_t jmin, size_t jmax) const;

  // This function add _this_ object to the input parameter object kriging_data_block.
  void                       AddDataToBlock(KrigingData2D & kriging_data_block,
                                            size_t imin, size_t imax,
                                            size_t jmin, size_t jmax) const;

  void                       WriteToFile(const std::string & name) const;

private:
  int                        GotBlock(int i, int j) const;

  std::vector<float>         data_;
  std::vector<int>           indexI_;
  std::vector<int>           indexJ_;
  std::vector<int>           count_;

};
}
#endif
