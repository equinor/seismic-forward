// $Id: grid2dgeneral.hpp 1308 2015-09-08 09:50:31Z perroe $

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

#ifndef NRFAULT_FAULTSEAL_GRID2DGENERAL_H
#define NRFAULT_FAULTSEAL_GRID2DGENERAL_H

#include <cassert>
#include <sstream>
#include <vector>

namespace NRLib {

template<class A>
class Grid2DGeneral {
  // This class implements a generalized Grid2D class. Grid2DGeneral takes as input global indices that defines a part of a large grid.
public:
  typedef typename std::vector<A>::iterator        iterator;
  typedef typename std::vector<A>::const_iterator  const_iterator;
  typedef typename std::vector<A>::reference       reference;
  typedef typename std::vector<A>::const_reference const_reference;

  // Constructor with no parameters: Initilizing size and data to 0.
  Grid2DGeneral();

  // Constructor with indices defining the grid
  Grid2DGeneral(size_t i_min, size_t i_max, size_t j_min, size_t j_max, const A& val = A());

  virtual ~Grid2DGeneral();

  virtual void           Resize(size_t i_min, size_t i_max, size_t j_min, size_t j_max, const A& val = A());

  inline reference       operator()(size_t i, size_t j);
  inline reference       operator()(size_t index);

  inline const_reference operator()(size_t i, size_t j) const;
  inline const_reference operator()(size_t index) const;

  iterator               begin()       { return data_.begin(); }
  iterator               end()         { return data_.end(); }

  const_iterator         begin() const { return data_.begin(); }
  const_iterator         end()   const { return data_.end(); }

  size_t GetIMin() const {return i_min_; }
  size_t GetIMax() const {return i_max_; }
  size_t GetJMin() const {return j_min_; }
  size_t GetJMax() const {return j_max_; }

  size_t                 GetNI() const { return ni_; }
  size_t                 GetNJ() const { return nj_; }
  size_t                 GetN()  const { return data_.size(); }

  inline size_t GetIndex(size_t i, size_t j) const;
  void          GetIJ(size_t index, size_t &i, size_t &j) const;
  bool          IsValidIndex(size_t i, size_t j) const;

  void                   Swap(Grid2DGeneral<A>& other);

  A                      FindMin(A missingValue) const;
  A                      FindMax(A missingValue) const;


private:
  size_t ni_;
  size_t nj_;
  std::vector<A> data_; /// The grid data, column-major ordering.

  // Start and endpoint (global) indices
  size_t i_min_;
  size_t i_max_;
  size_t j_min_;
  size_t j_max_;
};

template<class A>
Grid2DGeneral<A>::Grid2DGeneral()
{
  ni_ = 0;
  nj_ = 0;
  data_.resize(0);
  i_min_ = 0;
  i_max_ = 0;
  j_min_ = 0;
  j_max_ = 0;
}

template<class A>
Grid2DGeneral<A>::Grid2DGeneral(size_t i_min, size_t i_max, size_t j_min, size_t j_max, const A& val)
: i_min_(i_min),
i_max_(i_max),
j_min_(j_min),
j_max_(j_max)
{
  ni_ = i_max_ - i_min_ + 1;
  nj_ = j_max_ - j_min_ + 1;
  data_.resize(ni_*nj_, val);
}

template<class A>
Grid2DGeneral<A>::~Grid2DGeneral()
{}

template<class A>
void Grid2DGeneral<A>::Resize(size_t i_min, size_t i_max, size_t j_min, size_t j_max, const A& val)
{
  i_min_ = i_min;
  i_max_ = i_max;
  j_min_ = j_min;
  j_max_ = j_max;

  ni_ = i_max_ - i_min_ + 1;
  nj_ = j_max_ - j_min_ + 1;

  data_.resize(0); //To avoid copying of elements
  data_.resize(ni_ * nj_, val);
}

template<class A>
typename Grid2DGeneral<A>::reference Grid2DGeneral<A>::operator()(size_t i, size_t j)
{
  return(data_[GetIndex(i, j)]);
}

template<class A>
typename Grid2DGeneral<A>::reference Grid2DGeneral<A>::operator()(size_t index)
{
  assert(index < GetN());

  return(data_[index]);
}

template<class A>
typename Grid2DGeneral<A>::const_reference Grid2DGeneral<A>::operator()(size_t i, size_t j) const
{
  return(data_[GetIndex(i, j)]);
}

template<class A>
typename Grid2DGeneral<A>::const_reference Grid2DGeneral<A>::operator()(size_t index) const
{
  assert(index < GetN());

  return(data_[index]);
}

template<class A>
size_t Grid2DGeneral<A>::GetIndex(size_t i, size_t j) const
{
  assert(i >= i_min_);
  assert(i <= i_max_);
  assert(j >= j_min_);
  assert(j <= j_max_);

  return((i-i_min_)+(j-j_min_)*ni_);
}

template<class A>
void Grid2DGeneral<A>::GetIJ(size_t index, size_t &i, size_t &j) const
{
  assert (index < GetN());

  i = (index % ni_);
  i = i + i_min_;
  j = static_cast<int>(index/ni_); // Implicit floor
  j = j + j_min_;
}

template<class A>
bool Grid2DGeneral<A>::IsValidIndex(size_t i, size_t j) const
{
  if (i >= i_min_ && i <= i_max_ && j >= j_min_ && j <= j_max_)
    return true;

  return false;
}

template<class A>
void Grid2DGeneral<A>::Swap(NRLib::Grid2DGeneral<A> &other)
{
  std::swap(ni_, other.ni_);
  std::swap(nj_, other.nj_);
  data_.swap(other.data_);
}

template<class A>
A Grid2DGeneral<A>::FindMin(A missingValue) const
{
  A minVal = (*this)(0);
  typename std::vector<A>::const_iterator i;
  for (i = this->begin(); i < this->end(); i++) {
    if ((minVal == missingValue || (*i) < minVal) && (*i) != missingValue)
      minVal = *i;
  }
  return minVal;
}

template<class A>
A Grid2DGeneral<A>::FindMax(A missingValue) const
{
  A maxVal = (*this)(0);
  typename std::vector<A>::const_iterator i;
  for (i = this->begin(); i < this->end(); i++) {
    if ((maxVal == missingValue || (*i) > maxVal) && (*i) != missingValue)
      maxVal = *i;
  }
  return maxVal;
}

} // namespace NRLib

#endif // NRFAULT_FAULTSEAL_GRID2DGENERAL_H
