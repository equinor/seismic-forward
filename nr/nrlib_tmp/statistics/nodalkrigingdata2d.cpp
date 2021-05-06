
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

#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>

//#include "nrlib/iotools/logkit.hpp"
#include "nrlib/iotools/fileio.hpp"

#include "nodalkrigingdata2d.hpp"

using namespace NRLib;

NodalKrigingData2D::NodalKrigingData2D(int nData)
  : data_(0),   // Do not reserve space here (gives trouble with .push_back())
    indexI_(0),
    indexJ_(0),
    count_(0)
{
  //
  // Using .reserve() we set aside space for vectors, but such that
  // .push_back() used in AddData() adds the 0'th element first.
  //
  data_.reserve(nData);
  indexI_.reserve(nData);
  indexJ_.reserve(nData);
  count_.reserve(nData);
}

//----------------------------------------------------------------------
NodalKrigingData2D::NodalKrigingData2D(const std::vector<Point* >   & data_points,
                                       const RegularSurface<double> & trend)
  : data_(0),
    indexI_(0),
    indexJ_(0),
    count_(0)
{
// map over values from std::vector<Point*> to the KrigingData2D-structure
  for(size_t p = 0; p < data_points.size(); p++) {
    size_t i;
    size_t j;
    float value = static_cast<float>(data_points[p]->z);
    trend.FindIndex(data_points[p]->x, data_points[p]->y, i, j);
    AddData(static_cast<int>(i), static_cast<int>(j), value);
  }
  FindMeanValues();
}

//---------------------------------------------------------------------
NodalKrigingData2D::~NodalKrigingData2D(void)
{
}

//---------------------------------------------------------------------
void
NodalKrigingData2D::AddData(int   i,
                       int   j,
                       float value)
{
  if (value != 99999.000) // Do not add missing values
  {
    int ij = GotBlock(i,j);
    if (ij == -1) // Data in new location
    {
      indexI_.push_back(i);
      indexJ_.push_back(j);
      count_.push_back(1);
      data_.push_back(value);
    }
    else // Add more data to same location
    {
      count_[ij] += 1;
      data_[ij] += value;
    }
  }
}

//---------------------------------------------------------------------
int
NodalKrigingData2D::GotBlock(int i, int j) const
{
  for (unsigned int k = 0 ; k < data_.size() ; k++)
    if (indexI_[k]==i && indexJ_[k]==j) // Do we already have data in block (i,j)?
      return int(k);
  return -1;
}

//---------------------------------------------------------------------
void
NodalKrigingData2D::FindMeanValues(void)
{
  for (unsigned int k = 0 ; k < data_.size() ; k++) {
    data_[k] /= count_[k];
    count_[k] = 1;
  }
}

//---------------------------------------------------------------------
void
NodalKrigingData2D::WriteToFile(const std::string & fileName) const
{
  std::ofstream file;
  NRLib::OpenWrite(file, fileName);

  file << "    i     j     value\n";
  file << "---------------------\n";
  for (unsigned int k = 0 ; k < data_.size() ; k++)
    file  <<std::fixed << std::setprecision(2)
          << std::setw(5)  << indexI_[k] << " "
          << std::setw(5)  << indexJ_[k] << " "
          << std::setw(10) << data_[k]   << std::endl;
  file.close();
}

//--------------------------------------------------------------
size_t
NodalKrigingData2D::CountDataInBlock(size_t imin, size_t imax,
                                     size_t jmin, size_t jmax) const
//--------------------------------------------------------------
{
  size_t count = 0;
  for (int k = 0; k < static_cast<int>(data_.size()); ++k){
    size_t index_i = indexI_[k];
    size_t index_j = indexJ_[k];
    if (index_i <= imax)
      if (index_i >= imin)
        if (index_j <= jmax)
          if (index_j >= jmin)
            count++;
  }
  return count;
}

//--------------------------------------------------------------
void
NodalKrigingData2D::AddDataToBlock(KrigingData2D & kriging_data_block,
                                   size_t imin, size_t imax,
                                   size_t jmin, size_t jmax) const
//--------------------------------------------------------------
{

  NodalKrigingData2D &kriging_data_block_tmp = dynamic_cast<NodalKrigingData2D&>(kriging_data_block);
  for (int k = 0; k < static_cast<int>(data_.size()); ++k){
    size_t index_i = static_cast<size_t>(indexI_[k]);
    size_t index_j = static_cast<size_t>(indexJ_[k]);
    if (index_i <= imax)
      if (index_i >= imin)
        if (index_j <= jmax)
          if (index_j >= jmin)
            kriging_data_block_tmp.AddData(static_cast<int>(index_i), static_cast<int>(index_j), data_[k]);
  }
}

