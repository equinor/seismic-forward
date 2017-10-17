#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>

#include "nrlib/iotools/fileio.hpp"

#include "nrlib/statistics/continuouskrigingdata2d.hpp"
#include "nrlib/geometry/point.hpp"
#include "nrlib/surface/regularsurface.hpp"

using namespace NRLib;

ContinuousKrigingData2D::ContinuousKrigingData2D(int nData)
  : data_(0),   // Do not reserve space here (gives trouble with .push_back())
    xvector_(0),
    yvector_(0),
    xmin_(0),
    ymin_(0),
    dx_(0),
    dy_(0)
{
  //
  // Using .reserve() we set aside space for vectors, but such that
  // .push_back() used in AddData() adds the 0'th element first.
  //
  data_.reserve(nData);
  xvector_.reserve(nData);
  yvector_.reserve(nData);
}

//----------------------------------------------------------------------
ContinuousKrigingData2D::ContinuousKrigingData2D(const std::vector<Point* >    & data_points,
                                                 const RegularSurface<double>  & trend_grid)
  : data_(0),
    xvector_(0),
    yvector_(0),
    xmin_(trend_grid.GetXMin()),
    ymin_(trend_grid.GetYMin()),
    dx_(trend_grid.GetDX()),
    dy_(trend_grid.GetDY())
{
// map over values from std::vector<Point*> to vectors or coordinates and data (as indices and data in IndexKrigingData2D).
  for(size_t p = 0; p < data_points.size(); p++) {
    float value = static_cast<float>(data_points[p]->z);
    double x = data_points[p]->x;
    double y = data_points[p]->y;

    AddData(x, y, value);
  }
}

//---------------------------------------------------------------------
ContinuousKrigingData2D::~ContinuousKrigingData2D(void)
{
}

//---------------------------------------------------------------------
void
ContinuousKrigingData2D::AddData(double   x,
                                 double   y,
                                 double value)
{
  if (value != 99999.000) // Do not add missing values
  {
      xvector_.push_back(x);
      yvector_.push_back(y);
         data_.push_back(value);
  }
}

//---------------------------------------------------------------------
void
ContinuousKrigingData2D::WriteToFile(const std::string & fileName) const
{
  std::ofstream file;
  NRLib::OpenWrite(file, fileName);

  file << "    i     j     value\n";
  file << "---------------------\n";
  for (unsigned int k = 0 ; k < data_.size() ; k++)
    file  <<std::fixed << std::setprecision(2)
          << std::setw(5)  << xvector_[k] << " "
          << std::setw(5)  << yvector_[k] << " "
          << std::setw(10) << data_[k]   << std::endl;
  file.close();
}

//----------------------------------------------------------------------
size_t
ContinuousKrigingData2D::CountDataInBlock(size_t imin, size_t imax,
                                          size_t jmin, size_t jmax) const
//----------------------------------------------------------------------
{
  double xmax = xmin_ + imax*dx_;
  double xmin = xmin_ + imin*dx_;
  double ymax = ymin_ + jmax*dy_;
  double ymin = ymin_ + jmin*dy_;

  size_t count = 0;
  for (int k = 0; k < static_cast<int>(data_.size()); ++k){
   double x = xvector_[k];
   double y = yvector_[k];
    if (x <= xmax)
      if (x >= xmin)
        if (y <= ymax)
          if (y >= ymin)
            count++;
  }
  return count;
}

//----------------------------------------------------------------------
void
ContinuousKrigingData2D::AddDataToBlock(KrigingData2D & kriging_data_block,
                                        size_t imin, size_t imax,
                                        size_t jmin, size_t jmax) const
//----------------------------------------------------------------------
{

  ContinuousKrigingData2D &kriging_data_block_tmp = dynamic_cast<ContinuousKrigingData2D&>(kriging_data_block);
  double xmax = xmin_ + imax*dx_;
  double xmin = xmin_ + imin*dx_;
  double ymax = ymin_ + jmax*dy_;
  double ymin = ymin_ + jmin*dy_;
  for (int k = 0; k < static_cast<int>(data_.size()); ++k){
    double x = xvector_[k];
    double y = yvector_[k];
    if (x < xmax)
      if (x >= xmin)
        if (y < ymax)
          if (y >= ymin)
            kriging_data_block_tmp.AddData(x, y, data_[k]);
   }
 }
