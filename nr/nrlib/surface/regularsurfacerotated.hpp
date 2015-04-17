// $Id: regularsurfacerotated.hpp 882 2011-09-23 13:10:16Z perroe $

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

#ifndef NRLIB_REGULARSURFACEROTATED_HPP
#define NRLIB_REGULARSURFACEROTATED_HPP

#include "surface.hpp"
#include "regularsurface.hpp"
#include "../grid/grid2d.hpp"
#include "../iotools/stringtools.hpp"

#include <cmath>
#include <algorithm>

namespace NRLib {
template <class A>
class RegularSurfaceRotated : public Surface<A> {
public:
  RegularSurfaceRotated();

  RegularSurfaceRotated(double x0,
                        double y0,
                        double lx,
                        double ly,
                        size_t nx,
                        size_t ny,
                        double angle,
                        const A& value = A());

  RegularSurfaceRotated(const RegularSurface<A>& surface, double angle);

  /// \brief Read surface file on given format.
  /// \param filename  File name.
  /// \param format    File format. If SURF_UNKNOWN we try to determine the format.
  /// \throws IOError         If we failed to open the file
  /// \throws FileFormatError If we can not determine the file format, or the contents
  ///                         of the file does not match the file format.
  RegularSurfaceRotated(const std::string& filename,
                        SurfaceFileFormat  format = SURF_UNKNOWN);

  Surface<A>* Clone() const
  { return new RegularSurfaceRotated<A>(*this); }

  typedef typename std::vector<A>::iterator        iterator;
  iterator               begin()       { return surface_.begin(); }
  iterator               end()         { return surface_.end(); }

  typedef typename std::vector<A>::const_iterator  const_iterator;
  const_iterator         begin() const { return surface_.begin(); }
  const_iterator         end()   const { return surface_.end(); }

  typedef typename std::vector<A>::reference       reference;
  inline reference       operator()(size_t i, size_t j) { return surface_(i,j); }
  inline reference       operator()(size_t index)       { return surface_(index); }

  typedef typename std::vector<A>::const_reference const_reference;
  inline const_reference operator()(size_t i, size_t j) const { return surface_(i,j); }
  inline const_reference operator()(size_t index) const       { return surface_(index); }

  /// Return z. Returns missing if we are outside the grid.
  /// Returns also values when not all grid cell corners are present,
  /// e.g. when (x,y) is just outside the grid.
  A GetZ(double x, double y) const;

  /// Return z. Returns missing if not all grid cell corners are present,
  /// e.g. when (x,y) is just outside the grid.
  A GetZInside(double x, double y) const;

  /// Checks if point is inside definition area for surface.
  bool IsInsideSurface(double x, double y) const;

  bool EnclosesRectangle(double x_min, double x_max,
                         double y_min, double y_max) const;

/// Sets all values on the surface to a constant value.
  void Assign(A c) {surface_.Assign(c);}

  void Add(A c) { surface_.Add(c);}

  void Subtract(A c) {surface_.Subtract(c);}

  void Multiply(A c) {surface_.Multiply(c);}

  void GetXY(size_t i, size_t j, double & x, double & y) const;

  void GetXY(size_t index, double & x, double & y) const;

  double GetDX()      const { return surface_.GetDX(); }
  double GetDY()      const { return surface_.GetDY(); }
  double GetLengthX() const { return surface_.GetLengthX(); }
  double GetLengthY() const { return surface_.GetLengthY(); }

  size_t GetNI() const { return surface_.GetNI(); }
  size_t GetNJ() const { return surface_.GetNJ(); }
  size_t GetN()  const { return surface_.GetN(); }

   /// Check if grid value is missing
  bool IsMissing(A val)                const { return surface_.IsMissing(val); }
  bool IsMissingAt(size_t i, size_t j) const { return surface_.IsMissingAt(i, j); }
  void SetMissing(size_t i, size_t j)        { surface_.SetMissing(i, j); }

  A Min()             const { return surface_.Min(); }
  A Max()             const { return surface_.Max(); }
  A Avg()             const { return surface_.Avg(); }
  A Avg(int& n_nodes) const { return surface_.Avg(n_nodes); }

  double GetXRef() const { return x_ref_; }
  double GetYRef() const { return y_ref_; }
  double GetXMin() const { return x_min_; }
  double GetYMin() const { return y_min_; }
  double GetXMax() const { return x_max_; }
  double GetYMax() const { return y_max_; }

  RegularSurface<A> ResampleSurface() const ;
  double GetAngle() const {return angle_;}

  A GetMissingValue() const { return surface_.GetMissingValue() ;}
  void SetMissingValue(A missing_val) { surface_.SetMissingValue(missing_val); }

// We do not want to use the following methods if possible
//  RegularSurface<A> GetSurface() const {return surface_;}
//  Grid2D<A> GetSurfaceGrid2D() const { return surface_;}

  const std::string& GetName() const                  { return surface_.GetName(); }
  void               SetName(const std::string& name) { surface_.SetName(name); }

  /// \brief Read surface file on given format.
  /// \param filename  File name.
  /// \param format    File format. If SURF_UNKNOWN we try to determine the format.
  /// \throws IOError         If we failed to open the file
  /// \throws FileFormatError If we can not determine the file format, or the contents
  ///                         of the file does not match the file format.
  void ReadFromFile(const std::string& filename,
                    SurfaceFileFormat  format = SURF_UNKNOWN);

  /// \brief Write surface to file on given format.
  /// If the file format does not support rotation, the resampled surface is written to file.
  /// \throws IOError  If we failed to open the file
  void WriteToFile(const std::string& filename,
                   SurfaceFileFormat  format = SURF_STORM_BINARY) const;

private:
 void GlobalToLocalCoord(double global_x,
                          double global_y,
                          double& local_x,
                          double& local_y) const;
  void LocalToGlobalCoord(double local_x,
                          double local_y,
                          double& global_x,
                          double& global_y) const;
  void CalculateMinMaxXY();
  RegularSurface<A> surface_;

  double angle_;
  double x_ref_;
  double y_ref_;
  double x_min_;
  double y_min_;
  double x_max_;
  double y_max_;

};

// ==================== TEMPLATE IMPLEMENTATIONS ====================

template <class A>
RegularSurfaceRotated<A>::RegularSurfaceRotated()
:angle_(0),
 x_ref_(0.0),
 y_ref_(0.0)
{
  surface_ = RegularSurface<A>();
  CalculateMinMaxXY();
}


template <class A>
RegularSurfaceRotated<A>::RegularSurfaceRotated(double x_min, double y_min,
                                  double lx, double ly,
                                  size_t nx, size_t ny, double angle,
                                  const A& value)
  : angle_(angle),
    x_ref_(x_min),
    y_ref_(y_min)
{
  surface_ = RegularSurface<A>(x_min,y_min,lx,ly,nx,ny,value);
  CalculateMinMaxXY();
}


template <class A>
RegularSurfaceRotated<A>::RegularSurfaceRotated(const std::string& filename,
                                                SurfaceFileFormat  format)
{
  ReadFromFile(filename, format);
}


template <class A>
void
RegularSurfaceRotated<A>::GlobalToLocalCoord(double global_x,
                                double global_y,
                                double& local_x,
                                double& local_y) const
{
  double x_rel = global_x-x_ref_;
  double y_rel = global_y - y_ref_;

  local_x =   std::cos(angle_)*x_rel + std::sin(angle_)*y_rel ;
  local_y = - std::sin(angle_)*x_rel + std::cos(angle_)*y_rel ;
}

template <class A>
void RegularSurfaceRotated<A>::LocalToGlobalCoord(double local_x,
                                     double local_y,
                                     double& global_x,
                                     double& global_y) const
{
  global_x = std::cos(angle_)*local_x - std::sin(angle_)*local_y + x_ref_;
  global_y = std::sin(angle_)*local_x + std::cos(angle_)*local_y + y_ref_;
}

template <class A>
A RegularSurfaceRotated<A>::GetZ(double x, double y) const
{
  double local_x, local_y;
  GlobalToLocalCoord(x, y, local_x, local_y);
  A z = surface_.GetZ(local_x+x_ref_, local_y+y_ref_);
  return z;
}
template <class A>
A RegularSurfaceRotated<A>::GetZInside(double x, double y) const
{
  double local_x, local_y;
  GlobalToLocalCoord(x, y, local_x, local_y);
  A z = surface_.GetZInside(local_x+x_ref_, local_y+y_ref_);
  return z;
}
template <class A>
bool RegularSurfaceRotated<A>::IsInsideSurface(double x, double y) const
{
  double local_x, local_y;
  GlobalToLocalCoord(x, y, local_x, local_y);
  bool z = surface_.IsInsideSurface(local_x+x_ref_, local_y+y_ref_);
  return z;

}
template <class A>
bool RegularSurfaceRotated<A>::EnclosesRectangle(double x_min, double x_max,
                         double y_min, double y_max) const
{
 if (x_min < GetXMin() || x_max > GetXMax() ||
      y_min < GetYMin() || y_max > GetYMax()) {
    return false;
  }

  return true;

}

template <class A>
void RegularSurfaceRotated<A>::CalculateMinMaxXY()
{
  std::vector<double> x(4);
  std::vector<double> y(4);

  LocalToGlobalCoord(0, 0, x[0], y[0]);
  LocalToGlobalCoord(surface_.GetLengthX(), 0, x[1], y[1]);
  LocalToGlobalCoord(0, surface_.GetLengthY(), x[2], y[2]);
  LocalToGlobalCoord(surface_.GetLengthX(), surface_.GetLengthY(), x[3], y[3]);

  x_min_ = *(std::min_element(x.begin(), x.end()));
  y_min_ = *(std::min_element(y.begin(), y.end()));
  x_max_ = *(std::max_element(x.begin(), x.end()));
  y_max_ = *(std::max_element(y.begin(), y.end()));

}

template <class A>
void RegularSurfaceRotated<A>::GetXY(size_t i, size_t j, double & x, double & y) const
{
  double xloc, yloc;
  surface_.GetXY(i,j,xloc, yloc);
  LocalToGlobalCoord(xloc,yloc,x,y);

}

template <class A>
void RegularSurfaceRotated<A>::GetXY(size_t index, double & x, double & y) const
{
  double xloc, yloc;
  surface_.GetXY(index,xloc, yloc);
  LocalToGlobalCoord(xloc,yloc,x,y);

}

template <class A>
RegularSurface<A> RegularSurfaceRotated<A>::ResampleSurface() const
{
  double dx = surface_.GetDX()*cos(angle_);
  double dy = surface_.GetDY()*cos(angle_);
  int nx = int((x_max_-x_min_)/dx) + 1;
  int ny = int((y_max_-y_min_)/dy) +1;
  RegularSurface<A> surf(x_min_, y_min_, x_max_-x_min_, y_max_-y_min_, nx, ny, 0.0);
  surf.SetMissingValue(GetMissingValue());
  A value;
  double x,y;
  int i,j;
  for(i=0;i<nx;i++)
    for(j=0;j<ny;j++)
    {
      x = x_min_ + i*dx;
      y = y_min_ + j*dy;
      value = GetZ(x,y);
      surf(i,j) = value;
    }

  return surf;
}


template <class A>
void RegularSurfaceRotated<A>::ReadFromFile(const std::string& filename,
                                            SurfaceFileFormat  format)
{
  if (format == SURF_UNKNOWN) {
    format = FindSurfaceFileType(filename);
    if (format == SURF_UNKNOWN) {
      throw FileFormatError("Failed to determine file format for surface file: " + filename);
    }
  }

  switch (format) {
    case SURF_IRAP_CLASSIC_ASCII:
      ReadIrapClassicAsciiSurf(filename, surface_, angle_);
      break;
    case SURF_STORM_BINARY:
      ReadStormBinarySurf(filename, surface_);
      angle_ = 0.0;
      break;
    case SURF_SGRI:
      ReadSgriSurf(filename, surface_, angle_);
      break;
    default:
      throw FileFormatError("Reading of file " + filename + " on format " + ToString(format)
                             + " as a rotated grid is currently not supported.");
  }

  x_ref_ = surface_.GetXMin();
  y_ref_ = surface_.GetYMin();

  CalculateMinMaxXY();
}


template <class A>
void RegularSurfaceRotated<A>::WriteToFile(const std::string& filename,
                                           SurfaceFileFormat  format) const
{
  switch (format) {
    case SURF_IRAP_CLASSIC_ASCII:
      WriteIrapClassicAsciiSurf(surface_, angle_, filename);
      break;
    case SURF_STORM_BINARY:
      {
        if (angle_ != 0.0) {
          RegularSurface<A> resampled_surf = ResampleSurface();
          WriteStormBinarySurf(resampled_surf, filename);
        }
        else {
          WriteStormBinarySurf(surface_, filename);
        }
      }
      break;
    default:
      throw FileFormatError("Writing of surface to file " + filename + " on format "
                            + ToString(format) + " is currently not supported.");
  }
}


}
#endif
