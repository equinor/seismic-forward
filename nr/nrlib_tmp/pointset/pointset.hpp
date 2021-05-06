// $Id: pointset.hpp 1068 2012-09-18 11:21:53Z perroe $

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

#ifndef NRLIB_POINT_SET_H
#define NRLIB_POINT_SET_H

#include <cassert>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../geometry/point.hpp"
#include "../surface/regularsurface.hpp"

namespace NRLib {

class Polygon;

class PointSet {

public:
  enum FileFormat{
    UnknownFormat,
    RmsInternalPoints,  ///< RMS internal points format for tagged points.
    RoxarText,           ///< 3 columns: X, Y and Z. Delimiter: 999.0 999.0, 999.0
    JasonXYZ,            ///< 5 columns: IL, XL, X, Y and Z
    GeneralXYZ,          ///< 3 columns: X, Y and Z. No header. No delimiter.
    OpenWorksXYZ         ///< 3 columns: X, Y and Z. Header.
  };

  //Constructor :
  PointSet();
  explicit PointSet(const std::vector<NRLib::Point> & points);
  explicit PointSet(size_t n, const Point & init_val = Point());
  explicit PointSet(const std::string & filename,
                    bool                filter_inactive_horizon_points = false);
  explicit PointSet(const NRLib::RegularSurface<double> & storm_surface);

  /// Access operator
  inline NRLib::Point       & operator[](size_t i);

  /// Constant access operator.
  inline const NRLib::Point & operator[](size_t i) const;

  /// Add point to end of list. All parameters get default value.
  void AddPoint(const NRLib::Point& point);

  /// Number of points.
  size_t GetSize() const {return(points_.size());}

  void AddDiscreteParameter(const std::string &name, std::vector<int> &values);
  void AddContParameter(const std::string & name, std::vector<double> &values);
  void AddStringParameter(const std::string &name, std::vector<std::string> &values);

  void WriteToFile(const std::string& filename, FileFormat format) const;

  const std::vector<int>         & GetDiscreteParameter(const std::string& name) const;
  const std::vector<double>      & GetContParameter(const std::string& name)     const;
  const std::vector<std::string> & GetStringParameter(const std::string& name)   const;

  const std::map<std::string, std::vector<int> > & GetDiscreteParameters() const {return integer_tag_;}
  const std::map<std::string, std::vector<double> > & GetContParameters() const {return float_tag_;}
  const std::map<std::string, std::vector<std::string> > & GetStringParameters() const {return string_tag_;}

  const std::vector<NRLib::Point>& GetPoints() const { return points_;}

  void SetDiscreteParameter(const std::string& name,
                            size_t             point_nr,
                            int                new_parameter);


  static PointSet::FileFormat FindFileType(const std::string & filename);

  bool HasDiscreteParameter(const std::string& name) const;
  bool HasContParameter(const std::string& name)     const;
  bool HasStringParameter(const std::string& name)   const;

  void Resize(size_t new_size);

  void Erase(size_t index);

  /// Smallest convex polygon P for which each point in the PointSetSurface is either on the boundary of or in the interior of P (neglecting z-component, 2D algorithm, Graham Scan)
  Polygon GetConvexHull();

private:
  // ----------------- PRIVATE FUNCTIONS ---------------------------

  static bool CheckFileIsRoxarText(const std::string& filename);
  static bool CheckFileIsInternalPoints(const std::string& filename);
  static bool CheckFileIsJasonXYZ(const std::string& filename);
  /// Checks general xyz-file. Three columns, only data, no header.
  /// Obs, Roxar text files will also be accepted as general xyz-file.
  static bool CheckFileIsGeneralXYZ(const std::string& filename);
  static bool CheckFileIsOpenWorksXYZ(const std::string& filename);

  void ReadRoxarText(const std::string& filename);
  void ReadRMSInternalPoints(const std::string &filename, bool filter_inactive_horizon_points = false);
  void ReadJasonXYZ(const std::string& filename);
  void ReadGeneralXYZ(const std::string& filename);
  void ReadOpenWorksXYZ(const std::string& filename);

  void WriteRoxarText(const std::string& filename) const;
  void WriteInternalPoints(const std::string& filename) const;
  void WriteGeneralXYZ(const std::string& filename) const;

  static bool IsRoxarFileEnder(const std::vector<double> &xyz);
  bool HasTags() const;

  /// Function used by GetConvexHull to sort indexes from start to end according to angle wrt point p which has index p_index in points_ (neglecting z-component, 2D algorithm, MergeSort).
  void AngleSort(std::vector<int> &indexes, int p_index, int start, int end);

  /// Function used by GetConvexHull to loop the sorted index list (AngleSort) and for each pair of points with equal angle, the one closest to the point p is removed (neglecting z-component, 2D algorithm)
  void RemoveEqualAngles(std::vector<int> &indexes, int p_index);

  // --------------------- MEMBER VARIABLES -------------------------
  std::vector<NRLib::Point> points_;
  std::map<std::string, std::vector<int> > integer_tag_;
  std::map<std::string, std::vector<double> > float_tag_;
  std::map<std::string, std::vector<std::string> > string_tag_;

}; // end PointSet


NRLib::Point& PointSet::operator[](size_t i)
{
  assert(i < points_.size());
  return points_[i];
}


const NRLib::Point& PointSet::operator[](size_t i) const
{
  assert(i < points_.size());
  return points_[i];
}


} //namespace NRLib

#endif
