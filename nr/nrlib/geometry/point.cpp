// $Id: point.cpp 1693 2017-09-15 10:33:19Z vegard $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// �  Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// �  Redistributions in binary form must reproduce the above copyright notice, this list of
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

#include "point.hpp"
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace NRLib;


Point::Point()
{
  x = y = z = 0.0;
}


Point::Point( double x_in, double y_in, double z_in )
{
  x = x_in;
  y = y_in;
  z = z_in;
}


Point Point::GetNormalized() const
{
  // If the norm is 0, return the point which should be (0, 0, 0)
  if ( Norm() == 0.0 )
    return (*this);
  else
    return ( *this ) / Norm();
}


void Point::Normalize()
{
  if ( Norm() != 0.0 )
    *this /= Norm();
}


double Point::GetAngle(const Point& p_in) const
{
  double loc;
  double norm = Norm()*p_in.Norm();
  if (norm > 0.0) {
    loc = this->Dot(p_in) / norm;
    return std::acos(loc);
  }
  else
    return 0.0;   // returning unvalid value
}


std::ostream & NRLib::operator<<( std::ostream& out, const Point& p )
{
  out  << std::fixed << std::setprecision(2) << std::right <<  "( "
       << std::setw(10) << p.x << ", "
       << std::setw(10) << p.y << ", "
       << std::setw(10) << p.z << " )";
  return out;
}


std::istream & NRLib::operator>>( std::istream& in, Point& p )
{
  in >> p.x >> p.y >> p.z;
  return in;
}


void NRLib::NRLibPrivate::Matrix3DInverse(const Point& v1, const Point& v2, const Point& v3,
                                          double& det, Point& s1, Point& s2, Point& s3)
{
  det = v3.Dot(v1.Cross(v2));  // det = v3 * (v1 x v2)

  assert(det != 0.0);

  double inverse_det = 1.0 / det;
  s1.x = ( v2.y*v3.z - v2.z*v3.y ) * inverse_det;
  s2.x = ( v2.z*v3.x - v2.x*v3.z ) * inverse_det;
  s3.x = ( v2.x*v3.y - v2.y*v3.x ) * inverse_det;
  s1.y = ( v1.z*v3.y - v1.y*v3.z ) * inverse_det;
  s2.y = ( v1.x*v3.z - v1.z*v3.x ) * inverse_det;
  s3.y = ( v1.y*v3.x - v1.x*v3.y ) * inverse_det;
  s1.z = ( v1.y*v2.z - v1.z*v2.y ) * inverse_det;
  s2.z = ( v1.z*v2.x - v1.x*v2.z ) * inverse_det;
  s3.z = ( v1.x*v2.y - v1.y*v2.x ) * inverse_det;
}





/// algorithm from http://geomalgorithms.com
double
Point::LeftFromLineXY(const Point & P0, const Point & P1) const
{
  return (P1.x - P0.x)*(y - P0.y) - (x - P0.x)*(P1.y - P0.y);
}


bool
Point::IsInList(const std::vector<Point> & point_list) const
{
  if (point_list.size() ==0){
    return false;
  }
  else{
    Point point = *this;
    size_t length = point_list.size();
    for (size_t i = 0; i < length; i++){
      if (point == point_list[i]) {
        return true;
      }
    }
    return false;
  }
}

double NRLib::Point::CalculatePointVectorLength(const std::vector<NRLib::Point>& points)
{
  if (points.size() < 2)
    return 0.0;
  else {
    double acc_length = 0.0;
    for (size_t i = 1; i < points.size(); i++) {
      acc_length += (points[i] - points[i - 1]).Norm();
    }
    return acc_length;
  }
}
