// $Id: point.hpp 1037 2012-06-28 12:55:48Z anner $

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

#ifndef NRLIB_GEOMETRY_POINT
#define NRLIB_GEOMETRY_POINT

#include <cassert>
#include <cmath>
#include <iostream>

namespace NRLib {

class Point {
public:
  double x, y, z;

  // constructors :
  Point();
  Point(double x_in, double y_in, double z_in);

  inline void Swap(Point& other);

  inline void SetPoint(double x_in, double y_in, double z_in)
  { x = x_in;  y = y_in;  z = z_in; }

  // Add/sub p_in to def.point, update/return def.point :
  inline Point operator+=(const Point& rhs);
  inline Point operator-=(const Point& rhs);

  /// Scale with t_in.
  inline Point operator*=(double t_in);

  inline Point operator-() const;

  /// Scale with 1/t_in.
  /// \throw Exception if t_in == 0.
  inline Point operator/=(double t_in);

  inline bool operator==(const Point& rhs) const;

  inline bool operator!=(const Point& rhs) const;

  inline bool IsEqual(const Point& rhs, double tol = 1.0e-6) const;

  /// Scalar product.
  inline double Dot(const Point& rhs) const;

  /// Cross product. (Vector product)
  inline Point Cross(const Point& rhs) const;

  /// Distance between p_in and def.point :
  inline double GetDistance(const Point& p_in) const;

  /// Distance between p_in and def.point in the xy-plane:
  inline double GetDistanceXY(const Point& p_in) const;

  /// Norm (length of vector).
  inline double Norm() const;

  /// Normalize vector. If *this is the zero-vector it will remain so.
  void Normalize();

  /// Return the normalized point.
  Point GetNormalized() const;

  // Calc angle between input point/vector (p_in) and def.point. Return 0.0 if p_in or def.pt equal 0
  double GetAngle(const Point& p_in) const;  // return radians
};  // end of Class Point


inline Point operator+(const Point& lhs, const Point& rhs);
inline Point operator-(const Point& lhs, const Point& rhs);

/// Multiply scalar t_in and p_in, return new point (t*p) (p*t)
inline Point operator*(double t_in, const Point& p_in);

inline Point operator*(const Point& p_in, double t_in);

/// Scalar product
inline double Dot(const Point& rhs, const Point& lhs);

/// Vector cross product
inline Point Cross(const Point& rhs, const Point& lhs);

/// Divide def.point with scalar (t_in), return new point (p/t)
inline Point operator/(const Point& p_in, double t_in);

/// Vector norm
inline double Norm(const Point& p);

/// Write point with std::cout, format: "(x, y, z)"
std::ostream & operator<<(std::ostream & out, const Point& p);

/// Read point with std::cin, format: "x y z"
std::istream & operator>>(std::istream & in, Point& p);


namespace NRLibPrivate {
  /// Calculate the invertmatrix (s1, s2, s3 and determinant) for the 3D input matrix (v-points)
  /// \throw Exception if the matrix V = (v1, v2, v3) is singular (Det(V) == 0).
  void Matrix3DInverse(const Point& v1, const Point& v2, const Point& v3,
                       double& det, Point& s1, Point& s2, Point& s3);
}

// ------------ INLINE FUNCTION DEFINITIONS -----------------------

// ------------ POINT CLASS ---------------

void Point::Swap(Point& other)
{
  std::swap(x, other.x);
  std::swap(y, other.y);
  std::swap(z, other.z);
}


Point Point::operator+=(const Point& rhs)
{
  x += rhs.x;
  y += rhs.y;
  z += rhs.z;
  return *this;
}


Point Point::operator-=(const Point& rhs)
{
  x -= rhs.x;
  y -= rhs.y;
  z -= rhs.z;
  return *this;
}


Point Point::operator-() const
{
  Point p;
  p.x = -this->x;
  p.y = -this->y;
  p.z = -this->z;
  return p;
}


Point Point::operator*=(double t_in)
{
  x *= t_in;
  y *= t_in;
  z *= t_in;
  return( *this );
}


Point Point::operator/=(double t_in)
{
  assert(t_in != 0.0);

  x /= t_in;
  y /= t_in;
  z /= t_in;
  return *this;
}


bool Point::operator==(const Point& rhs) const
{
  return (x == rhs.x && y == rhs.y && z == rhs.z);
}

bool Point::IsEqual(const Point& rhs, double tol) const
{
  double abs_rhs_x = std::fabs(rhs.x);
  double nTol_x = abs_rhs_x < 1.0 ? tol : abs_rhs_x * tol;
  double abs_rhs_y = std::fabs(rhs.y);
  double nTol_y = abs_rhs_y < 1.0 ? tol : abs_rhs_y * tol;
  double abs_rhs_z = std::fabs(rhs.z);
  double nTol_z = abs_rhs_z < 1.0 ? tol : abs_rhs_z * tol;
  return (std::fabs(x-rhs.x) <= nTol_x && std::fabs(y-rhs.y) <= nTol_y && std::fabs(z-rhs.z) <= nTol_z);
}

bool Point::operator!=(const Point& rhs) const
{
  return !(operator==(rhs));
}


double Point::Dot(const Point& rhs) const
{
  return(  x*rhs.x + y*rhs.y + z*rhs.z );
}


Point Point::Cross(const Point& rhs) const
{
  Point p_out;
  p_out.x = y*rhs.z - z*rhs.y;
  p_out.y = z*rhs.x - x*rhs.z;
  p_out.z = x*rhs.y - y*rhs.x;
  return p_out;
}


double Point::GetDistance( const Point& p_in ) const
{
  return (p_in - *this).Norm();
}


double Point::GetDistanceXY( const Point& p_in ) const
{
  double dx = x-p_in.x;
  double dy = y-p_in.y;
  return std::sqrt(dx*dx+dy*dy);
}


double Point::Norm() const
{
  return std::sqrt(x*x + y*y + z*z);
}

// ------------ POINT HELPER FUNCTIONS ---------------

Point operator+(const Point& lhs, const Point& rhs)
{
  Point p_out = lhs;
  p_out += rhs;
  return p_out;
}


Point operator-(const Point& lhs, const Point& rhs)
{
  Point p_out = lhs;
  p_out -= rhs;
  return p_out;
}


Point operator*(double t_in, const Point& p_in)
{
  Point p_out = p_in;
  p_out      *= t_in;
  return p_out;
}


Point operator*(const Point& p_in, double t_in)
{
  return t_in * p_in;
}


double Dot(const Point& rhs, const Point& lhs)
{
  return rhs.Dot(lhs);
}


Point Cross(const Point& rhs, const Point& lhs)
{
  return rhs.Cross(lhs);
}


Point operator/(const Point& p_in, double t_in)
{
  Point p_out = p_in;
  p_out /= t_in;
  return p_out;
}


double Norm(const Point& p)
{
  return p.Norm();
}


}  // end of namespace NRLib


#endif
