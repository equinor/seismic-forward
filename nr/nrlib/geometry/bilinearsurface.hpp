// $Id: bilinearsurface.hpp 921 2011-11-02 12:10:23Z perroe $

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

#ifndef NRLIB_GEOMETRY_BILINEARSURFACE
#define NRLIB_GEOMETRY_BILINEARSURFACE

#include <iostream>

#include "line.hpp"
#include "point.hpp"

namespace NRLib {

class BilinearSurface {

public:
  /// Default constructor.
  BilinearSurface();

  /// Constructor
  /// default: using four cornerpoints for definition. Start 00 -> 10 -> 11 -> 01. Is also setting the
  /// defining vectors from given input points.
  BilinearSurface( const Point& pt00_in, const Point& pt10_in, const Point& pt11_in, const Point& pt01_in );

  /// set one local startpoint (pt00_) and three vectors (avec0,avec1,avec2) eq to given input-points. Also calculating the
  /// other cornerpoints (pt10,pt11,pt01)
  void SetPtVector( const Point& pt00_in, const Point& avec0_in, const Point& avec1_in, const Point& avec2_in );

  /// set local variables (pt00_) equal to given input points. Also calculating the defining vectors (avec_).
  void SetCornerPoints( const Point& pt00_in, const Point& pt10_in, const Point& pt11_in, const Point& pt01_in );

  /// returns point for given input index (u,v): p(u,v) = u*v*a + u*b + v*c + d
  Point GetPoint( double u_in, double v_in ) const;

  /// Returns true if pt_in is "above" the bilinearsurface, false if "below"
  bool IsAbove( const Point& pt_in ) const;

  /// Return-value is num of intersections, either (0,1,2), between bs and line_in. num=2: then intersec1_out is point
  /// corresponding to smallest t-value>=0 (first t-point inside), the other intersection is intersec2_out.
  /// num=1: intersec1_out is only intersection point of interest
  /// \note This function is unstable, and can give wrong results on edges of bilinear surfaces. It also fails for
  ///       surfaces with the wrong orientation (flat patches facing in x direction, I think...).
  ///       Should be fixed so that we use the optimal 2D projection of the surface before finding the intersection.
  int FindIntersections(const Line& line_in, Point& intersec1_out, Point& intersec2_out);

  /// Print the 4 defining cornerpoints
  void PrintPoints() const;

  /// Print the defining input-point and the other 3 def-vectors
  void PrintVectors() const;

private:
  Point pt00_, pt10_, pt11_, pt01_;  // the four cornerpoints defining the bilinearsurface
  Point avec0_, avec1_, avec2_;      // three vectors defining the bilinearsurface together with pt00_ (startpoint)

  static const double min_;   // used instead of zero when doubles are compared.

  // calculates the solution(s) of the given quadratic-eq (sol1,sol2). returns number of solutions found: (0,1,2)
  int QuadraticSol( double A_in, double B_in, double C_in, double& sol1_out, double& sol2_out );
  // calculate t-value (double) in: p_in = r_in + t*q_in
  double Compute_T( const Point& r_in, const Point& q_in, const Point& p_in );
};  // class BilinearSurface

}  // end of using namespace NRLib

#endif
