// $Id: bilinearsurface.cpp 882 2011-09-23 13:10:16Z perroe $

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

#include "bilinearsurface.hpp"
#include <cmath>

using namespace NRLib;

const double BilinearSurface::min_ = 0.001; // 0.00000001; Decreased instability...

BilinearSurface::BilinearSurface()
{
  // no default values defined here, use SetPtVector
}


// default constructor: four startpoints, direction pt00 -> pt10 -> pt11 -> pt01
BilinearSurface::BilinearSurface( const Point& pt00_in, const Point& pt10_in, const Point& pt11_in, const Point& pt01_in )
{
  pt00_ = pt00_in;
  pt10_ = pt10_in;
  pt11_ = pt11_in;
  pt01_ = pt01_in;

  // calculation of the three vectors defining the bilinear surface, togetgher with pt00_
  avec0_ = pt10_ - pt00_;
  avec1_ = pt01_ - pt00_;
  avec2_ = pt11_ - pt01_ - pt10_ + pt00_;
}


// one startpoint and three vectors defining the bilinearsurface. Also calculating the other cornerpoints (pt10_)
void BilinearSurface::SetPtVector( const Point& pt00_in, const Point& avec0_in, const Point& avec1_in, const Point& avec2_in )
{
  pt00_  = pt00_in;
  avec0_ = avec0_in;
  avec1_ = avec1_in;
  avec2_ = avec2_in;

  // calculation of the other cornerpoints from the input-vectors. (pt00_ is defined).
  pt10_ = pt00_ + avec0_;
  pt11_ = pt00_ + avec0_ + avec1_ + avec2_;
  pt01_ = pt00_ + avec1_;

  return;
}


// set local variables (pt00_) equal to given input points. Also calculating the local defining vectors (avec_).
void BilinearSurface::SetCornerPoints( const Point& pt00_in, const Point& pt10_in, const Point& pt11_in, const Point& pt01_in )
{
  pt00_ = pt00_in;
  pt10_ = pt10_in;
  pt11_ = pt11_in;
  pt01_ = pt01_in;

  // calculation of the three vectors defining the bilinear surface, togetgher with pt00_
  avec0_ = pt10_ - pt00_;
  avec1_ = pt01_ - pt00_;
  avec2_ = pt11_ - pt01_ - pt10_ + pt00_;

  return;
}


// returns point for given input index (u,v): p(u,v) = u*v*a + u*b + v*c + d
Point BilinearSurface::GetPoint( double u_in, double v_in ) const
{
// returns point corresponding to any (u,v)-value, not only inside [0,1]
  return ( pt00_ + u_in*v_in*avec2_ + v_in*avec1_ + u_in*avec0_ );
}


// printing the 4 defining cornerpoints
void BilinearSurface::PrintPoints() const
{
  std::cout << "Four cornerpoints defining the bilinear surface:\n" << "pt00 = " << pt00_
            << "    pt10 = " << pt10_ << "\npt11 = " << pt11_ << "    pt01 = " << pt01_ << "\n";
  return;
}


// printing the defining input-point and the other 3 def-vectors
void BilinearSurface::PrintVectors() const
{
  std::cout << "One startpoint and three vectors defining the bilinear surface:\n" << "pt00 = "
            << pt00_ << "\navec0 = " << avec0_ << "    avec1 = " << avec1_ << "    avec2 = "
            << avec2_ << "\n";
  return;
}


// returns true if pt_in is "above" bilinearsurface, false if "below"
bool BilinearSurface::IsAbove( const Point& pt_in ) const
{
  Point b = pt_in - pt00_;

  if ( avec2_.Norm() < min_ )   // then bilinearSurface is planar
    return( b.Dot( avec0_.Cross( avec1_ ) ) > min_ );
  else {
    Point  inv_avec0, inv_avec1, inv_avec2;
    double det_avec;
    NRLib::NRLibPrivate::Matrix3DInverse(avec0_, avec1_, avec2_,
                                         det_avec, inv_avec0, inv_avec1, inv_avec2);
    Point c0(inv_avec0.x, inv_avec1.x, inv_avec2.x);
    Point c1(inv_avec0.y, inv_avec1.y, inv_avec2.y);
    Point c2(inv_avec0.z, inv_avec1.z, inv_avec2.z);

    double inv_u = b.Dot(c0);
    double inv_v = b.Dot(c1);

    if ( ( inv_u >= 0.0 && inv_u <= 1.0 ) || ( inv_v >= 0.0 && inv_v <= 1.0 ) ) {
    //   std::cout << "No.1\t";
       return ((b.Dot(c2) - inv_u*inv_v )*det_avec > 0.0);
    }
    else if (inv_u < 0.0 && inv_v < 0.0) {
    //   std::cout << "No.2\t";
       return (b.Dot(c2) * det_avec > 0.0);
    }
    else if (inv_u < 0.0 && inv_v > 1.0) {
    //   std::cout << "No.3\t";
       return ((b.Dot(c2) - inv_u) * det_avec > 0.0);
    }
    else if (inv_u > 1.0 && inv_v < 0.0) {
    //   std::cout << "No.4\t";
       return ((b.Dot(c2) - inv_v) * det_avec > 0.0);
    }
    else {      // inv_u > 1 && inv_v > 1
    //   std::cout << "No.5\t";
       return ((b.Dot(c2) - inv_u - inv_v + 1.0) * det_avec > 0.0);
    }
  }   // end of if(inv_u/inv_v)-test
}


// calculate the solutions from the quad-eq:  (-b +/- sqrt(b^2-4*a*c)) / 2*a
int BilinearSurface::QuadraticSol( double A_in, double B_in, double C_in, double& sol1_out, double& sol2_out )
{
  double tmp = B_in*B_in - 4*A_in*C_in;
  if ( A_in != 0.0 && tmp >= 0.0 ) {   // then two solutions
    sol1_out = ( -B_in + std::sqrt(tmp) ) / ( 2*A_in );
    sol2_out = ( -B_in - std::sqrt(tmp) ) / ( 2*A_in );
    // std::cout << "QuadraticSol: v1 = " << sol1_out << "\tv2 = " << sol2_out << "\n";
    return 2;
  }
  else if ( ( A_in == 0.0 ) && ( B_in != 0.0 ) ) {  // linear equation and one solution (set in sol1_out)
    sol1_out = -C_in / B_in;
    sol2_out = 0.0;     // no valid value for sol2_out
    return 1;
  }
  else if ( ( B_in == 0.0 ) && ( -C_in/A_in > 0.0 ) ) {
    sol1_out = std::sqrt( -C_in/A_in );
    sol2_out = 0.0;     // no valid value for sol2_out
    return 1;
  }
  else { // no valid solutions.
    sol1_out = sol2_out = 0.0;
    return 0;
  }
}


// calculate t-value in: p_in = r_in + t*q_in
double BilinearSurface::Compute_T( const Point& r_in, const Point& q_in, const Point& p_in )
{
  double t;
  if ( ( std::abs( q_in.x ) >= std::abs( q_in.y ) ) &&  ( std::abs( q_in.x ) >= std::abs( q_in.z ) ) )
    t = ( p_in.x - r_in.x ) / q_in.x;
  else if ( std::abs( q_in.y ) >= std::abs( q_in.z ) )
    t = ( p_in.y - r_in.y ) / q_in.y;
  else
    t = ( p_in.z - r_in.z ) / q_in.z;
  return t;
}


// return-value is num of intersections, either (0,1,2). If num=2 then intersec1_out is point corresponding
// to smallest t-value>=0 (first t-point inside), the other intersec is intersec2_out.
// num=1: intersec1_out is only intersection point
int BilinearSurface::FindIntersections(const Line& line_in, Point& intersec1_out, Point& intersec2_out)
{
  intersec1_out.x = intersec1_out.y = intersec1_out.z = 0.0;   // define unvalid default-value
  intersec2_out.x = intersec2_out.y = intersec2_out.z = 0.0;   // define unvalid default-value

  // Any pt on surface is def by p = u*v*a + u*b + v*c + d. (u,v) in [0,1]:
  Point a( pt00_ - pt10_ - pt01_ + pt11_ );
  Point b( pt10_ - pt00_ );
  Point c( pt01_ - pt00_ );
  Point d( pt00_ );

  // Make sure that we do not handle planar surfaces wrong due to numerical instability>
  if (fabs(a.x/pt00_.x) < min_)
    a.x = 0;
  if (fabs(a.y/pt00_.y) < min_)
    a.y = 0;
  if (fabs(a.z/pt00_.z) < min_)
    a.z = 0;

  // Any pt along the ray: p = r + t*q. t >= 0. line_in given by two points: pt1=r (t=0) and pt2=r+q (t=1)
  Point r = line_in.GetPt1();
  Point q = line_in.GetPt2() - r;
  // std::cout << "r = " << r << "     q = " << q << "\n";
  Point p1, p2;

  double A1_loc = a.x*q.z - a.z*q.x;               //  std::cout << "A1 = " << A1_loc << "\n";
  double B1_loc = b.x*q.z - b.z*q.x;               //  std::cout << "B1 = " << B1_loc << "\n";
  double C1_loc = c.x*q.z - c.z*q.x;               //  std::cout << "C1 = " << C1_loc << "\n";
  double D1_loc = (d.x-r.x)*q.z - (d.z-r.z)*q.x;   //  std::cout << "D1 = " << D1_loc << "\n";
  double A2_loc = a.y*q.z - a.z*q.y;               //  std::cout << "A2 = " << A2_loc << "\n";
  double B2_loc = b.y*q.z - b.z*q.y;               //  std::cout << "B2 = " << B2_loc << "\n";
  double C2_loc = c.y*q.z - c.z*q.y;               //  std::cout << "C2 = " << C2_loc << "\n";
  double D2_loc = (d.y-r.y)*q.z - (d.z-r.z)*q.y;   //  std::cout << "D2 = " << D2_loc << "\n";

  // v-equation def by: A*v*v + B*v + C = 0:
  double A_loc = A2_loc*C1_loc - A1_loc*C2_loc;
  double B_loc = A2_loc*D1_loc - A1_loc*D2_loc + B2_loc*C1_loc - B1_loc*C2_loc;
  double C_loc = B2_loc*D1_loc - B1_loc*D2_loc;

  double v1, v2, u1, u2, t1, t2, aa, bb;
  t1 = t2 = -1.0;

  // std::cout << "QuadraticEq: A = " << A_loc << "    B = " << B_loc << "    C = " << C_loc << "\n";
  int num = QuadraticSol( A_loc, B_loc, C_loc, v1, v2 );

  // std::cout << "QuadraticSol: v1 = " << v1 << "    v2 = " << v2 << "\n";
  // std::cout << "Eq1-calc: A*v1*v1 + B*v1 + C = " << A_loc*v1*v1 + B_loc*v1 + C_loc << "\n";
  // std::cout << "Eq2-calc: A*v2*v2 + B*v2 + C = " << A_loc*v2*v2 + B_loc*v2 + C_loc << "\n";

  if ( num == 0 )  // no solutions from QuadraticEq. No intersections, return def-value intersec_out=0
    return 0;

  else if ( num == 1 ) { // 1 solution from QuadraticEq. v1 is calculated, but might be outside [0,1] (v2 is not calculated)
    aa = v1*A2_loc + B2_loc;
    bb = v1*( A2_loc - A1_loc ) + B2_loc - B1_loc;
    if ( std::abs(bb) >= std::abs(aa) ) // finding what eq to use
      u1 = ( v1*( C1_loc - C2_loc ) + D1_loc - D2_loc ) / bb;  // eq.13 from paper
    else
      u1 = ( -v1*C2_loc - D2_loc ) / aa;  // eq.11 from paper
    p1 = u1*v1*a + u1*b + v1*c + d;
    t1 = Compute_T( r, q, p1 ); // t-value corresponding to eq: p1 = r + t1*q
//    std::cout << "num = 1 : (u,v,t)     = ( " << u1 << ", " << v1 << ", " << t1 << " )\n"
//              << "          Intersec_pt = " << p1 << "\n";

    // checking if (u,v,t) are valid (i.e. on surface). (u,v) in [0,1] and t>=0:
    if ( v1 >= -0.001 && v1 <= 1.001 && u1 >= -0.05 && u1 <= 1.05 ) {
  //    std::cout << "Valid (u,v,t)-values : One intersection Line / BilinearSurface.\n";
      intersec1_out = p1;     // def value for intersec2_out=0
      return 1;
    }

    return 0;
   }  // end of if (num=1)

  else if ( num == 2 ) {  // both v1 and v2 been calculated, but might be outside [0,1]

    // calc corresponding to v1:
    aa = v1*A2_loc + B2_loc;
    bb = v1*( A2_loc - A1_loc ) + B2_loc - B1_loc;
    if ( std::abs(bb) >= std::abs(aa) ) // finding what eq to use
      u1 = ( v1*( C1_loc - C2_loc ) + D1_loc - D2_loc ) / bb;  // eq.13 from paper
    else
      u1 = ( -v1*C2_loc - D2_loc ) / aa;  // eq.11 from paper
    p1 = u1*v1*a + u1*b + v1*c + d;
    t1 = Compute_T( r, q, p1 ); // t-value corresponding to eq: p1 = r + t1*q

    // calc corresponding to v2:
    aa = v2*A2_loc + B2_loc;
    bb = v2*( A2_loc - A1_loc ) + B2_loc - B1_loc;
    if ( std::abs(bb) >= std::abs(aa) ) // finding what eq to use
      u2 = ( v2*( C1_loc - C2_loc ) + D1_loc - D2_loc ) / bb;  // eq.13 from paper
    else
      u2 = ( -v2*C2_loc - D2_loc ) / aa;  // eq.11 from paper
    p2 = u2*v2*a + u2*b + v2*c + d;
    t2 = Compute_T( r, q, p2 ); // t-value corresponding to eq: p2 = r + t2*q

//    std::cout << "num = 2 : (u1,v1,t1)   = ( " << u1 << ", " << v1 << ", " << t1 << " )\n"
//              << "          Intersec_pt1 = " << p1 << "\n";
//    std::cout << "          (u2,v2,t2)   = ( " << u2 << ", " << v2 << ", " << t2 << " )\n"
//              << "          Intersec_pt2 = " << p2 << "\n";

    // finding what intersec_pt (intersec1_out) to use (smallest t-value with t>=0)

    bool valid_1 = false;
    bool valid_2 = false;
    // Lower limits are extended a little bit to account for numerical instability.
    if ( v1 >= -0.001 && v1 <= 1.001 && u1 >= -0.05 && u1 <= 1.05 )
      valid_1 = true;
    if ( v2 >= -0.001 && v2 <= 1.001 && u2 >= -0.05 && u2 <= 1.05 )
      valid_2 = true;

    if ( valid_1 && valid_2 ) { // use intersec_pt corresponding to smallest t-value.
      if ( t1 <= t2 ) {
        intersec1_out = p1;   intersec2_out = p2;
      }
      else  {
        intersec1_out = p2;   intersec2_out = p1;
      }
      return 2;
    }
    else if ( valid_1 ) {
      intersec1_out = p1;   intersec2_out = p2;   // unvalid p2-value given to intersec_out
      return 1;
    }
    else if ( valid_2 ) {
      intersec1_out = p2;   intersec2_out = p1;   // unvalid p1-value given to intersec_out
      return 1;
    }
    // else unvalid (valid_1 and valid_2) : intersec1_out/intersec2_out is init-value =0

    return 0;
   }  // end of (num==2)

  // must have been unvalid num-value from QuadraticSol(), not inside (0,1,2):
   return (-1);

}
