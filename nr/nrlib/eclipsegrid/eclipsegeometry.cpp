// $Id: eclipsegeometry.cpp 1136 2013-01-21 13:13:56Z georgsen $

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

#include "eclipsegeometry.hpp"
#include "eclipsetools.hpp"

#include "../exception/exception.hpp"
#include "../geometry/polygon.hpp"
#include "../geometry/triangle.hpp"
#include "../geometry/bilinearsurface.hpp"
#include "../iotools/fileio.hpp"
#include "../iotools/stringtools.hpp"
#include "../math/constants.hpp"

#include <stdexcept>
#include <fstream>
#include <iostream>
#include <ostream>
#include <istream>
#include <vector>

using namespace NRLib;

EclipseGeometry::EclipseGeometry()
 : ni_(0),
   nj_(0),
   nk_(0),
   pillars_(),
   active_pillars_(),
   z_corners_(),
   active_()
{ }


EclipseGeometry::EclipseGeometry(size_t ni, size_t nj, size_t nk)
 : ni_(ni),
   nj_(nj),
   nk_(nk),
   pillars_(ni + 1, nj + 1),
   active_pillars_(ni + 1, nj + 1),
   z_corners_(8 * ni * nj * nk),
   active_(ni, nj, nk)
{}


void EclipseGeometry::Resize(size_t ni, size_t nj, size_t nk)
{
  ni_ = ni;
  nj_ = nj;
  nk_ = nk;
  pillars_.Resize(ni + 1, nj + 1);
  active_pillars_.Resize(ni + 1, nj + 1);
  z_corners_.resize(8 * ni * nj * nk);
  active_.Resize(ni, nj, nk);
}


void EclipseGeometry::WriteSpecGrid( std::ofstream& out_file ) const
{
  out_file << "SPECGRID\n";
  out_file << "  " << ni_ << "  " << nj_ << "  " << nk_ << "   1  F  /\n\n";
}


void EclipseGeometry::WriteCoord(std::ofstream& out_file) const
{
  out_file.precision(3);      // this value might be changed
  out_file.setf(std::ios_base::fixed);

  // (nx+1)*(ny+1) lines for each corner, two points: p1(upper), p2(lower)
  out_file << "COORD\n";
  for (size_t j=0; j <= nj_; j++)
    for (size_t i=0; i <= ni_; i++)
      out_file << "   " << pillars_(i,j).GetPt1().x << "   "
                        << pillars_(i,j).GetPt1().y << "   "
                        << pillars_(i,j).GetPt1().z << "   "
                        << pillars_(i,j).GetPt2().x << "   "
                        << pillars_(i,j).GetPt2().y << "   "
                        << pillars_(i,j).GetPt2().z << "\n";
  out_file << "  /\n\n";
}


void EclipseGeometry::WriteZCorn( std::ofstream& out_file ) const
{
  out_file.precision(3);      // this value might be changed
  out_file.setf(std::ios_base::fixed);

  out_file << "ZCORN\n";

  WriteAsciiArray(out_file, z_corners_.begin(), z_corners_.end());

  out_file << "\n  /\n\n";
}


void EclipseGeometry::WriteActNum(std::ofstream& out_file) const
{
  out_file << "ACTNUM\n";
  WriteAsciiArray(out_file, active_.begin(), active_.end(), 40);
  out_file << "  /\n\n";
}


void EclipseGeometry::WriteGeometry(std::ofstream& out_file) const
{
  WriteSpecGrid(out_file);
  WriteCoord(out_file);
  WriteZCorn(out_file);
  WriteActNum(out_file);
}


Point EclipseGeometry::FindCellCenterPoint(size_t i, size_t j, size_t k) const
{
  Point p = ( FindCornerPoint(i, j, k, 0,0,0) + FindCornerPoint(i, j, k, 1,0,0) +
              FindCornerPoint(i, j, k, 0,1,0) + FindCornerPoint(i, j, k, 1,1,0) +
              FindCornerPoint(i, j, k, 0,0,1) + FindCornerPoint(i, j, k, 1,0,1) +
              FindCornerPoint(i, j, k, 0,1,1) + FindCornerPoint(i, j, k, 1,1,1) ) / 8.0;
  return p;
}

Point EclipseGeometry::FindCellCenterPointTopBot(size_t i, size_t j, size_t k, bool top) const
{
  int c;
  if(top == true)
    c = 0;
  else c = 1;

  Point p = ( FindCornerPoint(i, j, k, 0,0,c) + FindCornerPoint(i, j, k, 1,0,c) +
              FindCornerPoint(i, j, k, 0,1,c) + FindCornerPoint(i, j, k, 1,1,c) ) / 4.0;
  return p;
}



// find the cell (i,j) the given point (x,y,z) is inside. Returns true if point is inside the grid, and has been found.
bool EclipseGeometry::GetCellIndex_ij(double x_in, double y_in, double z_in, size_t& i_out, size_t& j_out) const
{

  Point  p_in(x_in, y_in, z_in);
 // bool inside = largepoly.IsInsidePolygonXY(p_in);
  bool inside1 = polymin_.IsInsidePolygonXY(p_in);
  bool inside2 = polymax_.IsInsidePolygonXY(p_in);
  if(inside1==false && inside2==false)
    return false;

 // Generate check_plane, defined by given input point and any normal vector, here (0,0,1)

  Point  nvec_in( 0,0,1 );
  Plane  checkPlane( p_in, nvec_in );

  Grid2D<Point> interSecTab(ni_+1, nj_+1);      // intersection-points: pillars_ and checkPlane
  Grid2D<int>   pillar_ok(ni_+1, nj_+1, true);
  size_t i, j;
  for (j=0; j <= nj_; j++)
    for (i=0; i <= ni_; i++) {
      if (active_pillars_(i, j)){
      try {
        interSecTab(i, j) = checkPlane.FindIntersection(pillars_(i, j));
      }
      catch (NRLib::Exception& ) {
        // Degenerated pillar, or pillar parallel to plane.
        pillar_ok(i, j) = false;
      }
      }
      else
        pillar_ok(i, j) = false;
    } // end i-loop

  Polygon poly(4);
  bool is_inside = false;
  for ( j=0; j < nj_; j++ ) {
    for ( i=0; i < ni_; i++ ) {
      if ( pillar_ok( i, j ) &&  pillar_ok( i+1, j ) && pillar_ok( i+1, j+1 ) && pillar_ok( i, j+1 ) ) {
        poly(0)   = interSecTab( i, j );      // = pt1. side1: pt1->pt2
        poly(1)   = interSecTab( i+1, j );    // = pt2. side2: pt2->pt3
        poly(2)   = interSecTab( i+1, j+1 );  // = pt3. side3: pt3->pt4
        poly(3)   = interSecTab( i, j+1 );    // = pt4. side4: pt4->pt1
        is_inside = poly.IsInsidePolygonXY( p_in );

        if ( is_inside ) {
          i_out = i;
          j_out = j;
          // cout << "GetCellIndex_ij : (i_out,j_out) = (" << i_out << "," << j_out << ")\n";
          return true;
         }
      }
    }
  }

  return false;
}


Polygon EclipseGeometry::FindPolygonAroundActivePillars(double z_in)
{
// printf(" find polygon around active pillars\n");
  std::vector<std::vector<Point>*> polygons;
  std::vector<Point>* p;
  bool continuing=true;
  size_t checked_to_i=0;
  size_t checked_to_j=0;
  size_t last_index;
  size_t choice,i,j;
  NRLib::Grid2D<bool> checked(ni_+1,nj_+1,false);
  while (continuing) {
    continuing=false;
    for (i=checked_to_i; i<=ni_; i++) {
      last_index=0;
      for (j=checked_to_j; j<=nj_; j++) {
        if (active_pillars_(i,j)!=0 && checked(i,j)==false) {
          if (j==0 || active_pillars_(i,last_index)==0) {
            checked_to_i=i;
            checked_to_j=j+1;
            continuing=true;
            i=ni_+1;
            j=nj_+1;
          }
        }
        last_index=j;
      }
    }
    if (continuing) {
      p=new std::vector<Point>;
      checked(checked_to_i, checked_to_j-1)=true;
      p->push_back(FindPointAtPillar(checked_to_i, checked_to_j-1,z_in));
      choice=1;
      i=checked_to_i;
      j=checked_to_j-2;
      while (choice!=0) { //choice holds information about where to search in next step: 1<=>Up, 2<=>Right, 3<=>Down, 4<=>Left
        if (choice==1) {
          j++;
          choice=SearchUp(i,j,*p,checked,z_in,checked_to_i,checked_to_j-1);
        }
        if (choice==2) {
          i++;
          choice=SearchRight(i,j,*p,checked,z_in,checked_to_i,checked_to_j-1);
        }
        if (choice==3) {
          j--;
          choice=SearchDown(i,j,*p,checked,z_in,checked_to_i,checked_to_j-1);
        }
        if (choice==4) {
          i--;
          choice=SearchLeft(i,j,*p,checked,z_in,checked_to_i,checked_to_j-1);
        }
      }
      polygons.push_back(p);
    }
  }
  Polygon active_polygon;
  Point vec1, vec2, start_point;
  for (j=0; j<polygons.size();j++) {
    p=polygons[j];
    active_polygon.AddPoint((*p)[0]);
    for (i=1; i<(p->size()-1); i++) {
      vec1=(*p)[i]-(*p)[i-1];
      vec2=(*p)[i+1]-(*p)[i];
      if (vec1.x*vec2.y-vec1.y*vec2.x!=0.0) //
        active_polygon.AddPoint((*p)[i]);
    }
    active_polygon.AddPoint((*p)[p->size()-1]);
    if (j>0) {
      active_polygon.AddPoint(start_point);
    }
    else {
      start_point=(*p)[0];
    }
  }
  return active_polygon;
}

size_t EclipseGeometry::SearchUp(size_t i,
                                     size_t j,
                     std::vector<Point> &p,
                     NRLib::Grid2D<bool> &checked,
                     double z_in,
                     size_t start_i,
                     size_t start_j)
{
  bool searched=false;
  if (i>0) {
    if (active_pillars_(i-1,j)!=0) {
      searched=true;
      checked(i-1,j)=true;
      p.push_back(FindPointAtPillar(i-1,j,z_in));
      if (!(start_i==i-1 && start_j==j)) {
        return 4;
      }
    }
  }
  if (j<nj_ && searched==false) {
    if (active_pillars_(i,j+1)!=0) {
      searched=true;
      checked(i,j+1)=true;
      p.push_back(FindPointAtPillar(i,j+1,z_in));
      if (!(start_i==i && start_j==j+1)) {
        return 1;
      }
    }
  }
  if (i<ni_ && searched==false) {
    if (active_pillars_(i+1,j)!=0) {
      searched=true;
      checked(i+1,j)=true;
      p.push_back(FindPointAtPillar(i+1,j,z_in));
      if (!(start_i==i+1 && start_j==j)) {
        return 2;
      }
    }
  }
  return 0;
}

size_t EclipseGeometry::SearchDown(size_t i,
                                       size_t j,
                     std::vector<Point> &p,
                     NRLib::Grid2D<bool> &checked,
                     double z_in,
                     size_t start_i,
                     size_t start_j)
{
  bool searched=false;
  if (i<ni_ && searched==false) {
    if (active_pillars_(i+1,j)!=0) {
      searched=true;
      checked(i+1,j)=true;
      p.push_back(FindPointAtPillar(i+1,j,z_in));
      if (!(start_i==i+1 && start_j==j)) {
        return 2;
      }
    }
  }
  if (j>0 && searched==false) {
    if (active_pillars_(i,j-1)!=0) {
      searched=true;
      checked(i,j-1)=true;
      p.push_back(FindPointAtPillar(i,j-1,z_in));
      if (!(start_i==i && start_j==j-1)) {
        return 3;
      }
    }
  }
  if (i>0 && searched==false) {
    if (active_pillars_(i-1,j)!=0) {
      searched=true;
      checked(i-1,j)=true;
      p.push_back(FindPointAtPillar(i-1,j,z_in));
      if (!(start_i==i-1 && start_j==j)) {
        return 4;
      }
    }
  }
  return 0;
}


size_t EclipseGeometry::SearchLeft(size_t i,
                                       size_t j,
                     std::vector<Point> &p,
                     NRLib::Grid2D<bool> &checked,
                     double z_in,
                     size_t start_i,
                     size_t start_j)
{
  bool searched=false;
  if (j>0) {
    if (active_pillars_(i,j-1)!=0) {
      searched=true;
      checked(i,j-1)=true;
      p.push_back(FindPointAtPillar(i,j-1,z_in));
      if (!(start_i==i && start_j==j-1)) {
        return 3;
      }
    }
  }
  if (i>0 && searched==false) {
    if (active_pillars_(i-1,j)!=0) {
      searched=true;
      checked(i-1,j)=true;
      p.push_back(FindPointAtPillar(i-1,j,z_in));
      if (!(start_i==i-1 && start_j==j)) {
        return 4;
      }
    }
  }
  if (j<nj_ && searched==false) {
    if (active_pillars_(i,j+1)!=0) {
      searched=true;
      checked(i,j+1)=true;
      p.push_back(FindPointAtPillar(i,j+1,z_in));
      if (!(start_i==i && start_j==j+1)) {
        return 1;
      }
    }
  }
  return 0;
}

size_t EclipseGeometry::SearchRight(size_t i,
                                        size_t j,
                      std::vector<Point> &p,
                      NRLib::Grid2D<bool> &checked,
                      double z_in,
                      size_t start_i,
                      size_t start_j)
{
  bool searched=false;
  if (j<nj_) {
    if (active_pillars_(i,j+1)!=0) {
      searched=true;
      checked(i,j+1)=true;
      p.push_back(FindPointAtPillar(i,j+1,z_in));
      if (!(start_i==i && start_j==j+1)) {
        return 1;
      }
    }
  }
  if (i<ni_ && searched==false) {
    if (active_pillars_(i+1,j)!=0) {
      searched=true;
      checked(i+1,j)=true;
      p.push_back(FindPointAtPillar(i+1,j,z_in));
      if (!(start_i==i+1 && start_j==j)) {
        return 2;
      }
    }
  }
  if (j>0 && searched==false) {
    if (active_pillars_(i,j-1)!=0) {
      searched=true;
      checked(i,j-1)=true;
      p.push_back(FindPointAtPillar(i,j-1,z_in));
      if (!(start_i==i && start_j==j-1)) {
        return 3;
      }
    }
  }
  return 0;
}



bool EclipseGeometry::FindIndex(double x_in, double y_in, double z_in,
                                size_t& i_out, size_t& j_out, size_t& k_out) const
{
  bool is_inside_ij = false;
  Point p_in( x_in, y_in, z_in );

  k_out        = nk_+1;  // setting default unvalid value
  is_inside_ij = GetCellIndex_ij( x_in, y_in, z_in, i_out, j_out );

  if ( is_inside_ij ) {
    double u, v;
    FindUVCoordinates(x_in, y_in, z_in, i_out, j_out, u, v);

    bool above_top =false;
    bool above_bot = false;
    size_t newbot;
    size_t top = 0;
    size_t bot = nk_-1;
    if (z_in < FindPointCellSurface(i_out, j_out, top, 0, u, v).z )
      above_top = true;
    if (z_in < FindPointCellSurface(i_out, j_out, bot, 1, u, v).z )
      above_bot = true;
    if(above_top == true || above_bot == false)
      return false;
    while(top < bot-1){
        newbot = static_cast<size_t>(0.5*(top+bot));
        if (z_in < FindPointCellSurface(i_out, j_out, newbot, 1, u, v).z ){
          bot = newbot;
        }
        else if (z_in > FindPointCellSurface(i_out, j_out, newbot, 0, u, v).z )
          top = newbot;
        else {
          k_out = newbot;
          return true;
        }
    }
    above_top =false;
    above_bot = false;
    if (z_in < FindPointCellSurface(i_out, j_out, top, 0, u, v).z )
      above_top = true;
    if (z_in < FindPointCellSurface(i_out, j_out, top, 1, u, v).z )
      above_bot = true;
    if (above_top != above_bot) {
      // then point is inside cell
      k_out       = top;
      return true;
    }
    above_top =false;
    above_bot = false;
     if (z_in < FindPointCellSurface(i_out, j_out, bot, 0, u, v).z )
      above_top = true;
    if (z_in < FindPointCellSurface(i_out, j_out, bot, 1, u, v).z )
      above_bot = true;
    if (above_top != above_bot) {
      // then point is inside cell
      k_out       = bot;
      return true;
    }
  }

 /*   for ( size_t k = 0; k < nk_; k++ )  {
      bool above_top = false, above_bot = false;

      // Points on the top surface belongs to the cell, while points on the bottom surface don't
      if (z_in < FindPointCellSurface(i_out, j_out, k, 0, u, v).z )
        above_top = true;
      if (z_in < FindPointCellSurface(i_out, j_out, k, 1, u, v).z )
        above_bot = true;

      if (above_top != above_bot) {
        // then point is inside cell
        k_out       = k;
        return true;
      }
    }   // end of k-loop
  }     // end of if (is_inside_ij)-test
 */
  return false;
}


Point EclipseGeometry::FindPointAtPillar(size_t i, size_t j, double z) const
{
  Point p1 = GetPillar(i, j).GetPt1();
  Point p2 = GetPillar(i, j).GetPt2();
  Point pillar_vector = p2 - p1;  //vector along the pillar

  if (pillar_vector.z == 0.0) {
    // Degenerated/horizontal pillar.
    assert(p1 == p2);

    Point p_out;
    p_out.x = p1.x;
    p_out.y = p1.y;
    p_out.z = z;
    return p_out;
  }

  double t = (z - p1.z) / pillar_vector.z;
  Point p_out;
  p_out.x = p1.x + t*pillar_vector.x;
  p_out.y = p1.y + t*pillar_vector.y;
  p_out.z = z;

  return p_out;
}

NRLib::Point
EclipseGeometry::FindPointAtPillarInsideGrid(size_t i, size_t j, double z, bool & found) const
{
  double z_pillar_top = FindZTopAtPillar(i,j, found);
  double new_z(z);
  if (z < z_pillar_top)
     new_z = z_pillar_top;
  else {
    double z_pillar_bot = FindZBotAtPillar(i,j,found);
    if (z > z_pillar_bot)
      new_z = z_pillar_bot;
  }
  Point pt;
  if (found)
    pt = FindPointAtPillar(i,j,new_z);

  return(pt);
}

double
EclipseGeometry::FindZTopAtPillar(size_t i, size_t j, bool & found) const
{
  double z_top = numeric_limits<double>::infinity();
  if (i > 0 && j > 0) {
    size_t k = FindTopCell(i-1, j-1);
    if (k < GetNK()) {
      double z = GetZCorner(i-1, j-1, k, 1, 1, 0);
      if (z < z_top)
        z_top = z;
    }
  }

  if (i > 0 && j < GetNJ()) {
    size_t k = FindTopCell(i-1, j);
    if (k < GetNK()) {
      double z = GetZCorner(i-1, j, k, 1, 0, 0);
      if (z < z_top)
        z_top = z;
    }
  }

  if (i < GetNI() && j > 0) {
    size_t k = FindTopCell(i, j-1);
    if (k < GetNK()) {
      double z = GetZCorner(i, j-1, k, 0, 1, 0);
      if (z < z_top)
        z_top = z;
    }
  }

  if (i < GetNI() && j < GetNJ()) {
    size_t k = FindTopCell(i, j);
    if (k < GetNK()) {
      double z = GetZCorner(i, j, k, 0, 0, 0);
      if (z < z_top)
        z_top = z;
    }
  }

  if (z_top < numeric_limits<double>::infinity())
    found = true;
  else
    found = false;

  return(z_top);
}

double
EclipseGeometry::FindZBotAtPillar(size_t i, size_t j, bool & found) const
{
  double z_bot = -numeric_limits<double>::infinity();
  if (i > 0 && j > 0) {
    size_t k = FindBottomCell(i-1, j-1);
    if (k < GetNK()) {
      double z = GetZCorner(i-1, j-1, k, 1, 1, 1);
      if (z > z_bot)
        z_bot = z;
    }
  }

  if (i > 0 && j < GetNJ()) {
    size_t k = FindBottomCell(i-1, j);
    if (k < GetNK()) {
      double z = GetZCorner(i-1, j, k, 1, 0, 1);
      if (z > z_bot)
        z_bot = z;
    }
  }

  if (i < GetNI() && j > 0) {
    size_t k = FindBottomCell(i, j-1);
    if (k < GetNK()) {
      double z = GetZCorner(i, j-1, k, 0, 1, 1);
      if (z > z_bot)
        z_bot = z;
    }
  }

  if (i < GetNI() && j < GetNJ()) {
    size_t k = FindBottomCell(i, j);
    if (k < GetNK()) {
      double z = GetZCorner(i, j, k, 0, 0, 1);
      if (z > z_bot)
        z_bot = z;
    }
  }

  if (z_bot > -numeric_limits<double>::infinity())
    found = true;
  else
    found = false;

  return(z_bot);
}

Point EclipseGeometry::FindPointCellSurface(size_t i_in, size_t j_in, size_t k_in, int lower_or_upper, double u_in, double v_in) const
{
  //The p values are the corner points on the cell, while q and r are points with the same z-value as the p-point
  // at respectively the neighbouring pillar in i-direction and neighbouring pillar in j-direction.
  Point p00 = FindPointAtPillar(i_in, j_in, GetZCorner(i_in, j_in, k_in, 0,0,lower_or_upper) );
  Point q00 = FindPointAtPillar(i_in, j_in, GetZCorner(i_in, j_in, k_in, 1,0,lower_or_upper) );
  Point r00 = FindPointAtPillar(i_in, j_in, GetZCorner(i_in, j_in, k_in, 0,1,lower_or_upper) );


  Point p10 = FindPointAtPillar(i_in + 1, j_in, GetZCorner(i_in, j_in, k_in, 1,0,lower_or_upper) );
  Point q10 = FindPointAtPillar(i_in + 1, j_in, GetZCorner(i_in, j_in, k_in, 0,0,lower_or_upper) );
  Point r10 = FindPointAtPillar(i_in + 1, j_in, GetZCorner(i_in, j_in, k_in, 1,1,lower_or_upper) );


  Point p01 = FindPointAtPillar(i_in, j_in + 1, GetZCorner(i_in, j_in, k_in, 0,1,lower_or_upper) );
  Point q01 = FindPointAtPillar(i_in, j_in + 1, GetZCorner(i_in, j_in, k_in, 1,1,lower_or_upper) );
  Point r01 = FindPointAtPillar(i_in, j_in + 1, GetZCorner(i_in, j_in, k_in, 0,0,lower_or_upper) );


  Point p11 = FindPointAtPillar(i_in + 1, j_in + 1, GetZCorner(i_in, j_in, k_in, 1,1,lower_or_upper) );
  Point q11 = FindPointAtPillar(i_in + 1, j_in + 1, GetZCorner(i_in, j_in, k_in, 0,1,lower_or_upper) );
  Point r11 = FindPointAtPillar(i_in + 1, j_in + 1, GetZCorner(i_in, j_in, k_in, 1,0,lower_or_upper) );


  //Coefficients for the model f(u,v) = (a2*v^2 + a1*v + a0)u^2 + (a1*v^2 + b1*v + c1)u + a0*v^2 + b0*v + c0
  Point a2 = (p01 + p11 - q01 - q11) -(p00 + p10 - q00 - q10) -(r10 + r11 - 2*p10);
  Point b2 = r10 + r11 - 2*p10;
  Point c2 = p00 + p10 - q00 - q10;

  Point a1 = (p00 + p10 - q00 - q10) + (p10 + p11 - r10 - r11) - (p01 + p11 - q01 - q11)
             - (p00 + p01 - r00 - r01) + (r10 + r11 - 2*p10);
  Point b1 = -1*(r00 + r01 - 2*p00);
  Point c1 = q00 + q10 - 2*p00;

  Point a0 = p00 + p01 - r00 - r01;
  Point b0 = r00 + r01 - 2*p00;
  Point c0 = p00;

  Point f = (a2 * v_in*v_in + b2*v_in + c2) * u_in*u_in
            + (a1 * v_in*v_in + b1*v_in + c1) * u_in
            + a0 * v_in*v_in + b0*v_in + c0;

  return f;
}

Point EclipseGeometry::FindPointInCell(size_t i_in, size_t j_in, size_t k_in, double u_in, double v_in, double w_in)
{
  Point cell_point = (1-w_in) * FindPointCellSurface(i_in, j_in, k_in, 0, u_in, v_in)
    + w_in * FindPointCellSurface(i_in, j_in, k_in, 1, u_in, v_in);

  return cell_point;
}

void EclipseGeometry::FindUVCoordinates(double x_in, double y_in, double z_in, size_t i_out, size_t j_out, double& u_out, double& v_out) const
{
  Point p00, p10, p01, p11;
  p00 = FindPointAtPillar(i_out, j_out, z_in);
  p10 = FindPointAtPillar(i_out + 1, j_out, z_in);
  p01 = FindPointAtPillar(i_out, j_out + 1, z_in);
  p11 = FindPointAtPillar(i_out + 1, j_out + 1, z_in);

  double ax, ay, bx, by, cx, cy, dx, dy;

  ax = p00.x + p11.x - p10.x - p01.x;
  ay = p00.y + p11.y - p10.y - p01.y;
  bx = p10.x - p00.x;
  by = p10.y - p00.y;
  cx = p01.x - p00.x;
  cy = p01.y - p00.y;
  dx = p00.x;
  dy = p00.y;

  double seca, secb, secc;


  if (ax == 0 && bx == 0){
    if (cx != 0){
      v_out = (x_in - dx)/cx;
    }
    else{
      throw Exception( "Degenerated cell" );
    }
    if (ay != 0 || by != 0) {
      if (ay*v_out + by != 0){
        u_out = (y_in - cy*v_out - dy)/(ay*v_out + by);
      }
      else{
        u_out = 0; //free parameter
      }
    }
    else {
      // ay == 0 && by == 0
      throw Exception( "Degenerated cell" );
    }
  }
  else{
    //seca*v^2 + secb*v + secc = 0
    seca = ax*cy - ay*cx;
    secb = dy*ax + bx*cy + x_in*ay - dx*ay - by*cx - ax*y_in;
    secc = dy*bx + x_in*by - by*dx - y_in*bx;

    if (seca != 0){
      // The expression with + in front of the square root seems to be the correct one
      v_out = (-secb + sqrt(secb*secb - 4*seca*secc))/(2*seca);
      if (v_out <0 || v_out>1){
        v_out = (-secb - sqrt(secb*secb - 4*seca*secc))/(2*seca);
      }
    }
    else{
      if (secb != 0){
        v_out = (-secc)/secb;
      }
      else{
        if (secc == 0){
          v_out = 0; //free parameter
        }
        else{
          throw Exception( "Degenerated cell" );
        }
      }
    }
    if (ax*v_out + bx != 0){
      u_out = (x_in - cx*v_out - dx)/(ax*v_out + bx);
    }
    else{
      u_out = 0; //free parameter
    }
  }
}

void EclipseGeometry::ReadSpecGrid(std::ifstream& in_file)
{
  //the line number is not used here, so it is just set as 0.
  int line = 0;
  ni_ = ReadNext<int>(in_file, line);
  nj_ = ReadNext<int>(in_file, line);
  nk_ = ReadNext<int>(in_file, line);

  pillars_.Resize(ni_+1, nj_+1);         // local 2D-grid of class Line
  z_corners_.resize(8 * ni_ * nj_ * nk_);   // local 3D-grid, each cell with local 3D-grid (loc_grid)
  active_.Resize(ni_, nj_, nk_);         // local 3D-grid of booleans (active cells)
  active_pillars_.Resize(ni_ + 1, nj_ + 1);

  ReadNext<int>( in_file, line );
  ReadNext<std::string>( in_file, line );
  ReadNext<std::string>( in_file, line );
}


void EclipseGeometry::ReadZCorn(std::ifstream& in_file)
{
  std::string buffer = ReadParameterBuffer(in_file);
  ParseAsciiArrayFast(buffer, z_corners_.begin(), 8*ni_*nj_*nk_);
}


void EclipseGeometry::ReadCoord(std::ifstream& in_file)
{
  std::string buffer = ReadParameterBuffer(in_file);
  std::vector<double> data( 6*(ni_ + 1)*(nj_ + 1) );
  ParseAsciiArrayFast( buffer, data.begin(),6*(ni_ + 1)*(nj_ + 1) );
  size_t i, j;
  Point pt1, pt2;
  for (j = 0; j <= nj_; j++){
    for (i = 0; i <= ni_; i++){
      pt1.x = data[6*i + 6*j*(ni_+1)];
      pt1.y = data[6*i + 1 + 6*j*(ni_+1)];
      pt1.z = data[6*i + 2 + 6*j*(ni_+1)];
      pt2.x = data[6*i + 3 + 6*j*(ni_+1)];
      pt2.y = data[6*i + 4 + 6*j*(ni_+1)];
      pt2.z = data[6*i + 5 + 6*j*(ni_+1)];
      pillars_(i,j).SetPt(pt1, pt2, true, true);
    }
  }
}


void EclipseGeometry::ReadActNum(std::ifstream& in_file)
{
  std::string buffer = ReadParameterBuffer(in_file);

  ParseAsciiArrayFast( buffer, active_.begin(), ni_*nj_*nk_ );

  InitializeActivePillars();
}


void EclipseGeometry::InitializeActivePillars()
{
  active_pillars_.Resize(GetNI() + 1, GetNJ() + 1, false);
  for(size_t i = 0; i < ni_; i++) {
    for(size_t j = 0; j < nj_; j++){
      if (IsColumnActive(i, j)) {
        active_pillars_(i  , j  ) = true;
        active_pillars_(i+1, j  ) = true;
        active_pillars_(i  , j+1) = true;
        active_pillars_(i+1, j+1) = true;
      }
    }
  }
}


void EclipseGeometry::FindMinAndMaxZValueAndSetPolygons()
{
  //printf("find min and max z\n");
  size_t i, j, k;
  Point pt, botpt, toppt;
  toppt.z = -999.0;
  botpt.z = -999.0;
  for(i = 0; i <= ni_; i++)
    for(j = 0; j <= nj_; j++){
      for(k = 0; k < nk_; k++){
        if(i < ni_ && j < nj_){
          if(IsActive(i,j,k)){
            pt = FindCornerPoint(i,j,k,0,0,0);
            if(botpt.z == -999.0 || pt.z > botpt.z)
              botpt = pt;
            if(toppt.z == -999.0 || pt.z < toppt.z)
              toppt = pt;
            pt = FindCornerPoint(i,j,k, 0, 0, 1);
            if(pt.z > botpt.z)
              botpt = pt;
            if(pt.z < toppt.z)
              toppt = pt;
          }
        }
        if(i > 0 && j < nj_){
          if(IsActive(i-1,j,k)){
            pt = FindCornerPoint(i-1, j, k, 1, 0, 0);
            if(botpt.z == -999.0 || pt.z > botpt.z)
              botpt = pt;
            if(toppt.z == -999.0 || pt.z < toppt.z)
              toppt = pt;
            pt = FindCornerPoint(i-1, j, k, 1, 0, 1);
            if(botpt.z == -999.0 || pt.z > botpt.z)
              botpt = pt;
            if(toppt.z == -999.0 || pt.z < toppt.z)
              toppt = pt;
          }
        }
        if(j > 0 && i < ni_){
          if(IsActive(i,j-1,k)){
            pt = FindCornerPoint(i, j - 1, k, 0, 1, 0);
            if(botpt.z == -999.0 || pt.z > botpt.z)
              botpt = pt;
            if(toppt.z == -999.0 || pt.z < toppt.z)
              toppt = pt;
            pt = FindCornerPoint(i, j - 1, k, 0, 1, 1);
            if(botpt.z == -999.0 || pt.z > botpt.z)
              botpt = pt;
            if(toppt.z == -999.0 || pt.z < toppt.z)
              toppt = pt;
          }
        }
        if(i > 0 && j > 0){
          if(IsActive(i -1,j - 1,k)){
            pt = FindCornerPoint(i - 1, j - 1, k, 1, 1, 0);
            if(botpt.z == -999.0 || pt.z > botpt.z)
              botpt = pt;
            if(toppt.z == -999.0 || pt.z < toppt.z)
              toppt = pt;
            pt = FindCornerPoint(i - 1, j - 1, k, 1, 1, 1);
            if(botpt.z == -999.0 || pt.z > botpt.z)
              botpt = pt;
            if(toppt.z == -999.0 || pt.z < toppt.z)
              toppt = pt;
          }
        }
      }
    }
  polymin_ = FindPolygonAroundActivePillars(toppt.z);
  polymax_ = FindPolygonAroundActivePillars(botpt.z);
}


size_t EclipseGeometry::FindTopCell(size_t i, size_t j) const
{
  size_t k = 0;
  while (k < GetNK()) {
    if (IsActive(i, j, k))
      break;
    ++k;
  }
  return k;
}


size_t EclipseGeometry::FindBottomCell(size_t i, size_t j) const
{
  int k = static_cast<int>(GetNK() - 1);
  while (k >= 0) {
    if (IsActive(i, j, k))
      break;
    --k;
  }

  if (k < 0)
    return GetNK();
  return static_cast<size_t>(k);
}


void EclipseGeometry::FindEnclosingVolume(double& x0,
                                          double& y0,
                                          double& lx,
                                          double& ly,
                                          double& angle) const
{
  std::vector<NRLib::Point> points;
  points.reserve(polymax_.GetSize()+polymin_.GetSize());
  std::vector<NRLib::Point> under;
  polymax_.GetPoints(points);
  polymin_.GetPoints(under);
  for (size_t i=0; i<under.size(); i++) {
    points.push_back(under[i]);
  }
  NRLib::PointSet surface(points);
  NRLib::Polygon convex_hull=surface.GetConvexHull();
  convex_hull.MinEnclosingRectangle(x0,y0,lx,ly,angle);
  //DEBUGGING
  //convex_hull.GetPoints(under);
  //std::cout<<"\n\nConvex hull points\n";
  //for (size_t i=0; i<under.size(); i++) {
  //  under[i].z=0.0;
  //  std::cout<<under[i].x<<" "<<under[i].y<<"\n";
  //}
 // NRLib::PointSet polygon(under);
 // polygon.WriteToFile("convexpolygonm.dat", NRLib::PointSet::RoxarText);
}


size_t EclipseGeometry::FindTopLayer() const
{
  size_t top_ij;
  size_t top=FindTopCell(0,0);
  for(size_t i = 0; i < ni_; i++) {
    for(size_t j = 0; j < nj_; j++){
      top_ij=FindTopCell(i,j);
      if (top>top_ij)
        top=top_ij;
    }
  }
  return top;
}

size_t EclipseGeometry::FindBottomLayer() const
{
  size_t bot_ij;
  size_t bot=FindBottomCell(0,0);
  for(size_t i = 0; i < ni_; i++) {
    for(size_t j = 0; j < nj_; j++){
      bot_ij=FindBottomCell(i,j);
      if (bot_ij<nk_) {
        if (bot<bot_ij || bot==nk_)
          bot=bot_ij;
      }
    }
  }
  return bot;
}

void EclipseGeometry::BilinearFillInZValuesInArea(NRLib::Grid2D<double> &z_surface,
                                                  NRLib::Grid2D<int> &is_set,
                                                  double x0,
                                                  double y0,
                                                  const std::vector<NRLib::Point>& corners,
                                                  double dx,
                                                  double dy) const
{
  size_t m=z_surface.GetNJ();
  size_t n=z_surface.GetNI();
  NRLib::BilinearSurface bilinear_corners(corners[0], corners[1], corners[2], corners[3]);
  double min_x,max_x,min_y,max_y;
  NRLib::Point point_xy(0.0,0.0,0.0);
  NRLib::Point z_dir(0.0,0.0,1.0);
  NRLib::Point intersec1, intersec2;
  size_t n1,n2,m1,m2;
  int nu_of_intersec;
  //Calculate the min and max in x- and y-direction (rotated, i.e. same direction as z_surface)
  min_x=corners[0].x;
  max_x=min_x;
  min_y=corners[0].y;
  max_y=min_y;
  for (size_t four=1; four<4; four++) {
    if (corners[four].x<min_x) {
      min_x=corners[four].x;
    }
    else if (corners[four].x>max_x) {
      max_x=corners[four].x;
    }
    if (corners[four].y<min_y) {
      min_y=corners[four].y;
    }
    else if (corners[four].y>max_y) {
      max_y=corners[four].y;
    }
  }
  //For loop over all points in z_surface inside the rectangle given by the mins and max' calulated above
  n1=static_cast<size_t>(max((min_x-x0)/dx-0.5,0.0)); //min_x=(x0+dx/2)+n1*dx, if n1 not integer, let it be the smallest integer greater than solution. Zero if negative number
  n2=static_cast<size_t>(max((max_x-x0)/dx+1.0,0.0)); //max_x=(x0+dx/2)+(n2-1)*dx, if n2 not integer, let it be the greatest integer smaller than solution. Zero if negative number
  if (n2>n)
    n2=n; //Should stop before grid ends
  m1=static_cast<size_t>(max((min_y-y0)/dy-0.5,0.0));
  m2=static_cast<size_t>(max((max_y-y0)/dy+1.0,0.0));
  if (m2>m)
    m2=m;
  for (size_t it1=n1; it1<n2; it1++) {
    for (size_t it2=m1; it2<m2; it2++) {
      point_xy.x=x0+dx/2+it1*dx;
      point_xy.y=y0+dy/2+it2*dy;
      NRLib::Line line_xy(point_xy,point_xy+z_dir,false,false);
      nu_of_intersec=bilinear_corners.FindIntersections(line_xy,intersec1,intersec2);
      if (nu_of_intersec>=1) {
        if (is_set(it1,it2)) {
          z_surface(it1,it2)+=intersec1.z;
          z_surface(it1,it2)=z_surface(it1,it2)/2;
        }
        else {
          z_surface(it1,it2)=intersec1.z;
          is_set(it1,it2)=true;
        }
      }
    }
  }
}

void EclipseGeometry::TriangularFillInZValuesInArea(NRLib::Grid2D<double> &z_surface,
                                                    NRLib::Grid2D<int> &is_set,
                                                    double x0,
                                                    double y0,
                                                    const std::vector<NRLib::Point>& corners_in,
                                                    double dx,
                                                    double dy) const
{
  std::vector<NRLib::Point> corners = corners_in;
  size_t m=z_surface.GetNJ();
  size_t n=z_surface.GetNI();
  double min_x,max_x,min_y,max_y,vec1_vec2_angle;
  NRLib::Triangle triangle1, triangle2;
  NRLib::Point point_xy(0.0,0.0,0.0);
  NRLib::Point vec1(0.0,0.0,0.0);
  NRLib::Point vec2(0.0,0.0,0.0);
  NRLib::Point z_dir(0.0,0.0,1.0);
  NRLib::Point intersec;
  size_t n2,m1,m2;
  bool two_triangles=true;
  //Check if two of the points are the same
  size_t n1=3;
  for (m1=0;m1<4;m1++) {
    if (corners[m1]==corners[n1]) { //NOTE: Could do better check on whether two points are equal !!!!!!!
      two_triangles=false;
      corners[n1]=corners[3];
      triangle1.SetCornerPoints(corners[0], corners[1], corners[2]);
    }
    n1=m1;
  }
  if (two_triangles) {
    // Calculate the sum of two opposite angles
    vec1.x=corners[1].x-corners[0].x;
    vec1.y=corners[1].y-corners[0].y;
    vec2.x=corners[3].x-corners[0].x;
    vec2.y=corners[3].y-corners[0].y;
    vec1_vec2_angle=vec1.GetAngle(vec2);
    vec1.x=corners[1].x-corners[2].x;
    vec1.y=corners[1].y-corners[2].y;
    vec2.x=corners[3].x-corners[2].x;
    vec2.y=corners[3].y-corners[2].y;
    vec1_vec2_angle+=vec1.GetAngle(vec2);
    // Make delunay triangles (according to the sum of the angles)
    if (vec1_vec2_angle<=NRLib::Pi) {
      triangle1.SetCornerPoints(corners[3], corners[0], corners[1]);
      triangle2.SetCornerPoints(corners[1], corners[2], corners[3]);
    }
    else {
      triangle1.SetCornerPoints(corners[0], corners[1], corners[2]);
      triangle2.SetCornerPoints(corners[2], corners[3], corners[0]);
    }
  }
  //Calculate the min and max in x- and y-direction (rotated, i.e. same direction as z_surface)
  min_x=corners[0].x;
  max_x=min_x;
  min_y=corners[0].y;
  max_y=min_y;
  for (size_t four=1; four<4; four++) {
    if (corners[four].x<min_x)
      min_x=corners[four].x;
    else if (corners[four].x>max_x)
      max_x=corners[four].x;
    if (corners[four].y<min_y)
      min_y=corners[four].y;
    else if (corners[four].y>max_y)
      max_y=corners[four].y;
  }
  //For loop over all points in z_surface inside the rectangle given by the mins and max' calulated above
  n1=static_cast<size_t>(max( ((min_x-x0)/dx-0.5),0.0)); //min_x=(x0+dx/2)+n1*dx, if n1 not integer, let it be the smallest integer greater than solution. Zero if negative number
  n2=static_cast<size_t>(max( ((max_x-x0)/dx+1.0),0.0)); //max_x=(x0+dx/2)+(n2-1)*dx, if n2 not integer, let it be the greatest integer smaller than solution. Zero if negative number
  if (n2>n)
    n2=n; //Should stop before grid ends
  m1=static_cast<size_t>(max( ((min_y-y0)/dy-0.5),0.0));
  m2=static_cast<size_t>(max( ((max_y-y0)/dy+1.0),0.0));
  if (m2>m)
    m2=m;
  for (size_t it1=n1; it1<n2; it1++) {
    for (size_t it2=m1; it2<m2; it2++) {
      point_xy.x=x0+dx/2+it1*dx;
      point_xy.y=y0+dy/2+it2*dy;
      NRLib::Line line_xy(point_xy,(point_xy+z_dir),false,false);
      if (triangle1.FindIntersection(line_xy,intersec,true)) {
        if (is_set(it1,it2)>0 ) {
          z_surface(it1,it2)*=static_cast<double>(1.0*is_set(it1,it2)/(1.0*is_set(it1,it2)+1.0));
          z_surface(it1,it2)+=intersec.z/static_cast<double>(is_set(it1,it2)+1.0);
          is_set(it1,it2)++;
          //z_surface(it1,it2)=z_surface(it1,it2)/2;
        }
       // else if(!is_set(it1,it2)){
        else {
          z_surface(it1,it2)=intersec.z;
          is_set(it1,it2) = 1;
        }
      }
      else if (two_triangles) {
        if (triangle2.FindIntersection(line_xy,intersec,true)) {
          if (is_set(it1,it2)>0 ) {
            z_surface(it1,it2)*=static_cast<double>(1.0*is_set(it1,it2)/(1.0*is_set(it1,it2)+1.0));
          z_surface(it1,it2)+=intersec.z/static_cast<double>(is_set(it1,it2)+1.0);
          is_set(it1,it2)++;
            //z_surface(it1,it2)+=intersec.z;
            //z_surface(it1,it2)=z_surface(it1,it2)/2;
          }
          //else if(!is_set(it1,it2)){
          else{
            z_surface(it1,it2)=intersec.z;
            is_set(it1,it2) =1;
          }
        }
      }
    }
  }
}

// Corner point interpolation.  This routine does not work with reverse faults.
void EclipseGeometry::FindLayerSurfaceCornerpoint(NRLib::Grid2D<double> &z_surface,
                                       size_t k,
                                       int lower_or_upper,
                                       double dx,
                                       double dy,
                                       double x0,
                                       double y0,
                                       double angle,
                                       bool bilinear_else_triangles) const
{
  size_t m=z_surface.GetNJ();
  size_t n=z_surface.GetNI();
  double rot_x0=cos(angle)*x0+sin(angle)*y0;
  double rot_y0=cos(angle)*y0-sin(angle)*x0;
  double test_if_equal;
  bool fault;
  NRLib::Point nonrotated_corner;
  NRLib::Point prev_upper_corner,prev_lower_corner, cell_under_right_corner, cell_under_left_corner;
  std::vector<NRLib::Point> corners(4);
  std::vector<NRLib::Point> fault_corners(4);
  NRLib::Grid2D<int> is_set(n,m,0);
  for(size_t j = 0; j < nj_; j++) { //Loops over each cell in the given layer
    for(size_t i = 0; i < ni_; i++){
      // Find rotated coordinates for the corners of the cell
      nonrotated_corner=FindCornerPoint(i,j,k,0,0,lower_or_upper);
      corners[0].x=cos(angle)*nonrotated_corner.x+sin(angle)*nonrotated_corner.y;
      corners[0].y=cos(angle)*nonrotated_corner.y-sin(angle)*nonrotated_corner.x;
      corners[0].z=nonrotated_corner.z;
      nonrotated_corner=FindCornerPoint(i,j,k,0,1,lower_or_upper);
      corners[1].x=cos(angle)*nonrotated_corner.x+sin(angle)*nonrotated_corner.y;
      corners[1].y=cos(angle)*nonrotated_corner.y-sin(angle)*nonrotated_corner.x;
      corners[1].z=nonrotated_corner.z;
      nonrotated_corner=FindCornerPoint(i,j,k,1,1,lower_or_upper);
      corners[2].x=cos(angle)*nonrotated_corner.x+sin(angle)*nonrotated_corner.y;
      corners[2].y=cos(angle)*nonrotated_corner.y-sin(angle)*nonrotated_corner.x;
      corners[2].z=nonrotated_corner.z;
      nonrotated_corner=FindCornerPoint(i,j,k,1,0,lower_or_upper);
      corners[3].x=cos(angle)*nonrotated_corner.x+sin(angle)*nonrotated_corner.y;
      corners[3].y=cos(angle)*nonrotated_corner.y-sin(angle)*nonrotated_corner.x;
      corners[3].z=nonrotated_corner.z;
      if (FindTopCell(i,j)!=nk_) {
        if (bilinear_else_triangles) {
          BilinearFillInZValuesInArea(z_surface,is_set,rot_x0,rot_y0,corners,dx,dy);
        }
        else {
          TriangularFillInZValuesInArea(z_surface,is_set,rot_x0,rot_y0,corners,dx,dy);
        }
      }
      if (j>0) {
        //Find rotated coordinates for the corners of the cell under
        nonrotated_corner=FindCornerPoint(i,j-1,k,0,1,lower_or_upper);
        cell_under_left_corner.x=cos(angle)*nonrotated_corner.x+sin(angle)*nonrotated_corner.y;
        cell_under_left_corner.y=cos(angle)*nonrotated_corner.y-sin(angle)*nonrotated_corner.x;
        cell_under_left_corner.z=nonrotated_corner.z;
        nonrotated_corner=FindCornerPoint(i,j-1,k,1,1,lower_or_upper);
        cell_under_right_corner.x=cos(angle)*nonrotated_corner.x+sin(angle)*nonrotated_corner.y;
        cell_under_right_corner.y=cos(angle)*nonrotated_corner.y-sin(angle)*nonrotated_corner.x;
        cell_under_right_corner.z=nonrotated_corner.z;
        fault=false;
        test_if_equal=(cell_under_left_corner.y/corners[0].y)-1.0;
        if (abs(test_if_equal)>0.00000001)
          fault=true;
        else {
          test_if_equal=(cell_under_right_corner.y/corners[3].y)-1.0;
          if (abs(test_if_equal)>0.00000001)
            fault=true;
          else {
            test_if_equal=(cell_under_right_corner.x/corners[3].x)-1.0;
            if (abs(test_if_equal)>0.00000001)
              fault=true;
            else {
              test_if_equal=(cell_under_left_corner.x/corners[0].x)-1.0;
              if (abs(test_if_equal)>0.00000001)
                fault=true;
            }
          }
        }
        if (fault){
          //Fault along the i-coordinate
          fault_corners[0]=cell_under_left_corner;
          fault_corners[1]=corners[0];
          fault_corners[2]=corners[3];
          fault_corners[3]=cell_under_right_corner;
          if(corners[0].y<cell_under_left_corner.y && corners[3].y < cell_under_right_corner.y){
           // double meanz = 0.125*(corners[0].z+cell_under_left_corner.z+cell_under_right_corner.z+corners[3].z);
            fault_corners[0].z = corners[0].z;
            fault_corners[1].z = cell_under_left_corner.z;
            fault_corners[2].z = cell_under_right_corner.z;
            fault_corners[3].z = corners[3].z;
          }

          if (bilinear_else_triangles)
            BilinearFillInZValuesInArea(z_surface,is_set,rot_x0,rot_y0,fault_corners,dx,dy);
          else
            TriangularFillInZValuesInArea(z_surface,is_set,rot_x0,rot_y0,fault_corners,dx,dy);
        }
      }
      if (i>0) {
        fault=false;
        test_if_equal=(prev_upper_corner.x/corners[1].x)-1.0;
        if (abs(test_if_equal)>0.00000001)
          fault=true;
        else {
          test_if_equal=(prev_lower_corner.x/corners[0].x)-1.0;
          if (abs(test_if_equal)>0.00000001)
            fault=true;
          else {
            test_if_equal=(prev_upper_corner.y/corners[1].y)-1.0;
            if (abs(test_if_equal)>0.00000001)
              fault=true;
            else {
              test_if_equal=(prev_lower_corner.y/corners[0].y)-1.0;
              if (abs(test_if_equal)>0.00000001)
                fault=true;
            }
          }
        }
        if (fault) {
          //Fault along the j-coordinate
          fault_corners[0]=prev_lower_corner;
          fault_corners[1]=prev_upper_corner;
          fault_corners[2]=corners[1];
          fault_corners[3]=corners[0];
          if(corners[1].x < prev_upper_corner.x && corners[0].x < prev_lower_corner.x){
            //double meanz = 0.125*(corners[1].z+corners[0].z+prev_upper_corner.z+prev_lower_corner.z);
            fault_corners[0].z = corners[1].z;
            fault_corners[1].z = corners[0].z;
            fault_corners[2].z = prev_upper_corner.z;
            fault_corners[3].z = prev_lower_corner.z;

          }
          if (bilinear_else_triangles)
            BilinearFillInZValuesInArea(z_surface,is_set,rot_x0,rot_y0,fault_corners,dx,dy);
          else
            TriangularFillInZValuesInArea(z_surface,is_set,rot_x0,rot_y0,fault_corners,dx,dy);
        }
      }
      prev_upper_corner=corners[2];
      prev_lower_corner=corners[3];
    }
  }
  FillInZValuesByAveraging(z_surface,is_set);
}


// center point interpolation
void EclipseGeometry::FindLayerSurface(NRLib::Grid2D<double> &z_surface,
                                       size_t k,
                                       int lower_or_upper,
                                       double dx,
                                       double dy,
                                       double x0,
                                       double y0,
                                       double angle,
                                       bool bilinear_else_triangles) const
{
  size_t m=z_surface.GetNJ();
  size_t n=z_surface.GetNI();
  double rot_x0=cos(angle)*x0+sin(angle)*y0;
  double rot_y0=cos(angle)*y0-sin(angle)*x0;

  NRLib::Point nonrotated_corner;
  NRLib::Point prev_upper_corner,prev_lower_corner, cell_under_right_corner, cell_under_left_corner;
  std::vector<NRLib::Point> corners(4);
  NRLib::Grid2D<int> is_set(n,m,0);

  for(size_t j = 0; j < nj_-1; j++) { //Loops over each cell in the given layer
    for(size_t i = 0; i < ni_-1; i++){
      if(IsPillarActive(i,j) && IsPillarActive(i+1,j) && IsPillarActive(i, j+1) && IsPillarActive(i+1, j+1) &&
        IsPillarActive(i+2,j) && IsPillarActive(i+2, j+1) && IsPillarActive(i, j+2) && IsPillarActive(i+1, j+2) &&
        IsPillarActive(i+2, j+2)){
          nonrotated_corner=FindPointCellSurface(i, j, k,lower_or_upper,0.5, 0.5);
          corners[0].x=cos(angle)*nonrotated_corner.x+sin(angle)*nonrotated_corner.y;
          corners[0].y=cos(angle)*nonrotated_corner.y-sin(angle)*nonrotated_corner.x;
          corners[0].z=nonrotated_corner.z;
          nonrotated_corner=FindPointCellSurface(i+1, j, k,lower_or_upper,0.5, 0.5);
          corners[1].x=cos(angle)*nonrotated_corner.x+sin(angle)*nonrotated_corner.y;
          corners[1].y=cos(angle)*nonrotated_corner.y-sin(angle)*nonrotated_corner.x;
          corners[1].z=nonrotated_corner.z;
          // Find rotated coordinates for the corners of the cell
          nonrotated_corner=FindPointCellSurface(i+1, j+1, k,lower_or_upper,0.5, 0.5);
          corners[2].x=cos(angle)*nonrotated_corner.x+sin(angle)*nonrotated_corner.y;
          corners[2].y=cos(angle)*nonrotated_corner.y-sin(angle)*nonrotated_corner.x;
          corners[2].z=nonrotated_corner.z;
          nonrotated_corner=FindPointCellSurface(i, j+1, k,lower_or_upper,0.5, 0.5);
          corners[3].x=cos(angle)*nonrotated_corner.x+sin(angle)*nonrotated_corner.y;
          corners[3].y=cos(angle)*nonrotated_corner.y-sin(angle)*nonrotated_corner.x;
          corners[3].z=nonrotated_corner.z;

          if (FindTopCell(i,j)!=nk_) {
            if (bilinear_else_triangles)
              BilinearFillInZValuesInArea(z_surface,is_set,rot_x0,rot_y0,corners,dx,dy);
            else
              TriangularFillInZValuesInArea(z_surface,is_set,rot_x0,rot_y0,corners,dx,dy);
          }
      }
    }
  }
  FillInZValuesByAveraging(z_surface,is_set);
}





void EclipseGeometry::FillInZValuesByAveraging(NRLib::Grid2D<double> &z_surface,
                                             NRLib::Grid2D<int> &is_set) const {
  size_t m=z_surface.GetNJ();
  size_t n=z_surface.GetNI();
  size_t average_i=0;
  size_t average_j=0;
  size_t count=0;
  size_t distance, iter_i, iter_j;
  for (size_t j=0; j<m;j++) {
    for (size_t i=0; i<n; i++) {
      if (is_set(i,j)) {
        average_i+=i;
        average_j+=j;
        count++;
      }
    }
  }
  if (count>0 && count<m*n) {
    average_i=average_i/count;
    average_j=average_j/count;
    distance=min(min(average_i,average_j),min(n-1-average_i,m-1-average_j));
    iter_i=average_i;
    iter_j=average_j-1;
    for (size_t it_dist=1; it_dist<=distance; it_dist++) {
      for (size_t r=1; r<=8*it_dist;r++) {
        if (!is_set(iter_i,iter_j)) {
          z_surface(iter_i,iter_j)=0.0;
          count=0;
          if (iter_i>0) {
            if (is_set(iter_i-1,iter_j)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i-1,iter_j);
              count++;
            }
          }
          if ((iter_i+1)<n) {
            if (is_set(iter_i+1,iter_j)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i+1,iter_j);
              count++;
            }
          }
          if (iter_j>0) {
            if (is_set(iter_i,iter_j-1)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i,iter_j-1);
              count++;
            }
          }
          if ((iter_j+1)<m) {
            if (is_set(iter_i,iter_j+1)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i,iter_j+1);
              count++;
            }
          }
          z_surface(iter_i,iter_j)=z_surface(iter_i,iter_j)/count;
          is_set(iter_i,iter_j)=true;
        }
        if (r<2*it_dist) {
          iter_i++;
        }
        else if (r<4*it_dist) {
          iter_j++;
        }
        else if (r<6*it_dist) {
          iter_i--;
        }
        else {
          iter_j--;
        }
      }
    }
    iter_j=average_j-distance;
    while (iter_j>0) {
      iter_j--;
      iter_i=average_i-distance;
      while (iter_i<=average_i+distance) {
        if (!is_set(iter_i,iter_j)) {
          z_surface(iter_i,iter_j)=0.0;
          count=0;
          if (iter_i>0) {
            if (is_set(iter_i-1,iter_j)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i-1,iter_j);
              count++;
            }
          }
          if ((iter_i+1)<n) {
            if (is_set(iter_i+1,iter_j)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i+1,iter_j);
              count++;
            }
          }
          if (iter_j>0) {
            if (is_set(iter_i,iter_j-1)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i,iter_j-1);
              count++;
            }
          }
          if ((iter_j+1)<m) {
            if (is_set(iter_i,iter_j+1)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i,iter_j+1);
              count++;
            }
          }
          z_surface(iter_i,iter_j)=z_surface(iter_i,iter_j)/count;
          is_set(iter_i,iter_j)=true;
        }
        iter_i++;
      }
    }
    iter_j=average_j+distance+1;
    while (iter_j<m) {
      iter_i=average_i-distance;
      while (iter_i<=average_i+distance) {
        if (!is_set(iter_i,iter_j)) {
          z_surface(iter_i,iter_j)=0.0;
          count=0;
          if (iter_i>0) {
            if (is_set(iter_i-1,iter_j)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i-1,iter_j);
              count++;
            }
          }
          if ((iter_i+1)<n) {
            if (is_set(iter_i+1,iter_j)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i+1,iter_j);
              count++;
            }
          }
          if (iter_j>0) {
            if (is_set(iter_i,iter_j-1)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i,iter_j-1);
              count++;
            }
          }
          if ((iter_j+1)<m) {
            if (is_set(iter_i,iter_j+1)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i,iter_j+1);
              count++;
            }
          }
          z_surface(iter_i,iter_j)=z_surface(iter_i,iter_j)/count;
          is_set(iter_i,iter_j)=true;
        }
        iter_i++;
      }
      iter_j++;
    }
    iter_i=average_i-distance;
    while (iter_i>0) {
      iter_i--;
      iter_j=0;
      while (iter_j<m) {
        if (!is_set(iter_i,iter_j)) {
          z_surface(iter_i,iter_j)=0.0;
          count=0;
          if (iter_i>0) {
            if (is_set(iter_i-1,iter_j)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i-1,iter_j);
              count++;
            }
          }
          if ((iter_i+1)<n) {
            if (is_set(iter_i+1,iter_j)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i+1,iter_j);
              count++;
            }
          }
          if (iter_j>0) {
            if (is_set(iter_i,iter_j-1)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i,iter_j-1);
              count++;
            }
          }
          if ((iter_j+1)<m) {
            if (is_set(iter_i,iter_j+1)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i,iter_j+1);
              count++;
            }
          }
          z_surface(iter_i,iter_j)=z_surface(iter_i,iter_j)/count;
          is_set(iter_i,iter_j)=true;
        }
        iter_j++;
      }
    }
    iter_i=average_i+distance+1;
    while (iter_i<n) {
      iter_j=0;
      while (iter_j<m) {
        if (!is_set(iter_i,iter_j)) {
          z_surface(iter_i,iter_j)=0.0;
          count=0;
          if (iter_i>0) {
            if (is_set(iter_i-1,iter_j)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i-1,iter_j);
              count++;
            }
          }
          if ((iter_i+1)<n) {
            if (is_set(iter_i+1,iter_j)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i+1,iter_j);
              count++;
            }
          }
          if (iter_j>0) {
            if (is_set(iter_i,iter_j-1)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i,iter_j-1);
              count++;
            }
          }
          if ((iter_j+1)<m) {
            if (is_set(iter_i,iter_j+1)) {
              z_surface(iter_i,iter_j)+=z_surface(iter_i,iter_j+1);
              count++;
            }
          }
          z_surface(iter_i,iter_j)=z_surface(iter_i,iter_j)/count;
          is_set(iter_i,iter_j)=true;
        }
        iter_j++;
      }
      iter_i++;
    }
  }
}

void EclipseGeometry::FindTopAndBotValuesOfGrid(std::vector<NRLib::Point> &toppoints, std::vector<NRLib::Point> &botpoints) const
{
  size_t i, j, k;
  Point pt, botpt, toppt;
  for(i = 0; i <= ni_; i++) {
    for(j = 0; j <= nj_; j++) {
      toppt.z = -999.0;
      botpt.z = -999.0;
      for(k = 0; k < nk_; k++){
        if(i < ni_ && j < nj_){
          if(IsActive(i,j,k)){
            pt = FindCornerPoint(i,j,k,0,0,0);
            if(botpt.z == -999.0 || pt.z > botpt.z)
              botpt = pt;
            if(toppt.z == -999.0 || pt.z < toppt.z)
              toppt = pt;
            pt = FindCornerPoint(i,j,k, 0, 0, 1);
            if(pt.z > botpt.z)
              botpt = pt;
            if(pt.z < toppt.z)
              toppt = pt;
          }
        }
        if(i > 0 && j < nj_){
          if(IsActive(i-1,j,k)){
            pt = FindCornerPoint(i-1, j, k, 1, 0, 0);
            if(botpt.z == -999.0 || pt.z > botpt.z)
              botpt = pt;
            if(toppt.z == -999.0 || pt.z < toppt.z)
              toppt = pt;
            pt = FindCornerPoint(i-1, j, k, 1, 0, 1);
            if(botpt.z == -999.0 || pt.z > botpt.z)
              botpt = pt;
            if(toppt.z == -999.0 || pt.z < toppt.z)
              toppt = pt;
          }
        }
        if(j > 0 && i < ni_){
          if(IsActive(i,j-1,k)){
            pt = FindCornerPoint(i, j - 1, k, 0, 1, 0);
            if(botpt.z == -999.0 || pt.z > botpt.z)
              botpt = pt;
            if(toppt.z == -999.0 || pt.z < toppt.z)
              toppt = pt;
            pt = FindCornerPoint(i, j - 1, k, 0, 1, 1);
            if(botpt.z == -999.0 || pt.z > botpt.z)
              botpt = pt;
            if(toppt.z == -999.0 || pt.z < toppt.z)
              toppt = pt;
          }
        }
        if(i > 0 && j > 0){
          if(IsActive(i -1,j - 1,k)){
            pt = FindCornerPoint(i - 1, j - 1, k, 1, 1, 0);
            if(botpt.z == -999.0 || pt.z > botpt.z)
              botpt = pt;
            if(toppt.z == -999.0 || pt.z < toppt.z)
              toppt = pt;
            pt = FindCornerPoint(i - 1, j - 1, k, 1, 1, 1);
            if(botpt.z == -999.0 || pt.z > botpt.z)
              botpt = pt;
            if(toppt.z == -999.0 || pt.z < toppt.z)
              toppt = pt;
          }
        }
      }
      if(toppt.z !=-999.0)
        toppoints.push_back(toppt);
      if(botpt.z !=-999.0)
        botpoints.push_back(botpt);
    }
  }
}

double
  EclipseGeometry::GetDZ(size_t i, size_t j, size_t k) const
{

  double dz = 0;
  dz = GetZCorner(i, j, k, 0, 0, 1) - GetZCorner(i, j, k, 0, 0, 0);
  dz += GetZCorner(i, j, k, 1, 0, 1) - GetZCorner(i, j, k, 1, 0, 0);
  dz += GetZCorner(i, j, k, 1, 1, 1) - GetZCorner(i, j, k, 1, 1, 0);
  dz += GetZCorner(i, j, k, 0, 1, 1) - GetZCorner(i, j, k, 0, 1, 0);
  return dz;

}
