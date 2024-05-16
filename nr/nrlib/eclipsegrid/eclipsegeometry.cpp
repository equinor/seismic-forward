// $Id: eclipsegeometry.cpp 1461 2017-04-03 15:18:25Z eyaker $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// o  Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// o  Redistributions in binary form must reproduce the above copyright notice, this list of
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
#include "../extrapolation/extrapolategrid2d.hpp"
#include "../geometry/polygon.hpp"
#include "../geometry/triangle.hpp"
#include "../geometry/bilinearsurface.hpp"
#include "../iotools/logkit.hpp"
#include "../iotools/fileio.hpp"
#include "../iotools/stringtools.hpp"
#include "../math/constants.hpp"

#include "../surface/regularsurfacerotated.hpp"
#include "../iotools/stringtools.hpp"


#include <stdexcept>
#include <fstream>
#include <iostream>
#include <ostream>
#include <istream>
#include <vector>

#include <ctime>

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
EclipseGeometry::FindMeanPillarDistance(double z) const
{
  double sum_d = 0;
  double nd = 0;
  for (size_t i = 1; i <= ni_; i++) {
    for (size_t j = 1; j <= nj_; j++) {
      if (IsPillarActive(i, j))
      {
        size_t im1 = i - 1;
        size_t jm1 = j - 1;
        NRLib::Point pillar_pt = FindPointAtPillar(i, j, z);
        if (IsPillarActive(im1, j))
        {
          double d = pillar_pt.GetDistance(FindPointAtPillar(im1, j, z));
          sum_d += d;
          nd += 1;
        }
        if (IsPillarActive(i, jm1))
        {
          double d = pillar_pt.GetDistance(FindPointAtPillar(i, jm1, z));
          sum_d += d;
          nd += 1;
        }
      }
    }
  }

  for (size_t j = 1; j <= nj_; j++)
  {
    if (IsPillarActive(0, j))
    {
      size_t jm1 = j - 1;
      NRLib::Point pillar_pt = FindPointAtPillar(0, j, z);
      if (IsPillarActive(0, jm1))
      {
        double d = pillar_pt.GetDistance(FindPointAtPillar(0, jm1, z));
        sum_d += d;
        nd += 1;
      }
    }
  }
  for (size_t i = 1; i <= ni_; i++)
  {
    if (IsPillarActive(i, 0))
    {
      size_t im1 = i - 1;
      NRLib::Point pillar_pt = FindPointAtPillar(i, 0, z);
      if (IsPillarActive(im1, 0))
      {
        double d = pillar_pt.GetDistance(FindPointAtPillar(im1, 0, z));
        sum_d += d;
        nd += 1;
      }
    }
  }
  return sum_d / nd;
}

double
EclipseGeometry::FindZTopInCellActiveColumn(size_t i, size_t j, size_t k, bool & found) const
{
  found = false;
  double z;
  double z_top = numeric_limits<double>::infinity();
  if (k < GetNK()) {
    if (IsColumnActive(i, j)) {
      z = GetZCorner(i, j, k, 1, 1, 0);
      if (z < z_top)
        z_top = z;
      z = GetZCorner(i, j, k, 1, 0, 0);
      if (z < z_top)
        z_top = z;
      z = GetZCorner(i, j, k, 0, 1, 0);
      if (z < z_top)
        z_top = z;
      z = GetZCorner(i, j, k, 0, 0, 0);
      if (z < z_top)
        z_top = z;
      found = true;
    }
  }
  return(z_top);
}

double
EclipseGeometry::FindZBotInCellActiveColumn(size_t i, size_t j, size_t k, bool & found) const
{
  found = false;
  double z;
  double z_bot = -numeric_limits<double>::infinity();
  if (k < GetNK()) {
    if (IsColumnActive(i, j)) {
      z = GetZCorner(i, j, k, 1, 1, 1);
        if (z > z_bot)
          z_bot = z;
      z = GetZCorner(i, j, k, 1, 0, 1);
        if (z > z_bot)
          z_bot = z;
      z = GetZCorner(i, j, k, 0, 1, 1);
        if (z > z_bot)
          z_bot = z;
      z = GetZCorner(i, j, k, 0, 0, 1);
        if (z > z_bot)
          z_bot = z;
        found = true;
    }
  }
  return(z_bot);
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

  Point a1 = (p00 + p10 - q00 - q10) + (p10 + p11 - r10 - r11) - (p01 + p11 - q01 - q11) - (p00 + p01 - r00 - r01) + (r10 + r11 - 2*p10);
  Point b1 = -1*(r00 + r01 - 2*p00);
  Point c1 = q00 + q10 - 2*p00;

  Point a0 = p00 + p01 - r00 - r01;
  Point b0 = r00 + r01 - 2*p00;
  Point c0 = p00;

  Point f = (a2 * v_in*v_in + b2*v_in + c2) * u_in*u_in
            + (a1 * v_in*v_in + b1*v_in + c1) * u_in
            + a0 * v_in*v_in + b0*v_in + c0;

  /*
  For debugging

  PointSet ps1;
  ps1.AddPoint(p00);
  ps1.AddPoint(p10);
  ps1.AddPoint(p01);
  ps1.AddPoint(p11);
  ps1.AddPoint(q00);
  ps1.AddPoint(q10);
  ps1.AddPoint(q01);
  ps1.AddPoint(q11);
  ps1.AddPoint(r00);
  ps1.AddPoint(r10);
  ps1.AddPoint(r01);
  ps1.AddPoint(r11);
  std::vector<std::string> labels1{"p00", "p10", "p01", "p11", "q00", "q10", "q01", "q11", "r00", "r10", "r01", "r11"};
  ps1.AddStringParameter("Point", labels1);
  ps1.WriteToFile("01_first_points.rxat", PointSet::RmsInternalPoints);

  PointSet ps2;
  ps2.AddPoint(a0);
  ps2.AddPoint(a1);
  ps2.AddPoint(a2);
  ps2.AddPoint(b0);
  ps2.AddPoint(b1);
  ps2.AddPoint(b2);
  ps2.AddPoint(c0);
  ps2.AddPoint(c1);
  ps2.AddPoint(c2);
  std::vector<std::string> labels2{"a0", "a1", "a2", "b0", "b1", "b2", "c0", "c1", "c2"};
  ps2.AddStringParameter("Point", labels2);
  ps2.WriteToFile("02_second_points.rxat", PointSet::RmsInternalPoints);

  PointSet ps3;
  ps3.AddPoint(f);
  std::vector<std::string> labels3{"f"};
  ps3.AddStringParameter("Point", labels3);
  ps3.WriteToFile("03_final_point.rxat", PointSet::RmsInternalPoints);
  */

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

void EclipseGeometry::BilinearFillInZValuesInArea(NRLib::Grid2D<double>           & z_grid,
                                                  NRLib::Grid2D<int>              & is_set,
                                                  const std::vector<NRLib::Point> & corners,
                                                  const double                      dx,
                                                  const double                      dy) const
{
  // Calculate the min and max in x- and y-direction (rotated, i.e. same direction as z_grid)
  double min_x = corners[0].x;
  double min_y = corners[0].y;
  double max_x = min_x;
  double max_y = min_y;

  for (size_t four = 1 ; four < 4 ; four++) {
    if      (corners[four].x < min_x) { min_x = corners[four].x ;}
    else if (corners[four].x > max_x) { max_x = corners[four].x ;}
    if      (corners[four].y < min_y) { min_y = corners[four].y ;}
    else if (corners[four].y > max_y) { max_y = corners[four].y ;}
  }

  // For loop over all points in z_grid inside the rectangle given by the mins and max' calulated above
  // min_x = dx/2 + dx*n1, if n1 not integer, let it be the smallest integer greater than solution. Zero if negative number
  // max_x = dx/2 + dx*(n2-1), if n2 not integer, let it be the greatest integer smaller than solution. Zero if negative number

  size_t n1 = static_cast<size_t>(max(min_x/dx - 0.5, 0.0));
  size_t n2 = static_cast<size_t>(max(max_x/dx + 1.0, 0.0));
  size_t m1 = static_cast<size_t>(max(min_y/dy - 0.5, 0.0));
  size_t m2 = static_cast<size_t>(max(max_y/dy + 1.0, 0.0));
  size_t m  = z_grid.GetNJ();
  size_t n  = z_grid.GetNI();

  if (n2 > n) n2 = n; //Should stop before grid ends
  if (m2 > m) m2 = m;

  NRLib::BilinearSurface bilinear_corners(corners[0], corners[1], corners[2], corners[3]);
  NRLib::Point point_xy(0.0, 0.0, 0.0);
  NRLib::Point z_dir(0.0, 0.0, 1.0);
  NRLib::Point intersec1, intersec2;

  for (size_t it1 = n1; it1 < n2 ; it1++) {
    for (size_t it2 = m1; it2 < m2 ; it2++) {

      point_xy.x = dx/2 + it1*dx;
      point_xy.y = dy/2 + it2*dy;

      NRLib::Line line_xy(point_xy, point_xy + z_dir, false, false);

      int nu_of_intersec = bilinear_corners.FindIntersections(line_xy,
                                                              intersec1,
                                                              intersec2);

      if (nu_of_intersec > 0) {
        double count = static_cast<double>(is_set(it1,it2));

        if (is_set(it1,it2)) {
          z_grid(it1, it2) *= count/(count + 1.0);
          z_grid(it1, it2) += intersec1.z/(count + 1.0);
          is_set(it1,it2)++;
        }
        else {
          z_grid(it1,it2) = intersec1.z;
          is_set(it1,it2) = 1;
        }
      }
    }
  }
}

void EclipseGeometry::TriangularFillInZValuesInArea(NRLib::Grid2D<double>           & z_grid,
                                                    NRLib::Grid2D<int>              & is_set,
                                                    const std::vector<NRLib::Point> & corners_in,
                                                    const double                      dx,
                                                    const double                      dy,
                                                    const bool                        fixed_triangularization) const
{
  std::vector<NRLib::Point> corners       = corners_in;
  bool                      two_triangles = true;
  NRLib::Triangle           triangle1, triangle2;
  NRLib::Point              vec1(0.0, 0.0, 0.0);
  NRLib::Point              vec2(0.0, 0.0, 0.0);

  //Check if two of the points are the same
  size_t nn = 3;
  for (size_t mm = 0 ; mm < 4 ; mm++) {
    if (corners[mm] == corners[nn]) { //NOTE: Could do better check on whether two points are equal !!!!!!!
      two_triangles = false;
      corners[nn ]  = corners[3];
      triangle1.SetCornerPoints(corners[0], corners[1], corners[2]);
    }
    nn = mm;
  }

  if (two_triangles) {
    // Calculate the sum of two opposite angles
    vec1.x = corners[1].x - corners[0].x;
    vec1.y = corners[1].y - corners[0].y;
    vec2.x = corners[3].x - corners[0].x;
    vec2.y = corners[3].y - corners[0].y;

    double vec1_vec2_angle = vec1.GetAngle(vec2);

    vec1.x = corners[1].x - corners[2].x;
    vec1.y = corners[1].y - corners[2].y;
    vec2.x = corners[3].x - corners[2].x;
    vec2.y = corners[3].y - corners[2].y;

    vec1_vec2_angle += vec1.GetAngle(vec2);

    // Make delunay triangles (according to the sum of the angles)
    //
    // NBNB! To avoid GEOS-29 bug we choose to always use the same triangulation. This is
    // crucial within a given column, and as we do the modelling layer by layer, this is
    // the only way to ensure consistency.
    //
    if (fixed_triangularization || vec1_vec2_angle <= NRLib::Pi) {
      triangle1.SetCornerPoints(corners[3], corners[0], corners[1]);
      triangle2.SetCornerPoints(corners[1], corners[2], corners[3]);
    }
    else {
      triangle1.SetCornerPoints(corners[0], corners[1], corners[2]);
      triangle2.SetCornerPoints(corners[2], corners[3], corners[0]);
    }
  }

  //Calculate the min and max in x- and y-direction (rotated, i.e. same direction as z_grid)
  double min_x = corners[0].x;
  double min_y = corners[0].y;
  double max_x = min_x;
  double max_y = min_y;

  for (size_t four = 1 ; four < 4 ; four++) {
    if      (corners[four].x < min_x)  { min_x = corners[four].x ;}
    else if (corners[four].x > max_x)  { max_x = corners[four].x ;}
    if      (corners[four].y < min_y)  { min_y = corners[four].y ;}
    else if (corners[four].y > max_y)  { max_y = corners[four].y ;}
  }

  // For loop over all points in z_grid inside the rectangle given by the mins and max' calulated above
  // min_x = dx/2 + dx*n1, if n1 not integer, let it be the smallest integer greater than solution. Zero if negative number
  // max_x = dx/2 + dx*(n2-1), if n2 not integer, let it be the greatest integer smaller than solution. Zero if negative number

  double dx_ecl  = max_x - min_x; // Use full length of cell as shift!
  double dy_ecl  = max_y - min_y;

  double delta_x = std::max(dx, 0.5*dx_ecl);
  double delta_y = std::max(dy, 0.5*dy_ecl);

  size_t n1      = static_cast<size_t>(std::max((min_x - delta_x)/dx, 0.0));
  size_t n2      = static_cast<size_t>(std::max((max_x + delta_x)/dx, 0.0));
  size_t m1      = static_cast<size_t>(std::max((min_y - delta_y)/dy, 0.0));
  size_t m2      = static_cast<size_t>(std::max((max_y + delta_y)/dy, 0.0));

  size_t n       = z_grid.GetNI();
  size_t m       = z_grid.GetNJ();

  if (n2 > n) n2 = n; //Should stop before grid ends
  if (m2 > m) m2 = m;

  NRLib::Point point_xy(0.0, 0.0, 0.0);
  NRLib::Point z_dir(0.0, 0.0, 1.0);
  NRLib::Point intersec;

  for (size_t it1 = n1 ; it1 < n2 ; it1++) {
    for (size_t it2 = m1 ; it2 < m2 ; it2++) {
      double count = static_cast<double>(is_set(it1,it2));
      double w1    = count/(count + 1.0);
      double w2    =   1.0/(count + 1.0);

      point_xy.x = dx/2 + it1*dx;
      point_xy.y = dy/2 + it2*dy;

      NRLib::Line line_xy(point_xy, point_xy + z_dir, false, false);

      if (triangle1.FindIntersection(line_xy, intersec, true)) {
        if (is_set(it1,it2) > 0) {
          z_grid(it1, it2) *= w1;
          z_grid(it1, it2) += w2*intersec.z;
          is_set(it1,it2)++;
        }
        else {
          z_grid(it1,it2) = intersec.z;
          is_set(it1,it2) = 1;
        }
      }
      else if (two_triangles) {
        if (triangle2.FindIntersection(line_xy, intersec, true)) {
          if (is_set(it1,it2) > 0) {
            z_grid(it1,it2) *= w1;
            z_grid(it1,it2) += w2*intersec.z;
            is_set(it1,it2)++;
          }
          else {
            z_grid(it1,it2) = intersec.z;
            is_set(it1,it2) = 1;
          }
        }
      }
    }
  }
}

void EclipseGeometry::TriangularFillInZValuesAtEdges(NRLib::Grid2D<double>           & z_grid,
                                                     NRLib::Grid2D<int>              & is_set,
                                                     const std::vector<NRLib::Point> & corners,
                                                     const double                      dx,
                                                     const double                      dy) const
{
  //Calculate the min and max in x- and y-direction (rotated, i.e. same direction as z_grid)
  double min_x = corners[0].x;
  double min_y = corners[0].y;
  double max_x = min_x;
  double max_y = min_y;

  for (size_t four = 1 ; four < 4 ; four++) {
    if      (corners[four].x < min_x)  { min_x = corners[four].x ;}
    else if (corners[four].x > max_x)  { max_x = corners[four].x ;}
    if      (corners[four].y < min_y)  { min_y = corners[four].y ;}
    else if (corners[four].y > max_y)  { max_y = corners[four].y ;}
  }

  // For loop over all points in z_grid inside the rectangle given by the mins and max' calulated above
  // min_x = dx/2 + dx*n1, if n1 not integer, let it be the smallest integer greater than solution. Zero if negative number
  // max_x = dx/2 + dx*(n2-1), if n2 not integer, let it be the greatest integer smaller than solution. Zero if negative number

  double dx_ecl  = max_x - min_x; // Use full length of cell as shift!
  double dy_ecl  = max_y - min_y;

  double delta_x = std::max(dx, 0.5*dx_ecl);
  double delta_y = std::max(dy, 0.5*dy_ecl);

  size_t n1      = static_cast<size_t>(std::max((min_x - delta_x)/dx, 0.0));
  size_t n2      = static_cast<size_t>(std::max((max_x + delta_x)/dx, 0.0));
  size_t m1      = static_cast<size_t>(std::max((min_y - delta_y)/dy, 0.0));
  size_t m2      = static_cast<size_t>(std::max((max_y + delta_y)/dy, 0.0));

  size_t n       = z_grid.GetNI();
  size_t m       = z_grid.GetNJ();

  if (n2 > n) n2 = n; //Should stop before grid ends
  if (m2 > m) m2 = m;

  //
  // Fill edge cells that are still unset
  //
  for (size_t it1 = n1 ; it1 < n2 ; it1++) {
    for (size_t it2 = m1 ; it2 < m2 ; it2++) {

      if (is_set(it1,it2) == 0) {

        //xxx
        //std::cout << "it1, it2 = " << it1 << "  " << it2 << std::endl;

        //
        // Take depth value from nearest corner
        //
        int    min_indx = -1;
        double min_dist = std::numeric_limits<double>::max();
        for (size_t i = 0 ; i < 4 ; i++) {
          double distx = corners[i].x - it1*dx;
          double disty = corners[i].y - it2*dy;
          double dist2 = distx*distx + disty*disty;

          //xxxx
          //printf("%4lu  corners[i]. (x,y,z)  = %7.2f, %7.2f,   %7.2f\n",i, corners[i].x, corners[i].y, corners[i].z);
         //printf("            it1*dx it2*dy  = %7.2f  %7.2f\n",it1*dx, it2*dy);
          //printf("               sqrt(dist2) = %7.2f\n", std::sqrt(dist2));

          if (dist2 < min_dist) {
            min_dist = dist2;
            min_indx = i;
          }
        }

        //xxxx
        //std::cout << "dist, min_indx , z = " << sqrt(min_dist) << " "<< min_indx << " " <<   corners[min_indx].z << std::endl;
        //printf("%7.2f ", corners[min_indx].z);

        z_grid(it1,it2) = corners[min_indx].z;
        is_set(it1,it2) = 1;
      }
    }
    //printf("\n");
  }
  //printf("\n");
}

//----------------------------------------------------------------------------------------------------------------
void EclipseGeometry::FindRegularGridOfZValues(NRLib::StormContGrid                & zgrid,
                                               const NRLib::RegularSurface<double> & topeclipse,
                                               const NRLib::RegularSurface<double> & boteclipse,
                                               const size_t                          top_k,
                                               const size_t                          n_threads,
                                               const bool                            cornerpoint_interpolation,
                                               const bool                            interpolation_at_faults,
                                               const bool                            bilinear_else_triangles,
                                               const bool                            fixed_triangularization,
                                               const bool                            horizontal_interpolation,
                                               const bool                            vertical_interpolation,
                                               const double                          missingValue) const
//----------------------------------------------------------------------------------------------------------------
{
  const double xmin  = zgrid.GetXMin();
  const double ymin  = zgrid.GetYMin();
  const double dx    = zgrid.GetDX();
  const double dy    = zgrid.GetDY();
  const double angle = zgrid.GetAngle();

  const size_t ni    = zgrid.GetNI();
  const size_t nj    = zgrid.GetNJ();
  const size_t nk    = zgrid.GetNK();

  float        monitor_size;
  float        next_monitor;
  float        count;
  std::string  carets;
  std::string  rest_carets;

  //
  // Setup layers to fill
  //
  std::vector<NRLib::Grid2D<double> > layer(nk);
  for (size_t k = 0; k < nk; k++) {
    layer[k] = NRLib::Grid2D<double>(ni, nj, missingValue);
  }

  //
  // Fill base layer
  //
  if (cornerpoint_interpolation)
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nExtracting z-values from Eclipse grid using corner-point interpolation.\n");
  else
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nExtracting z-values from Eclipse grid using center-point interpolation.\n");

  FindLayer(layer[nk - 1],
            nk - 2 + top_k,
            1,              // Use lower corner
            dx,
            dy,
            xmin,
            ymin,
            angle,
            cornerpoint_interpolation,
            interpolation_at_faults,   // For corner point interpolation only
            bilinear_else_triangles,
            fixed_triangularization,
            horizontal_interpolation,
            missingValue);

  //
  // Fill rest of the layers
  //
  MonitorInitialize(nk - 1, monitor_size, next_monitor, count, carets, rest_carets);

#ifdef WITH_OMP
  int chunk_size = 1;
#pragma omp parallel for schedule(dynamic, chunk_size) num_threads(n_threads)
#endif
  for (int k = static_cast<int>(nk - 2) ; k >= 0 ; --k) {
    FindLayer(layer[k],
              k + top_k,
              0,          // Use upper corner
              dx,
              dy,
              xmin,
              ymin,
              angle,
              cornerpoint_interpolation,
              interpolation_at_faults,   // For corner point interpolation only
              bilinear_else_triangles,
              fixed_triangularization,
              horizontal_interpolation,
              missingValue);
    Monitor(nk - 1, nk - 2 - k,  monitor_size, next_monitor, count, carets);
  }
  MonitorFinish(rest_carets);

  if (vertical_interpolation) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nInterpolating z-values that are undefined using vertical interpolation.\n");
    VerticalInterpolation(layer,
                          zgrid,
                          topeclipse,
                          boteclipse,
                          ni, nj, nk,
                          missingValue);
  }
  else {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nVertical interpolation of undefined z-values is not applied!\n");
  }

  //
  // Copy layer to grid
  //
  for (size_t k = 0 ; k < nk ; k++) {
    for (size_t i = 0 ; i < ni ; i++) {
      for (size_t j = 0 ; j < nj ; j++) {
        zgrid(i, j, k) = static_cast<float>(layer[k](i, j));
      }
    }
  }
}

//------------------------------------------------------------------------------------------------
void EclipseGeometry::VerticalInterpolation(std::vector<NRLib::Grid2D<double> > & layer,
                                            const NRLib::StormContGrid          & zgrid,
                                            const NRLib::RegularSurface<double> & topeclipse,
                                            const NRLib::RegularSurface<double> & boteclipse,
                                            const size_t                          ni,
                                            const size_t                          nj,
                                            const size_t                          nk,
                                            const double                          missing) const
//------------------------------------------------------------------------------------------------
{
  std::vector<size_t> indices(nk);
  std::vector<std::pair<size_t, size_t> > range;

  size_t count = 0;

  for (size_t i = 0 ; i < ni ; i++) {
    for (size_t j = 0 ; j < nj ; j++) {
      //
      // If first/last grid cell is unset/missing set from surface
      //
      double x,y,z;
      zgrid.FindCenterOfCell(i, j, 0, x, y, z); // Using k = 0

      double ztop = topeclipse.GetZ(x, y);
      double zbot = boteclipse.GetZ(x, y);

      assert(ztop != missing);
      assert(zbot != missing);

      int i1 = -1;
      int i2 = -1;

      for (size_t k = 0 ; k < nk ; k++) {
        if (layer[k](i, j) != missing) {
          if (i1 == -1)
            i1 = k;
          i2 = k; // Always set this
        }
      }

      if (layer[0](i, j) == missing) {
        if (i1 != -1 && layer[i1](i, j) < ztop)
          layer[0](i, j) = layer[i1](i, j);
        else
          layer[0](i, j) = ztop;
      }
      if (layer[nk-1](i, j) == missing) {
        if (i2 != -1 && layer[i2](i, j) > zbot)
          layer[nk-1](i, j) = layer[i2](i, j);
        else
          layer[nk-1](i, j) = zbot;
      }


      //
      // Find all unset/missing values
      //
      indices.clear();
      for (size_t k = 0 ; k < nk ; k++) {
        if (layer[k](i, j) == missing) {
          indices.push_back(k);
        }
      }

      if (!indices.empty()) {
        //
        // Find range of unset/missing values
        //
        range.clear();

        if (indices.size() == 1) { // There is only one unset grid cell
          range.push_back(std::make_pair(indices[0], indices[0])); // Single layer with unset value
        }
        else {
          size_t kstart = indices[0];
          size_t kprev  = indices[0];

          for (size_t k = 1; k < indices.size(); k++) {
            if (indices[k] != kprev + 1) {                           // Next undefined values is not in series / last series has only one undef
              range.push_back(std::make_pair(kstart, kprev));        // Layers with unset values
              kstart = indices[k];                                   // Start of next undefined series
              if (k == indices.size() - 1) {                         // The last undef is single
                range.push_back(std::make_pair(kstart, kstart));     // Layers with unset values
              }
            }
            else if (k == indices.size() - 1) {                      // Special treatment for the last undef series
              range.push_back(std::make_pair(kstart, indices[k]));   // Layers with unset values
            }
            kprev = indices[k];
          }
        }
        count += indices.size();

        //
        // Find values used for interpolation
        //
        for (size_t k = 0 ; k < range.size() ; k++) {
          size_t k1 = range[k].first;
          size_t k2 = range[k].second;
          double z1 = layer[k1 - 1](i, j);            // Interpolation top value is grid cell value layer above
          double z2 = layer[k2 + 1](i, j);            // Interpolation top value is grid cell value layer below
          double dz = (z2 - z1)/static_cast<double>(k2 - k1 + 2);

          //
          // Do interpolation
          //
          for (size_t kk = k1 ; kk <= k2 ; kk++) {
            layer[kk](i, j) = z1 + dz*(kk - k1 + 1);
          }
        }
      }
    }
  }
  if (count > 0) {
    double ntot = static_cast<double>(ni*nj*nk);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n  %5.2f%% of the cells were interpolated.\n", 100.0*count/ntot);
  }
  else {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n  No interpolation needed.\n");
  }
}

//----------------------------------------------------------------------
void EclipseGeometry::MonitorInitialize(size_t        n,
                                        float       & monitor_size,
                                        float       & next_monitor,
                                        float       & count,
                                        std::string & carets,
                                        std::string & rest_carets) const
//----------------------------------------------------------------------
{
  count        = 1.0f;
  monitor_size = static_cast<float>(n)*0.02f;
  carets       = "^";
  rest_carets  =  "";

  if (monitor_size < 1.0f) {
    int n_carets = static_cast<int>(floor(1.0f/monitor_size));
    int m_carets = 50 - n*n_carets;
    carets       = std::string(n_carets,'^');
    rest_carets  = std::string(m_carets,'^');
    monitor_size = 1.0f;
  }
  next_monitor = monitor_size;

  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |"
    << "\n  ^"
    << std::flush;
}

//-------------------------------------------------------
void EclipseGeometry::Monitor(size_t        n,
                              size_t        k,
                              float         monitor_size,
                              float       & next_monitor,
                              float       & count,
                              std::string & carets) const
//-------------------------------------------------------
{
#ifdef PARALLEL
#pragma omp critical
#endif
  {
    if (count >= next_monitor || k == n - 1) {
      next_monitor += monitor_size;
      std::cout << carets << std::flush;
    }
    count += 1.0f;
  }
}

//------------------------------------------------------------------------
void EclipseGeometry::MonitorFinish(const std::string & rest_carets) const
//------------------------------------------------------------------------
{
  std::cout << rest_carets << std::endl;;
}

//
// We can use the method below when we extrapolate top and base surfaces of Eclipsegrid ...
//
void EclipseGeometry::FindLayer(NRLib::Grid2D<double> & z_grid,
                                const size_t            k,
                                const int               lower_or_upper,
                                const double            dx,
                                const double            dy,
                                const double            x0,
                                const double            y0,
                                const double            angle,
                                const bool              cornerpoint_interpolation,
                                const bool              interpolation_at_faults,
                                const bool              bilinear_else_triangles,
                                const bool              fixed_triangularization,
                                const bool              horizontal_interpolation,
                                const double            missingValue) const
{
  if (cornerpoint_interpolation)
    FindLayerCornerPointInterpolation(z_grid,
                                      k,
                                      lower_or_upper,
                                      dx,
                                      dy,
                                      x0,
                                      y0,
                                      angle,
                                      interpolation_at_faults,
                                      bilinear_else_triangles,
                                      fixed_triangularization,
                                      horizontal_interpolation,
                                      missingValue); // Currently not in use
  else
    FindLayerCenterPointInterpolation(z_grid,
                                      k,
                                      lower_or_upper,
                                      dx,
                                      dy,
                                      x0,
                                      y0,
                                      angle,
                                      bilinear_else_triangles,
                                      fixed_triangularization,
                                      horizontal_interpolation,
                                      missingValue);
}

// Corner point interpolation.  This routine does not work with reverse faults.
void EclipseGeometry::FindLayerCornerPointInterpolation(NRLib::Grid2D<double> & z_grid,
                                                        const size_t            k,
                                                        const int               lower_or_upper,
                                                        const double            dx,
                                                        const double            dy,
                                                        const double            x0,
                                                        const double            y0,
                                                        const double            angle,
                                                        const bool              interpolation_at_faults,
                                                        const bool              bilinear_else_triangles,
                                                        const bool              fixed_triangularization,
                                                        const bool              horizontal_interpolation,
                                                        const double            missingValue) const
{
  double                    cosA = cos(angle);
  double                    sinA = sin(angle);
  double                    thr  = 0.00000001;
  NRLib::Point              C0, C1, C2, C3;
  NRLib::Point              prev_upper_corner,prev_lower_corner, cell_under_right_corner, cell_under_left_corner;
  NRLib::Grid2D<int>        is_set(z_grid.GetNI(), z_grid.GetNJ(), 0);
  std::vector<NRLib::Point> fault_corners(4);
  std::vector<NRLib::Point> Crot(4);

  for (size_t j = 0 ; j < nj_ ; j++) { // Loops over each cell in the given layer
    for (size_t i = 0 ; i < ni_ ; i++) {

      bool surface_edge = horizontal_interpolation && (i == 0 || j == 0 || i == ni_ - 2 || j == nj_ - 2);

      C0 = FindCornerPoint(i, j, k, 0, 0, lower_or_upper); // Find rotated coordinates for the corners of the cell
      C1 = FindCornerPoint(i, j, k, 0, 1, lower_or_upper);
      C2 = FindCornerPoint(i, j, k, 1, 1, lower_or_upper);
      C3 = FindCornerPoint(i, j, k, 1, 0, lower_or_upper);

      TranslateAndRotate(Crot[0], C0, x0, y0, cosA, sinA);
      TranslateAndRotate(Crot[1], C1, x0, y0, cosA, sinA);
      TranslateAndRotate(Crot[2], C2, x0, y0, cosA, sinA);
      TranslateAndRotate(Crot[3], C3, x0, y0, cosA, sinA);

      if (FindTopCell(i,j) != nk_) {
        FillInZValuesInArea(z_grid, is_set, Crot, dx, dy, bilinear_else_triangles, surface_edge, fixed_triangularization, true);
      }

      if (interpolation_at_faults) {
        if (j > 0) {
          C0 = FindCornerPoint(i, j - 1, k, 0, 1, lower_or_upper); // Find rotated coordinates for the corners of the cell under
          C1 = FindCornerPoint(i, j - 1, k, 1, 1, lower_or_upper);

          TranslateAndRotate(cell_under_left_corner , C0, x0, y0, cosA, sinA);
          TranslateAndRotate(cell_under_right_corner, C1, x0, y0, cosA, sinA);

          bool diff1 = abs(cell_under_left_corner.y /Crot[0].y - 1.0) > thr;
          bool diff2 = abs(cell_under_right_corner.y/Crot[3].y - 1.0) > thr;
          bool diff3 = abs(cell_under_right_corner.x/Crot[3].x - 1.0) > thr;
          bool diff4 = abs(cell_under_left_corner.x /Crot[0].x - 1.0) > thr;

          if (diff1 || diff2 || diff3 || diff4) {  // Fault along the i-coordinate
            fault_corners[0] = cell_under_left_corner;
            fault_corners[1] = Crot[0];
            fault_corners[2] = Crot[3];
            fault_corners[3] = cell_under_right_corner;

            if (Crot[0].y < cell_under_left_corner.y && Crot[3].y < cell_under_right_corner.y){
              fault_corners[0].z = Crot[0].z;
              fault_corners[1].z = cell_under_left_corner.z;
              fault_corners[2].z = cell_under_right_corner.z;
              fault_corners[3].z = Crot[3].z;
            }
            FillInZValuesInArea(z_grid, is_set, fault_corners, dx, dy, bilinear_else_triangles, surface_edge, fixed_triangularization, false);
          }
        }

        if (i > 0) {
          bool diff1 = abs(prev_upper_corner.x/Crot[1].x - 1.0) > thr;
          bool diff2 = abs(prev_lower_corner.x/Crot[0].x - 1.0) > thr;
          bool diff3 = abs(prev_upper_corner.y/Crot[1].y - 1.0) > thr;
          bool diff4 = abs(prev_lower_corner.y/Crot[0].y - 1.0) > thr;

          if (diff1 || diff2 || diff3 || diff4) { // Fault along the j-coordinate
            fault_corners[0] = prev_lower_corner;
            fault_corners[1] = prev_upper_corner;
            fault_corners[2] = Crot[1];
            fault_corners[3] = Crot[0];

            if (Crot[1].x < prev_upper_corner.x && Crot[0].x < prev_lower_corner.x) {
              fault_corners[0].z = Crot[1].z;
              fault_corners[1].z = Crot[0].z;
              fault_corners[2].z = prev_upper_corner.z;
              fault_corners[3].z = prev_lower_corner.z;

            }
            FillInZValuesInArea(z_grid, is_set, fault_corners, dx, dy, bilinear_else_triangles, surface_edge, fixed_triangularization, false);
          }
        }
        prev_upper_corner = Crot[2];
        prev_lower_corner = Crot[3];
      }
    }
  }
  //
  // We do a horizontal extrapolation/interpolation for the top and base surfaces, but not
  // for grid layers. These will be interpolated vertically to ensure vertical consistency
  //
  if (horizontal_interpolation) {
    bool iterate = true;
    while (iterate) {
      FillInZValuesByAveraging(z_grid, is_set, iterate);
    }
  }
}

void EclipseGeometry::FillInZValuesInArea(NRLib::Grid2D<double>           & z_grid,
                                          NRLib::Grid2D<int>              & is_set,
                                          const std::vector<NRLib::Point> & corners,
                                          const double                      dx,
                                          const double                      dy,
                                          const bool                        bilinear_else_triangles,
                                          const bool                        surface_edge,
                                          const bool                        fixed_triangularization,
                                          const bool                        write_warning) const
{
  if (bilinear_else_triangles) {
    if (write_warning) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "WARNING: You have asked for bilinear interpolation to be used. Please note that this\n");
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "         option may show edge effects. If so, use the default triangular interpolation.\n");
    }
    BilinearFillInZValuesInArea(z_grid, is_set, corners, dx, dy);
  }
  else {
    TriangularFillInZValuesInArea(z_grid, is_set, corners, dx, dy, fixed_triangularization);
    if (surface_edge) {
      TriangularFillInZValuesAtEdges(z_grid, is_set, corners, dx, dy);
    }
  }
}

// center point interpolation
void EclipseGeometry::FindLayerCenterPointInterpolation(NRLib::Grid2D<double> & z_grid,
                                                        const size_t            k,
                                                        const int               lower_or_upper,
                                                        const double            dx,
                                                        const double            dy,
                                                        const double            x0,
                                                        const double            y0,
                                                        const double            angle,
                                                        const bool              bilinear_else_triangles,
                                                        const bool              fixed_triangularization,
                                                        const bool              horizontal_interpolation,
                                                        const double            missingValue) const
{
  double                    cosA = cos(angle);
  double                    sinA = sin(angle);
  size_t                    m    = z_grid.GetNJ();
  size_t                    n    = z_grid.GetNI();
  NRLib::Grid2D<int>        is_set(n, m, 0);
  NRLib::Point              C0, C1, C2, C3;
  std::vector<NRLib::Point> Crot(4);

  for (size_t j = 0 ; j < nj_ - 1 ; j++) { // Loops over each i,j trace in Eclipse grid
    for (size_t i = 0 ; i < ni_ - 1 ; i++) {
      bool active = IsColumnActive(i  , j  ) &&
                    IsColumnActive(i+1, j  ) &&
                    IsColumnActive(i  , j+1) &&
                    IsColumnActive(i+1, j+1);

      if (active) {
        bool surface_edge = horizontal_interpolation && (i == 0 || j == 0 || i == ni_ - 2 || j == nj_ - 2);

        //NRLib::LogKit::LogFormatted(NRLib:a:LogKit::Low, "i,j =  %d, %d\n",i,j);

        C0 = FindPointCellSurface(i  , j  , k, lower_or_upper, 0.5, 0.5); // Find centre of eclipse grid cell (i  , j  )
        C1 = FindPointCellSurface(i+1, j  , k, lower_or_upper, 0.5, 0.5); // Find centre of eclipse grid cell (i+1, j  )
        C2 = FindPointCellSurface(i+1, j+1, k, lower_or_upper, 0.5, 0.5); // Find centre of eclipse grid cell (i+1, j+1)
        C3 = FindPointCellSurface(i  , j+1, k, lower_or_upper, 0.5, 0.5); // Find centre of eclipse grid cell (i  , j+1)

        //NRLib::Triangle triangle1(C3, C0, C1);
        //NRLib::Triangle triangle2(C1, C2, C3);
        //triangle1.WriteToFile("triangle1.irap");
        //triangle2.WriteToFile("triangle2.irap");

        TranslateAndRotate(Crot[0], C0, x0, y0, cosA, sinA);
        TranslateAndRotate(Crot[1], C1, x0, y0, cosA, sinA);
        TranslateAndRotate(Crot[2], C2, x0, y0, cosA, sinA);
        TranslateAndRotate(Crot[3], C3, x0, y0, cosA, sinA);

        //NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "i,j,k =  %d, %d, %d\n",i,j,k);
        //NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "UnrotC 1:  %7.2f,  %7.2f, %7.2f\n",k,C0.x,C0.y,C0.z);
        //NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "UnrotC 2:  %7.2f,  %7.2f, %7.2f\n",k,C1.x,C1.y,C1.z);
        //NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "UnrotC 3:  %7.2f,  %7.2f, %7.2f\n",k,C2.x,C2.y,C2.z);
        //NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "UnrotC 4:  %7.2f,  %7.2f, %7.2f\n",k,C3.x,C3.y,C3.z);

        if (FindTopCell(i,j) != nk_) {
          //NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Corner 0:  %7.2f, %7.2f, %7.2f\n", k, Crot[0].x, Crot[0].y, Crot[0].z);
          //NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Corner 1:  %7.2f, %7.2f, %7.2f\n", k, Crot[1].x, Crot[1].y, Crot[1].z);
          //NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Corner 2:  %7.2f, %7.2f, %7.2f\n", k, Crot[2].x, Crot[2].y, Crot[2].z);
          //NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Corner 3:  %7.2f, %7.2f, %7.2f\n", k, Crot[3].x, Crot[3].y, Crot[3].z);

          FillInZValuesInArea(z_grid, is_set, Crot, dx, dy, bilinear_else_triangles, surface_edge, fixed_triangularization, true);
        }
      }
    }
  }

  //
  // We do a horizontal extrapolation/interpolation for the top and base surfaces, but not
  // for grid layers. These will be interpolated vertically to ensure vertical consistency
  //
  if (horizontal_interpolation) {
    bool iterate = true;
    while (iterate) {
      FillInZValuesByAveraging(z_grid, is_set, iterate);
    }
  }
}

void EclipseGeometry::TranslateAndRotate(NRLib::Point       & Crot,
                                         const NRLib::Point & C,
                                         const double         x0,
                                         const double         y0,
                                         const double         cosA,
                                         const double         sinA) const
{
  // Translate and rotate eclipse grid nodes to surface grid
  Crot.x = cosA*(C.x - x0) + sinA*(C.y - y0);
  Crot.y = cosA*(C.y - y0) - sinA*(C.x - x0);
  Crot.z = C.z;
}


void EclipseGeometry::FillInZValuesByAveraging(NRLib::Grid2D<double> & z_grid,
                                               NRLib::Grid2D<int>    & is_set,
                                               bool                  & iterate) const
{
  size_t m         = z_grid.GetNJ();
  size_t n         = z_grid.GetNI();
  size_t average_i = 0;
  size_t average_j = 0;
  size_t count1    = 0;
  for (size_t j = 0; j < m ; j++) {
    for (size_t i = 0; i < n ; i++) {
      if (is_set(i, j)) {
        average_i += i;
        average_j += j;
        count1++;
      }
    }
  }

  if (count1 > 0 && count1 < m*n) {
    average_i = average_i/count1;
    average_j = average_j/count1;

    size_t distance = min(min(average_i, average_j), min(n - 1 - average_i, m - 1 - average_j));
    size_t i    = average_i;
    size_t j    = average_j - 1;

    for (size_t it_dist = 1 ; it_dist <= distance ; it_dist++) {
      for (size_t r = 1 ; r <= 8*it_dist ; r++) {

        if (!is_set(i, j)) {
          double zij   = 0.0;
          size_t count = 0;
          if (i     > 0 && is_set(i - 1, j    )) { zij += z_grid(i - 1, j    ); count++; }
          if (i + 1 < n && is_set(i + 1, j    )) { zij += z_grid(i + 1, j    ); count++; }
          if (j     > 0 && is_set(i    , j - 1)) { zij += z_grid(i    , j - 1); count++; }
          if (j + 1 < m && is_set(i    , j + 1)) { zij += z_grid(i    , j + 1); count++; }
          if (count > 0) {
            z_grid(i, j) = zij/count;
            is_set(i, j) = true;
          }
        }
        if      (r < 2*it_dist) { i++ ;}
        else if (r < 4*it_dist) { j++ ;}
        else if (r < 6*it_dist) { i-- ;}
        else                    { j-- ;}
      }
    }

    j = average_j - distance;

    while (j > 0) {
      i = average_i - distance;
      j--;
      while (i <= average_i + distance) {
        if (!is_set(i, j)) {
          double zij   = 0.0;
          size_t count = 0;
          if (i     > 0 && is_set(i - 1, j    )) { zij += z_grid(i - 1, j    ); count++; }
          if (i + 1 < n && is_set(i + 1, j    )) { zij += z_grid(i + 1, j    ); count++; }
          if (j     > 0 && is_set(i    , j - 1)) { zij += z_grid(i    , j - 1); count++; }
          if (j + 1 < m && is_set(i    , j + 1)) { zij += z_grid(i    , j + 1); count++; }
          if (count > 0) {
            z_grid(i, j) = zij/count;
            is_set(i, j) = true;
          }
        }
        i++;
      }
    }

    j = average_j + distance + 1;

    while (j < m) {
      i = average_i - distance;
      while (i <= average_i + distance) {

        if (!is_set(i, j)) {
          double zij   = 0.0;
          size_t count = 0;
          if (i     > 0 && is_set(i - 1,j    )) { zij += z_grid(i - 1, j    ); count++; }
          if (i + 1 < n && is_set(i + 1,j    )) { zij += z_grid(i + 1, j    ); count++; }
          if (j     > 0 && is_set(i    ,j - 1)) { zij += z_grid(i    , j - 1); count++; }
          if (j + 1 < m && is_set(i    ,j + 1)) { zij += z_grid(i    , j + 1); count++; }
          if (count > 0) {
            z_grid(i,j) = zij/count;
            is_set(i,j) = true;
          }
        }
        i++;
      }
      j++;
    }

    i = average_i - distance;

    while (i > 0) {
      i--;
      j = 0;
      while (j < m) {
        if (!is_set(i,j)) {
          double zij   = 0.0;
          size_t count = 0;
          if (i     > 0 && is_set(i - 1, j    )) { zij += z_grid(i - 1, j    ); count++; }
          if (i + 1 < n && is_set(i + 1, j    )) { zij += z_grid(i + 1, j    ); count++; }
          if (j     > 0 && is_set(i    , j - 1)) { zij += z_grid(i    , j - 1); count++; }
          if (j + 1 < m && is_set(i    , j + 1)) { zij += z_grid(i    , j + 1); count++; }
          if (count > 0) {
            z_grid(i, j) = zij/count;
            is_set(i,j)  = true;
          }
        }
        j++;
      }
    }

    i = average_i + distance + 1;

    while (i < n) {
      j = 0;
      while (j < m) {
        if (!is_set(i, j)) {
          double zij   = 0.0;
          size_t count = 0;
          if (i     > 0 && is_set(i - 1, j    )) { zij += z_grid(i - 1, j    ); count++; }
          if (i + 1 < n && is_set(i + 1, j    )) { zij += z_grid(i + 1, j    ); count++; }
          if (j     > 0 && is_set(i    , j - 1)) { zij += z_grid(i    , j - 1); count++; }
          if (j +1  < m && is_set(i    , j + 1)) { zij += z_grid(i    , j + 1); count++; }
          if (count > 0) {
            z_grid(i, j) = zij/count;
            is_set(i, j) = true;
          }
        }
        j++;
      }
      i++;
    }

    bool empty_nodes = false;
    for (size_t j = 0; j < m ; j++) {
      for (size_t i = 0; i < n ; i++) {
        if (!is_set(i,j)) {
          empty_nodes = true;
          break;
        }
      }
    }
    if (!empty_nodes) {
      iterate = false;
    }
  }
  else {
    iterate = false;
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
