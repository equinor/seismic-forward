// $Id: polygon.cpp 1233 2014-01-21 09:41:46Z anner $

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

#include "polygon.hpp"
#include "../math/constants.hpp"
#include <vector>

using namespace NRLib;
using std::cout;

// Empty constructor
Polygon::Polygon()
{  }


//  Making a Polygon object of given size
Polygon::Polygon( size_t size_in )
 : polygon_points_( size_in )
{
}


// Return number of polygon-points
size_t Polygon::GetSize()  const
{
  return static_cast<int>( polygon_points_.size() );
}


// Add input-point to end of polygon, and increase size of polygon.
void Polygon::AddPoint( const Point& pt_in )
{
  polygon_points_.push_back( pt_in );
}


// Check if given pt_in is inside the polygon by counting the number of edges (on
// the polygon) that a line between (x, y) and (pt_in.x, infinity)
// would intersect. An odd number implies TRUE.
bool Polygon::IsInsidePolygonXY(double x, double y)  const
{
  bool inside = false;
  size_t i;
  size_t n = GetSize();
  Point p1, p2;

  if (GetSize() <= 1)
    return false;

  p1 = polygon_points_[n-1];
  for ( i=0; i < n; ++i )  {
    p2 = polygon_points_[i];

    if ( y <= std::max( p1.y, p2.y ) )  {
      if ( y > std::min( p1.y, p2.y ) )  {
        if ( x < std::max( p1.x, p2.x ) )  {

          if ( x < ( ( y - p1.y )/( p2.y - p1.y )*( p2.x - p1.x ) + p1.x ) )
            inside = !inside;

        }
      }
    }
    p1 = p2;
  }  // end of i-for-loop

  return inside;
}


bool Polygon::IsCounterclockwise()
{
  double doublearea=0;
  for (size_t i=0; i<(polygon_points_.size()-1) ;i++) {
    doublearea+=polygon_points_[i].x*polygon_points_[i+1].y-polygon_points_[i].y*polygon_points_[i+1].x;
  }
  return (doublearea>=0);
}

void Polygon::XYExtent(double& xmin, double& xmax, double& ymin, double& ymax)
{
  xmin = polygon_points_[0].x;
  xmax = xmin;
  ymin = polygon_points_[0].y;
  ymax = ymin;
  for (size_t i = 1; i < polygon_points_.size(); i++){
    if (polygon_points_[i].x < xmin)
      xmin = polygon_points_[i].x;
    else if (polygon_points_[i].x > xmax)
      xmax = polygon_points_[i].x;
    if (polygon_points_[i].y < ymin)
      ymin = polygon_points_[i].y;
    else if (polygon_points_[i].y > ymax)
      ymax = polygon_points_[i].y;
  }
}

bool Polygon::IsConvex()
{
  int clockwisefactor=1; //equal to +1 if the polygon points are ordered counterclockwise, equal to -1 if the polygon points are ordered clockwise
  size_t length=polygon_points_.size();
  if (!IsCounterclockwise())
    clockwisefactor=-1;
  Point vec1, vec2;
  double cross;
  for (size_t i=0; i<length-2; i++){
    vec1=polygon_points_[i+1]-polygon_points_[i];
    vec2=polygon_points_[i+2]-polygon_points_[i+1];
    cross=vec1.x*vec2.y-vec1.y*vec2.x;
    if (cross*clockwisefactor<0.0) //ALLOWS THREE POINTS TO BE ON A LINE
      return false;
  }
  vec1=polygon_points_[length-1]-polygon_points_[length-2];
  vec2=polygon_points_[0]-polygon_points_[length-1];
  cross=vec1.x*vec2.y-vec1.y*vec2.x;
  if (cross*clockwisefactor<0.0)
    return false;
  vec1=polygon_points_[0]-polygon_points_[length-1];
  vec2=polygon_points_[1]-polygon_points_[0];
  cross=vec1.x*vec2.y-vec1.y*vec2.x;
  if (cross*clockwisefactor<0.0)
    return false;
  return true;
}

void Polygon::MakeConvex()
{
  int clockwisefactor=1;
  if (!IsCounterclockwise())
    clockwisefactor=-1;
  size_t first_index, second_index, third_index;
  Point vec1, vec2;
  double cross;
  while (!IsConvex()) {
    third_index=0;
    second_index=polygon_points_.size()-1;
    first_index=second_index-1;
    while (third_index<polygon_points_.size()) {
      vec1=polygon_points_[second_index]-polygon_points_[first_index];
      vec2=polygon_points_[third_index]-polygon_points_[second_index];
      cross=vec1.x*vec2.y-vec1.y*vec2.x;
      if (cross*clockwisefactor<=0.0) { //deletes the point at index second_index and then ends the inner while loop
        size_t it1=second_index;
        while (it1<polygon_points_.size()-1) {
          polygon_points_[it1]=polygon_points_[it1+1];
          it1++;
        }
        polygon_points_.resize(polygon_points_.size()-1);
        third_index=polygon_points_.size();
      }
      else { //Checks the next three indexes
        first_index=second_index;
        second_index=third_index;
        third_index++;
      }
    }
  }
}

void Polygon::MinEnclosingRectangle(double &x,
                                    double &y,
                                    double &length1,
                                    double &length2,
                                    double &angle)
{
  std::vector<Point> points;
  points.reserve(polygon_points_.size());
  for (size_t i=0; i<polygon_points_.size(); i++) {
    points.push_back(Point(polygon_points_[i].x,polygon_points_[i].y,0.0));
  }
  std::vector<Point> calipers; //the rotating calipers
  std::vector<Point> edges; //the next edges
  std::vector<size_t> pointindexes; //the indexes of corner points the calipers intersect at
  std::vector<double> angles; //the angles between the calipers and the edges
  double min_angle=0;
  double rotated_angle=0; //the total angle the calipers have rotated
  double min_area; //the minimum area so far
  double curr_area; //area between current calipers
  Point vec1, vec2;
  double l1, l2;
  for (size_t i=0; i<4; i++) { //initializing s.t. the vectors get the right sizes
    pointindexes.push_back(0);
    edges.push_back(Point());
    calipers.push_back(Point(0.0,0.0,0.0));
    angles.push_back(0);
  }
  for (size_t i=1; i<points.size(); i++) { //computing indexes of points corresponding to min_x, min_y, max_x, max_y
    if (points[i].x<=points[pointindexes[3]].x) //min_x
      pointindexes[3]=i;
    if (points[i].x>=points[pointindexes[1]].x) //max_x
      pointindexes[1]=i;
    if (points[i].y<=points[pointindexes[0]].y) //min_y
      pointindexes[0]=i;
    if (points[i].y>=points[pointindexes[2]].y) //max_y
      pointindexes[2]=i;
  }
  min_area=10*(points[pointindexes[1]].x-points[pointindexes[3]].x)*(points[pointindexes[2]].y-points[pointindexes[0]].y);
  length1 = points[pointindexes[1]].x-points[pointindexes[3]].x;
  length2 = points[pointindexes[2]].y-points[pointindexes[0]].y;
  x=points[pointindexes[3]].x+calipers[3].x*(calipers[0].x*(points[pointindexes[3]].y-points[pointindexes[0]].y)-calipers[0].y*(points[pointindexes[3]].x-points[pointindexes[0]].x));
  y=points[pointindexes[3]].y+calipers[3].y*(calipers[0].x*(points[pointindexes[3]].y-points[pointindexes[0]].y)-calipers[0].y*(points[pointindexes[3]].x-points[pointindexes[0]].x));
  while ( (rotated_angle+min_angle)<NRLib::PiHalf) { //rotate edges, compute area, update if improvement
    rotated_angle+=min_angle;
    for (size_t i=0; i<4; i++) { //rotate the calipers
      calipers[i].x=cos(rotated_angle+i*NRLib::PiHalf);
      calipers[i].y=sin(rotated_angle+i*NRLib::PiHalf);
    }
    size_t j = 3; // j = (i-1) mod 4
    for (size_t i = 0; i < 4; i++) { //calculating the edges and the angles wrt the calipers
      size_t next_index = pointindexes[i];
      while (next_index == pointindexes[i]) {
        next_index = pointindexes[i] + 1;
        if (next_index == points.size())
          next_index = 0; //next_index=pointindexes[i]+1 mod points.size()
        edges[i] = points[next_index]-points[pointindexes[i]];
        angles[i] = calipers[i].GetAngle(edges[i]);
        if ( (Dot(points[pointindexes[i]], calipers[j])) <= Dot(points[next_index], calipers[j]) ) {
          pointindexes[i] = next_index;
        }
        /// \todo Replace 0.00001 with constant.
        if (angles[i] <= 0.00001)
          pointindexes[i] = next_index;
      }
      j = i;
    }
    min_angle=angles[0]; //Calculating the angle of the next rotation
    for (size_t i=1; i<4; i++) {
      if (angles[i]<min_angle) {
        min_angle=angles[i];
      }
    }
    vec1=points[pointindexes[2]]-points[pointindexes[0]];
    l2=Dot(vec1,calipers[1]);
    vec2=points[pointindexes[3]]-points[pointindexes[1]];
    l1=Dot(vec2, calipers[2]);
    curr_area=l1*l2;
    if (curr_area<min_area) {
      min_area=curr_area;
      length1=l1;
      length2=l2;
      x=points[pointindexes[3]].x+calipers[3].x*(calipers[0].x*(points[pointindexes[3]].y-points[pointindexes[0]].y)-calipers[0].y*(points[pointindexes[3]].x-points[pointindexes[0]].x));
      y=points[pointindexes[3]].y+calipers[3].y*(calipers[0].x*(points[pointindexes[3]].y-points[pointindexes[0]].y)-calipers[0].y*(points[pointindexes[3]].x-points[pointindexes[0]].x));
      angle=rotated_angle;
    }
  }
}
