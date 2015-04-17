// $Id: segygeometry.cpp 1199 2013-10-02 08:24:02Z anner $

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

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include "segy.hpp"
#include "commonheaders.hpp"
#include "traceheader.hpp"

#include "../exception/exception.hpp"
#include "../iotools/fileio.hpp"
#include "../iotools/logkit.hpp"
#include "../math/constants.hpp"
#include "../surface/surface.hpp"

#include "../iotools/stringtools.hpp"

const int IMISSING = -99999;

using namespace NRLib;


SegyGeometry::SegyGeometry(std::vector<SegYTrace *> &traces)
{
  if(traces.size() > 1) {//More than one trace allows computation of geometry.
    IL0_         = IMISSING;
    XL0_         = IMISSING;
    minIL_       = IMISSING;
    maxIL_       = IMISSING;
    ILStep_      = IMISSING;
    minXL_       = IMISSING;
    maxXL_       = IMISSING;
    XLStep_      = IMISSING;
    firstAxisIL_ = false;

    size_t ntraces = traces.size();

    size_t ii = 0;
    while(traces[ii]==0)
      ii++;

    // il0,xl0,xl1,il1,xl2,il2 are the three points used to find the geometry.

    int ilMin = traces[ii]->GetInline();
    int ilMax = ilMin;
    int xlMin = traces[ii]->GetCrossline();
    int xlMax = xlMin;
    int il0   = ilMin;
    int xl0   = xlMax;
    int il1   = ilMin;
    int xl1   = xlMax;
    int il2   = ilMin;
    int xl2   = xlMax;

    double dist1sq = 0.0; //Distance from 0 to 1
    size_t ind1    = ii;

    for (size_t i = ii + 1; i < ntraces; i++)
    {
      if (traces[i]!=0)
      {
        int il = traces[i]->GetInline();
        int xl = traces[i]->GetCrossline();
        if (il != IMISSING && xl != IMISSING)
        {
          if (il < ilMin)
            ilMin = il;
          if (il > ilMax)
            ilMax = il;
          if (xl < xlMin)
            xlMin = xl;
          if (xl > xlMax)
            xlMax = xl;
          double distsq = (xl-xl0)*(xl-xl0)+(il-il0)*(il-il0);
          if(distsq > dist1sq) {
            dist1sq = distsq;
            ind1    = i;
          }
        }
      }
    }
    il1 = traces[ind1]->GetInline();
    xl1 = traces[ind1]->GetCrossline();

    int  deltaIL = ilMax-ilMin;
    int  deltaXL = xlMax-xlMin;
    double dist2 = 0.0;
    size_t ind2 = ii;
    int   prevIL = il0;
    int   prevXL = xl0;
    for(size_t i = ii + 1 ; i < ntraces ; i++)
    {
      if (traces[i]!=0)
      {
        int il = traces[i]->GetInline();
        int xl = traces[i]->GetCrossline();
        if (il != IMISSING && xl != IMISSING)
        {
          double dist = sqrt(static_cast<double>((xl-xl0)*(xl-xl0) + (il-il0)*(il-il0)))
                      + sqrt(static_cast<double>((xl-xl1)*(xl-xl1) + (il-il1)*(il-il1)));
          if(dist > dist2) {
            dist2 = dist;
            ind2  = i;
          }
          int delta = abs(static_cast<int>(il)-prevIL);
          if(delta > 0 && delta < deltaIL)
            deltaIL = delta;

          delta = abs(static_cast<int>(xl)-prevXL);
          if(delta > 0 && delta < deltaXL)
            deltaXL = delta;
        }
      }
    }
    il2 = traces[ind2]->GetInline();
    xl2 = traces[ind2]->GetCrossline();

    double dxIL, dyIL, dxXL, dyXL;
    double x0 = traces[ii]->GetCoord1();
    double y0 = traces[ii]->GetCoord2();
    double x1 = traces[ind1]->GetCoord1();
    double y1 = traces[ind1]->GetCoord2();
    double x2 = traces[ind2]->GetCoord1();
    double y2 = traces[ind2]->GetCoord2();

    double dIL1 = static_cast<double>(il1-il0);
    double dXL1 = static_cast<double>(xl1-xl0);
    double dIL2 = static_cast<double>(il2-il0);
    double dXL2 = static_cast<double>(xl2-xl0);
    if(dist2 > sqrt(dist1sq)+0.005) { //Standard case, not a simple straight line along XL, IL or diagonally..
      double denominator;
      if(dXL2 > 0) { //Rarely fails, but may happen with diagonal line.
        denominator = dIL1-(dXL1/dXL2)*dIL2;
        dxIL = (x1 - x0 - dXL1/dXL2*(x2-x0))/denominator;
        dxXL = (x2 - x0 - dxIL*dIL2        )/dXL2;
        dyIL = (y1 - y0 - dXL1/dXL2*(y2-y0))/denominator;
        dyXL = (y2 - y0 - dyIL*dIL2        )/dXL2;
      }
      else {
        denominator = dIL2-(dXL2/dXL1)*dIL1;
        dxIL = (x2 - x0 - dXL2/dXL1*(x1-x0))/denominator;
        dxXL = (x1 - x0 - dxIL*dIL1        )/dXL1;
        dyIL = (y2 - y0 - dXL2/dXL1*(y1-y0))/denominator;
        dyXL = (y1 - y0 - dyIL*dIL1        )/dXL1;
      }
      //Regularize
      if(fabs(dxIL*dxXL+dyIL*dyXL) > 0.01)
        //Do not regularize if orthogonal already.
        Regularize(x0, y0, x1, y1, dIL1, dXL1, dxIL, dyIL, dxXL, dyXL);
    }
    else {
      if(dXL1 == 0 && dXL2 == 0) { //Single crossline or single trace.
        if(dIL1 > 0) {
          dxIL = (x1-x0)/dIL1;
          dyIL = (y1-y0)/dIL1;
        }
        else {
          dxIL = 0;
          dyIL = 0;
        }
        dxXL = 0;
        dyXL = 0;
      }
      else if(dIL1 == 0 && dIL2 == 0) { //Single inline
        dxXL = (x1-x0)/dXL1;
        dyXL = (y1-y0)/dXL1;
        dxIL = 0;
        dyIL = 0;
      }
      else { //Only care about endpoints. Diagonal line/two traces
        //Keep it simple. May make a bizarre geometry, but ok along the line.
        //Associate x with XL, y with IL.
        dxIL = 0;
        dyIL = (y1-y0)/dIL1;
        dxXL = (x1-x0)/dXL1;
        dyXL = 0;
      }
    }

    SetupGeometry(traces[ii]->GetCoord1(), traces[ii]->GetCoord2(), traces[ii]->GetInline(), traces[ii]->GetCrossline(),
                  ilMin, ilMax, xlMin, xlMax, static_cast<float>(dxIL), static_cast<float>(dyIL),
                  static_cast<float>(dxXL), static_cast<float>(dyXL), deltaIL, deltaXL);

    size_t nTraces = nx_ * ny_;
    std::vector<SegYTrace *> tracestmp;
    tracestmp.resize(nTraces);

    for (size_t k = 0; k < nTraces; k++)
      tracestmp[k] = 0;

    float distIL = static_cast<float>(dxIL*dxIL+dyIL*dyIL);
    float distXL = static_cast<float>(dxXL*dxXL+dyXL*dyXL);
    float error_dist = distIL;
    if(error_dist == 0 || (error_dist > distXL && distXL > 0))
      error_dist = distXL;
    error_dist = 0.25f*error_dist; //Compare square distances, use half distance as limit.

    double max_dist = 0;
    for (size_t k = 0; k < ntraces; k++)
    {
      if (traces[k] != 0)
      {
        int til = traces[k]->GetInline();
        int txl = traces[k]->GetCrossline();

        if (til != IMISSING && txl != IMISSING)
        {
          size_t i, j;
          double x = traces[k]->GetX();
          double y = traces[k]->GetY();

          float cx, cy;
          this->FindXYFromILXL(til, txl, cx, cy);
          double dtx  = x - cx;
          double dty  = y - cy;
          double dist = dtx*dtx + dty*dty;
          if(dist> max_dist)
            max_dist = dist;

          if(dist > error_dist) {
            WriteGeometry();
            //WriteILXL(); ILXL information is not available as class members at this point
            std::string txt = "In trace " + NRLib::ToString(k) + ", IL=" + NRLib::ToString(til) + " XL="+NRLib::ToString(txl) +
              " the trace has coordinates ("     + NRLib::ToString(x ,2) + ", " + NRLib::ToString( y,2)
              +") while the geometry predicts (" + NRLib::ToString(cx,2) + ", " + NRLib::ToString(cy,2)+").\n"
              + std::string("\nDetected Segy inline (IL) and cross-line (XL) values of corners: \n")
              + "  (il0 ,xl0) : (" + NRLib::ToString(il0, 2) + ", " + NRLib::ToString(NRLib::ToString(xl0, 2)) + ")\n"
              + "  (il1 ,xl1) : (" + NRLib::ToString(il1, 2) + ", " + NRLib::ToString(NRLib::ToString(xl1, 2)) + ")\n"
              + "  (il2 ,xl2) : (" + NRLib::ToString(il2, 2) + ", " + NRLib::ToString(NRLib::ToString(xl2, 2)) + ")\n"
              + std::string("\nDetected Segy XY corners: \n")
              + "  (x0 ,y0) :   (" + NRLib::ToString(x0, 2) + ", " + NRLib::ToString(NRLib::ToString(y0, 2)) + ")\n"
              + "  (x1 ,y1) :   (" + NRLib::ToString(x1, 2) + ", " + NRLib::ToString(NRLib::ToString(y1, 2)) + ")\n"
              + "  (x2 ,y2) :   (" + NRLib::ToString(x2, 2) + ", " + NRLib::ToString(NRLib::ToString(y2, 2)) + ")\n";
            throw(NRLib::Exception(txt));
          }

          try {
            FindIndex(x, y, i, j);
          }
          catch(NRLib::Exception & ) {
            WriteGeometry();
            //WriteILXL(); ILXL information is not avaialable ass class members at this point
            std::string txt = "In trace "+NRLib::ToString(k)+", IL="+NRLib::ToString(til)+" XL="+NRLib::ToString(txl)+
              " the trace has coordinates ("+NRLib::ToString(x,2)+","+NRLib::ToString(y,2)+") which is outside the geometry.\n"
              + std::string("\nDetected Segy inline (IL) and cross-line (XL) ranges: \n")
              + "  IL : (" + NRLib::ToString(ilMin) + ", " + NRLib::ToString(ilMax) + ")\n"
              + "  XL : (" + NRLib::ToString(xlMin) + ", " + NRLib::ToString(xlMax) + ")\n"
              + std::string("\nDetected Segy XY corners: \n")
              + "  (x0 ,y0):  (" + NRLib::ToString(x0, 2) + " , " + NRLib::ToString(NRLib::ToString(y0, 2)) + ")\n"
              + "  (x1 ,y1):  (" + NRLib::ToString(x1, 2) + " , " + NRLib::ToString(NRLib::ToString(y1, 2)) + ")\n"
              + "  (x2 ,y2):  (" + NRLib::ToString(x2, 2) + " , " + NRLib::ToString(NRLib::ToString(y2, 2)) + ")\n";
            throw(NRLib::Exception(txt));
          }
          tracestmp[i+nx_*j] = traces[k];
        }
      }
    }
    LogKit::LogFormatted(LogKit::High,"\nLargest distance between trace location and centre of assigned grid node: %.2fm\n",max_dist);

    traces.resize(nTraces);
    for (size_t k = 0 ; k < nTraces ; k++) {
      traces[k] = tracestmp[k];
      if (traces[k] != NULL)
        traces[k]->SetTableIndex(k);
    }
    // If (min,max,step)(IL,XL) is directly available in this constructor, we can use this instead.
    FindILXLGeometry();
  }
  else {
    SetupGeometrySingleTrace(traces[0]);
  }
}


void
SegyGeometry::Regularize(double   x0,   double   y0,
                         double   x1,   double   y1,
                         double   dIL1, double   dXL1,
                         double & dxIL, double & dyIL,
                         double & dxXL, double & dyXL)
{
  //Set the lengths, and orthogonalize.
  double lengthXL = sqrt(dxXL*dxXL+dyXL*dyXL);
  double factorXL = lengthXL/6.25;
  double scaleXL  = floor(factorXL+0.5)/factorXL;
  double lengthIL = sqrt(dxIL*dxIL+dyIL*dyIL);
  double factorIL = lengthIL/6.25;
  double scaleIL  = floor(factorIL+0.5)/factorIL;
  double factorILfromXL = floor(factorIL+0.5)/floor(factorXL+0.5);
  double factorXLfromIL = 1.0/factorILfromXL;
  double sign1, sign2, length;
  if(fabs(1-scaleXL) < 0.1 || fabs(1-scaleIL) < 0.1) { //Only regularize if the scale is reasonable
    if(fabs(1-scaleXL) < fabs(1-scaleIL)) { //Choose the least scaling.
      dxXL  *= scaleXL;
      dyXL  *= scaleXL;
      length = floor(factorXL+0.5)*6.25;
      if(dyXL > 0) {
        sign1 = dxIL*dyXL/fabs(dxIL*dyXL);
        sign2 = -sign1;
      }
      else {
        sign2 = dyIL*dxXL/fabs(dyIL*dxXL);
        sign1 = -sign2;
      }
      dxIL = sign1*factorILfromXL*dyXL;
      dyIL = sign2*factorILfromXL*dxXL;
    }
    else {
      dxIL  *= scaleIL;
      dyIL  *= scaleIL;
      length = floor(factorIL+0.5)*6.25*factorXLfromIL; //length is always measured along XL.
      if(dxIL > 0) {
        sign1 = dxIL*dyXL/fabs(dxIL*dyXL);
        sign2 = -sign1;
      }
      else  {
        sign2 = dyIL*dxXL/fabs(dyIL*dxXL);
        sign1 = -sign2;
      }
      dxXL = sign2*factorXLfromIL*dyIL;
      dyXL = sign1*factorXLfromIL*dxIL;
    }
    double r0;
    if(dxXL != 0) {
      r0 = atan(dyXL/dxXL);
      if(dxXL < 0)
        r0 = NRLib::Pi + r0;
    }
    else if(dyXL > 0)
      r0 = 0.5*NRLib::Pi;
    else
      r0 = -0.5*NRLib::Pi;

    double dx = x1 - (x0 + dIL1*dxIL + dXL1*dxXL);
    double dy = y1 - (y0 + dIL1*dyIL + dXL1*dyXL);
    double dist0 = dx*dx+dy*dy;

    double small_delta = 0.00001;
    double r1 = r0+small_delta;
    dxXL = length*cos(r1);
    dyXL = length*sin(r1);
    dxIL = sign1*factorILfromXL*dyXL;
    dyIL = sign2*factorILfromXL*dxXL;
    dx = x1 - (x0 + dIL1*dxIL + dXL1*dxXL);
    dy = y1 - (y0 + dIL1*dyIL + dXL1*dyXL);
    double dist1 = dx*dx+dy*dy;

    double delta_r;
    if(dist1<dist0)
      delta_r = 0.1;
    else {
      r1 = r0-small_delta;
      dxXL = length*cos(r1);
      dyXL = length*sin(r1);
      dxIL = sign1*factorILfromXL*dyXL;
      dyIL = sign2*factorILfromXL*dxXL;
      dx = x1 - (x0 + dIL1*dxIL + dXL1*dxXL);
      dy = y1 - (y0 + dIL1*dyIL + dXL1*dyXL);
      dist1 = dx*dx+dy*dy;
      if(dist1>dist0) //No improvement this way either, give up.
        return;
      else
        delta_r = -0.1;
    }
    double dist2 = dist1;
    double r2    = r1;


    while(dist2 <= dist1) {
      dist1 = dist2;
      r1    = r2;

      r2 = r1 + delta_r;
      dxXL = length*cos(r2);
      dyXL = length*sin(r2);
      dxIL = sign1*factorILfromXL*dyXL;
      dyIL = sign2*factorILfromXL*dxXL;
      dx = x1 - (x0 + dIL1*dxIL + dXL1*dxXL);
      dy = y1 - (y0 + dIL1*dyIL + dXL1*dyXL);
      dist2 = dx*dx+dy*dy;
    }
    dist1 = dist0;
    if(dist2<dist1)
      dist1 = dist2;

    while(dist1 > 1 && r2-r0 > 1e-7) {
      double r3 = 0.75*r0+0.25*r2;
      dxXL = length*cos(r3);
      dyXL = length*sin(r3);
      dxIL = sign1*factorILfromXL*dyXL;
      dyIL = sign2*factorILfromXL*dxXL;
      dx = x1 - (x0 + dIL1*dxIL + dXL1*dxXL);
      dy = y1 - (y0 + dIL1*dyIL + dXL1*dyXL);
      double dist3 = sqrt(dx*dx+dy*dy);

      r1 = 0.25*r0+0.75*r2;
      dxXL = length*cos(r1);
      dyXL = length*sin(r1);
      dxIL = sign1*factorILfromXL*dyXL;
      dyIL = sign2*factorILfromXL*dxXL;
      dx = x1 - (x0 + dIL1*dxIL + dXL1*dxXL);
      dy = y1 - (y0 + dIL1*dyIL + dXL1*dyXL);
      dist1 = sqrt(dx*dx+dy*dy);
      if(dist1 < dist3)
        r0 = 0.5*(r0+r2);
      else {
        dist1 = dist3;
        r1 = r3;
        r2 = 0.5*(r0+r2);
      }
    }
    dxXL = length*cos(r1);
    dyXL = length*sin(r1);
    dxIL = sign1*factorILfromXL*dyXL;
    dyIL = sign2*factorILfromXL*dxXL;
  }
}

void
SegyGeometry::SetupGeometry(double xRef, double yRef, int ilRef, int xlRef,
                            int minIL, int maxIL, int minXL, int maxXL,
                            float dxIL, float dyIL, float dxXL, float dyXL,
                            int deltaIL, int deltaXL)
{
  double lx0 = xRef-xlRef*dxXL-ilRef*dxIL;
  double ly0 = yRef-xlRef*dyXL-ilRef*dyIL;

  std::vector<double> cornerx(4);
  std::vector<double> cornery(4);
  std::vector<double> cornerxl(4);
  std::vector<double> corneril(4);

  cornerx[0] = lx0+(minXL-0.5*deltaXL)*dxXL+(minIL-0.5*deltaIL)*dxIL;
  cornery[0] = ly0+(minXL-0.5*deltaXL)*dyXL+(minIL-0.5*deltaIL)*dyIL;
  cornerx[1] = lx0+(minXL-0.5*deltaXL)*dxXL+(maxIL+0.5*deltaIL)*dxIL;
  cornery[1] = ly0+(minXL-0.5*deltaXL)*dyXL+(maxIL+0.5*deltaIL)*dyIL;
  cornerx[2] = lx0+(maxXL+0.5*deltaXL)*dxXL+(minIL-0.5*deltaIL)*dxIL;
  cornery[2] = ly0+(maxXL+0.5*deltaXL)*dyXL+(minIL-0.5*deltaIL)*dyIL;
  cornerx[3] = lx0+(maxXL+0.5*deltaXL)*dxXL+(maxIL+0.5*deltaIL)*dxIL;
  cornery[3] = ly0+(maxXL+0.5*deltaXL)*dyXL+(maxIL+0.5*deltaIL)*dyIL;

  cornerxl[0] = minXL-0.5*deltaXL;
  corneril[0] = minIL-0.5*deltaIL;
  cornerxl[1] = minXL-0.5*deltaXL;
  corneril[1] = maxIL+0.5*deltaIL;
  cornerxl[2] = maxXL+0.5*deltaXL;
  corneril[2] = minIL-0.5*deltaIL;
  cornerxl[3] = maxXL+0.5*deltaXL;
  corneril[3] = maxIL+0.5*deltaIL;
  double miny, maxx, minx, maxy;
  miny = cornery[0];
  maxx = cornerx[0];
  minx = cornerx[0];
  maxy = cornery[0];
  size_t index = 0;
  size_t index2 = 0;
  size_t index3 = 0;
  size_t index4 = 0;
  size_t i;
  for (i = 1; i < 4; i++)
  {
    if (cornery[i] < miny)
    {
      miny = cornery[i];
      index = i;
    }
    if (cornerx[i] > maxx)
    {
      maxx = cornerx[i];
      index2 = i;
    }
    if (cornerx[i] < minx)
    {
      minx = cornerx[i];
      index3 = i;
    }
    if (cornery[i] > maxy)
    {
      maxy = cornery[i];
      index4 = i;
    }
  }

  //If the grid is oriented along UTMX/UTMY, we may have a problem now,
  //since we have non-unique max and min. Check and fix.
  if(index2 == index || index3 == index || index2 == index4 || index3 == index4) {
    //Is this a single line?
    if((dxXL+dyXL) == 0.0 || (dxIL+dyIL) == 0.0) {
      if(maxx == minx) {
        //For line: If equal x, index and and index2 are lowermost, rest uppermost.
        index2 = index;
        index3 = index4;
      }
      else {
        //For line: If equal y, index and index3 are leftmost, rest rightmost.
        index  = index3;
        index4 = index2;
      }
      assert(index != index4); //Can not compute geometry for a single trace, should be handled earlier.
    }
    else {
      //Square oriented with UTMX/UTMY.
      //index  = lower left
      //index2 = lower right
      //index3 = upper left
      //index4 = upper right
      float dXL = sqrt(dxXL*dxXL+dyXL*dyXL);
      float dl  = sqrt(dxIL*dxIL+dyIL*dyIL);
      if(dl == 0 || (dXL > 0 && dXL < dl))
        dl = dXL;

      //Utlilize fact that sum of index in diagonally opposite corners is 3.
      if(fabs(cornerx[0]-cornerx[1]) < 0.2*dl) {
        //0 and 1 are along same constant x.
        if(cornerx[2] < cornerx[0]) { //2 and 3 are leftmost edge.
          if(cornery[2] < cornery[3])
            index = 2;
          else
            index = 3;
          index2 = index-2;
        }
        else {
          if(cornery[0] < cornery[1])
            index = 0;
          else
            index = 1;
          index2 = index+2;
        }
      }
      else {
        //0 and 1 are along same constant y.
        if(cornery[2] < cornery[0]) { //2 and 3 are bottom edge.
          if(cornerx[2] < cornerx[3])
            index = 2;
          else
            index = 3;
          index2 = 5-index;
        }
        else {
          if(cornerx[0] < cornerx[1])
            index = 0;
          else
            index = 1;
          index2 = 1-index;
        }
      }
      index3 = 3-index2;
      index4 = 3-index;
    }
  }

  LogKit::LogMessage(LogKit::DebugLow,"Corner 0:  ("+NRLib::ToString(cornerx[index]) +" , "+NRLib::ToString(cornery[index]) +").\n");
  LogKit::LogMessage(LogKit::DebugLow,"Corner 1:  ("+NRLib::ToString(cornerx[index2])+" , "+NRLib::ToString(cornery[index2])+").\n");
  LogKit::LogMessage(LogKit::DebugLow,"Corner 2:  ("+NRLib::ToString(cornerx[index3])+" , "+NRLib::ToString(cornery[index3])+").\n");
  LogKit::LogMessage(LogKit::DebugLow,"Corner 3:  ("+NRLib::ToString(cornerx[index4])+" , "+NRLib::ToString(cornery[index4])+").\n");

  double deltax = cornerx[index]-cornerx[index2];
  double deltay = cornery[index]-cornery[index2];
  if(deltax == 0)
    rot_ = 0.0;
  else
    rot_ = atan(deltay /deltax);

  double lx,ly;
  if (std::abs(rot_)<=0.25*NRLib::Pi)
  {
    x0_ = cornerx[index];
    y0_= cornery[index];
    lx = sqrt(deltay*deltay+deltax*deltax);
    //lx = sqrt((cornery[index]-cornery[index2])*(cornery[index]-cornery[index2])+(cornerx[index]-cornerx[index2])*(cornerx[index]-cornerx[index2]));
    if(deltaXL!=0)
      nx_ = static_cast<size_t>(std::abs(cornerxl[index]-cornerxl[index2])/deltaXL);
    if (nx_==0 || deltaXL ==0)
      nx_ = static_cast<size_t>(std::abs(corneril[index]-corneril[index2])/deltaIL);
    if (lx==0)
      nx_ = 1;
    dx_ = static_cast<double>(lx/nx_);

    deltax = cornerx[index]-cornerx[index3];
    deltay = cornery[index]-cornery[index3];
    ly = sqrt(deltay*deltay+deltax*deltax);
    //ly = sqrt((cornery[index]-cornery[index3])*(cornery[index]-cornery[index3])+(cornerx[index]-cornerx[index3])*(cornerx[index]-cornerx[index3]));
    if(deltaXL!=0)
      ny_ = static_cast<size_t>(std::abs(cornerxl[index]-cornerxl[index3])/deltaXL);
    if (ny_==0 || deltaXL==0)
      ny_ = static_cast<size_t>(std::abs(corneril[index]-corneril[index3])/deltaIL);
    if (ly==0)
    {
      ny_ = 1;
      ly = dx_;
    }
    dy_ = static_cast<double>(ly/ny_);
    if (index==2 || index==3)
      deltaXL = -deltaXL;
    if (index==1 || index==3)
      deltaIL = -deltaIL;
  }
  else
  {
    x0_ = cornerx[index3];
    y0_= cornery[index3];
    if (index3 == index)
      index = index2;
    deltax = cornerx[index]-cornerx[index3];
    deltay = cornery[index]-cornery[index3];
   // rot_ = atan((cornery[index]-cornery[index3])/(cornerx[index]-cornerx[index3]));
    rot_ = atan(deltay/deltax);
    //lx = sqrt((cornery[index]-cornery[index3])*(cornery[index]-cornery[index3])+(cornerx[index]-cornerx[index3])*(cornerx[index]-cornerx[index3]));
    lx = sqrt(deltay*deltay+deltax*deltax);
    nx_ = static_cast<size_t>(std::abs(cornerxl[index]-cornerxl[index3])/deltaXL);
    if (nx_ == 0 || deltaXL == 0)
      nx_ = static_cast<size_t>(std::abs(corneril[index]-corneril[index3])/deltaIL);
    if (lx == 0)
      nx_ = 1;
    dx_ = static_cast<double>(lx / nx_);
    if (index3 == index4)
      index3 = index;
    deltax = cornerx[index4]-cornerx[index3];
    deltay = cornery[index4]-cornery[index3];
    ly = sqrt(deltay*deltay+deltax*deltax);
    //ly = sqrt((cornery[index4]-cornery[index3])*(cornery[index4]-cornery[index3])+(cornerx[index4]-cornerx[index3])*(cornerx[index4]-cornerx[index3]));
    ny_ = static_cast<size_t>(std::abs(cornerxl[index4]-cornerxl[index3])/deltaXL);
    if (ny_ == 0 || deltaXL == 0)
      ny_ = static_cast<size_t>(std::abs(corneril[index4]-corneril[index3])/deltaIL);
    if (ly == 0)
    {
      ny_ = 1;
      ly = dx_; // avoid 0
    }
    dy_ = static_cast<double>(ly / ny_);
    if (index3 == 2 || index3 == 3)
      deltaXL = -deltaXL;
    if (index3 == 1 || index3 == 3)
      deltaIL = -deltaIL;
  }
  sin_rot_ = sin(rot_);
  cos_rot_ = cos(rot_);
  if ((nx_==1 || ny_==1) && std::abs(rot_)<=0.25*NRLib::Pi)
  {
    //Single line, either inline or crossline.
    if (deltaXL == 0) { //Inline
      xl_stepX_ = 0;
      xl_stepY_ = 0;
      //Straight line, underdetermined, set only one coefficient.
      if (dxIL > dyIL) {
        il_stepX_ = 1/dxIL;
        il_stepY_ = 0;
      }
      else {
        il_stepY_ = 1/dyIL;
        il_stepX_ = 0;
      }
      in_line0_    = ilRef-il_stepX_*(xRef-x0_)-il_stepY_*(yRef-y0_);
      cross_line0_ = xlRef; //Constant, so use one seen.
    }
    else { //Crossline
      il_stepX_ = 0;
      il_stepY_ = 0;
      //Straight line, underdetermined, set only one coefficient.
      if (dxXL > dyXL) {
        xl_stepX_ = 1/dxXL;
        xl_stepY_ = 0;
      }
      else {
        xl_stepY_ = 1/dyXL;
        xl_stepX_ = 0;
      }
      in_line0_    = ilRef; //Constant, so use one seen.
      cross_line0_ = xlRef-xl_stepX_*(xRef-x0_)-xl_stepY_*(yRef-y0_);
    }
  }
  else
  {
    if (dyIL != 0) {
      xl_stepX_ = 1/(dxXL-dxIL*dyXL/dyIL);
      il_stepX_ = -xl_stepX_*dyXL/dyIL;
    }
    else{  //implies that dxXL = 0, dxIL > 0)
      xl_stepX_ = 0;
      il_stepX_ = 1/dxIL;
    }
    if (dxIL != 0) {
      xl_stepY_ = 1/(dyXL-dyIL*dxXL/dxIL);
      il_stepY_ = -xl_stepY_*dxXL/dxIL;
    }
    else {
      xl_stepY_ = 0;
      il_stepY_ = 1/dyIL;
    }
    //IL(x,y) = inLine0_ + (x-x0_)*ilStepX_ + (y-y0)*ilStepY_
    //IL(lx0,ly0) = 0;
    in_line0_    = ilRef-il_stepX_*(xRef-x0_)-il_stepY_*(yRef-y0_);
    cross_line0_ = xlRef-xl_stepX_*(xRef-x0_)-xl_stepY_*(yRef-y0_);
  }
}

void
SegyGeometry::SetupGeometrySingleTrace(const SegYTrace * trace)
{
  dx_          = 0.1; //To ensure that transfroms between indexes and coordinates work.
  dy_          = 0.1;
  x0_          = trace->GetX()-0.5*dx_;
  y0_          = trace->GetY()-0.5*dy_;
  nx_          = 1;
  ny_          = 1;
  il_stepX_    = 1.0/dx_;
  xl_stepX_    = 0.0;
  il_stepY_    = 0.0;
  xl_stepY_    = 1.0/dy_;
  in_line0_    = static_cast<double>(trace->GetInline())-0.5*dx_*il_stepX_;
  cross_line0_ = static_cast<double>(trace->GetCrossline())-0.5*dy_*xl_stepY_;
  sin_rot_     = 0.0;
  cos_rot_     = 1.0;
  rot_         = 0.0;
  IL0_         = trace->GetInline();
  XL0_         = trace->GetCrossline();
  minIL_       = IL0_;
  maxIL_       = IL0_;
  ILStep_      = 1;
  minXL_       = XL0_;
  maxXL_       = XL0_;
  XLStep_      = 1;
  firstAxisIL_ = false;
}

SegyGeometry::SegyGeometry(double x0, double y0, double dx, double dy, size_t nx, size_t ny,
                           double rot)
{
  x0_          = x0;
  y0_          = y0;
  dx_          = dx;
  dy_          = dy;
  nx_          = nx;
  ny_          = ny;
  cos_rot_     = cos(rot);
  sin_rot_     = sin(rot);
  rot_         = rot;
  in_line0_    = -0.5; //Value in corner, gives 0 in center.
  cross_line0_ = -0.5;
  xl_stepX_    = cos_rot_/dx_;
  xl_stepY_    = sin_rot_/dx_;
  il_stepX_    = -sin_rot_/dy_;
  il_stepY_    = cos_rot_/dy_;
  IL0_         = IMISSING;
  XL0_         = IMISSING;
  minIL_       = IMISSING;
  maxIL_       = IMISSING;
  ILStep_      = IMISSING;
  minXL_       = IMISSING;
  maxXL_       = IMISSING;
  XLStep_      = IMISSING;
  firstAxisIL_ = false;
}

SegyGeometry::SegyGeometry(double x0, double y0, double dx, double dy, size_t nx, size_t ny,
                           double IL0, double XL0, double ilStepX, double ilStepY,
                           double xlStepX, double xlStepY, double rot)
{
  x0_          = x0;
  y0_          = y0;
  dx_          = dx;
  dy_          = dy;
  nx_          = nx;
  ny_          = ny;
  in_line0_    = IL0;
  cross_line0_ = XL0;
  il_stepX_    = ilStepX;
  il_stepY_    = ilStepY;
  xl_stepX_    = xlStepX;
  xl_stepY_    = xlStepY;
  cos_rot_     = cos(rot);
  sin_rot_     = sin(rot);
  rot_         = rot;
  IL0_         = IMISSING;
  XL0_         = IMISSING;
  minIL_       = IMISSING;
  maxIL_       = IMISSING;
  ILStep_      = IMISSING;
  minXL_       = IMISSING;
  maxXL_       = IMISSING;
  XLStep_      = IMISSING;
  firstAxisIL_ = false;
}

SegyGeometry::SegyGeometry(const RegularSurface<double> surf)
{
  x0_          = surf.GetXMin();
  y0_          = surf.GetYMin();
  dx_          = surf.GetDX();
  dy_          = surf.GetDY();
  nx_          = surf.GetNI()-1;
  ny_          = surf.GetNJ()-1;
  rot_         = 0;
  cos_rot_     = cos(rot_);
  sin_rot_     = sin(rot_);
  in_line0_    = -0.5; //Value in corner, gives 0 in center.
  cross_line0_ = -0.5;
  xl_stepX_    = cos_rot_/dx_;
  xl_stepY_    = sin_rot_/dx_;
  il_stepX_    = -sin_rot_/dy_;
  il_stepY_    = cos_rot_/dy_;
  IL0_         = IMISSING;
  XL0_         = IMISSING;
  minIL_       = IMISSING;
  maxIL_       = IMISSING;
  ILStep_      = IMISSING;
  minXL_       = IMISSING;
  maxXL_       = IMISSING;
  XLStep_      = IMISSING;
  firstAxisIL_ = false;
}

SegyGeometry::SegyGeometry(const RegularSurfaceRotated<double> surf)
{
  x0_          = surf.GetXRef();  // Use rotation origin and not minimum x-value
  y0_          = surf.GetYRef();  // Use rotation origin and not minimum y-value
  dx_          = surf.GetDX();
  dy_          = surf.GetDY();
  nx_          = surf.GetNI() - 1;
  ny_          = surf.GetNJ() - 1;
  rot_         = surf.GetAngle();
  cos_rot_     = cos(rot_);
  sin_rot_     = sin(rot_);
  in_line0_    = -0.5; //Value in corner, gives 0 in center.
  cross_line0_ = -0.5;
  xl_stepX_    = cos_rot_/dx_;
  xl_stepY_    = sin_rot_/dx_;
  il_stepX_    = -sin_rot_/dy_;
  il_stepY_    = cos_rot_/dy_;
  IL0_         = IMISSING;
  XL0_         = IMISSING;
  minIL_       = IMISSING;
  maxIL_       = IMISSING;
  ILStep_      = IMISSING;
  minXL_       = IMISSING;
  maxXL_       = IMISSING;
  XLStep_      = IMISSING;
  firstAxisIL_ = false;
}

SegyGeometry::SegyGeometry(const SegyGeometry *geometry)
{
  x0_          = geometry->x0_;
  y0_          = geometry->y0_;
  dx_          = geometry->dx_;
  dy_          = geometry->dy_;
  nx_          = geometry->nx_;
  ny_          = geometry->ny_;
  in_line0_    = geometry->in_line0_;
  cross_line0_ = geometry->cross_line0_;
  il_stepX_    = geometry->il_stepX_;
  il_stepY_    = geometry->il_stepY_;
  xl_stepX_    = geometry->xl_stepX_;
  xl_stepY_    = geometry->xl_stepY_;
  sin_rot_     = geometry->sin_rot_;
  cos_rot_     = geometry->cos_rot_;
  rot_         = geometry->rot_;
  IL0_         = geometry->IL0_;
  XL0_         = geometry->XL0_;
  minIL_       = geometry->minIL_;
  maxIL_       = geometry->maxIL_;
  ILStep_      = geometry->ILStep_;
  minXL_       = geometry->minXL_;
  maxXL_       = geometry->maxXL_;
  XLStep_      = geometry->XLStep_;
  firstAxisIL_ = geometry->firstAxisIL_;
}

SegyGeometry::~SegyGeometry()
{
}

bool
SegyGeometry::IsInside(double x, double y) const
{
  float xind, yind;
  return(FindContIndex(x,y,xind,yind));
}

size_t
SegyGeometry::FindIndex(float x, float y) const
{
  size_t i, j;
  FindIndex(x, y, i, j);
  return (j*nx_ + i);
}

void
SegyGeometry::FindIndex(double x, double y, size_t &i, size_t &j) const
{
  double teller1 = (x-x0_)*cos_rot_+(y-y0_)*sin_rot_;
  if (dx_ == 0)
    i = 0;
  else
    i = static_cast<int>(teller1/dx_);
  double teller2 = -(x-x0_)*sin_rot_+(y-y0_)*cos_rot_;
  if (dy_==0)
    j = 0;
  else
    j = static_cast<int>(teller2/dy_);
  if(i>=nx_ || j >= ny_ || i<0 || j<0) {
    throw Exception("Invalid index\n");
  }
}

void
SegyGeometry::FindIndex(int IL, int XL, size_t &i, size_t &j) const
{
  float x,y;
  FindXYFromILXL(IL, XL, x, y);
  FindIndex(x,y,i,j);
}

bool
SegyGeometry::FindContIndex(double x, double y, float &xind, float  &yind) const
{
  xind = static_cast<float>(( (x-x0_)*cos_rot_+(y-y0_)*sin_rot_)/dx_);
  if (dx_==0)
    xind = 0.0;
  yind = static_cast<float>((-(x-x0_)*sin_rot_+(y-y0_)*cos_rot_)/dy_);
  if (dy_==0)
    yind = 0.0;
  if (xind>=0.0f && xind<nx_ && yind>=0.0f && yind<ny_)
    return true;
  else
    return false;
}

void
SegyGeometry::FindILXL(float x, float y, int &IL, int &XL) const
{
  IL = static_cast<int>(0.5+in_line0_+(x-x0_)*il_stepX_+(y-y0_)*il_stepY_);
  XL = static_cast<int>(0.5+cross_line0_+(x-x0_)*xl_stepX_+(y-y0_)*xl_stepY_);
}

void
SegyGeometry::FindILXL(double x, double y, int &IL, int &XL) const
{
  IL = static_cast<int>(0.5+in_line0_+(x-x0_)*il_stepX_+(y-y0_)*il_stepY_);
  XL = static_cast<int>(0.5+cross_line0_+(x-x0_)*xl_stepX_+(y-y0_)*xl_stepY_);
}

void
SegyGeometry::FindILXL(size_t i, size_t j, int &IL, int &XL) const
{
  float x,y;
  FindXYFromIJ(i, j, x, y);
  FindILXL(x,y,IL,XL);
}

void
SegyGeometry::FindContILXL(float x, float y, double &IL, double &XL) const
{
  IL = in_line0_    + (x-x0_)*il_stepX_ + (y-y0_)*il_stepY_;
  XL = cross_line0_ + (x-x0_)*xl_stepX_ + (y-y0_)*xl_stepY_;
}

void
SegyGeometry::FindXYFromIJ(size_t i, size_t j, float &x, float &y) const
{
  double xind = 0.5+static_cast<double>(i);
  double yind = 0.5+static_cast<double>(j);
  x = static_cast<float>(x0_+xind*dx_*cos_rot_-yind*dy_*sin_rot_);
  y = static_cast<float>(y0_+xind*dx_*sin_rot_+yind*dy_*cos_rot_);
}

void
SegyGeometry::FindXYFromILXL(int IL, int XL, float &x, float &y) const
{
  FindXYFromContILXL(static_cast<float>(IL), static_cast<float>(XL), x, y);
}

void
SegyGeometry::FindXYFromContILXL(float IL, float XL, float &x, float &y) const
{
  double xd, yd;
  if (xl_stepX_ == 0.0 && xl_stepY_ == 0.0) {//inline only
    if (il_stepX_ > 0.0) {
      xd = x0_+(IL-in_line0_)/il_stepX_;
      yd = y0_+(xd-x0_)*sin_rot_/cos_rot_; //If cos_rot_ is 0, il_stepX_ should be 0 here.
    }
    else {
      yd = y0_+(IL-in_line0_)/il_stepY_;
      xd = x0_+(yd-y0_)*sin_rot_/cos_rot_;
    }
  }
  else if (il_stepX_ == 0.0 && il_stepY_ == 0.0) {// crossline only
    if (xl_stepX_ > 0.0) {
      xd = x0_+(XL-cross_line0_)/xl_stepX_;
      yd = y0_+(xd-x0_)*sin_rot_/cos_rot_; //If cos_rot_ is 0, xl_stepX_ should be 0 here.
    }
    else {
      yd = y0_+(XL-cross_line0_)/xl_stepY_;
      xd = x0_+(yd-y0_)*sin_rot_/cos_rot_;
    }
  }
  else if (il_stepX_ == 0.0) { //Implies xl_stepY_ == 0)
    xd = x0_+(XL-cross_line0_)/xl_stepX_;
    yd = y0_+(IL-in_line0_)/il_stepY_;
  }
  else {
    yd = y0_+(XL-cross_line0_-(IL-in_line0_)*xl_stepX_/il_stepX_)/(xl_stepY_-il_stepY_*xl_stepX_/il_stepX_);
    xd = x0_+(IL-in_line0_-il_stepY_*(yd-y0_))/il_stepX_;
  }
  x = static_cast<float>(xd);
  y = static_cast<float>(yd);
}



void SegyGeometry::WriteGeometry() const
{
  double geoangle = -rot_*180/(NRLib::Pi);
  if (geoangle < 0.0)
    geoangle += 360.0;

  LogKit::LogFormatted(LogKit::Low,"\n             ReferencePoint      Length Increment\n");
  LogKit::LogFormatted(LogKit::Low,"-------------------------------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"X-coordinate   %12.2f %11.2f %9.2f\n", x0_, nx_*dx_, dx_);
  LogKit::LogFormatted(LogKit::Low,"Y-coordinate   %12.2f %11.2f %9.2f\n", y0_, ny_*dy_, dy_);
  LogKit::LogFormatted(LogKit::Low,"Azimuth             %7.3f\n",geoangle);
}


void SegyGeometry::WriteILXL(bool errorMode) const
{
  LogKit::LogFormatted(LogKit::Low,"\n                 Start      End    Step\n");
  LogKit::LogFormatted(LogKit::Low,"---------------------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"In-line       %8d %8d  %6d\n", minIL_, maxIL_, ILStep_);
  LogKit::LogFormatted(LogKit::Low,"Cross-line    %8d %8d  %6d\n", minXL_, maxXL_, XLStep_);

  if (errorMode) {
    LogKit::LogFormatted(LogKit::Low,"\n                      Start      StepX     StepY\n");
    LogKit::LogFormatted(LogKit::Low,"------------------------------------------------\n");
    LogKit::LogFormatted(LogKit::Low,"In-line[ref]       %8.2f   %8.4f  %8.4f\n", in_line0_   , il_stepX_, il_stepY_);
    LogKit::LogFormatted(LogKit::Low,"Cross-line[ref]    %8.2f   %8.4f  %8.4f\n", cross_line0_, xl_stepX_, xl_stepY_);
  }
}


void
SegyGeometry::FindILXLGeometry(void)
{
  int IL0, XL0, IL1, XL1, IL2, XL2;

  size_t nix = 0; //Just to let the complier know what function to use on next line.

  FindILXL(nix  ,   nix, IL0, XL0);
  FindILXL(nx_-1,     0, IL1, XL1);
  FindILXL(0    , ny_-1, IL2, XL2);

  IL0_ = IL0;
  XL0_ = XL0;

  if (XL1 - XL0 != 0) {
    firstAxisIL_ = false;
    if (nx_ > 1) {  // In case we have read a single trace or a line
      XLStep_ = abs(XL1 - XL0)/(static_cast<int>(nx_) - 1);
    }
    if(XL1 < XL0) {
      minXL_ = XL1;
      maxXL_ = XL0;
    }
    else {
      minXL_ = XL0;
      maxXL_ = XL1;
    }
    if (ny_ > 1) {  // In case we have read a single trace or a line
      ILStep_ = abs(IL2 - IL0)/(static_cast<int>(ny_) - 1);
    }
    if(IL2 < IL0) {
      minIL_ = IL2;
      maxIL_ = IL0;
    }
    else {
      minIL_ = IL0;
      maxIL_ = IL2;
    }
  }
  else {
    firstAxisIL_ = true;
    if (ny_ > 1) {  // In case we have read a single trace or a line
      XLStep_ = abs(XL2 - XL0)/(static_cast<int>(ny_) - 1);
    }
    if(XL2 < XL0) {
      minXL_ = XL2;
      maxXL_ = XL0;
    }
    else {
      minXL_ = XL0;
      maxXL_ = XL2;
    }
    if (nx_ > 1) {  // In case we have read a single trace or a line
      ILStep_ = abs(IL1 - IL0)/(static_cast<int>(nx_) - 1);
    }
    if(IL1 < IL0) {
      minIL_ = IL1;
      maxIL_ = IL0;
    }
    else {
      minIL_ = IL0;
      maxIL_ = IL1;
    }
  }
}


SegyGeometry *
SegyGeometry::GetILXLSubGeometry(const std::vector<int> & ilxl,
                                 bool                   & interpolation,
                                 bool                   & aligned)
{
  int IL0      = IL0_;
  int XL0      = XL0_;
  int minIL    = minIL_;
  int maxIL    = maxIL_;
  int minXL    = minXL_;
  int maxXL    = maxXL_;
  int ILStep   = ILStep_;
  int XLStep   = XLStep_;

  int minSubIL = ilxl[0];
  int maxSubIL = ilxl[1];
  int minSubXL = ilxl[2];
  int maxSubXL = ilxl[3];

  if(minSubIL < minIL)
    throw Exception("Requested IL-start ("+ToString(minSubIL)+") is smaller than minimum IL in seismic cube ("+ToString(minIL)+").\n");
  else if(minSubIL > maxIL)
    throw Exception("Requested IL-start ("+ToString(minSubIL)+") is larger than maximum IL in seismic cube ("+ToString(maxIL)+").\n");

  if(maxSubIL < minIL)
    throw Exception("Requested IL-end ("+ToString(maxSubIL)+") is smaller than minimum IL in seismic cube ("+ToString(minIL)+").\n");
  else if(maxSubIL > maxIL)
    throw Exception("Requested IL-end ("+ToString(maxSubIL)+") is larger than maximum IL in seismic cube ("+ToString(maxIL)+").\n");

  if(minSubXL < minXL)
    throw Exception("Requested XL-start ("+ToString(minSubXL)+") is smaller than minimum XL in seismic cube ("+ToString(minXL)+").\n");
  else if(minSubXL > maxXL)
    throw Exception("Requested XL-start ("+ToString(minSubXL)+") is larger than maximum XL in seismic cube ("+ToString(maxXL)+").\n");

  if (maxSubXL < minXL)
    throw Exception("Requested XL-end ("+ToString(maxSubXL)+") is smaller than minimum XL in seismic cube ("+ToString(minXL)+").\n");
  else if (maxSubXL > maxXL)
    throw Exception("Requested XL-end ("+ToString(maxSubXL)+") is larger than maximum XL in seismic cube ("+ToString(maxXL)+").\n");

  interpolation = false;
  int subILStep = ilxl[4];
  if(subILStep < 0)
    subILStep = ILStep;
  else {
    if(subILStep < ILStep) {
      LogKit::LogFormatted(LogKit::Warning, "\nWARNING: Requested IL-step ("+ToString(subILStep)+") is smaller than IL-step in seismic data ("+ToString(ILStep)+"). Data vill be interpolated.\n");
      interpolation = true;
    }
    else if((subILStep % ILStep) > 0) {
      LogKit::LogFormatted(LogKit::Warning, "\nWARNING: Requested IL-step ("+ToString(subILStep)+") is not a multiple of IL-step in seismic data ("+ToString(ILStep)+"). Data vill be interpolated.\n");
      interpolation = true;
    }
  }
  int subXLStep = ilxl[5];
  if(subXLStep < 0)
    subXLStep = XLStep;
  else {
    if(subXLStep < XLStep) {
      LogKit::LogFormatted(LogKit::Warning, "\nWARNING: Requested XL-step ("+ToString(subXLStep)+") is smaller than XL-step in seismic data ("+ToString(XLStep)+"). Data vill be interpolated.\n");
      interpolation = true;
    }
    else if((subXLStep % XLStep) > 0) {
      LogKit::LogFormatted(LogKit::Warning, "\nWARNING: Requested XL-step ("+ToString(subXLStep)+") is not a multiple of XL-step in seismic data ("+ToString(XLStep)+"). Data vill be interpolated.\n");
      interpolation = true;
    }
  }

  bool modified = false;
  int residual;
  aligned = false;
  if(interpolation == false) { //Make sure initial trace fits existing trace to avoid resampling.
    residual = ((minSubIL-minIL) % ILStep);
    if(residual != 0) { //Match minimum IL to existing line
      modified = true;
      aligned  = true;
      int oldIL = minSubIL;
      if(ILStep-residual-subILStep/2 <= 0)
        minSubIL += ILStep-residual;
      else
        minSubIL -= residual; //Note that this will never be smaller than minIL, due to computation of residual.

      LogKit::LogFormatted(LogKit::Warning, "\nWARNING: Initial minimum IL ("+ToString(oldIL)+
        ") does not match a sampled trace in the input data. Moved to "+ToString(minSubIL)+".\n");
    }
    residual = ((minSubXL-minXL) % XLStep);
    if(residual != 0) { //Match minimum XL to existing line
      modified = true;
      aligned  = true;
      int oldXL = minSubXL;
      if(XLStep-residual-subXLStep/2 <= 0)
        minSubXL += XLStep-residual;
      else
        minSubXL -= residual; //Note that this will never be smaller than minIL, due to computation of residual.
      LogKit::LogFormatted(LogKit::Warning, "\nWARNING: Initial minimum XL ("+ToString(oldXL)+
        ") does not match a sampled trace in the input data. Moved to "+ToString(minSubXL)+".\n");
    }
  }

  //Let last trace fit the inital trace and step.
  residual = ((maxSubIL-minSubIL) % subILStep);
  if(residual != 0) { //Move max to fit.
    modified = true;
    aligned  = true;
    int oldIL = maxSubIL;
    maxSubIL += subILStep-residual; //Try increasing
    if(maxSubIL > maxIL)
      maxSubIL -= subILStep;        //Have to decrease.
    LogKit::LogFormatted(LogKit::Warning, "\nWARNING: Initial maximum IL ("+ToString(oldIL)+
      ") is not an integer number of steps away from minimum IL ("+ToString(minSubIL)+"). Moved to "+ToString(maxSubIL)+".\n");
  }
  residual = ((maxSubXL-minSubXL) % subXLStep);
  if(residual != 0) { //Move max to fit.
    modified = true;
    aligned  = true;
    int oldXL = maxSubXL;
    maxSubXL += subXLStep-residual; //Try increasing
    if(maxSubXL > maxXL)
      maxSubXL -= subXLStep;        //Have to decrease.
    LogKit::LogFormatted(LogKit::Warning, "\nWARNING: Initial maximum IL ("+ToString(oldXL)+
      ") is not an integer number of steps away from minimum IL ("+ToString(minSubXL)+"). Moved to "+ToString(maxSubXL)+".\n");
  }

  int subNx, subNy;
  double subDx, subDy;
  if(firstAxisIL_ == true) {
    subNx = (maxSubIL - minSubIL)/subILStep+1;
    subNy = (maxSubXL - minSubXL)/subXLStep+1;
    subDx = dx_*static_cast<double>(subILStep)/static_cast<double>(ILStep);
    subDy = dy_*static_cast<double>(subXLStep)/static_cast<double>(XLStep);
  }
  else {
    subNx = (maxSubXL - minSubXL)/subXLStep+1;
    subNy = (maxSubIL - minSubIL)/subILStep+1;
    subDx = dx_*static_cast<double>(subXLStep)/static_cast<double>(XLStep);
    subDy = dy_*static_cast<double>(subILStep)/static_cast<double>(ILStep);
  }

  float subIL0, subXL0;
  //Determine reference corner, move from center and out to corner.
  if(IL0 == minIL) { //Reference corner is minimum inline.
    subIL0 = minSubIL-0.5f*subILStep;
  }
  else {             //Reference corner is maximum inline.
    subIL0 = maxSubIL+0.5f*subILStep;
  }
  if(XL0 == minXL) { //Reference corner is minimum crossline.
    subXL0 = minSubXL-0.5f*subXLStep;
  }
  else {             //Reference corner is maximum crossline.
    subXL0 = maxSubXL+0.5f*subXLStep;
  }

  float subX0, subY0;
  FindXYFromContILXL(subIL0, subXL0, subX0, subY0);

  LogKit::LogFormatted(LogKit::Low, "\n                  Inline                Crossline    ");
  LogKit::LogFormatted(LogKit::Low, "\n             Start   End  Step      Start   End  Step");
  LogKit::LogFormatted(LogKit::Low, "\n-----------------------------------------------------\n");
  LogKit::LogFormatted(LogKit::Low, "SEGY-grid    %5d %5d %5d      %5d %5d %5d\n",minIL   , maxIL   , ILStep   , minXL   , maxXL   , XLStep);
  LogKit::LogFormatted(LogKit::Low, "Requested    %5d %5d %5d      %5d %5d %5d\n",ilxl[0] , ilxl[1] , subILStep, ilxl[2] , ilxl[3] , subXLStep);
  if (modified)
  LogKit::LogFormatted(LogKit::Low, "Obtained     %5d %5d %5d      %5d %5d %5d\n",minSubIL, maxSubIL, subILStep, minSubXL, maxSubXL, subXLStep);
  SegyGeometry * result = new SegyGeometry(subX0, subY0, subDx, subDy, subNx, subNy, subIL0, subXL0, il_stepX_, il_stepY_, xl_stepX_, xl_stepY_, rot_);
  return(result);
}

std::vector<int>
SegyGeometry::findAreaILXL(SegyGeometry * tempGeometry)
{
  double x0   = tempGeometry->GetX0();
  double y0   = tempGeometry->GetY0();
  double lx   = tempGeometry->Getlx();
  double ly   = tempGeometry->Getly();
  double rot  = tempGeometry->GetAngle();

  double ax   =  lx*cos(rot);
  double ay   =  lx*sin(rot);
  double bx   = -ly*sin(rot);
  double by   =  ly*cos(rot);

  //
  // Find coordinates of four corners of the "template" geometry
  //
  double x1   = x0 + ax;
  double y1   = y0 + ay;
  double x2   = x0 + bx;
  double y2   = y0 + by;
  double x3   = x0 + ax + bx;  // = (x1 + bx)  = (x2 + ax)
  double y3   = y0 + ay + by;  // = (x1 + by)  = (y2 + ay)

  //
  // Find "continuous" IL and XL values for each corner.
  //
  double il0, il1, il2, il3;
  double xl0, xl1, xl2, xl3;
  FindContILXL(static_cast<float> (x0), static_cast<float> (y0), il0, xl0);
  FindContILXL(static_cast<float> (x1), static_cast<float> (y1), il1, xl1);
  FindContILXL(static_cast<float> (x2), static_cast<float> (y2), il2, xl2);
  FindContILXL(static_cast<float> (x3), static_cast<float> (y3), il3, xl3);

  //
  // Find smallest and largest IL and XL values
  //
  double minIL, maxIL;
  double minXL, maxXL;

  minIL = std::min(il0, il1);
  minIL = std::min(minIL, il2);
  minIL = std::min(minIL, il3);

  maxIL = std::max(il0, il1);
  maxIL = std::max(maxIL, il2);
  maxIL = std::max(maxIL, il3);

  minXL = std::min(xl0, xl1);
  minXL = std::min(minXL, xl2);
  minXL = std::min(minXL, xl3);

  maxXL = std::max(xl0, xl1);
  maxXL = std::max(maxXL, xl2);
  maxXL = std::max(maxXL, xl3);

  int intMinIL = static_cast<int>(floor(minIL));
  int intMaxIL = static_cast<int>(ceil (maxIL));
  int intMinXL = static_cast<int>(floor(minXL));
  int intMaxXL = static_cast<int>(ceil (maxXL));

  //
  // To ensure that the IL XL we find are existing traces
  //
  intMinIL -= (intMinIL - minIL_   ) % ILStep_;
  intMinXL -= (intMinXL - minXL_   ) % XLStep_;
  intMaxIL += (maxIL_   - intMaxIL ) % ILStep_;
  intMaxXL += (maxXL_   - intMaxXL ) % XLStep_;

  //
  // To ensure that we do not step outside seismic data
  //
  intMinIL = std::max(intMinIL, minIL_);
  intMaxIL = std::min(intMaxIL, maxIL_);
  intMinXL = std::max(intMinXL, minXL_);
  intMaxXL = std::min(intMaxXL, maxXL_);

  std::vector<int> areaILXL(6);
  areaILXL[0] = intMinIL;
  areaILXL[1] = intMaxIL;
  areaILXL[2] = intMinXL;
  areaILXL[3] = intMaxXL;
  areaILXL[4] = ILStep_;   // Enforce same step as original seismic data
  areaILXL[5] = XLStep_;   // Enforce same step as original seismic data

  return areaILXL;
}
