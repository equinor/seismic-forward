// $Id: segygeometry.cpp 1678 2017-08-30 12:06:06Z perroe $

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

#include "segygeometry.hpp"
#include "segy.hpp"
#include "commonheaders.hpp"
#include "traceheader.hpp"

#include "../exception/exception.hpp"
#include "../iotools/fileio.hpp"
#include "../iotools/logkit.hpp"
#include "../math/mathutility.hpp"
#include "../math/constants.hpp"
#include "../surface/surface.hpp"

#include "../iotools/stringtools.hpp"

const int IMISSING = -99999;

using namespace NRLib;


SegyGeometry::SegyGeometry(std::vector<SegYTrace *> & traces,
                           bool                       regularize_if_needed)
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
    double x0 = traces[ii]->GetCoord1();   // First defined trace
    double y0 = traces[ii]->GetCoord2();   // First defined trace
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

      LogKit::LogFormatted(LogKit::High,"\n                      x            y         IL       XL");
      LogKit::LogFormatted(LogKit::High,"\n--------------------------------------------------------\n");
      LogKit::LogFormatted(LogKit::High,"Origin     %12.2f %12.2f   %8d %8d\n", x0, y0, il0, xl0);
      LogKit::LogFormatted(LogKit::High,"Corner 1   %12.2f %12.2f   %8d %8d\n", x1, y1, il1, xl1);
      LogKit::LogFormatted(LogKit::High,"Corner 2   %12.2f %12.2f   %8d %8d\n", x2, y2, il2, xl2);

      double dIL   = std::sqrt(dxIL*dxIL + dyIL*dyIL);
      double dXL   = std::sqrt(dxXL*dxXL + dyXL*dyXL);
      double lenIL = (ilMax - ilMin + 1)*dIL;
      double lenXL = (xlMax - xlMin + 1)*dXL;
      double rotIL =               atan(std::abs(dyIL/dxIL));
      double rotXL = NRLib::Pi/2 - atan(std::abs(dyXL/dxXL));

      LogKit::LogFormatted(LogKit::High,"\n               dx           dy    increment       length        angle");
      LogKit::LogFormatted(LogKit::High,"\n---------------------------------------------------------------------");
      LogKit::LogFormatted(LogKit::High,"\nIL   %12.4f %12.4f %12.4f %12.4f %12.4f"  , dxIL, dyIL, dIL, lenIL, NRLib::RadToDeg(rotIL));
      LogKit::LogFormatted(LogKit::High,"\nXL   %12.4f %12.4f %12.4f %12.4f %12.4f\n", dxXL, dyXL, dXL, lenXL, NRLib::RadToDeg(rotXL));

      if (regularize_if_needed) {
        bool regularized = Regularize(dxIL, dyIL, dxXL, dyXL,
                                      dIL1, dIL2, dXL1, dXL2);

        if (regularized) {
          rotIL =               atan(std::abs(dyIL/dxIL));
          rotXL = NRLib::Pi/2 - atan(std::abs(dyXL/dxXL));

          LogKit::LogFormatted(LogKit::High,"\nAfter regularization:");
          LogKit::LogFormatted(LogKit::High,"\n               dx           dy    increment       length        angle");
          LogKit::LogFormatted(LogKit::High,"\n---------------------------------------------------------------------");
          LogKit::LogFormatted(LogKit::High,"\nIL   %12.4f %12.4f %12.4f %12.4f %12.4f"  , dxIL, dyIL, dIL, lenIL, NRLib::RadToDeg(rotIL));
          LogKit::LogFormatted(LogKit::High,"\nXL   %12.4f %12.4f %12.4f %12.4f %12.4f\n", dxXL, dyXL, dXL, lenXL, NRLib::RadToDeg(rotXL));
        }
      }
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

    // See pptx in doc folder of nrlib repository for explanation of (some of) the geometry variables
    SetupGeometry(traces[ii]->GetCoord1(), traces[ii]->GetCoord2(), traces[ii]->GetInline(), traces[ii]->GetCrossline(),
                  ilMin, ilMax, xlMin, xlMax, dxIL, dyIL,
                  dxXL, dyXL, deltaIL, deltaXL);

    // If (min,max,step)(IL,XL) is directly available in this constructor, we can use this instead.
    FindILXLGeometry(ilMin, ilMax, deltaIL, xlMin, xlMax, deltaXL);

    size_t nTraces = nx_ * ny_;
    std::vector<SegYTrace *> tracestmp;
    tracestmp.resize(nTraces);

    for (size_t k = 0; k < nTraces; k++)
      tracestmp[k] = 0;

    double distIL = dxIL*dxIL+dyIL*dyIL;
    double distXL = dxXL*dxXL+dyXL*dyXL;
    double error_dist = distIL;
    if(error_dist == 0 || (error_dist > distXL && distXL > 0))
      error_dist = distXL;
    error_dist = 0.25f*error_dist; //Compare square distances, use half distance as limit.

    int    mxd_IL   = 0;
    int    mxd_XL   = 0;
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

          double cx, cy;
          this->FindXYFromILXL(til, txl, cx, cy);
          double dtx  = x - cx;
          double dty  = y - cy;
          double dist = dtx*dtx + dty*dty;
          if(dist> max_dist) {
            mxd_IL   = til;
            mxd_XL   = txl;
            max_dist = dist;
          }

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
    //
    // Note that the trace with the largest mismatch may be in the middle
    // of the grid. The grid we setup is based on two of the SEGY corners,
    // and if there are internal irregularities, we may end up with
    // mismatches like this.
    //
    LogKit::LogFormatted(LogKit::High,"\nLargest distance between trace location and centre of assigned grid node: %.2fm",max_dist);
    LogKit::LogFormatted(LogKit::High,"\nThis distance is observed for IL=%d and XL=%d\n",mxd_XL,mxd_IL);

    traces.resize(nTraces);
    for (size_t k = 0 ; k < nTraces ; k++) {
      traces[k] = tracestmp[k];
      if (traces[k] != NULL)
        traces[k]->SetTableIndex(k);
    }
  }
  else if (traces.size() == 1) {
    SetupGeometrySingleTrace(traces[0]);
  }
  else {
    // No traces
    throw NRLib::Exception("Can not set up geometry: No traces in file.");
  }
}

bool
SegyGeometry::Regularize(double & dxIL, double & dyIL,
                         double & dxXL, double & dyXL,
                         double   dIL1, double   dIL2,
                         double   dXL1, double   dXL2)
{
  bool regularized = false;

  // Assume angle is correct for inline direction. With the same
  // angle enforced for the XL-direction, how does this affect the
  // positions of corners 1 and 2

  double a = std::sqrt(dxIL*dxIL + dyIL*dyIL);
  double b = std::sqrt(dxXL*dxXL + dyXL*dyXL);
  double R = acos((dxIL*dxXL + dyIL*dyXL)/(a*b));

  LogKit::LogFormatted(LogKit::High,"\nAngle between IL and XL axes : %.4f\n",RadToDeg(R));

  double dxILn = 0.0;
  double dyILn = 0.0;
  double dxXLn = 0.0;
  double dyXLn = 0.0;

  Orthogonalize(dxIL , dyIL,  dxXL , dyXL,  // Original values
                dxILn, dyILn, dxXLn, dyXLn, // New values
                dyIL/dxIL);                 // IL axis is not rotated in this case as original ratio is used.

  LogKit::LogFormatted(LogKit::DebugLow,"\nChanging angle of XL axis so that it becomes orthogonal to IL axis.\n");
  LogKit::LogFormatted(LogKit::DebugLow,"\n               dx           dy");
  LogKit::LogFormatted(LogKit::DebugLow,"\n------------------------------");
  LogKit::LogFormatted(LogKit::DebugLow,"\nIL   %12.4f %12.4f"  , dxILn, dyILn);
  LogKit::LogFormatted(LogKit::DebugLow,"\nXL   %12.4f %12.4f\n", dxXLn, dyXLn);


  double dx1  = (dxILn - dxIL)*dIL1 + (dxXLn - dxXL)*dXL1;  // Corner corrections
  double dy1  = (dyILn - dyIL)*dIL1 + (dyXLn - dyXL)*dXL1;
  double dx2  = (dxILn - dxIL)*dIL2 + (dxXLn - dxXL)*dXL2;
  double dy2  = (dyILn - dyIL)*dIL2 + (dyXLn - dyXL)*dXL2;

  double dxy1 = std::sqrt(dx1*dx1 + dy1*dy1);
  double dxy2 = std::sqrt(dx2*dx2 + dy2*dy2);

  LogKit::LogFormatted(LogKit::DebugLow,"\nFinding maximum mismatch between corner trace locations in SEGY volume and in grid.\n");
  LogKit::LogFormatted(LogKit::DebugLow,"\n              Delta_x    Delta_y   Delta_xy");
  LogKit::LogFormatted(LogKit::DebugLow,"\n-------------------------------------------");
  LogKit::LogFormatted(LogKit::DebugLow,"\nCorner 1   %10.2f %10.2f %10.2f"  , dx1, dy1, dxy1);
  LogKit::LogFormatted(LogKit::DebugLow,"\nCorner 2   %10.2f %10.2f %10.2f\n", dx2, dy2, dxy2);

  double max_dx = std::max(std::abs(dx1), std::abs(dx2));
  double max_dy = std::max(std::abs(dy1), std::abs(dy2));

  // Allow a mismatch of 1/3 of the distance between traces.
  double tol    = 1.0/3.0;
  double tol_x  = tol*std::min(std::abs(dxIL), std::abs(dxXL));
  double tol_y  = tol*std::min(std::abs(dyIL), std::abs(dyXL));

  LogKit::LogFormatted(LogKit::DebugLow,"\nAllowing a mismatch of %.2f of IL and XL increments.\n",tol);
  LogKit::LogFormatted(LogKit::DebugLow,"\n                           x          y");
  LogKit::LogFormatted(LogKit::DebugLow,"\n---------------------------------------");
  LogKit::LogFormatted(LogKit::DebugLow,"\nAllowed mismatch     %7.2f %10.2f\n", tol_x, tol_y);

  if (max_dx > tol_x || max_dy > tol_y) {
    regularized = true;

    LogKit::LogFormatted(LogKit::High,"\nMismatch too large. Regularization activated.\n");

    double fracIL = std::abs(dyIL/dxIL);
    double fracXL = std::abs(dxXL/dyXL);
    double minf   = std::min(fracIL, fracXL);
    double maxf   = std::max(fracIL, fracXL);
    double step   = (maxf - minf)/10;

    double dmin   = std::numeric_limits<double>::infinity();
    double fmin   = 0.0;

    LogKit::LogFormatted(LogKit::DebugLow,"\n                Angle    Delta_x    Delta_y   Delta_xy");
    LogKit::LogFormatted(LogKit::DebugLow,"\n------------------------------------------------------");

    for (int i = 0 ; i < 11 ; i++) {
      double frac = minf + step*static_cast<double>(i);

      Orthogonalize(dxIL , dyIL,  dxXL , dyXL,  // Original values
                    dxILn, dyILn, dxXLn, dyXLn, // New values
                    frac);                      //

      dx1  = (dxILn - dxIL)*dIL1 + (dxXLn - dxXL)*dXL1;  // Corner corrections
      dy1  = (dyILn - dyIL)*dIL1 + (dyXLn - dyXL)*dXL1;
      dx2  = (dxILn - dxIL)*dIL2 + (dxXLn - dxXL)*dXL2;
      dy2  = (dyILn - dyIL)*dIL2 + (dyXLn - dyXL)*dXL2;

      dxy1 = std::sqrt(dx1*dx1 + dy1*dy1);
      dxy2 = std::sqrt(dx2*dx2 + dy2*dy2);

      if (std::max(dxy1, dxy2) < dmin) {
        dmin = std::max(dxy1, dxy2);
        fmin = frac;
      }

      double angle = NRLib::RadToDeg(atan(frac));

      LogKit::LogFormatted(LogKit::DebugLow,"\nCorner 1   %10.4f %10.2f %10.2f %10.2f", angle, dx1, dy1, dxy1);
      LogKit::LogFormatted(LogKit::DebugLow,"\nCorner 2   %10.4f %10.2f %10.2f %10.2f", angle, dx2, dy2, dxy2);
    }

    LogKit::LogFormatted(LogKit::High,"\n\nOptimal angle     :  %8.4f", NRLib::RadToDeg(atan(fmin)));
    LogKit::LogFormatted(LogKit::High,"\nSmallest mismatch :  %8.2f\n", dmin);

    // Use fraction that minimizes error to find new dxIL, dyIL, ....
    Orthogonalize(dxIL , dyIL,  dxXL , dyXL,
                  dxILn, dyILn, dxXLn, dyXLn,
                  fmin);

    dxIL  = dxILn;
    dyIL  = dyILn;
    dxXL  = dxXLn;
    dyXL  = dyXLn;

    a = std::sqrt(dxIL*dxIL + dyIL*dyIL);
    b = std::sqrt(dxXL*dxXL + dyXL*dyXL);
    R = acos((dxIL*dxXL + dyIL*dyXL)/(a*b));

    LogKit::LogFormatted(LogKit::High,"\nAngle between IL and XL axes : %.4f\n",RadToDeg(R));
  }
  else {
    LogKit::LogFormatted(LogKit::High,"\nRegularization not needed.\n");
  }
  return regularized;
}

void
SegyGeometry::Orthogonalize(double   dxIL , double   dyIL,
                            double   dxXL , double   dyXL,
                            double & dxILn, double & dyILn,
                            double & dxXLn, double & dyXLn,
                            double   frac)
{
  double dIL, dXL, dILn, dXLn;
  // Rotate IL axis by changing dyIL/dxIL. Then rescale to original length.
  dxILn  = dxIL;
  dyILn  = dxIL*frac;
  dIL    = std::sqrt(dxIL*dxIL   + dyIL*dyIL);
  dILn   = std::sqrt(dxILn*dxILn + dyILn*dyILn);
  dxILn  = dxILn*(dIL/dILn);
  dyILn  = dyILn*(dIL/dILn);

  // Rotate XL axis by changing dyIL/dxIL. Then rescale to original length.
  dxXLn  =  dxXL;
  dyXLn  = -dxXL*(1.0/frac);
  dXL    = std::sqrt(dxXL*dxXL   + dyXL*dyXL);
  dXLn   = std::sqrt(dxXLn*dxXLn + dyXLn*dyXLn);
  dxXLn  = dxXLn*(dXL/dXLn);
  dyXLn  = dyXLn*(dXL/dXLn);
}

void
SegyGeometry::SetupGeometry(double xRef, double yRef, int ilRef, int xlRef,
                            int minIL, int maxIL, int minXL, int maxXL,
                            double dxIL, double dyIL, double dxXL, double dyXL,
                            int deltaIL, int deltaXL)
{
  if (minXL == maxXL) {
    // Line going in XL direction. (Only IL coordinates)
    dxXL = dyIL;
    dyXL = -dxIL;
    deltaXL = deltaIL;
  }
  else if (minIL == maxIL) {
    // Line going in IL direction. (Only XL coordinates)
    dxIL = -dyXL;
    dyIL = dxXL;
    deltaIL = deltaXL;
  }

  double lx0 = xRef - xlRef*dxXL - ilRef*dxIL;
  double ly0 = yRef - xlRef*dyXL - ilRef*dyIL;

  std::vector<double> cornerx(4);
  std::vector<double> cornery(4);
  std::vector<double> cornerxl(4);
  std::vector<double> corneril(4);

  cornerx[0] = lx0 + (minXL - 0.5*deltaXL)*dxXL + (minIL - 0.5*deltaIL)*dxIL;
  cornery[0] = ly0 + (minXL - 0.5*deltaXL)*dyXL + (minIL - 0.5*deltaIL)*dyIL;
  cornerx[1] = lx0 + (minXL - 0.5*deltaXL)*dxXL + (maxIL + 0.5*deltaIL)*dxIL;
  cornery[1] = ly0 + (minXL - 0.5*deltaXL)*dyXL + (maxIL + 0.5*deltaIL)*dyIL;
  cornerx[2] = lx0 + (maxXL + 0.5*deltaXL)*dxXL + (minIL - 0.5*deltaIL)*dxIL;
  cornery[2] = ly0 + (maxXL + 0.5*deltaXL)*dyXL + (minIL - 0.5*deltaIL)*dyIL;
  cornerx[3] = lx0 + (maxXL + 0.5*deltaXL)*dxXL + (maxIL + 0.5*deltaIL)*dxIL;
  cornery[3] = ly0 + (maxXL + 0.5*deltaXL)*dyXL + (maxIL + 0.5*deltaIL)*dyIL;

  cornerxl[0] = minXL - 0.5*deltaXL;
  corneril[0] = minIL - 0.5*deltaIL;
  cornerxl[1] = minXL - 0.5*deltaXL;
  corneril[1] = maxIL + 0.5*deltaIL;
  cornerxl[2] = maxXL + 0.5*deltaXL;
  corneril[2] = minIL - 0.5*deltaIL;
  cornerxl[3] = maxXL + 0.5*deltaXL;
  corneril[3] = maxIL + 0.5*deltaIL;
  double miny, maxx, minx, maxy;
  miny = cornery[0];
  maxx = cornerx[0];
  minx = cornerx[0];
  maxy = cornery[0];
  size_t index_min_y = 0;
  size_t index_max_x = 0;
  size_t index_min_x = 0;
  size_t index_max_y = 0;
  size_t i;
  for (i = 1; i < 4; i++)
  {
    if (cornery[i] < miny)
    {
      miny = cornery[i];
      index_min_y = i;
    }
    if (cornerx[i] > maxx)
    {
      maxx = cornerx[i];
      index_max_x = i;
    }
    if (cornerx[i] < minx)
    {
      minx = cornerx[i];
      index_min_x = i;
    }
    if (cornery[i] > maxy)
    {
      maxy = cornery[i];
      index_max_y = i;
    }
  }

  //If the grid is oriented along UTMX/UTMY, we may have a problem now,
  //since we have non-unique max and min. Check and fix.
  if (index_max_x == index_min_y || index_min_x == index_min_y || index_max_x == index_max_y || index_min_x == index_max_y) {
    //Square oriented with UTMX/UTMY.
    //index_min_y  = lower left
    //index_max_x = lower right
    //index_min_x = upper left
    //index_max_y = upper right
    double dXL = sqrt(dxXL*dxXL + dyXL*dyXL);
    double dl = sqrt(dxIL*dxIL + dyIL*dyIL);
    if (dl == 0 || (dXL > 0 && dXL < dl))
      dl = dXL;

    //Utlilize fact that sum of index in diagonally opposite corners is 3.
    if (fabs(cornerx[0] - cornerx[1]) < 0.2*dl) {
      //0 and 1 are along same constant x.
      if (cornerx[2] < cornerx[0]) { //2 and 3 are leftmost edge.
        if (cornery[2] < cornery[3])
          index_min_y = 2;
        else
          index_min_y = 3;
        index_max_x = index_min_y - 2;
      }
      else {
        if (cornery[0] < cornery[1])
          index_min_y = 0;
        else
          index_min_y = 1;
        index_max_x = index_min_y + 2;
      }
    }
    else {
      //0 and 1 are along same constant y.
      if (cornery[2] < cornery[0]) { //2 and 3 are bottom edge.
        if (cornerx[2] < cornerx[3])
          index_min_y = 2;
        else
          index_min_y = 3;
        index_max_x = 5 - index_min_y;
      }
      else {
        if (cornerx[0] < cornerx[1])
          index_min_y = 0;
        else
          index_min_y = 1;
        index_max_x = 1 - index_min_y;
      }
    }
    index_min_x = 3 - index_max_x;
    index_max_y = 3 - index_min_y;
  }

  LogKit::LogMessage(LogKit::DebugLow, "Corner 0:  (" + NRLib::ToString(cornerx[index_min_y]) + " , " + NRLib::ToString(cornery[index_min_y]) + ").\n");
  LogKit::LogMessage(LogKit::DebugLow, "Corner 1:  (" + NRLib::ToString(cornerx[index_max_x]) + " , " + NRLib::ToString(cornery[index_max_x]) + ").\n");
  LogKit::LogMessage(LogKit::DebugLow, "Corner 2:  (" + NRLib::ToString(cornerx[index_min_x]) + " , " + NRLib::ToString(cornery[index_min_x]) + ").\n");
  LogKit::LogMessage(LogKit::DebugLow, "Corner 3:  (" + NRLib::ToString(cornerx[index_max_y]) + " , " + NRLib::ToString(cornery[index_max_y]) + ").\n");

  double deltax = cornerx[index_min_y] - cornerx[index_max_x];
  double deltay = cornery[index_min_y] - cornery[index_max_x];
  if (deltax == 0.0) {
    assert(deltay != 0.0);
    rot_ = 0.0;
  }
  else
    rot_ = atan(deltay / deltax);

  double lx, ly;
  if (std::abs(rot_) <= 0.25*NRLib::Pi)
  {
    x0_ = cornerx[index_min_y];
    y0_ = cornery[index_min_y];
    lx = sqrt(deltay*deltay + deltax*deltax);
    nx_ = static_cast<size_t>(std::abs(cornerxl[index_min_y] - cornerxl[index_max_x]) / deltaXL);
    if (nx_ == 0)
      nx_ = static_cast<size_t>(std::abs(corneril[index_min_y] - corneril[index_max_x]) / deltaIL);

    deltax = cornerx[index_min_y] - cornerx[index_min_x];
    deltay = cornery[index_min_y] - cornery[index_min_x];
    ly = sqrt(deltay*deltay + deltax*deltax);
    ny_ = static_cast<size_t>(std::abs(cornerxl[index_min_y] - cornerxl[index_min_x]) / deltaXL);
    if (ny_ == 0)
      ny_ = static_cast<size_t>(std::abs(corneril[index_min_y] - corneril[index_min_x]) / deltaIL);
  }
  else
  {
    x0_ = cornerx[index_min_x];
    y0_ = cornery[index_min_x];
    deltax = cornerx[index_min_y] - cornerx[index_min_x];
    deltay = cornery[index_min_y] - cornery[index_min_x];

    if (deltax > 0)
      rot_ = atan(deltay / deltax);
    else {
      rot_ -= 0.5*NRLib::Pi;
    }

    lx = sqrt(deltay*deltay + deltax*deltax);
    nx_ = static_cast<size_t>(std::abs(cornerxl[index_min_y] - cornerxl[index_min_x]) / deltaXL);
    if (nx_ == 0)
      nx_ = static_cast<size_t>(std::abs(corneril[index_min_y] - corneril[index_min_x]) / deltaIL);

    deltax = cornerx[index_max_y] - cornerx[index_min_x];
    deltay = cornery[index_max_y] - cornery[index_min_x];
    ly = sqrt(deltay*deltay + deltax*deltax);
    ny_ = static_cast<size_t>(std::abs(cornerxl[index_max_y] - cornerxl[index_min_x]) / deltaXL);
    if (ny_ == 0)
      ny_ = static_cast<size_t>(std::abs(corneril[index_max_y] - corneril[index_min_x]) / deltaIL);
  }

  assert(lx != 0.0 && ly != 0.0);

  dx_ = static_cast<double>(lx / nx_);
  dy_ = static_cast<double>(ly / ny_);
  sin_rot_ = sin(rot_);
  cos_rot_ = cos(rot_);

  if (dyIL != 0) {
    xl_stepX_ = 1 / (dxXL - dxIL*dyXL / dyIL);
    il_stepX_ = -xl_stepX_*dyXL / dyIL;
  }
  else {  //implies that dxXL = 0, dxIL > 0)
    xl_stepX_ = 0;
    il_stepX_ = 1 / dxIL;
  }
  if (dxIL != 0) {
    xl_stepY_ = 1 / (dyXL - dyIL*dxXL / dxIL);
    il_stepY_ = -xl_stepY_*dxXL / dxIL;
  }
  else {
    xl_stepY_ = 0;
    il_stepY_ = 1 / dyIL;
  }
  in_line0_ = ilRef - il_stepX_*(xRef - x0_) - il_stepY_*(yRef - y0_);
  cross_line0_ = xlRef - xl_stepX_*(xRef - x0_) - xl_stepY_*(yRef - y0_);
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

SegyGeometry::SegyGeometry()
{}

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
  IL0_         = 0;
  XL0_         = 0;
  minIL_       = 0;
  maxIL_       = static_cast<int>(nx - 1);
  ILStep_      = 1;
  minXL_       = 0;
  maxXL_       = static_cast<int>(ny - 1);
  XLStep_      = 1;
  firstAxisIL_ = true;
}

SegyGeometry::SegyGeometry(double x0, double y0, double dx, double dy, size_t nx, size_t ny,
                           double ILStart, double XLStart, double ilStepX, double ilStepY,
                           double xlStepX, double xlStepY, double rot)
{
  x0_          = x0;
  y0_          = y0;
  dx_          = dx;
  dy_          = dy;
  nx_          = nx;
  ny_          = ny;
  in_line0_    = ILStart;
  cross_line0_ = XLStart;
  il_stepX_    = ilStepX;
  il_stepY_    = ilStepY;
  xl_stepX_    = xlStepX;
  xl_stepY_    = xlStepY;
  cos_rot_     = cos(rot);
  sin_rot_     = sin(rot);
  rot_         = rot;

  FindILXLGeometry();
}

SegyGeometry::SegyGeometry(const RegularSurface<double> & surf)
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
  IL0_         = 0;
  XL0_         = 0;
  minIL_       = 0;
  maxIL_       = static_cast<int>(nx_ - 1);
  ILStep_      = 1;
  minXL_       = 0;
  maxXL_       = static_cast<int>(ny_ - 1);
  XLStep_      = 1;
  firstAxisIL_ = true;
}

SegyGeometry::SegyGeometry(const RegularSurfaceRotated<double> & surf)
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
  IL0_         = 0;
  XL0_         = 0;
  minIL_       = 0;
  maxIL_       = static_cast<int>(nx_ - 1);
  ILStep_      = 1;
  minXL_       = 0;
  maxXL_       = static_cast<int>(ny_ - 1);
  XLStep_      = 1;
  firstAxisIL_ = true;
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
  double xind, yind;
  return(FindContIndex(x,y,xind,yind));
}

size_t
SegyGeometry::FindIndex(double x, double y) const
{
  size_t i, j;
  FindIndex(x, y, i, j);
  return (j*nx_ + i);
}

void
SegyGeometry::FindIndex(double x, double y, size_t &i, size_t &j) const
{
  double ci, cj;
  bool is_inside = FindContIndex(x, y, ci, cj);
  i = static_cast<int>(ci);
  j = static_cast<int>(cj);

  if (!is_inside)
    throw Exception("Coordinates (UTM-X, UTM-Y) = (" + NRLib::ToString(x) + ", " +
                    NRLib::ToString(y) + ") is outside the SEGY cube.");
}

void
SegyGeometry::FindIndex(int IL, int XL, size_t &i, size_t &j) const
{
  int ii, jj;

  if (IL0_ == minIL_) {
    ii = (IL - IL0_) / ILStep_;
  }
  else {
    ii = (IL0_ - IL) / ILStep_;
  }

  if (XL0_ == minXL_) {
    jj = (XL - XL0_) / XLStep_;
  }
  else {
    jj = (XL0_ - XL) / XLStep_;
  }

  if (!firstAxisIL_) {
    std::swap(ii, jj);
  }

  if (ii < 0 || jj < 0 || ii > static_cast<int>(nx_) || jj > static_cast<int>(ny_))
    throw Exception("Indexes (i, j) = (" + NRLib::ToString(ii) + ", " + NRLib::ToString(jj) + ") outside grid "
                    "found for (IL, XL) = (" + NRLib::ToString(IL) + ", " + NRLib::ToString(XL) + ")");

  i = static_cast<size_t>(ii);
  j = static_cast<size_t>(jj);
}

bool
SegyGeometry::FindContIndex(double x, double y, double &xind, double &yind) const
{
  double IL, XL;
  FindContILXL(x, y, IL, XL);
  FindContIndexFromContILXL(IL, XL, xind, yind);

  if (xind>=0.0 && xind<static_cast<double>(nx_) && yind>=0.0 && yind<static_cast<double>(ny_))
    return true;
  else
    return false;
}


void
SegyGeometry::FindContIndexFromContILXL(double il, double xl, double &xind, double &yind) const
{
  double xStrides, yStrides;
  int xStep, yStep;

  if (firstAxisIL_)
  {
    xStep = ILStep_;
    yStep = XLStep_;
    if (IL0_ == minIL_)
      xStrides = il - minIL_;
    else
      xStrides = maxIL_ - il;

    if (XL0_ == minXL_)
      yStrides = xl - minXL_;
    else
      yStrides = maxXL_ - xl;
  }
  else
  {
    xStep = XLStep_;
    yStep = ILStep_;
    if (IL0_ == minIL_)
      yStrides = il - minIL_;
    else
      yStrides = maxIL_ - il;

    if (XL0_ == minXL_)
      xStrides = xl - minXL_;
    else
      xStrides = maxXL_ - xl;
  }
  xind = 0.5 + xStrides / xStep;
  yind = 0.5 + yStrides / yStep;
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
  int ii = static_cast<int>(i);
  int jj = static_cast<int>(j);
  if (!firstAxisIL_)
    std::swap(ii, jj);

  if (IL0_ == minIL_) {
    IL = IL0_  + ii * ILStep_;
  }
  else {
    IL = IL0_ - ii * ILStep_;
  }

  if (XL0_ == minXL_) {
    XL = XL0_ + jj * XLStep_;
  }
  else {
    XL = XL0_ - jj * XLStep_;
  }
}

void
SegyGeometry::FindContILXL(double x, double y, double &IL, double &XL) const
{
  IL = in_line0_    + (x-x0_)*il_stepX_ + (y-y0_)*il_stepY_;
  XL = cross_line0_ + (x-x0_)*xl_stepX_ + (y-y0_)*xl_stepY_;
}

void
SegyGeometry::FindXYFromIJ(size_t i, size_t j, double &x, double &y) const
{
  int IL, XL;
  FindILXL(i, j, IL, XL);
  FindXYFromILXL(IL, XL, x, y);
}

void
SegyGeometry::FindXYFromILXL(int IL, int XL, double &x, double &y) const
{
  FindXYFromContILXL(static_cast<double>(IL), static_cast<double>(XL), x, y);
}

void
SegyGeometry::FindXYFromContILXL(double IL, double XL, double &x, double &y) const
{
  double xd, yd;
  if (xl_stepX_ == 0.0 && xl_stepY_ == 0.0) {//inline only
    if (sin_rot_ > 0.0) {
      xd = x0_+(IL-in_line0_)/il_stepX_;
      yd = y0_+(xd-x0_)*sin_rot_/cos_rot_; //If cos_rot_ is 0, il_stepX_ should be 0 here.
    }
    else {
      yd = y0_+(IL-in_line0_)/il_stepY_;
      xd = x0_-(yd-y0_)*sin_rot_/cos_rot_;
    }
  }
  else if (il_stepX_ == 0.0 && il_stepY_ == 0.0) {// crossline only
    if (sin_rot_ > 0.0) {
      xd = x0_+(XL-cross_line0_)/xl_stepX_;
      yd = y0_+(xd-x0_)*sin_rot_/cos_rot_; //If cos_rot_ is 0, xl_stepX_ should be 0 here.
    }
    else {
      yd = y0_+(XL-cross_line0_)/xl_stepY_;
      xd = x0_-(yd-y0_)*sin_rot_/cos_rot_;
    }
  }
  else if (il_stepX_ == 0.0 || xl_stepY_ == 0.0) {
    xd = x0_+(XL-cross_line0_)/xl_stepX_;
    yd = y0_+(IL-in_line0_)/il_stepY_;
  }
  else {
    yd = y0_+(XL-cross_line0_-(IL-in_line0_)*xl_stepX_/il_stepX_)/(xl_stepY_-il_stepY_*xl_stepX_/il_stepX_);
    xd = x0_+(IL-in_line0_-il_stepY_*(yd-y0_))/il_stepX_;
  }
  x = xd;
  y = yd;
}



void SegyGeometry::WriteGeometry() const
{
  double geoangle = -rot_*180/(NRLib::Pi);
  if (geoangle < 0.0)
    geoangle += 360.0;

  LogKit::LogFormatted(LogKit::Low,"\n             ReferencePoint      Length Increment\n");
  LogKit::LogFormatted(LogKit::Low,"-------------------------------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"X-coordinate   %12.2f %11.2f %9.5f\n", x0_, nx_*dx_, dx_);
  LogKit::LogFormatted(LogKit::Low,"Y-coordinate   %12.2f %11.2f %9.5f\n", y0_, ny_*dy_, dy_);
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
SegyGeometry::FindILXLGeometry(int minIL, int maxIL, int ILStep, int minXL, int maxXL, int XLStep)
// Assign values to geometry variables and identify if IL is along the firstAxis or not by moving
// a small step from the reference position (IL0_,XL0_)
// See pptx in doc folder of nrlib repository
{
  double x = x0_ + 0.5*dx_*cos_rot_ - 0.5*dy_*sin_rot_;
  double y = y0_ + 0.5*dx_*sin_rot_ + 0.5*dy_*cos_rot_;
  FindILXL(x, y, IL0_, XL0_);

  minIL_  = minIL;
  maxIL_  = maxIL;
  ILStep_ = ILStep;
  minXL_  = minXL;
  maxXL_  = maxXL;
  XLStep_ = XLStep;

  if (ILStep_ == 0)
    ILStep_ = 1;
  if (XLStep_ == 0)
    XLStep_ = 1;

  // This could maybe be set in SetupGeometry?
  // Find position of (x,y) for (i, j) = (1, 0)
  int IL1, XL1;

  bool second_axis = false;
  if (nx_ > 1) {
    x = x0_ + 1.5*dx_*cos_rot_ - 0.5*dy_*sin_rot_;
    y = y0_ + 1.5*dx_*sin_rot_ + 0.5*dy_*cos_rot_;
    FindILXL(x, y, IL1, XL1); // Moved along i (x), so change in XL now means that firstAxisIL_ should be false.
  }
  else {
    x = x0_ + 0.5*dx_*cos_rot_ - 1.5*dy_*sin_rot_;
    y = y0_ + 0.5*dx_*sin_rot_ + 1.5*dy_*cos_rot_;
    FindILXL(x, y, IL1, XL1);
    second_axis = true; //Moved along j (y), so change in XL now means that firstAxisIL_ should be true.
  }

  if (XL1 - XL0_ != 0) {
    firstAxisIL_ = second_axis;
  }
  else {
    firstAxisIL_ = !second_axis;
  }
}


void
SegyGeometry::FindILXLGeometry()
{
  int nx = static_cast<int>(nx_);
  int ny = static_cast<int>(ny_);

  double x = x0_ + 0.5*dx_*cos_rot_ - 0.5*dy_*sin_rot_;
  double y = y0_ + 0.5*dx_*sin_rot_ + 0.5*dy_*cos_rot_;
  FindILXL(x, y, IL0_, XL0_);

  int IL1, XL1;
  x = x0_ + 1.5*dx_*cos_rot_ - 0.5*dy_*sin_rot_;
  y = y0_ + 1.5*dx_*sin_rot_ + 0.5*dy_*cos_rot_;
  FindILXL(x, y, IL1, XL1);

  int IL2, XL2;
  x = x0_ + 0.5*dx_*cos_rot_ - 1.5*dy_*sin_rot_;
  y = y0_ + 0.5*dx_*sin_rot_ + 1.5*dy_*cos_rot_;
  FindILXL(x, y, IL2, XL2);

  if (IL1 != IL0_) {
    assert(XL1 == XL0_);
    assert(XL2 != XL0_);

    firstAxisIL_ = true;

    ILStep_ = std::abs(IL1 - IL0_);
    XLStep_ = std::abs(XL2 - XL0_);

    if (IL1 > IL0_) {
      minIL_ = IL0_;
      maxIL_ = IL0_ + (nx - 1) * ILStep_;
    }
    else {
      maxIL_ = IL0_;
      minIL_ = IL0_ - (nx - 1) * ILStep_;
    }
    if (XL2 > XL0_) {
      minXL_ = XL0_;
      maxXL_ = XL0_ + (ny - 1) * XLStep_;
    }
    else {
      maxXL_ = XL0_;
      minXL_ = XL0_ - (ny - 1) * XLStep_;
    }
  }
  else {
    assert(XL1 != XL0_);
    assert(IL2 != IL0_);

    firstAxisIL_ = false;

    ILStep_ = std::abs(IL2 - IL0_);
    XLStep_ = std::abs(XL1 - XL0_);

    if (IL2 > IL0_) {
      minIL_ = IL0_;
      maxIL_ = IL0_ + (ny - 1) * ILStep_;
    }
    else {
      maxIL_ = IL0_;
      minIL_ = IL0_ - (ny - 1) * ILStep_;
    }
    if (XL1 > XL0_) {
      minXL_ = XL0_;
      maxXL_ = XL0_ + (nx - 1) * XLStep_;
    }
    else {
      maxXL_ = XL0_;
      minXL_ = XL0_ - (nx - 1) * XLStep_;
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
      LogKit::LogFormatted(LogKit::Warning, "\nWARNING: Requested IL-step (%d) is smaller than IL-step in seismic data (%d). Data vill be interpolated.\n",
                           subILStep, ILStep);
      interpolation = true;
    }
    else if((subILStep % ILStep) > 0) {
      LogKit::LogFormatted(LogKit::Warning, "\nWARNING: Requested IL-step (%d) is not a multiple of IL-step in seismic data (%d). Data vill be interpolated.\n",
                           subILStep, ILStep);
      interpolation = true;
    }
  }
  int subXLStep = ilxl[5];
  if(subXLStep < 0)
    subXLStep = XLStep;
  else {
    if(subXLStep < XLStep) {
      LogKit::LogFormatted(LogKit::Warning, "\nWARNING: Requested XL-step (%d) is smaller than XL-step in seismic data (%d). Data vill be interpolated.\n",
                           subXLStep, XLStep);
      interpolation = true;
    }
    else if((subXLStep % XLStep) > 0) {
      LogKit::LogFormatted(LogKit::Warning, "\nWARNING: Requested XL-step (%d) is not a multiple of XL-step in seismic data (%d). Data vill be interpolated.\n",
                           subXLStep, XLStep);
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

      LogKit::LogFormatted(LogKit::Warning, "\nWARNING: Initial minimum IL (%d) does not match a sampled trace in the input data. Moved to %d.\n", oldIL, minSubIL);
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
      LogKit::LogFormatted(LogKit::Warning, "\nWARNING: Initial minimum XL (%d) does not match a sampled trace in the input data. Moved to %d.\n", oldXL, minSubXL);
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
    LogKit::LogFormatted(LogKit::Warning, "\nWARNING: Initial maximum IL (%d) is not an integer number of steps away from minimum IL (%d)."
                         " Moved to %d.\n", oldIL, minSubIL, maxSubIL);
  }
  residual = ((maxSubXL-minSubXL) % subXLStep);
  if(residual != 0) { //Move max to fit.
    modified = true;
    aligned  = true;
    int oldXL = maxSubXL;
    maxSubXL += subXLStep-residual; //Try increasing
    if(maxSubXL > maxXL)
      maxSubXL -= subXLStep;        //Have to decrease.
    LogKit::LogFormatted(LogKit::Warning, "\nWARNING: Initial maximum IL (%d) is not an integer number of steps away from minimum IL (%d)."
                         " Moved to %d.\n", oldXL, minSubXL, maxSubXL);
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

  double subIL0, subXL0;
  //Determine reference corner, move from center and out to corner.
  if(IL0 == minIL) { //Reference corner is minimum inline.
    subIL0 = minSubIL-0.5*subILStep;
  }
  else {             //Reference corner is maximum inline.
    subIL0 = maxSubIL+0.5*subILStep;
  }
  if(XL0 == minXL) { //Reference corner is minimum crossline.
    subXL0 = minSubXL-0.5*subXLStep;
  }
  else {             //Reference corner is maximum crossline.
    subXL0 = maxSubXL+0.5*subXLStep;
  }

  double subX0, subY0;
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
  FindContILXL(x0, y0, il0, xl0);
  FindContILXL(x1, y1, il1, xl1);
  FindContILXL(x2, y2, il2, xl2);
  FindContILXL(x3, y3, il3, xl3);

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
