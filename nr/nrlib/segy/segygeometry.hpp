// $Id: segygeometry.hpp 1678 2017-08-30 12:06:06Z perroe $

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

#ifndef SEGYGEOMETRY_HPP
#define SEGYGEOMETRY_HPP

#include <fstream>
#include <string>
#include <vector>

#include "traceheader.hpp"
#include "commonheaders.hpp"
#include "../volume/volume.hpp"
#include "../surface/regularsurface.hpp"
#include "../surface/regularsurfacerotated.hpp"

namespace NRLib {


class SegYTrace;


class SegyGeometry{
public:
  /**
  Constructor
  @param[in] traces  Traces from segY cube
  */
  explicit SegyGeometry(std::vector<SegYTrace *> & traces,
                        bool                       regularize_if_needed = false);
  SegyGeometry(double x0,double y0,double dx,double dy, size_t nx,size_t ny, ///< When XL, IL information not available, generates internally.
               double rot);
  SegyGeometry(double x0,double y0,double dx,double dy, size_t nx,size_t ny, ///< When XL, IL is available.
               double ILStart, double XLStart, double ilStepX, double ilStepY,
               double xlStepX, double xlStepY, double rot);
  explicit SegyGeometry(const RegularSurface<double> &surf);
  explicit SegyGeometry(const RegularSurfaceRotated<double> &surf);
  explicit SegyGeometry(const SegyGeometry *geometry);   ///< Copy constructor
  SegyGeometry(); //Default constructor for use with SetupGeometry.

  ~SegyGeometry();

  bool   IsInside(double x, double y) const;

  size_t FindIndex(double x, double y) const;
  void   FindIndex(double x, double y, size_t &i, size_t  &j) const;               ///< Return grid index for i and j
  void   FindIndex(int IL, int XL, size_t &i, size_t &j) const;
  bool   FindContIndex(double x, double y, double &xind, double &yind) const; ///< Compute decimal index, to give location in cell. Return true if inside grid.
  void   FindILXL(double x, double y, int &IL, int &XL) const;
  void   FindILXL(size_t i, size_t j, int &IL, int &XL) const;
  void   FindContILXL(double x, double y, double &IL, double &XL) const;
  void   FindXYFromIJ(size_t i, size_t j, double &x, double &y) const;
  void   FindXYFromILXL(int IL, int XL, double &x, double &y) const;
  void   FindXYFromContILXL(double IL, double XL, double &x, double &y) const;
  void   FindContIndexFromContILXL(double il, double xl, double &xind, double &yind) const;

  size_t GetNx()           const { return nx_        ;}            ///< return nx
  size_t GetNy()           const { return ny_        ;}            ///< return ny
  // Fyll på med get-funksjonar. etter behov.
  double GetDx()           const { return dx_        ;}
  double GetDy()           const { return dy_        ;}
  double GetX0()           const { return x0_        ;}
  double GetY0()           const { return y0_        ;}
  double Getlx()           const { return nx_*dx_    ;}
  double Getly()           const { return ny_*dy_    ;}
  double GetAngle()        const { return rot_       ;}
  double GetCosRot()       const { return cos_rot_    ;}
  double GetSinRot()       const { return sin_rot_    ;}
  double GetInLine0()      const { return in_line0_   ;}
  double GetCrossLine0()   const { return cross_line0_;}
  double GetILStepX()      const { return il_stepX_   ;}
  double GetILStepY()      const { return il_stepY_   ;}
  double GetXLStepX()      const { return xl_stepX_   ;}
  double GetXLStepY()      const { return xl_stepY_   ;}

  int    GetMinIL()        const { return minIL_      ;}
  int    GetMaxIL()        const { return maxIL_      ;}
  int    GetILStep()       const { return ILStep_     ;}
  int    GetMinXL()        const { return minXL_      ;}
  int    GetMaxXL()        const { return maxXL_      ;}
  int    GetXLStep()       const { return XLStep_     ;}
  int    GetIL0()          const { return IL0_        ;}
  int    GetXL0()          const { return XL0_        ;}

  bool   GetFirstAxisIL()  const { return firstAxisIL_ ;}

  void   WriteGeometry() const;
  void   WriteILXL(bool errorMode = false) const;

  void   SetupGeometry(double xRef, double yRef, int ilRef, int xlRef,
                       int minIL, int maxIL, int minXL, int maxXL,
                       double dxIL, double dyIL, double dxXL, double dyXL,
                       int deltaXL, int deltaIL);

  SegyGeometry * GetILXLSubGeometry(const std::vector<int> & ilxl, bool & interpolation, bool & aligned);

  std::vector<int> findAreaILXL(SegyGeometry * templateGeometry);

  void   FindILXLGeometry();
  void   FindILXLGeometry(int minIL, int maxIL, int ILStep, int minXL, int maxXL, int XLStep);

private:
  bool   Regularize(double & dxIL, double & dyIL,
                    double & dxXL, double & dyXL,
                    double   dIL1, double   dIL2,
                    double   dXL1, double   dXL2);
  void   Orthogonalize(double   dxIL , double   dyIL,
                       double   dxXL , double   dyXL,
                       double & dxILn, double & dyILn,
                       double & dxIXn, double & dyXLn,
                       double   frac);

  void   SetupGeometrySingleTrace(const SegYTrace * trace);

  double x0_, y0_;                  ///<
  double dx_, dy_;                  ///< Cell increments
  size_t nx_, ny_;                  ///< Grid dimensions
  double in_line0_, cross_line0_;   ///< Start value for inline and crossline
  double il_stepX_, il_stepY_;      ///< Change in IL when moving one unit along x/y (global)
  double xl_stepX_, xl_stepY_;      ///< Change in XL when moving one unit along x/y (global)
  double sin_rot_, cos_rot_;        ///<
  double rot_;                      ///<

  int IL0_, XL0_;                   ///<
  int minIL_, maxIL_, ILStep_;      ///< inline start,stop and step
  int minXL_, maxXL_, XLStep_;      ///< crossline start, stop and step
  bool firstAxisIL_;                ///<
};

} // namespace NRLib

#endif

