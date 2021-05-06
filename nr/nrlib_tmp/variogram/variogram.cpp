  // $Id: variogram.cpp 1360 2016-09-07 09:44:25Z vegard $

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

#include <string>
#include <stdio.h>
#include <assert.h>
#include "../exception/exception.hpp"
#include "variogram.hpp"
#include "variogramtypes.hpp"


using namespace NRLib;

Variogram::Variogram(double range_x,
                     double range_y,
                     double range_z,
                     double azimuth_angle,
                     double dip_angle,
                     double std_dev)
  : range_x_(range_x),
    range_y_(range_y),
    range_z_(range_z),
    azimuth_angle_(azimuth_angle),
    dip_angle_(dip_angle)
{
  var_ = std_dev * std_dev;
  EstimateFactors();
}


Variogram::Variogram()
  : range_x_(1.0),
    range_y_(1.0),
    range_z_(1.0),
    azimuth_angle_(0.0),
    dip_angle_(0.0),
    var_(1.0)
{
  EstimateFactors();
}


Variogram* Variogram::Create(Type type, double power, double range_x,
                             double range_y, double range_z, double azimuth_angle,
                             double dip_angle, double std_dev)
{
  if (type == EXPONENTIAL)
    return new ExpVario(range_x, range_y, range_z, azimuth_angle, dip_angle, std_dev);

  if (type == SPHERICAL)
    return new SphVario(range_x, range_y, range_z, azimuth_angle, dip_angle, std_dev);

  if (type == GAUSSIAN)
    return new GauVario(range_x, range_y, range_z, azimuth_angle, dip_angle, std_dev);

  if (type == GENERAL_EXPONENTIAL)
    return new GenExpVario(power, range_x, range_y, range_z, azimuth_angle, dip_angle, std_dev);

  if (type == CONSTANT)
    return new ConstVario(range_x, range_y, range_z, azimuth_angle, dip_angle, std_dev);

  if (type == MATERN32)
    return new Matern32Vario(range_x, range_y, range_z, azimuth_angle, dip_angle, std_dev);

  if (type == MATERN52)
    return new Matern52Vario(range_x, range_y, range_z, azimuth_angle, dip_angle, std_dev);

  if (type == MATERN72)
    return new Matern72Vario(range_x, range_y, range_z, azimuth_angle, dip_angle, std_dev);

  throw Exception("Unknown variogram type.");
}


Variogram::~Variogram()
{}


//
// FUNCTION: Vario::SetAngles
//
// PURPOSE:
// Set the rotation and dip angles.
//
void Variogram::SetAngles(double azimuth_angle, double dip_angle) {

  azimuth_angle_ = azimuth_angle;
  dip_angle_ = dip_angle;
  EstimateFactors();
}

//
// FUNCTION: Vario::RotateAngles
//
// PURPOSE:
// Manipulate rotation and dip angles.
//
void Variogram::RotateAngles(double rot_angle, double dip_angle) {
  azimuth_angle_     += rot_angle;
  dip_angle_     += dip_angle;
  EstimateFactors();
}

//
// FUNCTION: Vario::SetStd
//
// PURPOSE:
// Set constant standard error: sqrt(var).
// Negative sign (error) will be ignored.
//
void Variogram::SetStdDev(double std_dev) {
  var_ = std_dev*std_dev;
}

//
// FUNCTION: Vario::SetRange, Vario::SetSubRange, Vario::SetSubSubRange
//
// PURPOSE:
// Set x, y and z ranges respectively.
//
// RETURN VALUE:
// void
//
// SIDE EFFECTS:
//
// SPECIAL INSTRUCTIONS & NOTES:
//

void Variogram::SetRanges(double range_x, double range_y, double range_z) {
   if (range_x <= 0.0 || range_y <= 0.0 || range_z <= 0.0 ) {
    throw Exception("Negative ranges are not allowed.\n");
   }
   range_x_ = range_x;
   range_y_ = range_y;
   range_z_ = range_z;
   EstimateFactors();
}

void Variogram::SetRangeX(double r)
{
  if (r <= 0.0)
    throw Exception("Negative ranges are not allowed.\n");

  range_x_ = r;
  EstimateFactors();
}

void Variogram::SetRangeY(double r)
{
  if (r <= 0.0)
    throw Exception("Negative ranges are not allowed.\n");

  range_y_ = r;
  EstimateFactors();
}

void Variogram::SetRangeZ(double r)
{
  if (r <= 0.0)
    throw Exception("Negative ranges are not allowed.\n");

  range_z_ = r;
  EstimateFactors();
}

void Variogram::LimitRangeX(double limit)
{
  if (range_x_ > limit)
    SetRangeX(limit);
}

void Variogram::LimitRangeY(double limit)
{
  if (range_y_ > limit)
    SetRangeY(limit);
}

void Variogram::LimitRangeZ(double limit)
{
  if (range_z_ > limit)
    SetRangeZ(limit);
}


double Variogram::GetCorr(double dx, double dy, double dz) const
{
  double dist = Distance(dx, dy, dz);
  double corr = Corr1D(dist);
  return corr;
}

double Variogram::GetCorr(double dx, double dy) const
{
  double dist = Distance(dx, dy);
  double corr = Corr1D(dist);
  return corr;
}
double Variogram::GetCorr(double dx) const
{
  double dist = Distance(dx);
  double corr = Corr1D(dist);
  return corr;
}

//- Variograms:
double Variogram::GetVariogram(double dx, double dy, double dz) const
{
  return var_ * (1 - GetCorr(dx, dy, dz));
}

double Variogram::GetVariogram(double dx, double dy) const
{
  return var_ * (1 - GetCorr(dx, dy));
}

double Variogram::GetVariogram(double dx) const
{
  return var_* (1 - GetCorr(dx));
}

void Variogram::EstimateFactors()
{
  // Calculate factors:
  double cosRot = cos(azimuth_angle_);
  double sinRot = sin(azimuth_angle_);
  double cosDip = cos(dip_angle_);
  double sinDip = sin(dip_angle_);
  double fac1 = 1.0/(range_x_ * range_x_);
  double fac2 = 1.0/(range_y_ * range_y_);
  double fac3 = 1.0/(range_z_ * range_z_);

  txx_ = cosRot*cosRot*cosDip*cosDip*fac1
    + sinRot*sinRot*fac2
    + cosRot*cosRot*sinDip*sinDip*fac3;

  tyy_ = sinRot*sinRot*cosDip*cosDip*fac1
    + cosRot*cosRot*fac2
    + sinRot*sinRot*sinDip*sinDip*fac3;

  tzz_ = sinDip*sinDip*fac1 + cosDip*cosDip*fac3;

  txy_ = 2*(cosDip*cosDip*cosRot*sinRot*fac1
     - sinRot*cosRot*fac2
     + sinDip*sinDip*cosRot*sinRot*fac3);

  txz_ = 2*(cosRot*cosDip*sinDip*fac1
     - cosRot*cosDip*sinDip*fac3);

  tyz_ = 2*(sinRot*cosDip*sinDip*fac1
     - sinRot*cosDip*sinDip*fac3);
}
