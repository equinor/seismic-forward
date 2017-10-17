// $Id: variogramtypes.cpp 1360 2016-09-07 09:44:25Z vegard $

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

#include "variogramtypes.hpp"

namespace NRLib {

ConstVario::ConstVario(const double range_x,
           const double range_y,
           const double range_z,
           const double azimuth_angle,
           const double dip_angle,
           const double std_dev)
  : Variogram(range_x, range_y, range_z, azimuth_angle, dip_angle, std_dev)
{}

ConstVario::ConstVario() : Variogram()
{}

ExpVario::ExpVario(const double range_x,
       const double range_y,
       const double range_z,
       const double azimuth_angle,
       const double dip_angle,
       const double std_dev)
  : Variogram(range_x, range_y, range_z, azimuth_angle, dip_angle, std_dev)
{}

ExpVario::ExpVario() : Variogram()
{}

SphVario::SphVario(const double range_x,
       const double range_y,
       const double range_z,
       const double azimuth_angle,
       const double dip_angle,
       const double std_dev)
  : Variogram(range_x, range_y, range_z, azimuth_angle, dip_angle, std_dev)
{}

SphVario::SphVario() : Variogram() {}

GauVario::GauVario(const double range_x,
       const double range_y,
       const double range_z,
       const double azimuth_angle,
       const double dip_angle,
       const double std_dev)
  : Variogram(range_x, range_y, range_z, azimuth_angle, dip_angle, std_dev)
{}

GauVario::GauVario() : Variogram()
{}

GenExpVario::GenExpVario(const double power,
       const double range_x,
       const double range_y,
       const double range_z,
       const double azimuth_angle,
       const double dip_angle,
       const double std_dev)
  : Variogram(range_x, range_y, range_z, azimuth_angle, dip_angle, std_dev)
{
  power_ = power;
}

GenExpVario::GenExpVario() : Variogram()
{
  power_ = 1.0;
}

Matern32Vario::Matern32Vario(const double range_x,
       const double range_y,
       const double range_z,
       const double azimuth_angle,
       const double dip_angle,
       const double std_dev)
  : Variogram(range_x, range_y, range_z, azimuth_angle, dip_angle, std_dev)
{}

Matern32Vario::Matern32Vario() : Variogram()
{}

Matern52Vario::Matern52Vario(const double range_x,
       const double range_y,
       const double range_z,
       const double azimuth_angle,
       const double dip_angle,
       const double std_dev)
  : Variogram(range_x, range_y, range_z, azimuth_angle, dip_angle, std_dev)
{}

Matern52Vario::Matern52Vario() : Variogram()
{}

Matern72Vario::Matern72Vario(const double range_x,
       const double range_y,
       const double range_z,
       const double azimuth_angle,
       const double dip_angle,
       const double std_dev)
  : Variogram(range_x, range_y, range_z, azimuth_angle, dip_angle, std_dev)
{}

Matern72Vario::Matern72Vario() : Variogram()
{}

} // namespace NRLib
