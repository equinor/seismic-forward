// $Id: variogramtypes.hpp 1360 2016-09-07 09:44:25Z vegard $

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

#ifndef NRLIB_VARIOGRAM_VARIOGRAMTYPES_HPP
#define NRLIB_VARIOGRAM_VARIOGRAMTYPES_HPP

#include <cmath>
#include <string>
#include "variogram.hpp"

namespace NRLib {
class ConstVario : public Variogram
{

public:
  /// \param[in] range_x
  /// \param[in] range_y
  /// \param[in] range_z
  /// \param[in] azimuth_angle
  /// \param[in] dip_angle
  /// \param[in] std_dev Standard deviation
  /// This class is very little tested.
  /// Can be used in 1D, 2D and 3D
  ConstVario(const double range_x,
       const double range_y = 1.0,
       const double range_z = 1.0,
       const double azimuth_angle = 0.0,
       const double dip_angle = 0.0,
       const double std_dev = 1.0);
  ConstVario();
  ConstVario(const ConstVario &vario) : Variogram(vario) {}
  // ~ConstVario();

  /// This is not tested
  virtual Variogram *Clone() const {return new ConstVario(*this);}
  virtual std::string    GetName()         const {return "constant";}

protected:
  double Corr1D(double) const {return 1.0;}
};

class ExpVario : public Variogram
{
public:
/// \param[in] range_x
  /// \param[in] range_y
  /// \param[in] range_z
  /// \param[in] azimuth_angle
  /// \param[in] dip_angle
  /// \param[in] std_dev Standard deviation
  /// This class is very little tested.
  /// Can be used in 1D, 2D and 3D
  ExpVario(const double range_x,
     const double range_y = 1.0,
     const double range_z = 1.0,
     const double azimuth_angle = 0.0,
     const double dip_angle = 0.0,
     const double std_dev = 1.0);
  ExpVario();
  ExpVario(const ExpVario &vario) : Variogram(vario) {}


  virtual Variogram *Clone() const {return new ExpVario(*this);}
  virtual std::string    GetName()         const {return "exponential";}


protected:
  virtual double Corr1D(double dist) const {
                         return exp(-3.0*dist);}
};

class SphVario : public Variogram
{
public:
  /// \param[in] range_x
  /// \param[in] range_y
  /// \param[in] range_z
  /// \param[in] azimuth_angle
  /// \param[in] dip_angle
  /// \param[in] std_dev Standard deviation
  /// This class is very little tested.
  /// Can be used in 1D, 2D and 3D
  SphVario(const double range_x,
     const double range_y = 1.0,
     const double range_z = 1.0,
     const double azimuth_angle = 0.0,
     const double dip_angle = 0.0,
     const double std_dev = 1.0);
  SphVario(const SphVario &vario) : Variogram(vario) {}
  SphVario();

  virtual Variogram *Clone() const {return new SphVario(*this);}
  virtual std::string    GetName()      const {return "spherical";}

protected:
  virtual double Corr1D(double dist) const
  {
    if(dist<1.0) return(1.0-dist*(1.5-(0.5*dist*dist)));
    else return 0.0;
  }
};

class GauVario : public Variogram
{
public:
  /// \param[in] range_x
  /// \param[in] range_y
  /// \param[in] range_z
  /// \param[in] azimuth_angle
  /// \param[in] dip_angle
  /// \param[in] std_dev Standard deviation
  /// This class is very little tested.
  /// Can be used in 1D, 2D and 3D
  GauVario(const double range_x,
     const double range_y = 1.0,
     const double range_z = 1.0,
     const double azimuth_angle = 0.0,
     const double dip_angle = 0.0,
     const double std_dev = 1.0);
  GauVario();
  GauVario(const GauVario &vario) : Variogram(vario) {}


  virtual Variogram *Clone() const {return new GauVario(*this);}
  virtual std::string    GetName()      const {return "Gaussian";}

  protected:
  virtual double Corr1D(double dist) const
  { return std::exp(-3.0*dist*dist); }
};

class GenExpVario : public Variogram
{
public:
  /// \param[in] range_x
  /// \param[in] range_y
  /// \param[in] range_z
  /// \param[in] azimuth_angle
  /// \param[in] dip_angle
  /// \param[in] std_dev Standard deviation
  /// This class is very little tested.
  /// Can be used in 1D, 2D and 3D
  GenExpVario(const double power,
    const double range_x,
    const double range_y = 1.0,
    const double range_z = 1.0,
    const double azimuth_angle = 0.0,
    const double dip_angle = 0.0,
    const double std_dev = 1.0);
  GenExpVario();
  GenExpVario(const GenExpVario &vario) : Variogram(vario) {power_ = vario.power_;}
  virtual std::string  GetName()      const {return "general exponential";}
  virtual Variogram *Clone() const {return new GenExpVario(*this);}


protected:
  virtual double Corr1D(double dist) const {
    return std::exp(-3.0 * std::pow(dist,power_));}

private:
  double      power_;
};

class Matern32Vario : public Variogram
{
public:
  /// \param[in] range_x
  /// \param[in] range_y
  /// \param[in] range_z
  /// \param[in] azimuth_angle
  /// \param[in] dip_angle
  /// \param[in] std_dev Standard deviation
  /// This class is not tested.
  Matern32Vario(const double range_x,
    const double range_y = 1.0,
    const double range_z = 1.0,
    const double azimuth_angle = 0.0,
    const double dip_angle = 0.0,
    const double std_dev = 1.0);
  Matern32Vario();
  Matern32Vario(const Matern32Vario &vario) : Variogram(vario) {}
  virtual std::string  GetName()      const { return "matern 3/2"; }
  virtual Variogram *Clone() const { return new Matern32Vario(*this); }

protected:
  virtual double Corr1D(double dist) const {
    const double sd = 4.744 * dist; // distance to ensure 0.05 correlation at range
    return std::exp(-sd) * (1.0 + sd);
  }

private:
};

class Matern52Vario : public Variogram
{
public:
  /// \param[in] range_x
  /// \param[in] range_y
  /// \param[in] range_z
  /// \param[in] azimuth_angle
  /// \param[in] dip_angle
  /// \param[in] std_dev Standard deviation
  /// This class is not tested.
  Matern52Vario(const double range_x,
    const double range_y = 1.0,
    const double range_z = 1.0,
    const double azimuth_angle = 0.0,
    const double dip_angle = 0.0,
    const double std_dev = 1.0);
  Matern52Vario();
  Matern52Vario(const Matern52Vario &vario) : Variogram(vario) {}
  virtual std::string  GetName()      const { return "matern 5/2"; }
  virtual Variogram *Clone() const { return new Matern52Vario(*this); }

protected:
  virtual double Corr1D(double dist) const {
    const double sd = 5.918 * dist; // distance to ensure 0.05 correlation at range
    return std::exp(-sd) * (1.0 + sd + sd * sd / 3.0);
  }

private:
};

class Matern72Vario : public Variogram
{
public:
  /// \param[in] range_x
  /// \param[in] range_y
  /// \param[in] range_z
  /// \param[in] azimuth_angle
  /// \param[in] dip_angle
  /// \param[in] std_dev Standard deviation
  /// This class is not tested.
  Matern72Vario(const double range_x,
    const double range_y = 1.0,
    const double range_z = 1.0,
    const double azimuth_angle = 0.0,
    const double dip_angle = 0.0,
    const double std_dev = 1.0);
  Matern72Vario();
  Matern72Vario(const Matern72Vario &vario) : Variogram(vario) {}
  virtual std::string  GetName()      const { return "matern 7/2"; }
  virtual Variogram *Clone() const { return new Matern72Vario(*this); }

protected:
  virtual double Corr1D(double dist) const {
    const double sd = 6.877 * dist; // distance to ensure 0.05 correlation at range
    return std::exp(-sd) * (1.0 + sd + 2.0/5.0 * std::pow(sd, 2) + std::pow(sd, 3) / 15.0);
  }

private:
};


}

#endif // NRLIB_VARIOGRAM_VARIOGRAMTYPES_HPP
