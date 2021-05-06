// $Id: variogram.hpp 1360 2016-09-07 09:44:25Z vegard $

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

#ifndef NRLIB_VARIOGRAM_VARIOGRAM_HPP
#define NRLIB_VARIOGRAM_VARIOGRAM_HPP

#include <string>
#include <cmath>

namespace NRLib {
class Variogram
{
public:
  /// \param[in] range_x
  /// \param[in] range_y
  /// \param[in] range_z
  /// \param[in] azimuth_angle
  /// \param[in] dip_angle
  /// \param[in] std_dev       Standard deviation
  /// Can be used in 1D, 2D and 3D
  enum Type {CONSTANT, EXPONENTIAL, SPHERICAL, GAUSSIAN, GENERAL_EXPONENTIAL, MATERN32, MATERN52, MATERN72, N_TYPES };
  Variogram(double range_x,
            double range_y       = 1.0,
            double range_z       = 1.0,
            double azimuth_angle = 0.0,
            double dip_angle     = 0.0,
            double std_dev       = 1.0);
  Variogram();
  virtual ~Variogram();

  static  Variogram*  Create(Type type, double power = 1.5, double range_x = 1.0,
                             double range_y = 1.0, double range_z = 1.0, double azimuth_angle = 0.0,
                             double dip_angle = 0.0, double std_dev = 1.0);
  virtual Variogram*  Clone()    const = 0;
  virtual std::string GetName()  const = 0;

  double              GetRangeX() const {return range_x_;}
  double              GetRangeY() const {return range_y_;}
  double              GetRangeZ() const {return range_z_;}
  double              GetAzimuthAngle() const {return azimuth_angle_;}
  double              GetDipAngle()     const {return dip_angle_;}
  double              GetStdDev()       const {return sqrt(var_);}

  void RotateAngles(double azimuth_angle, double dip_angle = 0.0);
  void SetAngles(double azimuth_angle, double dip_angle = 0.0);
  void SetRangeX(double r);
  void SetRangeY(double r);
  void SetRangeZ(double r);
  /// If range is larger than limit, set it to limit
  void LimitRangeX(double limit);
  /// If subrange is larger than limit, set it to limit
  void LimitRangeY(double limit);
  /// If subrange is larger than limit, set it to limit
  void LimitRangeZ(double limit);

  void SetStdDev(double std_dev);
  void SetRanges(double range_x, double range_y, double range_z);
  /// Correlation function with distance as input for 3D.
  virtual double GetCorr(double dx, double dy, double dz) const;
  /// Correlation function with distance as input for 2D.
  virtual double GetCorr(double dx, double dy) const;
  /// Correlation function with distance as input for 1D.
  virtual double GetCorr(double dx) const;
  /// Correlation functions with points as input for 3D.
  double GetCorrpoint(double x1, double y1, double z1, double x2, double y2, double z2) const {return GetCorr(x2-x1, y2-y1, z2-z1);}
  /// Correlation functions with points as input for 2D.
  double GetCorrpoint(double x1, double y1, double x2, double y2) const {return GetCorr(x2-x1, y2-y1);}
  /// Correlation functions with points as input for 1D.
  double GetCorrpoint(double x1, double x2) const {return GetCorr(x2-x1);}
  /// Variogram function in 3D
  double GetVariogram(double dx, double dy, double dz) const;
  /// Variogram function in 2D.
  double GetVariogram(double dx, double dy) const;
  /// Variogram function in 1D.
  double GetVariogram(double dx) const;

  /// Covariance function with distance as input in 3D.
  double GetCov(double dx, double dy, double dz) const { return var_ * GetCorr(dx, dy, dz);}
  /// Covariance function with distance as input in 2D.
  double GetCov(double dx, double dy) const { return var_ * GetCorr(dx, dy);}
  /// Covariance function with distance as input in 1D.
  double GetCov(double dx) const {return var_*GetCorr(dx);}
  /// Covariance function with points as input in 3D.
  double GetCovpoint(double x1, double y1, double z1, double x2, double y2, double z2) const { return var_*GetCorr(x2-x1, y2-y1, z2-z1);}
  /// Covariance function with points as input in 2D.
  double GetCovpoint(double x1, double y1, double x2, double y2) const { return var_*GetCorr(x2-x1, y2-y1);}
  /// Covariance function with points as input in 1D.
  double GetCovpoint(double x1, double x2) const { return var_*GetCorr(x2-x1);}

protected:
  /// Defines 1D correlation function. Different for every variogram.
  virtual double Corr1D(double dist) const = 0 ;

private:
  void EstimateFactors();
  /// Normalized distance considering anisotropy: (3D)
  double Distance(double dx, double dy, double dz) const
  { return std::sqrt(txx_*dx*dx + tyy_*dy*dy + tzz_*dz*dz + txy_*dx*dy
                     + txz_*dx*dz + tyz_*dy*dz); }
  /// Normalized distance considering anisotropy: (2D)
  double Distance(double dx, double dy) const
  { return std::sqrt(txx_*dx*dx + tyy_*dy*dy + txy_*dx*dy); }
  /// Normalized distance considering anisotropy: (1D)
  double Distance(double dx) const
  { return std::sqrt(txx_*dx*dx);}

  double range_x_;
  double range_y_;
  double range_z_;
  double azimuth_angle_;
  double dip_angle_;
  /// Variance
  double var_;

  // Variables precalculated for anisotropy:
  double      txx_;
  double      tyy_;
  double      tzz_;
  double      txy_;
  double      txz_;
  double      tyz_;
};

} // namespace NRLib

#endif //NRLIB_VARIOGRAM_VARIOGRAM_HPP
