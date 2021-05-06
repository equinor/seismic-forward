// $Id: mathutility.hpp 1191 2013-08-07 10:02:00Z perroe $

#ifndef NRLIB_MATH_FLOATUTIL_HPP
#define NRLIB_MATH_FLOATUTIL_HPP

#include <cmath>

#include "../math/constants.hpp"

namespace NRLib {

  // Using macro to allow the definition of the constant in the header file.
  #define NRLIB_MACRO_STANDARD_TOLERANCE 1.0e-6

  /// Equality.
  /// Function checking if a equals b (numerically).
  inline bool IsEqual(double lhs, double rhs, double tol = NRLIB_MACRO_STANDARD_TOLERANCE)
  {
    double abs_lhs = std::fabs(lhs);
    double nTol = abs_lhs < 1.0 ? tol : abs_lhs * tol;
    return (std::fabs(lhs-rhs) <= nTol);
  }


  /// Equals zero.
  /// Function checking if a equals zero (numerically).
  inline bool IsZero(double a,  double tol = NRLIB_MACRO_STANDARD_TOLERANCE) {
    return (IsEqual(a, 0.0, tol));
  }

  /// Function checking if number is significantly positive.
  inline bool IsPositive(double a,  double tol = NRLIB_MACRO_STANDARD_TOLERANCE)
  {
    if (a < 0.0 || IsZero(a, tol))
      return false;
    return true;
  }

  /// Function checking if number is significantly positive.
  inline bool IsNegative(double a,  double tol = NRLIB_MACRO_STANDARD_TOLERANCE)
  {
    if (a > 0.0 || IsZero(a, tol))
      return false;
    return true;
  }

  inline bool InInterval(double x, double a,
                         double b, double tol = NRLIB_MACRO_STANDARD_TOLERANCE) {

   double fa = std::fabs(a);
   double fb = std::fabs(b);

   double tolL = (fa < 1.0 ? tol : fa*tol);
   double tolR = (fb < 1.0 ? tol : fb*tol);

   return (x >= a - tolL && x <= b + tolR);
  }


  inline bool IsSameSign(double a, double b) {
    if ((a >= 0.0 && b >= 0.0) || (a < 0.0 && b < 0.0))
      return (true);
    else
      return (false);
  }

  inline bool ZeroComparedTo(double a, double b)
  {
   if (IsZero(a))
     return true;
   else if (b!=0)
     return IsZero(a/b);
   else
     return false;
  }



  /// Convert.
  /// Converts from degrees to radians.
  /// \param deg angle in degrees.
  /// \return angle in radians.
  inline double DegToRad(double deg) {
    return (deg * NRLib::Degree);
  }


  /// Convert.
  /// Convert from radians to degrees.
  /// \param deg angle in radians.
  /// \return angle in degrees.
  inline double RadToDeg(double rad) {
    return (rad * NRLib::Radian);
  }

  inline double UnwindRad(double rad) {
   return rad - NRLib::Pi * floor(rad / NRLib::Pi);
  }

  inline double Ricker(double t, double peakF)
  {
    double c = NRLib::Pi * NRLib::Pi * peakF * peakF * t * t * 1e-6;
    return (1 - 2*c) * exp(-c);
  }


  inline double ExponentialCorr(double t, double range)
  {
    return exp(-3 * std::abs(t) / range);
  }


} // namespace NRLib

#endif  // NRLIB_MATH_FLOATUTIL_HPP
