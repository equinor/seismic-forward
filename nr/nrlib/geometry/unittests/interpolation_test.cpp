// $Id: interpolation_test.cpp 1302 2015-06-03 09:24:11Z veralh $

/// \file Unit tests for the Spline function in geometry

#include <nrlib/geometry/interpolation.hpp>

#include <boost/test/unit_test.hpp>
#include <math.h>
#include <algorithm>

using namespace NRLib;


BOOST_AUTO_TEST_CASE( Linear1D )
{
  size_t n = 10;
  double x_start = -1;
  double x_end   = 3.25;
  double dx = (x_end-x_start)/(n-1);

  double const_value = -5.3;
  std::vector<double> x_in(n, -1.0);
  std::vector<double> y_in_const (n, const_value);
  std::vector<double> y_in_linear(n, -1.0);

  for (size_t i = 0; i < x_in.size(); ++i){
    x_in[i] = x_start + i*dx;
    double x = x_in[i];

    y_in_linear[i] = 10.1 - 0.4*x;
  }

  n = 37;
  dx = (x_end - x_start)/(n-1);
  std::vector<double> x_out(n, -1.0);
  for (size_t i = 0; i < x_out.size(); ++i){
    x_out[i] = x_start + i*dx;
  }

  std::vector<double> y_out_const  = Interpolation::Interpolate1D(x_in, y_in_const,  x_out, "linear");
  std::vector<double> y_out_linear = Interpolation::Interpolate1D(x_in, y_in_linear, x_out, "linear");


  for(size_t i = 0; i < x_out.size(); ++i )
  {
    double x = x_out[i];
    BOOST_CHECK_CLOSE(y_out_const [i], const_value,  0.01);
    BOOST_CHECK_CLOSE(y_out_linear[i], 10.1 - 0.4*x, 0.01);
  }
}


BOOST_AUTO_TEST_CASE( Linear1DWithExtrapolation )
{
  size_t n = 10;
  double x_start = -1;
  double x_end   = 3.25;
  double dx = (x_end-x_start)/(n-1);

  double const_value = -5.3;
  std::vector<double> x_in(n, -1.0);
  std::vector<double> y_in_const (n, const_value);
  std::vector<double> y_in_linear(n, -1.0);

  for (size_t i = 0; i < x_in.size(); ++i){
    x_in[i] = x_start + i*dx;
    double x = x_in[i];

    y_in_linear[i] = 10.1 - 0.4*x;
  }

  n = 37;
  x_start = -10;
  x_end   = 5.5;
  dx = (x_end - x_start)/(n-1);
  std::vector<double> x_out(n, -1.0);
  for (size_t i = 0; i < x_out.size(); ++i){
    x_out[i] = x_start + i*dx;
  }

  std::vector<double> y_out_const  = Interpolation::Interpolate1D(x_in, y_in_const,  x_out, "linear");
  std::vector<double> y_out_linear = Interpolation::Interpolate1D(x_in, y_in_linear, x_out, "linear");


  for(size_t i = 0; i < x_out.size(); ++i )
  {
    double x = x_out[i];
    BOOST_CHECK_CLOSE(y_out_const [i], const_value,  0.01);
    BOOST_CHECK_CLOSE(y_out_linear[i], 10.1 - 0.4*x, 0.01);
  }
}



BOOST_AUTO_TEST_CASE( Spline1D )
{
  size_t n = 10;
  double x_start = 0;
  double x_end   = 2.25;
  double dx = (x_end-x_start)/(n-1);

  double const_value = 2.4;
  std::vector<double> x_in(n, -1.0);
  std::vector<double> y_in_const (n, const_value);
  std::vector<double> y_in_linear(n, -1.0);
  std::vector<double> y_in_square(n, -1.0);
  std::vector<double> y_in_cubic (n, -1.0);
  std::vector<double> y_in_sin   (n, -1.0);
  std::vector<double> y_in_sqrt  (n, -1.0);

  for (size_t i = 0; i < x_in.size(); ++i){
    x_in[i] = x_start + i*dx;
    double x = x_in[i];

    y_in_linear[i] = -4.4 + 2.3*x;
    y_in_square[i] = x*x;
    y_in_cubic [i] = -3 + x + 5.5*x*x - 4.4*x*x*x;
    y_in_sin   [i] = sin(x);
    y_in_sqrt  [i] = sqrt(x);
  }

  n = 37;
  dx = (x_end - x_start)/(n-1);
  std::vector<double> x_out(n, -1.0);
  for (size_t i = 0; i < x_out.size(); ++i){
    x_out[i] = x_start + i*dx;
  }

  std::vector<double> y_out_const  = Interpolation::Interpolate1D(x_in, y_in_const,  x_out, "spline");
  std::vector<double> y_out_linear = Interpolation::Interpolate1D(x_in, y_in_linear, x_out, "spline");
  std::vector<double> y_out_square = Interpolation::Interpolate1D(x_in, y_in_square, x_out, "spline");
  std::vector<double> y_out_cubic  = Interpolation::Interpolate1D(x_in, y_in_cubic,  x_out, "spline");
  std::vector<double> y_out_sin    = Interpolation::Interpolate1D(x_in, y_in_sin,    x_out, "spline");
  std::vector<double> y_out_sqrt   = Interpolation::Interpolate1D(x_in, y_in_sqrt,   x_out, "spline");


  for(size_t i = 0; i < x_out.size(); ++i )
  {
    double x = x_out[i];
    if(std::find(x_in.begin(), x_in.end(), x) != x_in.end())
    {
      // Exact at knots
      BOOST_CHECK_CLOSE(y_out_const [i], const_value,  1.0e-7);
      BOOST_CHECK_CLOSE(y_out_linear[i], -4.4 + 2.3*x, 1.0e-7);
      BOOST_CHECK_CLOSE(y_out_square[i], x*x,          1.0e-7);
      BOOST_CHECK_CLOSE(y_out_cubic [i], -3 + x + 5.5*x*x - 4.4*x*x*x, 1.0e-7);
      BOOST_CHECK_CLOSE(y_out_sin   [i], sin(x),       1.0e-7);
      BOOST_CHECK_CLOSE(y_out_sqrt  [i], sqrt(x),      1.0e-7);
    }
    else if(x > x_in[1])
    {
      // Not comparing in first interval
      BOOST_CHECK_CLOSE(y_out_const [i], const_value,  0.1);
      BOOST_CHECK_CLOSE(y_out_linear[i], -4.4 + 2.3*x, 0.1);
      BOOST_CHECK_CLOSE(y_out_square[i], x*x,          2.0);
      BOOST_CHECK_CLOSE(y_out_cubic [i], -3 + x + 5.5*x*x - 4.4*x*x*x, 1.0);
      BOOST_CHECK_CLOSE(y_out_sin   [i], sin(x),       2.0);
      BOOST_CHECK_CLOSE(y_out_sqrt  [i], sqrt(x),      2.8);
    }
  }
}


BOOST_AUTO_TEST_CASE( Spline1DExtrapolation)
{

  size_t n = 11;
  double x_start = -1;
  double x_end   = 4;
  double dx = (x_end-x_start)/(n-1);

  std::vector<double> x_in(n, -1.0);
  std::vector<double> y_in_cubic (n, -1.0);
  std::vector<double> y_in_sin   (n, -1.0);

  for (size_t i = 0; i < x_in.size(); ++i){
    x_in[i] = x_start + i*dx;
    double x = x_in[i];

    y_in_cubic[i] = -3 + x + 5.5*x*x - 4.4*x*x*x;;
    y_in_sin[i]   = sin(x);
  }

  n = 61;
  x_start = -6;
  x_end   = 9;
  dx = (x_end - x_start)/(n-1);
  std::vector<double> x_out(n, -1.0);
  for (size_t i = 0; i < x_out.size(); ++i){
    x_out[i] = x_start + i*dx;
  }

  std::vector<double> y_out_cubic_extrap = Interpolation::Interpolate1D(x_in, y_in_cubic,  x_out, "spline", 0.0);
  std::vector<double> y_out_sin_extrap   = Interpolation::Interpolate1D(x_in, y_in_sin,    x_out, "spline", 0.0);

  for(size_t i = 0; i < x_out.size(); ++i ){
    double x = x_out[i];

    // Check for equality to pre-set extrapapolation value outside defined intervals
    if(x < x_in[0] || x > x_in[x_in.size()-1])
    {
      BOOST_CHECK_EQUAL(y_out_cubic_extrap [i], 0.0);
      BOOST_CHECK_EQUAL(y_out_sin_extrap   [i], 0.0);
    }

    // Check for equality at the knots
    else if(std::find(x_in.begin(), x_in.end(), x) != x_in.end())
    {
      // Exact at knots
      BOOST_CHECK_CLOSE(y_out_cubic_extrap [i], -3 + x + 5.5*x*x - 4.4*x*x*x, 1.0e-7);
      BOOST_CHECK_CLOSE(y_out_sin_extrap   [i], sin(x),                       1.0e-7);

    }

    // Check for close value inside defined intervals
    else if(x > x_in[0] && x < x_in[2]) {
      // In first two interval we allow larger mismatch
      BOOST_CHECK_CLOSE(y_out_cubic_extrap [i], -3 + x + 5.5*x*x - 4.4*x*x*x, 36.0);
      BOOST_CHECK_CLOSE(y_out_sin_extrap   [i], sin(x),                       20.0);
    }

    else {
      BOOST_CHECK_CLOSE(y_out_cubic_extrap [i], -3 + x + 5.5*x*x - 4.4*x*x*x, 1.6);
      BOOST_CHECK_CLOSE(y_out_sin_extrap   [i], sin(x),                       2.2);
    }

  }
}
