#include <vector>
#include <string>
#include <iostream>


//#include <cstring>
//#include <typeinfo>
#include <assert.h>     /* assert */

#include "interpolation.hpp"

using namespace NRLib;



std::vector<double> Interpolation::Interpolate1D(const std::vector<double> &x_in,
                                                 const std::vector<double> &y_in,
                                                 const std::vector<double> &x_out,
                                                 const std::string          method)
{
  if(method == "linear")
  {
    return Linear1D(x_in, y_in, x_out);
  }
   else if(method == "spline")
  {
    return Spline1D(x_in, y_in, x_out);
  }
   else
   {
     std::cout << "Interpolation method " << method  << " is not implemented.\n"
               << "Valid methods are: linear, spline\n";
   }
}

std::vector<double> Interpolation::Linear1D(const std::vector<double> &x_in,
                                            const std::vector<double> &y_in,
                                            const std::vector<double> &x_out){

  std::vector<double> dx(x_in.size()), dy(x_in.size()), slope(x_in.size()), intercept(x_in.size());
  std::vector<double> y_out(x_out.size());

  for (size_t i = 0; i < x_in.size() - 1; i++) {
    if (x_in[i+1] == x_in[i] && y_in[i+1] == y_in[i] && i > 0) {
      dx[i] = dx[i-1];
      dy[i] = dy[i-1];
    }
    else {
      dx[i] = x_in[i+1] - x_in[i];
      dy[i] = y_in[i+1] - y_in[i];
    }
    slope[i] = dy[i] / dx[i];
    intercept[i] = y_in[i] - x_in[i] * slope[i];
  }
  dx[x_in.size()-1]        = dx[x_in.size()-2];
  dy[x_in.size()-1]        = dy[x_in.size()-2];
  slope[x_in.size()-1]     = slope[x_in.size()-2];
  intercept[x_in.size()-1] = intercept[x_in.size()-2];
  for (size_t i = 0; i < y_out.size(); i++) {
    int idx = FindNearestNeighborIndex(x_out[i], x_in);
    if (idx == -1)
      idx = 0;
    y_out[i] = slope[idx] * x_out[i] + intercept[idx];
  }

  return y_out;
}



std::vector<double> Interpolation::Spline1D(const std::vector<double> &x_in,
                                            const std::vector<double> &y_in,
                                            const std::vector<double> &x_out)
{
  // Cubic spline algorithm using natural boundary conditions from:
  // Numerical Analysis, David Kincaid and Ward Cheney, Second Ed. (1996)

  assert(x_in.size()  == y_in.size());

  size_t n_in  = x_in.size();
  size_t n_out = x_out.size();
  std::vector<double> y_out(n_out, -1.0);

  std::vector<double> h(n_in-1, 0);
  std::vector<double> b(n_in-1, 0);
  for(size_t i = 0; i < n_in - 1; ++i) {
    h[i] = x_in[i+1] - x_in[i];
    b[i] = 6*(y_in[i+1] - y_in[i])/h[i];
  }

  std::vector<double> u(n_in, 0);
  std::vector<double> v(n_in, 0);
  u[1] = 2*(h[0] + h[1]);
  v[1] = b[1] - b[0];
  for (size_t i = 2; i < n_in - 1; ++i) {
    u[i] = 2*(h[i] + h[i-1]) - h[i-1]*h[i-1]/u[i-1];
    v[i] = b[i] - b[i-1] - h[i-1]*v[i-1]/u[i-1];
  }

  std::vector<double> z(n_in, 0);
  z[n_in-1] = 0;                          // Natural boundary conditions
  for (size_t i = n_in-2; i > 0; --i) {
    z[i] = (v[i] - h[i]*z[i+1])/u[i];
  }
  z[0] = 0;                              // Natural boundary conditions

  std::vector<double> A(n_in-1, 0);
  std::vector<double> B(n_in-1, 0);
  std::vector<double> C(n_in-1, 0);
  for(size_t i = 0; i < n_in-1; ++i) {
    A[i] = 1/(6*h[i])*(z[i+1] - z[i]);
    B[i] = z[i]/2;
    C[i] = -h[i]/6*z[i+1] - h[i]/3*z[i] + 1/h[i]*(y_in[i+1]-y_in[i]);
  }
  for(size_t i = 0; i < n_out; ++i) {
    double x = x_out[i];
    int idx = FindNearestNeighborIndex(x, x_in);

    if(idx == -1) {
      idx = 0;
    }
    y_out[i] = y_in[idx] + (x - x_in[idx])*(C[idx] + (x - x_in[idx])*(B[idx] + (x - x_in[idx])*A[idx]));
  }

return y_out;
}

//--------------------------------------------------------------
size_t Interpolation::FindNearestNeighborIndex(const double x, const std::vector<double> &x_in)
//--------------------------------------------------------------
{
  double dist = 1e10;
  int idx = -1;
  for (size_t i = 0; i < x_in.size(); ++i) {
    float newDist = x - x_in[i];   // Er det en grunn til at er float her?
    if ( newDist > 0 && newDist < dist ) {
      dist = newDist;
      idx = i;
    }
  }
  return idx;
}

