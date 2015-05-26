#ifndef NRLIB_GEOMETRY_INTERPOLATION_HPP
#define NRLIB_GEOMETRY_INTERPOLATION_HPP

#include <cstdlib>
#include <vector>


namespace NRLib {
class Interpolation
{
public:
  static std::vector<double> Interpolate1D(const std::vector<double> &x_in,
                                           const std::vector<double> &y_in,
                                           const std::vector<double> &x_out,
                                           const std::string          method);



private:
  static std::vector<double> Spline1D(const std::vector<double> &x_in,
                                      const std::vector<double> &y_in,
                                      const std::vector<double> &x_out);

  static std::vector<double> Linear1D(const std::vector<double> &x_in,
                                      const std::vector<double> &y_in,
                                      const std::vector<double> &x_out);

  static size_t FindNearestNeighborIndex(const double x, const std::vector<double> &x_in);

};
}
#endif

