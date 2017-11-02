#ifndef NRLIB_STATISTICS_CONTINUSOUSKRIGINGDATA2D_HPP
#define NRLIB_STATISTICS_CONTINUSOUSKRIGINGDATA2D_HPP

#include <vector>
#include "nrlib/statistics/krigingdata2d.hpp"
#include "nrlib/geometry/point.hpp"
#include "nrlib/surface/regularsurface.hpp"

namespace NRLib {
class ContinuousKrigingData2D : public KrigingData2D
{
public:
  ContinuousKrigingData2D(int n_data = 0);
  ContinuousKrigingData2D(const std::vector<Point* >   & data_points,
                          const RegularSurface<double> & trend_grid);
  ~ContinuousKrigingData2D(void);

  KrigingData2D* Clone()                              const { return new ContinuousKrigingData2D(*this); }
  KrigingData2D* NewInstance()                        const { return new ContinuousKrigingData2D(); }

  void                         AddData(double   x, double   y, double   value);

  const std::vector<double>  & GetData(void)          const { return data_                         ;}
  const std::vector<double>  & GetXVector(void)       const { return xvector_                      ;}
  const std::vector<double>  & GetYVector(void)       const { return yvector_                      ;}

  inline int                   GetNumberOfData(void) const { return static_cast<int>(data_.size()) ;}

  void                         WriteToFile(const std::string & name) const;

  size_t                       CountDataInBlock(size_t imin, size_t imax,
                                                size_t jmin, size_t jmax) const;
  void                         AddDataToBlock  (KrigingData2D & kriging_data_block,
                                                size_t imin, size_t imax,
                                                size_t jmin, size_t jmax) const;

private:
  std::vector<double>         data_;
  std::vector<double>         xvector_;
  std::vector<double>         yvector_;

  // To convert from indices to coordinates
  double                      xmin_;
  double                      ymin_;
  double                      dx_;
  double                      dy_;
};
}
#endif
