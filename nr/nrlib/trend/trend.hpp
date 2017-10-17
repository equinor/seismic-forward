// $Id: trend.hpp 1422 2017-02-03 15:21:38Z aarnes $
#ifndef NRLIB_TREND_HPP
#define NRLIB_TREND_HPP

#include <vector>
#include "../grid/grid.hpp"
#include "../grid/grid2d.hpp"
#include "../volume/volume.hpp"
#include "../geometry/point.hpp"
#include "../iotools/fileio.hpp"

namespace NRLib {
class Trend {
public:
  Trend();

  virtual ~Trend();

  virtual Trend           * Clone()                                                                                const = 0;
  virtual void              AddConstant(double /*c*/)                                                                    = 0;
  virtual double            GetValue(double /*s1*/, double /*s2*/, double /*s3*/)                                  const = 0;
  virtual double            GetValue(double /*s1*/, double /*s2*/, double /*s3*/, const NRLib::Volume &/*volume*/) const = 0;
  virtual double            GetValue2D(double s1, double s2, const std::vector<double> &trend, size_t ns1, size_t ns2, size_t start_ind) const;
  virtual std::vector<double> GetIncrement()                                                                       const = 0;
  virtual std::vector<int>  GetTrendSize(void)                                                                     const = 0;
  virtual double            GetTrendElement(int /*i*/, int /*j*/, int /*k*/)                                       const = 0;
  virtual int               GetTrendDimension(void)                                                                const = 0;
  virtual int               GetReference(void)                                                                     const = 0;
  virtual NRLib::Point      DrawPoint(const NRLib::Volume & /*volume*/)                                            const = 0;
  virtual double            GetMaxValue(void)                                                                      const = 0;
  virtual double            GetMinValue(void)                                                                      const = 0;
  virtual double            GetMeanValue(void)                                                                     const = 0;
  virtual void              ScaleToAverage(double factor)                                                                = 0;
  virtual bool              CheckForNegative()                                                                     const = 0;
};

class TrendConstant : public Trend {
public:
  TrendConstant(double trend);
  TrendConstant(const TrendConstant & trend);
  TrendConstant();

  virtual ~TrendConstant();

  virtual Trend           * Clone()                                                                                const { return new TrendConstant(*this) ;}
  virtual void              AddConstant(double c)                                                                  {trend_+=c; }
  virtual double            GetValue(double /*s1*/, double /*s2*/, double /*s3*/)                                  const { return trend_                   ;}
  virtual double            GetValue(double /*s1*/, double /*s2*/, double /*s3*/, const NRLib::Volume &/*volume*/) const { return trend_                   ;}
  virtual double            GetValue()                                                                             const { return trend_                   ;}
  virtual double            GetTrendElement(int /*i*/, int /*j*/, int /*k*/)                                       const { return trend_                   ;}
  virtual int               GetReference(void)                                                                     const { return(0)                       ;}

  virtual NRLib::Point      DrawPoint(const NRLib::Volume & volume)                                                const;
  virtual std::vector<double> GetIncrement()                                                                       const;
  virtual std::vector<int>  GetTrendSize(void)                                                                     const;
  virtual int               GetTrendDimension(void)                                                                const;
  virtual double            GetMaxValue(void)                                                                      const { return trend_                   ;}
  virtual double            GetMinValue(void)                                                                      const { return trend_                   ;}
  virtual double            GetMeanValue(void)                                                                     const { return trend_                   ;}
  virtual void              ScaleToAverage(double factor)                                                          { trend_ = factor ;}
  virtual bool              CheckForNegative()                                                                     const;
private:
  double trend_;

};

class Trend1D : public Trend {
public:
  Trend1D(const std::vector<double> & trend, int reference = 1);
  Trend1D(const std::vector<double> & trend, int reference, double dz);
  Trend1D(const std::string &filename, int reference = 1);
  Trend1D(const Trend1D & trend);
  Trend1D();

  virtual ~Trend1D();

  virtual Trend           * Clone()                                                                               const { return new Trend1D(*this) ;}
  virtual void              AddConstant(double c);
  virtual double            GetTrendElement(int i, int j, int k)                                                  const;
  virtual int               GetReference(void)                                                                    const { return reference_         ;}
  virtual NRLib::Point      DrawPoint(const NRLib::Volume & /*volume*/)                                           const { return Point(0,0,0)       ;}

  virtual std::vector<double> GetIncrement()                                                                      const;
  virtual std::vector<int>  GetTrendSize(void)                                                                    const;
  virtual int               GetTrendDimension(void)                                                               const;
  virtual double            GetValue(double s1,
                                     double s2,
                                     double s3)                                                                   const;
  virtual double            GetValue(double s1,
                                     double s2,
                                     double s3,
                                     const NRLib::Volume &volume)                                                 const;
  virtual double            GetValue(double s1)                                                                   const { return GetValue(s1,0,0)   ;}
  virtual double            GetValue(double s1, double s2)                                                        const { return GetValue(s1,s2,0)   ;}
  virtual double            GetMaxValue(void)                                                                     const;
  virtual double            GetMinValue(void)                                                                     const;
  virtual double            GetMeanValue(void)                                                                    const;
  virtual void              ScaleToAverage(double factor);
  virtual bool              CheckForNegative()                                                                    const;

private:
  std::vector<double>         trend_;
  int                         reference_;  // 1 means s1 direction, 2 means s2 direction, 3 means s3 direction
  double                      inv_s1_inc_;

};

class Trend2D : public Trend {
public:
  Trend2D(const NRLib::Grid2D<double> & trend, int reference1 = 2, int reference2 = 3);
  Trend2D(const NRLib::Grid2D<double> & trend, int reference1, int reference2, double dz1, double dz2);
  Trend2D(const std::string &filename, int reference1 = 2, int reference2 = 3);
  Trend2D(const Trend2D & trend);
  Trend2D();

  virtual ~Trend2D();

  virtual Trend          * Clone()                                              const { return new Trend2D(*this) ;}
  virtual void             AddConstant(double c);
  virtual double           GetTrendElement(int i, int j, int k)                 const;
  virtual int              GetReference(void)                                   const { return reference1_        ;}
  virtual NRLib::Point     DrawPoint(const NRLib::Volume & /*volume*/)          const { return Point(0,0,0)       ;}

  virtual std::vector<double> GetIncrement()                                    const;
  virtual std::vector<int> GetTrendSize(void)                                   const;
  virtual int              GetTrendDimension(void)                              const;
  virtual double           GetValue(double s1,
                                    double s2,
                                    double /*s3*/)                              const;
  virtual double           GetValue(double s1,
                                    double s2,
                                    double s3,
                                    const NRLib::Volume &volume)                const;
  virtual double           GetValue(double s1,double s2)                        const { return GetValue(s1, s2, 0.0);}
  virtual double           GetMaxValue(void)                                    const;
  virtual double           GetMinValue(void)                                    const;
  virtual double           GetMeanValue(void)                                   const;
  virtual void             ScaleToAverage(double factor);
  virtual bool             CheckForNegative()                                   const;

private:
  int                           reference1_; // 1 means s1 direction, 2 means s2 direction, 3 means s3 direction
  int                           reference2_; // 1 means s1 direction, 2 means s2 direction, 3 means s3 direction
  NRLib::Grid2D<double>         trend_;
  int                           ns1_, ns2_;
  double                        inv_s1_inc_;
  double                        inv_s2_inc_;

};

class Trend3D : public Trend {
public:
  Trend3D(NRLib::Grid<double> & values);
  Trend3D(const std::string   & file_name,
          bool                  binary = false,
          Endianess             file_format = END_LITTLE_ENDIAN);
  Trend3D(const Trend3D & trend);
  Trend3D();

  virtual ~Trend3D();

  virtual Trend          * Clone()                                              const { return new Trend3D(*this) ;}
  virtual void             AddConstant(double c);
  virtual double           GetTrendElement(int i, int j, int k)                 const { return trend_(i,j,k)      ;}
  virtual int              GetReference(void)                                   const { return(0)                 ;}

  virtual NRLib::Point     DrawPoint(const NRLib::Volume & /*volume*/)          const { return Point(0,0,0)       ;}

  virtual int              GetTrendDimension(void)                              const;
  virtual std::vector<double> GetIncrement()                                    const;
  virtual std::vector<int> GetTrendSize(void)                                   const;

  virtual double           GetValue(double s1,
                                    double s2,
                                    double s3)                                  const;

  virtual double           GetValue(double s1,
                                    double s2,
                                    double s3,
                                    const NRLib::Volume &volume)                const;

  virtual double           GetMaxValue(void)                                    const;
  virtual double           GetMinValue(void)                                    const;
  virtual double           GetMeanValue(void)                                   const;
  virtual void             ScaleToAverage(double factor);
  virtual bool             CheckForNegative()                                   const;
private:
  NRLib::Grid<double> trend_;


};

}

#endif
