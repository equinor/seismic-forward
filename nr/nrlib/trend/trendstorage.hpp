// $Id: trendstorage.hpp 1499 2017-05-29 19:01:53Z perroe $
#ifndef NRLIB_TRENDSTORAGE_HPP
#define NRLIB_TRENDSTORAGE_HPP

#include "trend.hpp"
#include <stdio.h>
#include <string>

namespace NRLib {
class TrendStorage {
public:
  TrendStorage();

  virtual ~TrendStorage();

  virtual TrendStorage * Clone()       const = 0;
  virtual bool           GetEstimate() const = 0;

  virtual Trend * GenerateTrend(const std::string                       & /*path*/,
                                const std::vector<std::string>          & /*trend_cube_parameters*/,
                                const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                const std::vector<std::vector<double> > & /*blocked_logs*/,
                                const std::vector<std::vector<double> > & /*s1*/,
                                const std::vector<std::vector<double> > & /*s2*/,
                                const int                               & /*type, mean or var */,
                                const NRLib::Trend                      * /*mean_trend*/,
                                std::string                             & /*errTxt*/) const = 0;
  enum type{MEAN,
            VAR,
            SD};

};

//-------------------------------------------------------------------//

class TrendConstantStorage : public TrendStorage {
public:
  TrendConstantStorage(const double & value,
                       const bool   & estimate);

  TrendConstantStorage(const TrendConstantStorage & trend_storage);

  TrendConstantStorage();

  virtual ~TrendConstantStorage();

  virtual TrendStorage * Clone()   const { return new TrendConstantStorage(*this) ;}
  virtual double         GetMean() const { return mean_value_                     ;}
  virtual bool           GetEstimate()                                                                 const { return estimate_                 ;}

  virtual Trend * GenerateTrend(const std::string                       & path,
                                const std::vector<std::string>          & trend_cube_parameters,
                                const std::vector<std::vector<double> > & trend_cube_sampling,
                                const std::vector<std::vector<double> > & blocked_logs,
                                const std::vector<std::vector<double> > & s1,
                                const std::vector<std::vector<double> > & s2,
                                const int                               & type,
                                const NRLib::Trend                      * mean_trend,
                                std::string                             & errTxt) const;

private:
  double mean_value_;
  bool   estimate_;
};

//-------------------------------------------------------------------//

class Trend1DStorage : public TrendStorage {
public:
  Trend1DStorage(const std::string & file_name,
                 const std::string & reference_parameter,
                 const bool        & estimate);

  Trend1DStorage(const Trend1DStorage & trend_storage);

  Trend1DStorage();

  virtual ~Trend1DStorage();

  virtual TrendStorage * Clone()                                                                       const { return new Trend1DStorage(*this) ;}
  virtual bool           GetEstimate()                                                                 const { return estimate_                 ;}

  virtual Trend * GenerateTrend(const std::string                       & path,
                                const std::vector<std::string>          & trend_cube_parameters,
                                const std::vector<std::vector<double> > & trend_cube_sampling,
                                const std::vector<std::vector<double> > & blocked_logs,
                                const std::vector<std::vector<double> > & s1,
                                const std::vector<std::vector<double> > & s2,
                                const int                               & type,
                                const NRLib::Trend                      * mean_trend,
                                std::string                             & errTxt) const;

private:
  std::string file_name_;
  std::string reference_parameter_;
  bool        estimate_;
};

//-------------------------------------------------------------------//

class Trend2DStorage : public TrendStorage {
public:
  Trend2DStorage(const std::string & file_name,
                 const std::string & reference_parameter1,
                 const std::string & reference_parameter2,
                 const bool        & estimate);

  Trend2DStorage(const Trend2DStorage & trend_storage);

  virtual ~Trend2DStorage();

  virtual TrendStorage * Clone()                                                                       const { return new Trend2DStorage(*this) ;}
  virtual bool           GetEstimate()                                                                 const { return estimate_                 ;}

  virtual Trend * GenerateTrend(const std::string                       & path,
                                const std::vector<std::string>          & trend_cube_parameters,
                                const std::vector<std::vector<double> > & trend_cube_sampling,
                                const std::vector<std::vector<double> > & blocked_logs,
                                const std::vector<std::vector<double> > & s1,
                                const std::vector<std::vector<double> > & s2,
                                const int                               & type,
                                const NRLib::Trend                      * mean_trend,
                                std::string                             & errTxt)   const;

private:
  std::string file_name_;
  std::string reference_parameter_one_;
  std::string reference_parameter_two_;
  bool        estimate_;
};

}
#endif
