// $Id: trendstorage.cpp 1499 2017-05-29 19:01:53Z perroe $
#include "trendstorage.hpp"
#include "trendkit.hpp"
#include "trend.hpp"
#include "../grid/grid2d.hpp"

#include "../surface/regularsurface.hpp"
#include "../iotools/fileio.hpp"

#include <fstream>

using namespace NRLib;

TrendStorage::TrendStorage()
{
}

TrendStorage::~TrendStorage()
{
}

//-------------------------------------------------------------------//

TrendConstantStorage::TrendConstantStorage(const double & value,
                                           const bool   & estimate)
: mean_value_(value),
  estimate_(estimate)
{
}

TrendConstantStorage::TrendConstantStorage()
{
}

TrendConstantStorage::TrendConstantStorage(const TrendConstantStorage & trend_storage)
: mean_value_(trend_storage.mean_value_),
  estimate_(trend_storage.estimate_)
{
}

TrendConstantStorage::~TrendConstantStorage()
{
}

Trend *
TrendConstantStorage::GenerateTrend(const std::string                       & /*path*/,
                                    const std::vector<std::string>          & /*trend_cube_parameters*/,
                                    const std::vector<std::vector<double> > & trend_cube_sampling,
                                    const std::vector<std::vector<double> > & blocked_logs,
                                    const std::vector<std::vector<double> > & s1,
                                    const std::vector<std::vector<double> > & s2,
                                    const int                               & type,
                                    const NRLib::Trend                      * mean_trend,
                                    std::string                             & errTxt) const
{
  Trend * trend = NULL;
  std::string tmp_err_txt;

  double estimated_constant = 0.0;

  if (estimate_ == false)
    estimated_constant = mean_value_;

  else {
    if (type == TrendStorage::MEAN) {
      double mean = 0.0;
      EstimateConstantTrend(blocked_logs, mean);
      estimated_constant = mean;
    }

    else if (type == TrendStorage::VAR) {
      double var = 0.0;

      if (typeid(*mean_trend) == typeid(NRLib::TrendConstant)) {

        /* --- trend --- */
        std::vector<double> z;
        PreprocessData0D(blocked_logs, z);

        int    i_dummy = 0;
        int    j_dummy = 0;
        int    k_dummy = 0;
        double z_mean  = mean_trend->GetTrendElement(i_dummy, j_dummy, k_dummy);

        /* --- variance --- */
        int z_n = 0;
        for (size_t i = 0; i < z.size(); i++) {
          if (z_mean != RMISSING && z[i] != RMISSING) {
            var = var + (z[i] - z_mean)*(z[i] - z_mean);
            z_n++;
          }
        }
        if (z_n > 0) {
          var = var/(z_n - 1);
        }
        else {
          var = RMISSING;
          tmp_err_txt +=  "Error: Unable to estimate variance as the wells only contain missing values.\n";
        }
      }
      else if (typeid(*mean_trend) == typeid(NRLib::Trend1D)) {

        int mean_reference = mean_trend->GetReference();

        std::vector<std::vector<double> >  s;
        if (mean_reference == 1)
          s = s1;
        else
          s = s2;

        std::vector<double> x;
        std::vector<double> z;
        PreprocessData1D(s,
                         blocked_logs,
                         trend_cube_sampling[mean_reference - 1],
                         x,
                         z);

        /* --- trend --- */
        std::vector<double> trend_cube_trend(trend_cube_sampling[mean_reference - 1].size(), 0.0);
        int k_dummy = 0;

        if (mean_reference == 1) {
          int j_dummy = 0;

          for (int i = 0; i < static_cast<int>(trend_cube_trend.size()); ++i)
            trend_cube_trend[i] = mean_trend->GetTrendElement(i, j_dummy, k_dummy);
        }
        else {
          int i_dummy = 0;

          for (int j = 0; j < static_cast<int>(trend_cube_trend.size()); ++j)
            trend_cube_trend[j] = mean_trend->GetTrendElement(i_dummy, j, k_dummy);
        }


        std::vector<double> z_mean;
        LinearInterpolation(trend_cube_sampling[mean_reference - 1], trend_cube_trend, x, z_mean);

        /* --- variance --- */
        int                 z_n   = 0;
        for (size_t i = 0; i < z.size(); i++) {
          if (z_mean[i] != RMISSING && z[i] != RMISSING) {
            var = var + (z[i] - z_mean[i])*(z[i] - z_mean[i]);
            z_n++;
          }
        }
        if (z_n > 0)
          var = var/(z_n - 1);
        else {
          var = RMISSING;
          tmp_err_txt +=  "Error: Unable to estimate variance as the wells only contain missing values.\n";
        }
      }

      else if (typeid(*mean_trend) == typeid(NRLib::Trend2D)) {
        size_t trend_cube_n1 = trend_cube_sampling[0].size();
        size_t trend_cube_n2 = trend_cube_sampling[1].size();

        /* -- obtain trend surface -- */
        std::vector<std::vector<double> > trend_cube_trend_2D(trend_cube_n1, std::vector<double>(trend_cube_n2, RMISSING));
        int k_dummy = 0;
        for (int i = 0; i < static_cast<int>(trend_cube_n1); ++i) {
          for (int j = 0; j < static_cast<int>(trend_cube_n2); ++j)
            trend_cube_trend_2D[i][j] = mean_trend->GetTrendElement(i, j, k_dummy);
        }

        /* preprocess and variance */
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;

        PreprocessData2D(s1,
                         s2,
                         blocked_logs,
                         trend_cube_sampling,
                         x,
                         y,
                         z);

        int z_n   = 0;
        for (size_t i = 0; i < z.size(); i++) {
          double z_mean_i = Interpolate(trend_cube_sampling[0], trend_cube_sampling[1], trend_cube_trend_2D, x[i], y[i]);
          if (z_mean_i != RMISSING && z[i] != RMISSING) {
            var = var + (z[i] - z_mean_i)*(z[i] - z_mean_i);
            z_n++;
          }
        }
        if (z_n > 0)
          var = var/(z_n - 1);
        else {
          var = RMISSING;
          tmp_err_txt +=  "Error: Unable to estimate variance as the wells only contain missing values.\n";
        }
      }
      estimated_constant = var;
    }
  }

  if (tmp_err_txt == "")
    trend = new TrendConstant(estimated_constant);
  else
    errTxt += tmp_err_txt;

  return trend;
}

//-------------------------------------------------------------------//

Trend1DStorage::Trend1DStorage(const std::string & file_name,
                               const std::string & reference_parameter,
                               const bool        & estimate)
: file_name_(file_name),
  reference_parameter_(reference_parameter),
  estimate_(estimate)
{
}

Trend1DStorage::Trend1DStorage()
{
}

Trend1DStorage::Trend1DStorage(const Trend1DStorage & trend_storage)
: file_name_(trend_storage.file_name_),
  reference_parameter_(trend_storage.reference_parameter_),
  estimate_(trend_storage.estimate_)
{
}

Trend1DStorage::~Trend1DStorage()
{
}

Trend *
Trend1DStorage::GenerateTrend(const std::string                       & path,
                              const std::vector<std::string>          & trend_cube_parameters,
                              const std::vector<std::vector<double> > & trend_cube_sampling,
                              const std::vector<std::vector<double> > & blocked_logs,
                              const std::vector<std::vector<double> > & s1,
                              const std::vector<std::vector<double> > & s2,
                              const int                               & type,
                              const NRLib::Trend                      * mean_trend,
                              std::string                             & errTxt) const
{
  Trend * trend = NULL;

  int reference = 0;
  for(int i=0; i<static_cast<int>(trend_cube_parameters.size()); i++) {
    if(reference_parameter_ == trend_cube_parameters[i])
      reference = i+1;
  }

  if(reference == 0) {
    errTxt += "The reference parameter of the 1D trend in "+file_name_+" \n"
              " is not the same as the parameter names of the trend cubes\n";
    return(0);
  }

  int                 n_cube_samples = static_cast<int>(trend_cube_sampling[reference-1].size());
  double              increment      = trend_cube_sampling[reference-1][1]-trend_cube_sampling[reference-1][0];

  std::vector<double> resampled_trend(static_cast<int>(n_cube_samples));;
  std::vector<double> trend_values;
  std::vector<double> trend_sampling;

  if(estimate_ == true) {

    double scale = 1.414214;

    if (type == TrendStorage::MEAN) {
      std::vector<std::vector<double> > s;
      std::vector<double>               trend_estimate;
      std::vector<double>               x;
      std::vector<double>               y;

      if (reference == 1) {
        s = s1;
      }
      else if (reference == 2) {
        s = s2;
      }
      if (reference == 1 || reference == 2) {
        PreprocessData1D(s,
                         blocked_logs,
                         trend_cube_sampling[reference - 1],
                         x,
                         y);

        double bandwidth = CalculateBandwidth(x, scale, 0.2);

        if (CheckConfigurations1D(x, y, trend_cube_sampling[reference - 1], bandwidth,errTxt)) {
          EstimateTrend1D(x,
                          y,
                          trend_cube_sampling[reference - 1],
                          trend_estimate,
                          bandwidth,
                          errTxt);
        }
      }
      if (errTxt == "") {
        trend = new Trend1D(trend_estimate, reference, increment);
      }
    }

    if (type == TrendStorage::VAR) {
      std::vector<double>               variance_estimate;

      if (typeid(*mean_trend) == typeid(NRLib::TrendConstant) || typeid(*mean_trend) == typeid(NRLib::Trend1D)) {
        std::vector<std::vector<double> > s;
        std::vector<double>               trend_cube_trend;

        if (reference == 1) {
          s = s1;
          /* -- obtain trend -- */
          trend_cube_trend.resize(trend_cube_sampling[reference - 1].size(), 0.0);
          int j_dummy = 0;
          int k_dummy = 0;

          for (int i = 0; i < static_cast<int>(trend_cube_trend.size()); ++i) {
            trend_cube_trend[i] = mean_trend->GetTrendElement(i, j_dummy, k_dummy);
          }
        }
        else if (reference == 2) {
          s = s2;
          /* -- obtain trend -- */
          trend_cube_trend.resize(trend_cube_sampling[reference - 1].size(), 0.0);
          int i_dummy = 0;
          int k_dummy = 0;

          for (int j = 0; j < static_cast<int>(trend_cube_trend.size()); ++j) {
            trend_cube_trend[j] = mean_trend->GetTrendElement(i_dummy, j, k_dummy);
          }
        }

        if (reference == 1 || reference == 2) {
          std::vector<double>               x;
          std::vector<double>               y;

          PreprocessData1D(s,
                           blocked_logs,
                           trend_cube_sampling[reference - 1],
                           x,
                           y);

          double bandwidth = CalculateBandwidth(x, scale, 0.2);

          if (CheckConfigurations1D(x, y, trend_cube_sampling[0], bandwidth,errTxt)) {
            EstimateVariance1D(x,
                               y,
                               trend_cube_sampling[reference - 1],
                               trend_cube_trend,
                               variance_estimate,
                               bandwidth,
                               errTxt);
          }
        }
      }
      else if (typeid(*mean_trend) == typeid(NRLib::Trend2D)) {
        size_t trend_cube_n1 = trend_cube_sampling[0].size();
        size_t trend_cube_n2 = trend_cube_sampling[1].size();

        /* -- obtain trend surface -- */
        std::vector<std::vector<double> > trend_cube_trend_2D(trend_cube_n1, std::vector<double>(trend_cube_n2, RMISSING));
        int k_dummy = 0;
        for (int i = 0; i < static_cast<int>(trend_cube_n1); ++i) {
          for (int j = 0; j < static_cast<int>(trend_cube_n2); ++j) {
            trend_cube_trend_2D[i][j] = mean_trend->GetTrendElement(i, j, k_dummy);
          }
        }
        /* -------------------------- */

        /* preprocess and mean-correct observations */
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;

        PreprocessData2D(s1,
                         s2,
                         blocked_logs,
                         trend_cube_sampling,
                         x,
                         y,
                         z);

        for (size_t i = 0; i < z.size(); i++) {
          double z_trend_i = Interpolate(trend_cube_sampling[0], trend_cube_sampling[1], trend_cube_trend_2D, x[i], y[i]);
          if (z_trend_i != RMISSING && z[i] != RMISSING) {
            z[i] = (z[i] - z_trend_i);
          } else {
            z[i] = RMISSING;
          }
        }
        /* ---------------------------------------- */

        // the observations are already mean-corrected.
        std::vector<double> trend_cube_trend_dummy(trend_cube_sampling[reference - 1].size(), 0.0);

        if (reference == 1) {
          double bandwidth = CalculateBandwidth(x, scale, 0.2);

          if (CheckConfigurations1D(x, z, trend_cube_sampling[reference - 1], bandwidth, errTxt)) {
            EstimateVariance1D(x,
                               z,
                               trend_cube_sampling[reference - 1],
                               trend_cube_trend_dummy,
                               variance_estimate,
                               bandwidth,
                               errTxt);
          }
        }
        else if (reference == 2) {
          double bandwidth = CalculateBandwidth(y, scale, 0.2);

          if (CheckConfigurations1D(y, z, trend_cube_sampling[reference - 1], bandwidth, errTxt)) {
            EstimateVariance1D(y,
                               z,
                               trend_cube_sampling[reference - 1],
                               trend_cube_trend_dummy,
                               variance_estimate,
                               bandwidth,
                               errTxt);
          }
        }
      }
      if (errTxt == "")
        trend = new Trend1D(variance_estimate, reference, increment);
    }
  }
  else {

    std::string file_name = path + file_name_;

    int file_format = GetTrend1DFileFormat(file_name, errTxt);

    if(file_format < 0) {
      errTxt += "Invalid 1D trend file\n";
      return(0);
    }
    else {
      double              s_min;
      double              dz;

      ReadTrend1D(file_name,errTxt,trend_values,s_min,dz);

      double s_max          = s_min + dz*static_cast<int>(trend_values.size());

      if(trend_cube_sampling[reference-1][0] < s_min) {
        errTxt += "The mimimum value of the 1D trend in "+file_name_+" \n"
                  " is lower than the minimum value of "+trend_cube_parameters[reference-1]+"\n";
      }
      else if(trend_cube_sampling[reference-1][n_cube_samples-1] > s_max) {
        errTxt += "The maximum value of the 1D trend in "+file_name_+" \n"
                  " is higher than the maximum value of "+trend_cube_parameters[reference-1]+"\n";
      }
      else {

        trend_sampling.resize(trend_values.size());

        for(int i=0; i<static_cast<int>(trend_values.size()); i++)
          trend_sampling[i] = s_min + i*dz;

      }
    }

    if(errTxt == "") {
      ResampleTrend1D(trend_sampling,
                      trend_values,
                      trend_cube_sampling[reference-1],
                      resampled_trend);

      trend = new Trend1D(resampled_trend, reference, increment);

    }
  }
  return trend;
}

//-------------------------------------------------------------------//

Trend2DStorage::Trend2DStorage(const std::string & file_name,
                               const std::string & reference_parameter1,
                               const std::string & reference_parameter2,
                               const bool        & estimate)
: file_name_(file_name),
  reference_parameter_one_(reference_parameter1),
  reference_parameter_two_(reference_parameter2),
  estimate_(estimate)
{
}

Trend2DStorage::Trend2DStorage(const Trend2DStorage & trend_storage)
: file_name_(trend_storage.file_name_),
  reference_parameter_one_(trend_storage.reference_parameter_one_),
  reference_parameter_two_(trend_storage.reference_parameter_two_),
  estimate_(trend_storage.estimate_)
{
}

Trend2DStorage::~Trend2DStorage()
{
}

Trend *
Trend2DStorage::GenerateTrend(const std::string                       & path,
                              const std::vector<std::string>          & trend_cube_parameters,
                              const std::vector<std::vector<double> > & trend_cube_sampling,
                              const std::vector<std::vector<double> > & blocked_logs,
                              const std::vector<std::vector<double> > & s1,
                              const std::vector<std::vector<double> > & s2,
                              const int                               & type,
                              const NRLib::Trend                      * mean_trend,
                              std::string                             & errTxt) const
{
  Trend * trend = NULL;

  if(estimate_ == true) {
    size_t                            trend_cube_n1 = trend_cube_sampling[0].size();
    size_t                            trend_cube_n2 = trend_cube_sampling[1].size();
    double                            increment1    = trend_cube_sampling[0][1] - trend_cube_sampling[0][0];
    double                            increment2    = trend_cube_sampling[1][1] - trend_cube_sampling[1][0];
    Grid2D<double>                    trend_grid(trend_cube_n1, trend_cube_n2, RMISSING);

    double                            scale        = 1.0;
    double                            bandwidth_x;
    double                            bandwidth_y;
    std::vector<double>               x;
    std::vector<double>               y;
    std::vector<double>               z;

    if (type == TrendStorage::MEAN) {
      PreprocessData2D(s1,
                       s2,
                       blocked_logs,
                       trend_cube_sampling,
                       x,
                       y,
                       z);

      /* -- estimate mean surface -- */
      std::vector<std::vector<double> > trend_surface(trend_cube_n1, std::vector<double>(trend_cube_n2, RMISSING));

      // Note that: n^{-0.167} \approx n^{-1/6}
      bandwidth_x = CalculateBandwidth(x, scale, 0.167);
      bandwidth_y = CalculateBandwidth(y, scale, 0.167);

      bool valid_dataset = CheckConfigurations2D(x,
                                                 y,
                                                 z,
                                                 trend_cube_sampling[0],
                                                 trend_cube_sampling[1],
                                                 bandwidth_x,
                                                 bandwidth_y,
                                                 errTxt);

      if (valid_dataset) {

        EstimateTrend2D(x,
                        y,
                        z,
                        trend_cube_sampling[0],
                        trend_cube_sampling[1],
                        trend_surface,
                        bandwidth_x,
                        bandwidth_y,
                        errTxt);
      }

      for (size_t i = 0; i < trend_cube_n1; i++) {
        for (size_t j = 0; j < trend_cube_n2; j++) {
          trend_grid(i, j) = trend_surface[i][j];
        }
      }
      trend = new Trend2D(trend_grid, 1, 2, increment1, increment2);
    }

    if (type == TrendStorage::VAR) {
      if (mean_trend != NULL) {

        /* -- obtain trend surface -- */
        std::vector<std::vector<double> > trend_cube_trend(trend_cube_n1, std::vector<double>(trend_cube_n2, RMISSING));
        int k_dummy = 0;
        for (int i = 0; i < static_cast<int>(trend_cube_n1); ++i) {
          for (int j = 0; j < static_cast<int>(trend_cube_n2); ++j) {
            trend_cube_trend[i][j] = mean_trend->GetTrendElement(i, j, k_dummy);
          }
        }
        /* -------------------------- */

        PreprocessData2D(s1,
                         s2,
                         blocked_logs,
                         trend_cube_sampling,
                         x,
                         y,
                         z);

        std::vector<std::vector<double> > trend_cube_variance(trend_cube_n1, std::vector<double>(trend_cube_n1, RMISSING));

        bandwidth_x = CalculateBandwidth(x, scale, 0.167);
        bandwidth_y = CalculateBandwidth(y, scale, 0.167);

        bool valid_dataset = CheckConfigurations2D(x,
                                                   y,
                                                   z,
                                                   trend_cube_sampling[0],
                                                   trend_cube_sampling[1],
                                                   bandwidth_x,
                                                   bandwidth_y,
                                                   errTxt);
        if (valid_dataset) {
          EstimateVariance2D(x,
                             y,
                             z,
                             trend_cube_sampling[0],
                             trend_cube_sampling[1],
                             trend_cube_trend,
                             trend_cube_variance,
                             bandwidth_x,
                             bandwidth_y,
                             errTxt);
        }

        for (size_t i = 0; i < trend_cube_n1; i++) {
          for (size_t j = 0; j < trend_cube_n2; j++) {
            trend_grid(i, j) = trend_cube_variance[i][j];
          }
        }
        trend = new Trend2D(trend_grid, 1, 2, increment1, increment2);
      }
      else {
        errTxt += "Error: Invalid trend surface for variance estimation.\n";
      }
    }
  } else {

    int reference1 = 0;
    int reference2 = 0;
    for(int i=0; i<static_cast<int>(trend_cube_parameters.size()); i++) {
      if(reference_parameter_one_ == trend_cube_parameters[i])
        reference1 = i+1;
      else if(reference_parameter_two_ == trend_cube_parameters[i])
        reference2 = i+1;
    }

    if(reference1 == 0 || reference2 == 0)
      errTxt += "The reference parameters of the 2D trend in "+file_name_+" \n"
                " is not the same as the parameter names of the trend cubes\n";

    std::string file_name = path + file_name_;
    RegularSurface<double> surface(file_name);

    RegularSurface<double> resampled_surface;
    bool transpose;
    if(reference1 == 1)
      transpose = false;
    else {
      transpose = true;
      reference1 = 1;
      reference2 = 2;
    }
    resampled_surface = ResampleTrend2D(surface, trend_cube_sampling[0], trend_cube_sampling[1], transpose);

    double increment1 = trend_cube_sampling[0][1] - trend_cube_sampling[0][0];
    double increment2 = trend_cube_sampling[1][1] - trend_cube_sampling[1][0];

    trend = new Trend2D(resampled_surface, reference1, reference2, increment1, increment2);
  }

  return trend;

}

