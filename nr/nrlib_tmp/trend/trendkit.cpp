// $Id: trendkit.cpp 1499 2017-05-29 19:01:53Z perroe $
#include "trendkit.hpp"
#include "../iotools/logkit.hpp"
#include "../iotools/fileio.hpp"
#include "../surface/regularsurface.hpp"
#include <fstream>

namespace NRLib {

void EstimateConstantTrend(const std::vector<std::vector<double> > & blocked_logs,
                           double                                  & trend)
{
  int nWells  = static_cast<int>(blocked_logs.size());

  if(nWells > 0) {
    int nBlocks = static_cast<int>(blocked_logs[0].size());

    trend = 0;

    int count = 0;

    for (int w = 0 ; w < nWells ; w++) {
      if(blocked_logs[w].size() > 0) {
        for (int k = 0 ; k < nBlocks ; k++) {
          if(blocked_logs[w][k] != RMISSING) {
            trend += blocked_logs[w][k];
            count++;
          }
        }
      }
    }


    if (count > 0)
      trend = trend/count;
  }
}
//-------------------------------------------------------------------------------
void PreprocessData0D(const std::vector<std::vector<double> > & blocked_logs,
                      std::vector<double>                     & y)
{
  y.clear();

  if (blocked_logs.size() > 0) {
    for (size_t i = 0; i < blocked_logs.size(); i++) {
      if (blocked_logs[i].size() > 0) {
        for (size_t j = 0; j < blocked_logs[i].size(); j++) {
          double y_j_log = blocked_logs[i][j];
          if (y_j_log != RMISSING) {
            y.push_back(y_j_log);
          }
        }
      }
    }
  }
}
//-------------------------------------------------------------------------------
void ReadTrend1D(const std::string        & file_name,
                      std::string         & errText,
                      std::vector<double> & trend1d,
                      double              & s_min,
                      double              & dz)
{
  std::ifstream file;
  OpenRead(file,file_name);
  std::string dummyStr;
  bool lineIsComment = true;
  int  line          = 0;
  int  thisLine      = 0;

  while( lineIsComment == true) {

    if(CheckEndOfFile(file)) {
      errText += "Error: End of file "+file_name+" premature.\n";
      return;
    }

    ReadNextToken(file,dummyStr,line);
    if (line == thisLine)
      DiscardRestOfLine(file,line,false);
    thisLine = line;
    if((dummyStr[0]!='*') &  (dummyStr[0]!='"')) {
      lineIsComment = false;
    }
  }

  s_min = ParseType<double>(dummyStr);

  if (CheckEndOfFile(file))  {
    errText += "Error: End of file "+file_name+" premature.\n";
    return;
  }
  ReadNextToken(file,dummyStr,line);
  if (line == thisLine)
    DiscardRestOfLine(file,line,false);
  thisLine = line;

  dz = ParseType<double>(dummyStr);

  if (CheckEndOfFile(file)) {
    errText += "Error: End of file "+file_name+" premature.\n";
    return;
  }
  ReadNextToken(file,dummyStr,line);
  if (line == thisLine)
    DiscardRestOfLine(file,line,false);

  int nz = ParseType<int>(dummyStr);

  trend1d.resize(nz);

  for(int i=0; i<nz; i++) {
    if (CheckEndOfFile(file)) {
      errText += "Error: End of file "+file_name+" premature.\n";
      return;
    }
    ReadNextToken(file,dummyStr,line);

    trend1d[i] = ParseType<double>(dummyStr);
  }
  file.close();
}

//----------------------------------------------------//

void WriteTrend1D(const std::string         & filename,
                  const std::vector<double> & s,
                  const std::vector<double> & trend)
{
  std::ofstream file;
  OpenWrite(file, filename);

  file << std::fixed
       << std::setprecision(6)
       << s.front()   << "\n"
       << s[1] - s[0] << "\n"
       << s.back()    << "\n";

  for (size_t i = 0; i < trend.size(); i++) {
    file << trend[i] << "\n";
  }
  file.close();
}

//----------------------------------------------------//


void ReadTrend1DPlainAscii(const std::string   & file_name,
                           std::string         & /*errText*/,
                           std::vector<double> & trend1d)
{
  std::ifstream file;
  OpenRead(file,file_name);
  std::string dummyStr;
  int line = 0;

  while(CheckEndOfFile(file)==false){
    ReadNextToken(file,dummyStr,line);
    trend1d.push_back(ParseType<double>(dummyStr));
  }
  file.close();

}

//----------------------------------------------------//

int
GetTrend1DFileFormat(const std::string & file_name,
                     std::string       & errText)
{
  std::string   dummyStr;
  std::ifstream file;

  // test for jason file format
  OpenRead(file,file_name);

  int  fileformat    = -1;
  int  line          = 0;
  int  thisLine      = 0;
  bool lineIsComment = true;
  bool commentFound  = false;

  while (lineIsComment == true) {
    ReadNextToken(file,dummyStr,line);
    if (CheckEndOfFile(file)) {
      errText += "End of trend file "+file_name+" is premature\n";
      return 0;
    }
    else {
      if (thisLine == line) {
        DiscardRestOfLine(file,line,false);
        thisLine = line;
      }
      if((dummyStr[0]!='*') &  (dummyStr[0]!='"'))
        lineIsComment = false;
      else
        commentFound = true;
    }
  }
  file.close();
  if (IsNumber(dummyStr) && commentFound == true) // not convertable number
    fileformat = 0; //Same file format as WAVELET::JASON
  else
    fileformat = 1; // plain ascii file
  return fileformat;
}

//----------------------------------------------------//

void InterpolateTrend1DValues(const double & xi,
                              const double & xi1,
                              const double & fxi,
                              const double & fxi1,
                              const double & yj,
                              double       & fyj)
{
  double t;

  t   = (yj-xi1)/(xi-xi1);

  fyj = fxi*t + fxi1*(1-t);

}

//----------------------------------------------------//

void ResampleTrend1D(const std::vector<double> & x,
                     const std::vector<double> & fx,
                     const std::vector<double> & y,
                     std::vector<double>       & fy)
{
  // Resample fx with sampling x into fy with sampling y

  int nx = static_cast<int>(x.size());
  int ny = static_cast<int>(y.size());

  int i=0;
  int j=0;

  while(i<nx-1) {

    while(j<ny && y[j]>=x[i] && y[j]<x[i+1]) {

      InterpolateTrend1DValues(x[i],
                               x[i+1],
                               fx[i],
                               fx[i+1],
                               y[j],
                               fy[j]);
      j++;
    }

    i++;
  }
}
//-------------------------------------------------------------------------------
void EstimateTrend1D(const std::vector<double> & x,
                     const std::vector<double> & y,
                     const std::vector<double> & x0,
                     std::vector<double>       & y0,
                     double                      bandwidth,
                     std::string               & errTxt,
                     const std::vector<double> * weights)
{
  std::vector<double> x_binned     = x;
  std::vector<double> y_binned     = y;
  std::vector<double> x0_regridded = x0;
  std::vector<double> y0_regridded = y0;

  std::vector<double> w_binned;

  if (weights == NULL)
    w_binned.resize(x_binned.size(), 1.0);
  else
    w_binned = *weights;

  // if the number of sample points is too large, we regrid to save computation time.
  if (x0.size() > 128) {
    MakeNewGridResolution(bandwidth, x0, x0_regridded);
  }

  y0_regridded.resize(x0_regridded.size());
  for (size_t i = 0; i < y0_regridded.size(); i++)
    y0_regridded[i] = RMISSING;

  // if the sample size is large (> 128^2 = 16384), linear binning is used to save computation time.
  if (x.size() > 16384) {
    MakeBinnedDataset1D(x,
                        y,
                        weights,
                        x0_regridded,
                        y_binned,
                        w_binned);

    x_binned = x0_regridded;
  }

  // finds the effective sample size (in ''the best of all worlds'') and use this to
  // make an adaptive bandwidth for the Kernel estimator.
  double effective_sample_size = CalculateEffectiveSampleSize1D(x_binned, bandwidth, w_binned);

  double w_binned_size = 0;
  for (size_t i = 0; i < w_binned.size(); i++)
    w_binned_size = w_binned_size + w_binned[i];
  // the effective weighted sample size has to be less than the actual weighted sample size
  if (w_binned_size < effective_sample_size) {
    effective_sample_size = w_binned_size;
    LogKit::LogFormatted(LogKit::High,"\nWARNING : The effective sample size used in estimating 1D trend is larger than the actual sample size.");
    LogKit::LogFormatted(LogKit::High,"\n          The local linear model is reduced to a standard linear model.");
  }

  // local linear regression with varying bandwidth
  bool   complete_line = false;
  size_t nr_iterations = 0;
  while (!complete_line) {
    LocalLinearRegression1D(x_binned,
                            y_binned,
                            w_binned,
                            bandwidth,
                            effective_sample_size,
                            x0_regridded,
                            y0_regridded,
                            complete_line);

    bandwidth     = (1.00 + 0.005)*bandwidth;
    nr_iterations = nr_iterations + 1;

    if (nr_iterations > 1e4) {
      errTxt +=  "Error: Unable to estimate complete 1D trend.\n";
      errTxt +=  "       The interpolated region is too large compared to the support of the data.\n";
      errTxt +=  "       Reduce the size of the estimated region or estimate a constant value.\n";
      break;
    }
  }

  // linear interpolation back to the original scale
  if (x0.size() != x0_regridded.size()) {
    LinearInterpolation(x0_regridded, y0_regridded, x0, y0);
  } else {
    y0 = y0_regridded;
  }
}
//-------------------------------------------------------------------------------
void PreprocessData1D(const std::vector<std::vector<double> > & s,
                      const std::vector<std::vector<double> > & blocked_logs,
                      const std::vector<double>               & trend_cube_sampling,
                      std::vector<double>                     & x,
                      std::vector<double>                     & y)
{
  size_t n_wells   = blocked_logs.size();

  x.clear();
  y.clear();

  // Note that s1 and s2 are corrected/centered according to trend_cube_sampling,
  // we have to add x_min or y_min to obtain the original s1 and s2 scale.
  double x_min = trend_cube_sampling[0];

  if (n_wells > 0) {
    for (size_t i = 0; i < n_wells; i++) {
      size_t n_samples_i = blocked_logs[i].size();

      if (n_samples_i > 0) {
        for (size_t j = 0; j < n_samples_i; j++) {
          double x_j     = s[i][j];
          double y_j_log = blocked_logs[i][j];

          if (x_j != RMISSING && y_j_log != RMISSING) {
            x.push_back(x_j + x_min);
            y.push_back(y_j_log);
          }
        }
      }
    }
  }
}
//-------------------------------------------------------------------------------
bool CheckConfigurations1D(const std::vector<double> & x,
                           const std::vector<double> & y,
                           const std::vector<double> & x0,
                           const double              & bandwidth,
                           std::string               & errTxt)
{
  bool   valid_dataset = true;

  double upp       = 3.0;
  double large     = 1.0;

  size_t n_samples = y.size();

  double x_min     = *std::min_element(x.begin(), x.end());
  double x_max     = *std::max_element(x.begin(), x.end());
  double x_delta   = std::abs(x_max - x_min);
  double x_upp     = upp*x_delta;
  double x_large   = large*x_delta;
  double x0_min    = *std::min_element(x0.begin(), x0.end());
  double x0_max    = *std::max_element(x0.begin(), x0.end());

  if (!CheckIfVectorIsSorted(x0)) {
    errTxt        += "Error: Unable to estimate 1D trend/variance.\n";
    errTxt        += "       The inputs vector is not sorted. \n";
    valid_dataset  = false;
  }
  else {
    if (x0_min < x_min - x_upp && x_max + x_upp < x0_max) {
      errTxt        += "Error: Unable to estimate 1D trend/variance.\n";
      errTxt        += "       The interpolated region is too large compared to the support of the data.\n";
      errTxt        += "       Reduce the size of the estimated region or estimate a constant value.\n";
      valid_dataset  = false;
    }
    else {
      if (x0_min < x_min - x_large || x_max + x_large < x0_max) {
        LogKit::LogFormatted(LogKit::High,"\nWARNING : The defined region for 1D trend/variance is large compared to the support of the data.");
        LogKit::LogFormatted(LogKit::High,"\n          This can result in unstable estimates. Consider using an alternative (low-dimensional) method.\n");
      }

      else {
        if (n_samples < 10) {
          errTxt        += "Error: Unable to estimate 1D trend/variance.\n";
          errTxt        += "       Too few observations.\n";
          errTxt        += "       The algorithm requires at least 11 observations to run.\n";
          valid_dataset  = false;
        }
        else {
          if (bandwidth < 1e-5) {
            errTxt        += "Error: Unable to estimate 1D trend/variance.\n";
            errTxt        += "       The spread in the data is too low to provide a valid bandwidth.\n ";
            valid_dataset  = false;
          }
          else {
            if (n_samples < 32) {
              LogKit::LogFormatted(LogKit::High,"\nWARNING : The sample size is relatively small for 1D trend/variance estiamtion.");
              LogKit::LogFormatted(LogKit::High,"\n          This can result in unstable estimates. Consider using an alternative (low-dimensional) method.\n");
            }
          }
        }
      }
    }
  }
  return(valid_dataset);
}
//-------------------------------------------------------------------------------
void MakeBinnedDataset1D(const std::vector<double>         & x,
                         const std::vector<double>         & y,
                         const std::vector<double>         * w,
                         const std::vector<double>         & x_new,
                         std::vector<double>               & y_new,
                         std::vector<double>               & w_new)
{
  double epsilon            = 1e-5;

  size_t x_new_n            = x_new.size();

  double x_new_low          = x_new.front();
  double x_new_upp          = x_new.back();
  double x_new_inc          = x_new[1] - x_new[0];
  double x_new_one_over_inc = 1/x_new_inc;

  //We also bin the weights if they are given as input, else we return the local weights
  bool w_input = w != NULL;

  y_new.clear();
  w_new.clear();
  std::vector<double> w_local;

  y_new.resize(x_new_n, 0.0);
  w_new.resize(x_new_n, 0.0);
  w_local.resize(x_new_n, 0.0);

  size_t j;
  double weight;
  for (size_t i = 0; i < y.size(); i++) {
    if (x[i] < x_new_low + epsilon) {
      j        = 0;
      weight   = std::max(1.0 - std::abs(x[i] - x_new[j])*x_new_one_over_inc, 0.0);

      w_local[j] = w_local[j] + weight;
      y_new[j]   = y_new[j] + weight*y[i];
      if (w_input)
        w_new[j] = w_new[j] + weight * (*w)[i];
    }
    else if (x[i] > x_new_upp - epsilon) {
      j        = x_new_n - 1;
      weight   = std::max(1.0 - std::abs(x[i] - x_new[j])*x_new_one_over_inc, 0.0);

      w_local[j] = w_local[j] + weight;
      y_new[j]   = y_new[j] + weight*y[i];
      if (w_input)
        w_new[j] = w_new[j] + weight * (*w)[i];
    }
    else {
      j            = static_cast<size_t>(std::floor((x[i] - x_new_low)*x_new_one_over_inc));
      weight       = std::max(1.0 - std::abs(x[i] - x_new[j])*x_new_one_over_inc, 0.0);

      w_local[j]   = w_local[j] + weight;
      y_new[j]     = y_new[j] + weight*y[i];
      if (w_input)
        w_new[j]   = w_new[j] + weight*(*w)[i];

      w_local[j + 1] = w_local[j + 1] + (1 - weight);
      y_new[j + 1]   = y_new[j + 1] + (1 - weight)*y[i];
      if (w_input)
        w_new[j + 1] = w_new[j + 1] + (1 - weight)*(*w)[i];
    }
  }
  for (size_t i = 0; i < y_new.size(); i++) {
    if (w_local[i] != 0.0) {
      y_new[i] = y_new[i]/w_local[i];
      if (w_input)
        w_new[i] = w_new[i]/w_local[i];
    }
    else {
      y_new[i] = RMISSING;
      w_new[i] = RMISSING;
    }
  }

  if (w_input == false)
    w_new = w_local;

}
//-------------------------------------------------------------------------------
void LocalLinearRegression1D(const std::vector<double> & x,
                             const std::vector<double> & y,
                             const std::vector<double> & w,
                             const double              & bandwidth,
                             const double              & effective_sample_size,
                             const std::vector<double> & x0,
                             std::vector<double>       & y0,
                             bool                      & complete_line)
{
  size_t              n_y                       = y.size();

  double              threshold_times_bandwidth = (std::pow(std::log(1e5), 0.5))*bandwidth;
  double              one_over_bandwidth_sq     = 1.0/(bandwidth*bandwidth);
  double              delta_limit               = 2*std::log(1e5);

  std::vector<double> x_order(n_y);
  std::vector<size_t> x_backward(n_y);
  std::vector<size_t> x_forward(n_y);
  SortOrderAndRank(x, x_order, x_backward, x_forward);

  std::vector<double> weights;
  weights.reserve(n_y);

  std::vector<size_t> index;
  index.reserve(n_y);

  complete_line = true;
  for (size_t i = 0; i < x0.size(); i++) {
    if (y0[i] == RMISSING) {
      double weight_i = 0.0;
      double x_low    = x0[i] - threshold_times_bandwidth;
      double x_upp    = x0[i] + threshold_times_bandwidth;

      size_t x_low_j  = FindLowerBoundInSortedVector(x_order, x_low);
      size_t x_upp_j  = FindUpperBoundInSortedVector(x_order, x_upp);

      for (size_t j = x_low_j; j < x_upp_j + 1; j++) {
        size_t k       = x_backward[j];
        double delta_j = (x0[i] - x[k])*(x0[i] - x[k])*one_over_bandwidth_sq;

        if (delta_j < delta_limit && y[k] != RMISSING) {
          double weight_j = w[k]*std::exp(-0.5*delta_j);
          weight_i        = weight_i + weight_j;

          weights.push_back(weight_j);
          index.push_back(k);
        }
      }

      if (weight_i > effective_sample_size) {
        // calcualte (X^{t}WX) = A and X^{t}Wy = c
        double a00 = 0.0;
        double a01 = 0.0;
        double a11 = 0.0;
        double c0  = 0.0;
        double c1  = 0.0;

        for (size_t j = 0; j < index.size(); j++) {
          size_t k      = index[j];
          double weight = weights[j];

          a00 = a00 +                weight;
          a01 = a01 + x[k]          *weight;
          a11 = a11 + x[k]*x[k]     *weight;
          c0  = c0  +           y[k]*weight;
          c1  = c1  + x[k]     *y[k]*weight;
        }

        // add a small constant epsilon = 1e-4 to make sure that A is invertible
        //a00                   = a00*(1 + 1e-4);
        a11                   = a11*(1 + 1e-4);

        // calcualte B = A^{-1}
        double b00            =  a11;
        double b01            = -a01;
        double b11            =  a00;
        double det_A          =  a00*a11 - a01*a01;

        // if det_A is to small, the esitamted curve becomes unstable
        if (det_A > 1e-6) {
          // calculate x_{0}A^{-1}c = x0(X^{t}WX)^{-1}X^{t}y = \hat{y}_{0}
          double tmp = 0.0;
          tmp        = tmp +       (b00*c0 + b01*c1)/det_A;
          tmp        = tmp + x0[i]*(b01*c0 + b11*c1)/det_A;
          y0[i]      = tmp;
        } else {
          y0[i]         = RMISSING;
          complete_line = false;
        }
      } else {
        y0[i]         = RMISSING;
        complete_line = false;
      }
      index.clear();
      weights.clear();
    }
  }
}
//-------------------------------------------------------------------------------
void LinearInterpolation(const std::vector<double> & x,
                         const std::vector<double> & y,
                         const std::vector<double> & x0,
                         std::vector<double>       & y0)
{
  // Note: x is assumed to be a sorted vector
  double delta          = 1e-5;
  double x_low          = x.front();
  double x_upp          = x.back();
  double x_inc          = x[1] - x[0];
  double x_one_over_inc = 1/x_inc;

  y0.resize(x0.size(), RMISSING);

  for (size_t i = 0; i < x0.size(); i++) {
    if (x_low - delta < x0[i] && x0[i] < x_upp + delta) {
      if (x0[i] < x_low + delta) {
        y0[i] = y.front();
      }
      else if (x0[i] > x_upp - delta) {
        y0[i] = y.back();
      }
      else {
        size_t j      = static_cast<size_t>(std::floor((x0[i] - x_low)*x_one_over_inc));
        double weight = std::max(1.0 - std::abs(x0[i] - x[j])*x_one_over_inc, 0.0);
        y0[i]         = weight*y[j] + (1 - weight)*y[j + 1];
      }
    } else {
      y0[i] = RMISSING;
    }
  }
}
//-------------------------------------------------------------------------------
void EstimateVariance1D(const std::vector<double> & x,
                        const std::vector<double> & y,
                        const std::vector<double> & x0,
                        const std::vector<double> & trend,
                        std::vector<double>       & variance,
                        double                      bandwidth,
                        std::string               & errTxt,
                        const std::vector<double> * weights)
{
  std::vector<double> x_binned;
  std::vector<double> w_binned;
  std::vector<double> x0_regridded = x0;
  std::vector<double> w0_regridded(x.size(), RMISSING);
  std::vector<double> variance_regridded(x.size(), RMISSING);

  /* -- estiamte global variance -- */
  std::vector<double> y_y_mean_squared_binned;
  std::vector<double> y_y_mean_squared;
  std::vector<double> y_mean;

  LinearInterpolation(x0, trend, x, y_mean);

  double              y_var = 0.0;
  for (size_t i = 0; i < y.size(); i++) {
    if (y_mean[i] != RMISSING && y[i] != RMISSING) {
      double y_y_mean_squared_i = (y[i] - y_mean[i])*(y[i] - y_mean[i]);
      y_y_mean_squared.push_back(y_y_mean_squared_i);
      x_binned.push_back(x[i]);

      if (weights == NULL)
        w_binned.push_back(1.0);
      else
        w_binned.push_back((*weights)[i]);

      y_var = y_var + y_y_mean_squared_i;
    }
  }

  if (y_y_mean_squared.size() > 0) {
    y_var = y_var/(y_y_mean_squared.size() - 1);
  }
  else {
    y_var = RMISSING;
  }
  /* -------------------------------*/

  if (y_var == RMISSING) {
    errTxt += "Error: Unable to estimate 1D variance.\n";
    errTxt += "       The data contains only missing values.\n";
  }
  else {

    // resampling
    if (x0.size() > 128) {
      MakeNewGridResolution(bandwidth, x0, x0_regridded);
    }

    // linear binning
    y_y_mean_squared_binned = y_y_mean_squared;
    if (x.size() > 16384) {
      MakeBinnedDataset1D(x,
                          y_y_mean_squared,
                          weights,
                          x0_regridded,
                          y_y_mean_squared_binned,
                          w_binned);
     x_binned = x0_regridded;
    }

    variance_regridded.resize(x0_regridded.size(), RMISSING);
    w0_regridded.resize(x0_regridded.size(), RMISSING);

    double y_var_weight = CalculateEffectiveSampleSize1D(x_binned, bandwidth, w_binned);

    double w_binned_size = 0;
    for (size_t i = 0; i < w_binned.size(); i++)
      w_binned_size = w_binned_size + w_binned[i];
    if (w_binned_size < y_var_weight) {
      y_var_weight = w_binned_size;
      LogKit::LogFormatted(LogKit::High,"\nWARNING : The effective sample size used in estimating 1D variance is larger than the actual sample size.");
      LogKit::LogFormatted(LogKit::High,"\n          The local linear model is reduced to a standard linear model.");
    }

    /* -- estiamte weighted variance -- */
    bool complete_line_dummy = false;
    KernelSmoother1DLine(x_binned,
                         y_y_mean_squared_binned,
                         w_binned,
                         bandwidth,
                         0.0,
                         x0_regridded,
                         variance_regridded,
                         w0_regridded,
                         complete_line_dummy);
    for (size_t i = 0; i < variance_regridded.size(); i++) {
      if (variance_regridded[i] != RMISSING) {
        double weight_i     = w0_regridded[i]/(w0_regridded[i] + y_var_weight);
        variance_regridded[i] = weight_i*variance_regridded[i] + (1 - weight_i)*y_var;
      }
      else {
        variance_regridded[i] = y_var;
      }
    }
    /* ---------------------------------*/

    /* -- linear interpolation -- */
    if (x0.size() != x0_regridded.size()) {
      LinearInterpolation(x0_regridded, variance_regridded, x0, variance);
    }
    else {
      variance = variance_regridded;
    }
    /* ---------------------------*/
  }
}
//-------------------------------------------------------------------------------
void KernelSmoother1DLine(const std::vector<double> & x,
                          const std::vector<double> & y,
                          const std::vector<double> & w,
                          const double              & bandwidth,
                          const double              & effective_sample_size,
                          const std::vector<double> & x0,
                          std::vector<double>       & y0,
                          std::vector<double>       & w0,
                          bool                      & complete_line)
{
  size_t               n_y                       = y.size();

  double               threshold_times_bandwidth = (std::pow(std::log(1e5), 0.5))*bandwidth;
  double               one_over_bandwidth_sq     = 1.0/(bandwidth*bandwidth);
  double               delta_limit               = 2*std::log(1e5);

  std::vector<double> x_order(n_y);
  std::vector<size_t> x_backward(n_y);
  std::vector<size_t> x_forward(n_y);
  SortOrderAndRank(x, x_order, x_backward, x_forward);

  complete_line = true;
  for (size_t i = 0; i < x0.size(); i++) {
    if (y0[i] == RMISSING) {
      double y0_hat   = 0.0;
      double weight_i = 0.0;
      double x_low    = x0[i] - threshold_times_bandwidth;
      double x_upp    = x0[i] + threshold_times_bandwidth;

      size_t x_low_j  = FindLowerBoundInSortedVector(x_order, x_low);
      size_t x_upp_j  = FindUpperBoundInSortedVector(x_order, x_upp);

      for (size_t j = x_low_j; j < x_upp_j + 1; j++) {
        size_t k       = x_backward[j];
        double delta_j = (x0[i] - x[k])*(x0[i] - x[k])*one_over_bandwidth_sq;

        if (delta_j < delta_limit && y[k] != RMISSING) {
          double weight = w[k]*std::exp(-0.5*delta_j);
          y0_hat        = y0_hat + weight*y[k];
          weight_i      = weight_i + weight;
        }
      }

      if (weight_i > effective_sample_size) {
        y0[i] = y0_hat/weight_i;
        w0[i] = weight_i;
      } else {
        y0[i]         = RMISSING;
        w0[i]         = 0.0;
        complete_line = false;
      }
    }
  }
}
//-------------------------------------------------------------------------------
void EstimateTrend2D(const std::vector<double>         & x,
                     const std::vector<double>         & y,
                     const std::vector<double>         & z,
                     const std::vector<double>         & x0,
                     const std::vector<double>         & y0,
                     std::vector<std::vector<double> > & z0,
                     double                              bandwidth_x,
                     double                              bandwidth_y,
                     std::string                       & errTxt)
{
  std::vector<double> x_binned     = x;
  std::vector<double> y_binned     = y;
  std::vector<double> z_binned     = z;
  std::vector<double> w_binned(x_binned.size(), 1.0);
  std::vector<double> x0_regridded = x0;
  std::vector<double> y0_regridded = y0;

  MakeNewGridResolution(bandwidth_x, x0, x0_regridded);
  MakeNewGridResolution(bandwidth_y, y0, y0_regridded);

  std::vector<std::vector<double> > z0_regridded(x0_regridded.size(), std::vector<double>(y0_regridded.size(), RMISSING));
  std::vector<std::vector<double> > w0_regridded_dummy(x0_regridded.size(), std::vector<double>(y0_regridded.size(), RMISSING));

  if (z_binned.size() > 16384) {
    MakeBinnedDataset2D(x0_regridded, y0_regridded, x_binned, y_binned, z_binned, w_binned);
  }

  double effective_sample_size = CalculateEffectiveSampleSize2D(x_binned, y_binned, bandwidth_x, bandwidth_y);

  double z_binned_size = static_cast<double>(z_binned.size());
  if (z_binned_size < effective_sample_size) {
    effective_sample_size = z_binned_size;
    LogKit::LogFormatted(LogKit::High,"\nWARNING : The effective sample size used in estimating 2D trend is larger than the actual sample size.");
    LogKit::LogFormatted(LogKit::High,"\n          The local linear is reduced to a standard linear model.");
  }

  bool   complete_surface = false;

  size_t l = 0;

  while (l < 1e4) {
    LocalLinearRegression2DSurface(x_binned,
                                   y_binned,
                                   z_binned,
                                   w_binned,
                                   bandwidth_x,
                                   bandwidth_y,
                                   effective_sample_size,
                                   x0_regridded,
                                   y0_regridded,
                                   z0_regridded,
                                   w0_regridded_dummy,
                                   complete_surface);

    bandwidth_x = (1 + 0.02)*bandwidth_x;
    bandwidth_y = (1 + 0.02)*bandwidth_y;

    if (complete_surface || l > 1e4) {
      break;
    }
    l = l + 1;
  }
  if (l > 1e4) {
    errTxt +=  "Error: Unable to compute complete 2D trend surface.\n";
    errTxt +=  "       The interpolated region is too large compared to the support of the data.\n";
    errTxt +=  "       Reduce the size of the estiamted region or use a low-dimensional method.\n";
  }
  BilinearInterpolation(x0_regridded, y0_regridded, z0_regridded, x0, y0, z0);
}
//-------------------------------------------------------------------------------
void PreprocessData2D(const std::vector<std::vector<double> > & s1,
                      const std::vector<std::vector<double> > & s2,
                      const std::vector<std::vector<double> > & blocked_logs,
                      const std::vector<std::vector<double> > & trend_cube_sampling,
                      std::vector<double>                     & x,
                      std::vector<double>                     & y,
                      std::vector<double>                     & z)
{
  size_t n_wells   = blocked_logs.size();

  x.clear();
  y.clear();
  z.clear();

  // Note that s1 and s2 are corrected/centered according to trend_cube_sampling,
  // we have to add x_min and y_min to obtain the original s1 and s2 scale
  double x_min = trend_cube_sampling[0][0];
  double y_min = trend_cube_sampling[1][0];

  if (n_wells > 0) {
    for (size_t i = 0; i < n_wells; i++) {
      size_t n_samples_i = blocked_logs[i].size();

      if (n_samples_i > 0) {
        for (size_t j = 0; j < n_samples_i; j++) {
          double x_j     = s1[i][j];
          double y_j     = s2[i][j];
          double z_j_log = blocked_logs[i][j];

          if (x_j != RMISSING && y_j != RMISSING && z_j_log != RMISSING) {
            x.push_back(x_j + x_min);
            y.push_back(y_j + y_min);
            z.push_back(z_j_log);
          }
        }
      }
    }
  }
}
//-------------------------------------------------------------------------------
bool CheckConfigurations2D(const std::vector<double> & x,
                           const std::vector<double> & y,
                           const std::vector<double> & z,
                           const std::vector<double> & x0,
                           const std::vector<double> & y0,
                           double                      bandwidth_x,
                           double                      bandwidth_y,
                           std::string               & errTxt)
{
  bool                valid_dataset = true;

  double              upp           = 3.0;
  double              large         = 1.0;

  size_t              n_samples     = z.size();

  double              x_min         = *std::min_element(x.begin(), x.end());
  double              x_max         = *std::max_element(x.begin(), x.end());
  double              x_delta       = std::abs(x_max - x_min);
  double              x_upp         = upp*x_delta;
  double              x_large       = large*x_delta;

  double              y_min         = *std::min_element(y.begin(), y.end());
  double              y_max         = *std::max_element(y.begin(), y.end());
  double              y_delta       = std::abs(y_max - y_min);
  double              y_upp         = upp*y_delta;
  double              y_large       = large*y_delta;

  if (!CheckIfVectorIsSorted(x0) || !CheckIfVectorIsSorted(y0)) {
    errTxt        += "Error: Unable to estimate 2D trend/variance.\n";
    errTxt        += "       At least one of the inputs are not sorted.\n";
    valid_dataset  = false;
  }
  else {
    if (x0.front() < x_min - x_upp &&  x_max + x_upp < x0.back() && y0.front() < y_min - y_upp &&  y_max + y_upp < y0.back()) {
      errTxt        += "Error: Unable to estimate 2D trend/variance.\n";
      errTxt        += "       The interpolated region is too large compared to the support of the data.\n";
      errTxt        += "       Reduce the size of the estiamted region or use a low dimensional method.\n";
      valid_dataset  = false;
    }
    else {
      if (x0.front() < x_min - x_large || x0.back() > x_max + x_large || y0.front() < y_min - y_large || y0.back() > y_max + y_large) {
        LogKit::LogFormatted(LogKit::High,"\nWARNING : The defined region for 2D trend/variance is large compared to the support of the data.");
        LogKit::LogFormatted(LogKit::High,"\n          This can result in unstable estimates. Consider using an alternative (low-dimensional) method.\n");
      }

      else {
        if (n_samples < 10) {
          errTxt        += "Error: Unable to estimate 2D trend/variance.\n";
          errTxt        += "       Too few observations. The algorithm requires at least 11 observations to run.\n";
          valid_dataset  = false;
        } else {
          if (bandwidth_x < 1e-5 || bandwidth_y < 1e-5) {
            errTxt        += "Error: Unable to estimate 2D trend/variance.\n";
            errTxt        += "       The spread in the data is too low to provide a valid bandwidth.\n ";
            valid_dataset  = false;
          }
          else {
            if (n_samples < 48) {
              LogKit::LogFormatted(LogKit::High,"\nWARNING : The sample size is relatively small for 2D trend/variance estiamtion.");
              LogKit::LogFormatted(LogKit::High,"\n          This can result in unstable estimates. Consider using an alternative (low-dimensional) method.\n");
            }
          }
        }
      }
    }
  }
  return(valid_dataset);
}

//-------------------------------------------------------------------------------
bool CheckIfVectorIsSorted(const std::vector<double> & x) {

  bool is_sorted = true;

  for (size_t i = 0; i < x.size() - 1; i++) {
    if (x[i] > x[i + 1]) {
      is_sorted = false;
    }
  }
  return(is_sorted);
}
//-------------------------------------------------------------------------------
double CalculateBandwidth(const std::vector<double> & x,
                          const double              & scale,
                          const double              & power)
{
  double n               = static_cast<double>(x.size());

  double one_over_root_2 = 0.7071068;
  double sd_x            = std::pow(CalculateVariance(x), 0.5);

  return(scale*one_over_root_2*sd_x*std::pow(static_cast<double>(n), -power));
}
//-------------------------------------------------------------------------------
double CalculateEffectiveSampleSize1D(const std::vector<double> & x,
                                      const double              & bandwidth,
                                      const std::vector<double> & weights)
{
  double w_sum = 0;
  for (size_t i = 0; i < weights.size(); i++)
    w_sum += weights[i];

  double n     = w_sum;
  double x_min = *std::min_element(x.begin(), x.end());
  double x_max = *std::max_element(x.begin(), x.end());

  return(2.506628*bandwidth*n/(x_max - x_min));
}
//-------------------------------------------------------------------------------
double CalculateEffectiveSampleSize2D(const std::vector<double> & x,
                                      const std::vector<double> & y,
                                      const double              & bandwidth_x,
                                      const double              & bandwidth_y)
{
  size_t n     = x.size();

  double x_min = *std::min_element(x.begin(), x.end());
  double x_max = *std::max_element(x.begin(), x.end());

  double y_min = *std::min_element(y.begin(), y.end());
  double y_max = *std::max_element(y.begin(), y.end());

  return(6.283185*bandwidth_x*bandwidth_y*n/((x_max - x_min)*(y_max - y_min)));
}
//-------------------------------------------------------------------------------
double CalculateVariance(const std::vector<double> & x)
{
  int n         = 0;

  double sum_x  = 0.0;
  double sum_x2 = 0.0;

  for (size_t i = 0; i < x.size(); i++) {
    if (x[i] != RMISSING) {
      sum_x  = sum_x  + x[i];
      sum_x2 = sum_x2 + x[i]*x[i];
      n      = n + 1;
    }
  }
  double var_x  = sum_x2/(n - 1) - sum_x*sum_x/(n*(n - 1));

  return(var_x);
}
//-------------------------------------------------------------------------------
void MakeNewGridResolution(const double              & bandwidth_x,
                           const std::vector<double> & x0,
                           std::vector<double>       & x0_regridded)
{
  double              x0_min           = *std::min_element(x0.begin(), x0.end());
  double              x0_max           = *std::max_element(x0.begin(), x0.end());

  size_t              n_max            = 128;
  size_t              x0_n_regridded   = std::min(static_cast<size_t>(std::ceil(7*(x0_max - x0_min)/(2*bandwidth_x))), n_max);

  double              x0_regridded_inc = (x0_max - x0_min)/(x0_n_regridded - 1);
  x0_regridded.resize(x0_n_regridded);

  for (size_t i = 0; i < x0_n_regridded; i++) {
    x0_regridded[i] = x0_min + x0_regridded_inc*i;
  }
}
//-------------------------------------------------------------------------------
void MakeBinnedDataset2D(const std::vector<double>         & x0,
                         const std::vector<double>         & y0,
                         std::vector<double>               & x,
                         std::vector<double>               & y,
                         std::vector<double>               & z,
                         std::vector<double>               & w)
{
  double delta = 1e-5;

  size_t n     = x0.size();
  size_t m     = y0.size();

  std::vector<std::vector<double> > z0(n, std::vector<double>(m, 0.0));
  std::vector<std::vector<double> > w0(n, std::vector<double>(m, 0.0));

  BilinearBinning(x, y, z, x0, y0, z0, w0);

  std::vector<double> x_new(n*m, RMISSING);
  std::vector<double> y_new(n*m, RMISSING);
  std::vector<double> z_new(n*m, RMISSING);
  std::vector<double> w_new(n*m, RMISSING);

  size_t k = 0;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      if (w0[i][j] > 0.0 + delta) {
        x_new[k] = x0[i];
        y_new[k] = y0[j];
        z_new[k] = z0[i][j]/w0[i][j];
        w_new[k] = w0[i][j];

        k = k + 1;
      }
    }
  }
  if (k > 0) {
    x_new.resize(k);
    y_new.resize(k);
    z_new.resize(k);
    w_new.resize(k);
  } else {
    x_new.resize(1, RMISSING);
    y_new.resize(1, RMISSING);
    z_new.resize(1, RMISSING);
    w_new.resize(1, RMISSING);
  }
  x = x_new;
  y = y_new;
  z = z_new;
  w = w_new;
}
//-------------------------------------------------------------------------------
void BilinearBinning(const std::vector<double>               & x,
                     const std::vector<double>               & y,
                     const std::vector<double>               & z,
                     const std::vector<double>               & x0,
                     const std::vector<double>               & y0,
                     std::vector<std::vector<double> >       & z0,
                     std::vector<std::vector<double> >       & w0)
{
  double delta    = 1e-5;

  double x0_first = x0.front();
  double x0_last  = x0.back();
  double x0_inc   = x0[1] - x0[0];
  double x0_scale = 1/x0_inc;

  double y0_first = y0.front();
  double y0_last  = y0.back();
  double y0_inc   = y0[1] - y0[0];
  double y0_scale = 1/y0_inc;

  for (size_t i = 0; i < z.size(); i++) {

    if ((x0_first - delta < x[i]) && (x[i] < x0_last + delta) && (y0_first - delta < y[i]) && (y[i] < y0_last + delta)) {

      size_t k;
      size_t l;
      double w_y;
      double w_x;

      if (x[i] < x0_first + delta) {
        k   = 0;
        w_x = 0.0;
      } else if (x[i] > x0_last - delta) {
        k   = x0.size() - 2;
        w_x = 1.0;
      } else {
        k   = static_cast<size_t>(std::floor((x[i] - x0_first)*x0_scale));
        w_x = (x[i] - x0_first - k*x0_inc)*x0_scale;
      }

      if (y[i] < y0_first + delta) {
        l   = 0;
        w_y = 0.0;
      } else if (y[i] > y0_last - delta) {
        l   = y0.size() - 2;
        w_y = 1.0;
      } else {
        l   = static_cast<size_t>(std::floor((y[i] - y0_first)*y0_scale));
        w_y = (y[i] - y0_first - l*y0_inc)*y0_scale;
      }

      double w00 = (1 - w_x)*(1 - w_y);
      double w10 = w_x*(1 - w_y);
      double w01 = (1 - w_x)*w_y;
      double w11 = w_x*w_y;

      z0[k    ][l    ] = z0[k    ][l    ] + z[i]*w00;
      z0[k + 1][l    ] = z0[k + 1][l    ] + z[i]*w10;
      z0[k    ][l + 1] = z0[k    ][l + 1] + z[i]*w01;
      z0[k + 1][l + 1] = z0[k + 1][l + 1] + z[i]*w11;

      w0[k    ][l    ] = w0[k    ][l    ] + w00;
      w0[k + 1][l    ] = w0[k + 1][l    ] + w10;
      w0[k    ][l + 1] = w0[k    ][l + 1] + w01;
      w0[k + 1][l + 1] = w0[k + 1][l + 1] + w11;
    }
  }
}
//-------------------------------------------------------------------------------
void LocalLinearRegression2DSurface(const std::vector<double>         & x,
                                    const std::vector<double>         & y,
                                    const std::vector<double>         & z,
                                    const std::vector<double>         & w,
                                    const double                      & bandwidth_x,
                                    const double                      & bandwidth_y,
                                    const double                      & effective_sample_size,
                                    const std::vector<double>         & x0,
                                    const std::vector<double>         & y0,
                                    std::vector<std::vector<double> > & z0,
                                    std::vector<std::vector<double> > & w0,
                                    bool                              & complete_surface)
{
  double               epsilon               = 1e-4;
  double               threshold             = std::pow(std::log(1e5), 0.5);
  double               one_over_bandwidth_x  = 1.0/bandwidth_x;
  double               one_over_bandwidth_y  = 1.0/bandwidth_y;
  double               weight_total_limit    = effective_sample_size;
  double               delta_limit           = 2*std::log(1e5);

  size_t               nSamples              = z.size();

  std::vector<double> weights(nSamples, RMISSING);
  std::vector<size_t> index(nSamples, nSamples + 1);

  std::vector<double> x_order;
  std::vector<size_t> x_backward;
  std::vector<size_t> x_forward;
  SortOrderAndRank(x, x_order, x_backward, x_forward);

  std::vector<double> y_order;
  std::vector<size_t> y_backward;
  std::vector<size_t> y_forward;
  SortOrderAndRank(y, y_order, y_backward, y_forward);


  complete_surface = true;

  for (size_t i = 0; i < x0.size(); i++) {

    double x_low        = x0[i] - threshold*bandwidth_x;
    double x_upp        = x0[i] + threshold*bandwidth_x;

    size_t x_low_k      = FindLowerBoundInSortedVector(x_order, x_low);
    size_t x_upp_k      = FindUpperBoundInSortedVector(x_order, x_upp);

    for (size_t j = 0; j < y0.size(); j++) {

      if (z0[i][j] == RMISSING) {

        double weight_total = 0.0;

        double y_low        = y0[j] - threshold*bandwidth_y;
        double y_upp        = y0[j] + threshold*bandwidth_y;

        size_t y_low_l      = FindLowerBoundInSortedVector(y_order, y_low);
        size_t y_upp_l      = FindUpperBoundInSortedVector(y_order, y_upp);

        size_t K = 0;
        for (size_t k = x_low_k; k < x_upp_k + 1; k++) {

          size_t l = x_backward[k];

          if (y_low_l <= y_forward[l] && y_forward[l] <= y_upp_l) {
            double x_delta = (x0[i] - x[l])*one_over_bandwidth_x;
            double y_delta = (y0[j] - y[l])*one_over_bandwidth_y;
            double delta   = std::pow(x_delta, 2) + std::pow(y_delta, 2);

            if (delta < delta_limit && z[l] != RMISSING) {
              double weight = w[l]*std::exp(-0.5*delta);

              weight_total  = weight_total + weight;
              weights[K]    = weight;
              index[K]      = l;
              K             = K + 1;
            }
          }
        }
        w0[i][j] = weight_total;

        if (K > 0 && weight_total > weight_total_limit) {
          // exact calculation of the inverse.
          // A_inv = B/det_A, where A = X^{t}WX and A_inv = (X^{t}WX)^{-1}
          double a11 = 0.0;
          double a12 = 0.0;
          double a13 = 0.0;
          double a22 = 0.0;
          double a23 = 0.0;
          double a33 = 0.0;

          // C = X^{t}Wz
          double c1  = 0.0;
          double c2  = 0.0;
          double c3  = 0.0;

          for (size_t k = 0; k < K; k++) {
            size_t m      = index[k];
            double weight = weights[k];

            a11 = a11 +                          weight;
            a12 = a12 + x[m]                    *weight;
            a13 = a13 +           y[m]          *weight;
            a22 = a22 + x[m]*x[m]               *weight;
            a23 = a23 + x[m]     *y[m]          *weight;
            a33 = a33 +           y[m]*y[m]     *weight;

            c1 = c1  +                     z[m]*weight;
            c2 = c2  + x[m]               *z[m]*weight;
            c3 = c3  +           y[m]     *z[m]*weight;
          }

          // add a smll number to the diagonal in case the matrix is singular.
          //a11        = a11*(1 + epsilon);
          a22          = a22*(1 + epsilon);
          a33          = a33*(1 + epsilon);

          double b11   =  (a33*a22 - a23*a23);
          double b12   = -(a33*a12 - a23*a13);
          double b13   =  (a23*a12 - a22*a13);
          double b22   =  (a33*a11 - a13*a13);
          double b23   = -(a23*a11 - a12*a13);
          double b33   =  (a22*a11 - a12*a12);

          double det_A = a11*b11 + a12*b12 + a13*b13;

          if (det_A > epsilon) {
            double tmp = 0.0;

            tmp        = tmp +       (b11*c1 + b12*c2 + b13*c3)/det_A;
            tmp        = tmp + x0[i]*(b12*c1 + b22*c2 + b23*c3)/det_A;
            tmp        = tmp + y0[j]*(b13*c1 + b23*c2 + b33*c3)/det_A;

            z0[i][j]   = tmp;

          }
          else {
            z0[i][j]         = RMISSING;
            complete_surface = false;
          }
        }
        else {
          z0[i][j]         = RMISSING;
          complete_surface = false;
        }
      }
    }
  }
}
//--------------------------------------------------------------------------------------------------------------
size_t FindLowerBoundInSortedVector(const std::vector<double> & x,
                                    const double              & x_0)
{
  size_t n   = x.size();
  size_t low = 0;
  size_t upp = n - 1;

  if (x_0 <= x[low]) {
    return(low);
  }
  else if (x_0 >= x[upp]) {
    return(upp);
  }
  else {
    size_t j = static_cast<size_t>(std::floor((upp - low)*0.5));
    while (upp - low > 1) {
      while (x[j] <= x_0 && upp - low > 1) {
        low = j;
        j   = low + static_cast<size_t>(std::floor((upp - low)*0.5));
      }
      while (x_0 < x[j] && upp - low > 1) {
        upp = j;
        j   = low + static_cast<size_t>(std::floor((upp - low)*0.5));
      }
    }
    return(j);
  }
}
//--------------------------------------------------------------------------------------------------------------
size_t FindUpperBoundInSortedVector(const std::vector<double> & x,
                                    const double              & x_0)
{
  size_t n   = x.size();
  size_t low = 0;
  size_t upp = n - 1;

  if (x_0 <= x[low]) {
    return(low);
  }
  else if (x_0 >= x[upp]) {
    return(upp);
  }
  else {
    size_t j = static_cast<size_t>(std::floor((upp - low)*0.5));
    while (upp - low > 1) {
      while (x[j] < x_0 && upp - low > 1) {
        low = j;
        j   = low + static_cast<size_t>(std::floor((upp - low)*0.5));
      }
      while (x_0 <= x[j] && upp - low > 1) {
        upp = j;
        j   = low + static_cast<size_t>(std::floor((upp - low)*0.5));
      }
    }
    return(std::min(j + 1, n - 1));
  }
}
//-------------------------------------------------------------------------------
void SortOrderAndRank(const std::vector<double> & x,
                      std::vector<double>       & x_sort,
                      std::vector<size_t>       & x_order,
                      std::vector<size_t>       & x_rank)
{
  size_t n = x.size();

  x_sort.resize(n, 0);
  x_order.resize(n, 0);
  x_rank.resize(n, 0);

  std::vector<std::pair<double, size_t> > x_order_tmp(n);
  std::vector<std::pair<size_t, size_t> > x_rank_tmp(n);

  for (size_t i = 0; i < n; i++) {
    x_order_tmp[i].first  = x[i];
    x_order_tmp[i].second = i;
  }

  std::sort(x_order_tmp.begin(), x_order_tmp.end());

  for (size_t i = 0; i < n; i++) {
    x_sort[i]            = x_order_tmp[i].first;
    x_order[i]           = x_order_tmp[i].second;

    x_rank_tmp[i].first  = x_order_tmp[i].second;
    x_rank_tmp[i].second = i;
  }

  std::sort(x_rank_tmp.begin(), x_rank_tmp.end());

  for (size_t i = 0; i < n; i++) {
    x_rank[i] = x_rank_tmp[i].second;
  }
}
//-------------------------------------------------------------------------------
void BilinearInterpolation(const std::vector<double>               & x,
                           const std::vector<double>               & y,
                           const std::vector<std::vector<double> > & z,
                           const std::vector<double>               & x0,
                           const std::vector<double>               & y0,
                           std::vector<std::vector<double> >       & z0)
{
  for (size_t i = 0; i < x0.size(); i++) {
    for (size_t j = 0; j < y0.size(); j++) {
      z0[i][j] = Interpolate(x, y, z, x0[i], y0[j]);
    }
  }
}
//-------------------------------------------------------------------------------
double Interpolate(const std::vector<double>               & x,
                   const std::vector<double>               & y,
                   const std::vector<std::vector<double> > & z,
                   const double                            & x0,
                   const double                            & y0)
{
  double delta   = 1e-5;

  double x_first = x.front();
  double x_last  = x.back();
  double x_inc   = x[1] - x[0];
  double x_scale = 1/x_inc;

  double y_first = y.front();
  double y_last  = y.back();
  double y_inc   = y[1] - y[0];
  double y_scale = 1/y_inc;

  if ((x_first - delta < x0) && (x0 < x_last + delta) && (y_first - delta < y0) && (y0 < y_last + delta)) {

    size_t k;
    size_t l;
    double w_y;
    double w_x;

    if (x0 < x_first + delta) {
      k   = 0;
      w_x = 0.0;
    }
    else if (x0 > x_last - delta) {
      k   = x.size() - 2;
      w_x = 1.0;
    }
    else {
      k   = static_cast<size_t>(std::floor((x0 - x_first)*x_scale));
      w_x = (x0 - x_first - k*x_inc)*x_scale;
    }

    if (y0 < y_first + delta) {
      l   = 0;
      w_y = 0.0;
    }
    else if (y0 > y_last - delta) {
      l   = y.size() - 2;
      w_y = 1.0;
    }
    else {
      l   = static_cast<size_t>(std::floor((y0 - y_first)*y_scale));
      w_y = (y0 - y_first - l*y_inc)*y_scale;
    }

    double v00 = z[k    ][l    ];
    double v10 = z[k + 1][l    ];
    double v01 = z[k    ][l + 1];
    double v11 = z[k + 1][l + 1];

    if (v00 != RMISSING && v10 != RMISSING && v01 != RMISSING && v11 != RMISSING) {
      return(v00*(1 - w_x)*(1 - w_y) + v10*w_x*(1 - w_y) + v01*(1 - w_x)*w_y + v11*w_x*w_y);
    }
    else {
      return(RMISSING);
    }
  }
  else {
    return(RMISSING);
  }
}
//-------------------------------------------------------------------------------
void EstimateVariance2D(const std::vector<double>               & x,
                        const std::vector<double>               & y,
                        const std::vector<double>               & z,
                        const std::vector<double>               & x0,
                        const std::vector<double>               & y0,
                        const std::vector<std::vector<double> > & trend,
                        std::vector<std::vector<double> >       & variance,
                        double                                    bandwidth_x,
                        double                                    bandwidth_y,
                        std::string                             & errTxt)
{

  std::vector<double> x0_regridded = x0;
  std::vector<double> y0_regridded = y0;

  /* -- estiamte global variance --*/
  size_t              nSamples     = z.size();


  std::vector<double> x_binned;
  std::vector<double> y_binned;
  std::vector<double> w_binned;
  std::vector<double> z_z_mean_squared;

  x_binned.reserve(nSamples);
  y_binned.reserve(nSamples);
  w_binned.reserve(nSamples);
  z_z_mean_squared.reserve(nSamples);

  double              z_var       = 0.0;
  for (size_t i = 0; i < nSamples; i++) {
    double z_mean_i = Interpolate(x0, y0, trend, x[i], y[i]);
    if (z_mean_i != RMISSING && z[i] != RMISSING) {

      double z_z_mean_squared_i = (z[i] - z_mean_i)*(z[i] - z_mean_i);
      z_var                     = z_var + z_z_mean_squared_i;

      z_z_mean_squared.push_back(z_z_mean_squared_i);
      x_binned.push_back(x[i]);
      y_binned.push_back(y[i]);
      w_binned.push_back(1.0);
    }
  }

  nSamples = z_z_mean_squared.size();
  if (nSamples > 1) {
    z_var = z_var/(nSamples - 1);
  }
  else {
    z_var = RMISSING;
  }
  /* ------------------------------*/


  if (z_var == RMISSING) {
    errTxt += "Error: Unable to estimate 2D variance.\n";
    errTxt += "       The data contains only missing values.\n";
    }
    else {
    MakeNewGridResolution(bandwidth_x, x0, x0_regridded);
    MakeNewGridResolution(bandwidth_y, y0, y0_regridded);

    size_t                            x0_regridded_n = x0_regridded.size();
    size_t                            y0_regridded_n = y0_regridded.size();

    std::vector<std::vector<double> > z0_weighted_var_regridded(x0_regridded_n, std::vector<double>(y0_regridded_n, RMISSING));
    std::vector<std::vector<double> > z0_var_regridded(x0_regridded_n, std::vector<double>(y0_regridded_n, RMISSING));
    std::vector<std::vector<double> > w0_regridded(x0_regridded_n, std::vector<double>(y0_regridded_n, RMISSING));

    std::vector<double> z_z_mean_squared_binned = z_z_mean_squared;
    if (z_z_mean_squared.size() > 16384) {
      MakeBinnedDataset2D(x0_regridded, y0_regridded, x_binned, y_binned, z_z_mean_squared_binned, w_binned);
    }

    double effective_sample_size = CalculateEffectiveSampleSize2D(x_binned, y_binned, bandwidth_x, bandwidth_y);

    if (x_binned.size() < effective_sample_size) {
      effective_sample_size = static_cast<double>(x_binned.size());
      LogKit::LogFormatted(LogKit::High,"\nWARNING : The effective sample size used in estimating 2D variance is larger than the actual sample size.");
      LogKit::LogFormatted(LogKit::High,"\n          The local linear is reduced to a standard linear model.");
    }

    double z_var_weight = effective_sample_size;


    bool complete_surface_dummy = false;
    KernelSmoother2DSurface(x_binned,
                            y_binned,
                            z_z_mean_squared_binned,
                            w_binned,
                            bandwidth_x,
                            bandwidth_y,
                            0.0,
                            x0_regridded,
                            y0_regridded,
                            z0_var_regridded,
                            w0_regridded,
                            complete_surface_dummy);

    for (size_t i = 0; i < x0_regridded_n; i++) {
      for (size_t j = 0; j < y0_regridded_n; j++) {
        double numerator   = z_var_weight*z_var;
        double denominator = z_var_weight;

        if (z0_var_regridded[i][j] != RMISSING) {
          numerator   = numerator   + w0_regridded[i][j]*z0_var_regridded[i][j];
          denominator = denominator + w0_regridded[i][j];
        }
        z0_weighted_var_regridded[i][j] = numerator/denominator;
      }
    }
    BilinearInterpolation(x0_regridded, y0_regridded, z0_weighted_var_regridded, x0, y0, variance);
  }
}
//-------------------------------------------------------------------------------
void KernelSmoother2DSurface(const std::vector<double>         & x,
                             const std::vector<double>         & y,
                             const std::vector<double>         & z,
                             const std::vector<double>         & w,
                             const double                      & bandwidth_x,
                             const double                      & bandwidth_y,
                             const double                      & effective_sample_size,
                             const std::vector<double>         & x0,
                             const std::vector<double>         & y0,
                             std::vector<std::vector<double> > & z0,
                             std::vector<std::vector<double> > & w0,
                             bool                              & complete_surface)
{
  double               threshold             = std::pow(std::log(1e5), 0.5);
  double               one_over_bandwidth_x  = 1.0/bandwidth_x;
  double               one_over_bandwidth_y  = 1.0/bandwidth_y;
  double               weight_total_limit    = effective_sample_size;
  double               delta_limit           = 2*std::log(1e5);

  size_t               nSamples              = z.size();

  std::vector<double> weights(nSamples, RMISSING);
  std::vector<size_t> index(nSamples, nSamples + 1);

  std::vector<double> x_order;
  std::vector<size_t> x_backward;
  std::vector<size_t> x_forward;
  SortOrderAndRank(x, x_order, x_backward, x_forward);

  std::vector<double> y_order;
  std::vector<size_t> y_backward;
  std::vector<size_t> y_forward;
  SortOrderAndRank(y, y_order, y_backward, y_forward);


  complete_surface = true;

  for (size_t i = 0; i < x0.size(); i++) {

    double x_low        = x0[i] - threshold*bandwidth_x;
    double x_upp        = x0[i] + threshold*bandwidth_x;

    size_t x_low_k      = FindLowerBoundInSortedVector(x_order, x_low);
    size_t x_upp_k      = FindUpperBoundInSortedVector(x_order, x_upp);

    for (size_t j = 0; j < y0.size(); j++) {

      if (z0[i][j] == RMISSING) {

        double weight_total = 0.0;
        double z0_hat       = 0.0;

        double y_low        = y0[j] - threshold*bandwidth_y;
        double y_upp        = y0[j] + threshold*bandwidth_y;

        size_t y_low_l      = FindLowerBoundInSortedVector(y_order, y_low);
        size_t y_upp_l      = FindUpperBoundInSortedVector(y_order, y_upp);

        for (size_t k = x_low_k; k < x_upp_k + 1; k++) {

          size_t l = x_backward[k];

          if (y_low_l <= y_forward[l] && y_forward[l] <= y_upp_l) {
            double x_delta = (x0[i] - x[l])*one_over_bandwidth_x;
            double y_delta = (y0[j] - y[l])*one_over_bandwidth_y;
            double delta   = std::pow(x_delta, 2) + std::pow(y_delta, 2);

            if (delta < delta_limit && z[l] != RMISSING) {
              double weight = w[l]*std::exp(-0.5*delta);

              weight_total  = weight_total + weight;
              z0_hat        = z0_hat       + weight*z[l];
            }
          }
        }
        w0[i][j] = weight_total;

        if (weight_total > weight_total_limit) {
          z0[i][j] = z0_hat/weight_total;
        }
        else {
          z0[i][j]         = RMISSING;
          complete_surface = false;
        }
      }
    }
  }
}

//----------------------------------------------------//

RegularSurface<double>
ResampleTrend2D(const RegularSurface<double> & surface,
                const std::vector<double>    & x,
                const std::vector<double>    & y,
                const bool                   & transpose)
{
  int length_x = static_cast<int>(x.size());
  int length_y = static_cast<int>(y.size());

  Grid2D<double> resampled_grid(length_x, length_y);

  for(int i=0; i<length_x; i++) {
    for(int j=0; j<length_y; j++) {
      if(transpose)
        resampled_grid(i,j) = surface.GetZ(y[j],x[i]);
      else
        resampled_grid(i,j) = surface.GetZ(x[i],y[j]);
    }
  }

  double x0 = x[0];
  double y0 = y[0];

  RegularSurface<double> resampled_surface(x0, y0, length_x, length_y, resampled_grid);

  return(resampled_surface);
}

//----------------------------------------------------//

void ReadTrend2DPlainAscii(const std::string     & file_name,
                           std::string           & err_txt,
                           NRLib::Grid2D<double> & trend2d)
{
  std::ifstream file;
  OpenRead(file,file_name);
  int line = 0;

  int ni = ReadNext<int>(file, line);
  int nj = ReadNext<int>(file, line);
  trend2d.Resize(ni, nj);

  for (int i = 0; i < ni; i++) {
    for (int j = 0; j < nj; j++) {
      if (!CheckEndOfFile(file))
        trend2d(i,j) = ReadNext<double>(file, line);
      else
        err_txt += "Premature end of file for 2D trend: " + file_name + ".\n";
    }
  }

  file.close();
}


void ReadTrend3DPlainAscii(const std::string   & file_name,
                           std::string         & /*err_txt*/,
                           NRLib::Grid<double> & trend3d)
{
  std::ifstream file;
  OpenRead(file, file_name);

  int line = 0;

  int ni = ReadNext<int>(file, line);
  int nj = ReadNext<int>(file, line);
  int nk = ReadNext<int>(file, line);
  trend3d.Resize(ni, nj, nk);
  ReadAsciiArrayFast(file, trend3d.begin(), trend3d.GetN());
  file.close();

}


void ReadTrend3DBinary(const std::string   & file_name,
                       Endianess             file_format,
                       std::string           /*err_txt*/,
                       NRLib::Grid<double> & trend3d)
{
  std::ifstream file;
  OpenRead(file, file_name, std::ios::in | std::ios::binary);
  std::vector<int> dim(3);
  ReadBinaryIntArray(file, dim.begin(), 3, file_format);
  trend3d.Resize(dim[0], dim[1], dim[2]);
  ReadBinaryFloatArray(file, trend3d.begin(), trend3d.GetN(), file_format);
  file.close();
}
}
