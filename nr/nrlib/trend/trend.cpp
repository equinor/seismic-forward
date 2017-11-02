// $Id: trend.cpp 1649 2017-07-24 07:44:22Z aarnes $
#include "trend.hpp"
#include "trendkit.hpp"
#include "../exception/exception.hpp"
#include "../random/random.hpp"
#include "../surface/regularsurface.hpp"
#include "../iotools/stringtools.hpp"
#include <cassert>
#include <cmath>

using namespace NRLib;

Trend::Trend()
{
}

Trend::~Trend()
{
}

double
Trend::GetValue2D(double s1, double s2,  const std::vector<double> &trend,
                  size_t ns1_in, size_t ns2_in, size_t start_ind) const {
  //Bilinear interpolation method - always assume 2D plane (can be used by 2D trend and 3D trend)

  int    ns1        = static_cast<int>(ns1_in);
  int    ns2        = static_cast<int>(ns2_in);
  double inv_s1_inc = ns1 - 1.0;
  double inv_s2_inc = ns2 - 1.0;
  double i_frac     = s1*inv_s1_inc;
  double j_frac     = s2*inv_s2_inc;
  int    i0         = static_cast<int>(floor(i_frac));  //left edge
  int    j0         = static_cast<int>(floor(j_frac));  //bottom edge
  int    i1         = i0 + 1;                      //right edge
  int    j1         = j0 + 1;                      //top edge
  double delta_s1   = i_frac - i0;       //proportion inside unit cell i, delta left for s1
  double delta_s2   = j_frac - j0;       //proportion inside unit cell j

  double val = 0;
  if (i0 < 0) {
    if (j0 < 0)
      val = trend[start_ind];
    else if (j1 > ns2 - 1)
      val = trend[start_ind + (ns2 - 1)*ns1];
    else
      val = delta_s2*trend[start_ind + j1*ns1] + (1 - delta_s2)*trend[start_ind  +  j0*ns1];
  }
  else if (i1 > ns1 - 1) {
    if (j0 < 0)
      val = trend[start_ind + ns1 - 1 +  0];
    else if (j1 > ns2 - 1)
      val = trend[start_ind + ns1 - 1 + (ns2 - 1)*ns1];
    else
      val = delta_s2*trend[start_ind + ns1 - 1 + j1*ns1] + (1 - delta_s2)*trend[start_ind + ns1 - 1 +  j0*ns1];
  }
  else if (j0 < 0)// i0 is inside
    val = delta_s1*trend[start_ind + i1 +  0] + (1 - delta_s1)*trend[start_ind + i0 +  0];
  else if (j1 > ns2 - 1)
    val = delta_s1*trend[start_ind + i1 +  (ns2 - 1)*ns1] + (1 - delta_s1)*trend[start_ind + i0 +  (ns2 - 1)*ns1];
  else { //i0 AND j0 are inside
    double val_s1_1 = delta_s1*trend[start_ind + i1 +  j0*ns1] + (1 - delta_s1)*trend[start_ind + i0 +  j0*ns1];
    double val_s1_2 = delta_s1*trend[start_ind + i1 +  j1*ns1] + (1 - delta_s1)*trend[start_ind + i0 +  j1*ns1];

    val = delta_s2*val_s1_2 + (1 - delta_s2)*val_s1_1;
  }

  return val;
}


//-------------------------------------------------------------------//

TrendConstant::TrendConstant(double trend)
: trend_(trend)
{
}

TrendConstant::TrendConstant(const TrendConstant & trend)
: Trend(),
  trend_(trend.trend_)
{
}

TrendConstant::TrendConstant()
{
}

TrendConstant::~TrendConstant()
{
}

std::vector<double>
TrendConstant::GetIncrement(void) const
{
  std::vector<double> increment(3);
  for(int i=0; i<static_cast<int>(increment.size()); i++)
    increment[i] = 0;

  return(increment);
}

std::vector<int>
TrendConstant::GetTrendSize(void) const
{
  std::vector<int> size(3);
  size[0] = 1;
  size[1] = 0;
  size[2] = 0;

  return(size);
}

int
TrendConstant::GetTrendDimension(void) const
{
  int dim = 0;
  return(dim);
}

NRLib::Point
TrendConstant::DrawPoint(const NRLib::Volume &volume) const
{
  double x_rel = NRLib::Random::Unif01();
  double y_rel = NRLib::Random::Unif01();
  double z_rel = NRLib::Random::Unif01();
  double x, y;
  volume.GetXYFromRelative(x_rel, y_rel, x, y);
  NRLib::Point pt;
  pt.x = x;
  pt.y = y;
  pt.z = z_rel*volume.GetLZ();
  return pt;

}

bool
TrendConstant::CheckForNegative() const
{
  if (trend_ < 0)
    return (false);
  else
    return (true);
}
//-------------------------------------------------------------------//

Trend1D::Trend1D(const std::vector<double> & trend, int reference)
: trend_(trend)
{
  if(trend_.empty() || reference > 3 || reference < 1)
    throw NRLib::Exception("Trend1D: Input parameters are not valid.");
  assert(trend_.size() > 1);

  inv_s1_inc_ = (trend_.size() - 1.0);
  reference_  = reference;
}

Trend1D::Trend1D(const std::vector<double> & trend, int reference, double dz)
: trend_(trend)
{
  //Use this constructor when the increment needs to be specified

  if(trend_.empty() || reference > 3 || reference < 1)
    throw NRLib::Exception("Trend1D: Input parameters are not valid.");

  assert(dz > 0);
  reference_  = reference;
  inv_s1_inc_ = 1/dz;
}

Trend1D::Trend1D(const Trend1D & trend)
: Trend(),
  trend_(trend.trend_),
  reference_(trend.reference_),
  inv_s1_inc_(trend.inv_s1_inc_)
{
}

Trend1D::Trend1D(const std::string &filename, int reference)
{
  std::string errTxt;
  int file_format = NRLib::GetTrend1DFileFormat(filename, errTxt);
  if(file_format == 0){
    double s_min, s_max;
    NRLib::ReadTrend1D(filename,errTxt, trend_, s_min, s_max);
    if(errTxt !="")
      throw NRLib::Exception(errTxt);
  }
  else if(file_format==1){
    // read plain ascii from file
    NRLib::ReadTrend1DPlainAscii(filename, errTxt, trend_);
    if(errTxt !="")
      throw NRLib::Exception(errTxt);
  }
  assert(trend_.size() > 1);
  inv_s1_inc_ = (trend_.size() - 1.0);
  if(reference > 3 || reference < 1){
    errTxt = "Wrong reference value in Trend1D. Reference must be either 1, 2 or 3. Value used is " + NRLib::ToString(reference) + "\n";
    throw NRLib::Exception(errTxt);
  }
  reference_ = reference;

}

Trend1D::Trend1D() {
}

Trend1D::~Trend1D() {

}

void
Trend1D::AddConstant(double c)
{
  for(size_t i = 0; i < trend_.size(); i++)
    trend_[i]+=c;

}

double
Trend1D::GetTrendElement(int i, int j, int k) const
{
  double value = 0;
  if (reference_ == 1)
    value = trend_[i];

  else if (reference_ == 2)
    value = trend_[j];

  else
    value = trend_[k];

  return value;
}

double
Trend1D::GetValue(double s1, double s2, double s3) const {
  //Linear interpolation method

  double i = 0;
  if(reference_ == 1)
    i = s1*inv_s1_inc_;
  else if(reference_ == 2)
    i = s2*inv_s1_inc_;
  else if(reference_ == 3)
    i = s3*inv_s1_inc_;

  int i1   = static_cast<int>(floor(i));

  double val = 0;
  if (i1 < 0)
    val = trend_.front();
  else if (i1+1 > static_cast<int>(trend_.size()) - 1)
    val = trend_.back();
  else {
    double t = i - i1;
    val = t*trend_[i1+1] + (1.0 - t)*trend_[i1];
  }

  return val;

}

double
Trend1D::GetValue(double s1, double s2, double s3, const NRLib::Volume &volume) const {
  //Linear interpolation method
   double s1rel, s2rel, s3rel;
   s1rel = 0.0;
   s2rel = 0.0;
   s3rel = 0.0;
   if(reference_ == 1){
    s1rel = s1/volume.GetLX();
   }
   else if(reference_ == 2){
    s2rel = s2/volume.GetLY();
   }
   else{
    s3rel = s3/volume.GetLZ();
   }
   double val = GetValue(s1rel, s2rel, s3rel);  //GetValue checks for values outside volume
   return val;


}

std::vector<double>
Trend1D::GetIncrement(void) const
{
  std::vector<double> increment(3);
  for(int i=0; i<static_cast<int>(increment.size()); i++)
    increment[i] = 0;
  increment[reference_-1] = 1/inv_s1_inc_;

  return(increment);
}

std::vector<int>
Trend1D::GetTrendSize(void) const
{
  std::vector<int> size(3);
  for(int i=0; i<3; i++)
    size[i] = 0;
  size[reference_-1] = static_cast<int>(trend_.size());

  return(size);
}

int
Trend1D::GetTrendDimension(void) const
{
  int dim = 1;
  return(dim);
}

double
Trend1D::GetMaxValue(void) const
{
  assert (trend_.size() > 0);
  double max = trend_[0];
  for (size_t i = 1; i < trend_.size(); i++)
    if (trend_[i] > max)
      max = trend_[i];

  return (max);
}

double
Trend1D::GetMinValue(void) const
{
  assert (trend_.size() > 0);
  double min = trend_[0];
  for (size_t i = 1; i < trend_.size(); i++)
    if (trend_[i] < min)
      min = trend_[i];

  return (min);
}
double
Trend1D::GetMeanValue(void) const
{
  assert (trend_.size() > 0);
  double mean = trend_[0];
  for (size_t i = 1; i < trend_.size(); i++)
    mean +=trend_[i];

  return (mean/trend_.size());
}

void
Trend1D::ScaleToAverage(double factor) {
  //normalize trend to have expectation = factor.
  double sum_trend = 0.0;
  for (size_t i = 0; i < trend_.size(); i++)
    sum_trend += trend_[i];

  double avg_trend = sum_trend / trend_.size();
  double scale = factor/avg_trend;
  for (size_t i = 0; i < trend_.size(); i++)
    trend_[i] *= scale;
}

bool
Trend1D::CheckForNegative() const
{
  bool ok = true;
  for (size_t i = 1; i < trend_.size(); i++)
    if (trend_[i] < 0)
      ok = false;

  return ok;
}

//-------------------------------------------------------------------//
Trend2D::Trend2D(const NRLib::Grid2D<double> & trend,
                 int                           reference1,
                 int                           reference2)
: trend_(trend)
{
  // Use this constructor when the total length of the trend is one

  ns1_ = static_cast<int>(trend_.GetNI());
  ns2_ = static_cast<int>(trend_.GetNJ());

  inv_s1_inc_ = (ns1_ - 1.0);
  inv_s2_inc_ = (ns2_ - 1.0);

  reference1_ = reference1;
  reference2_ = reference2;
}


Trend2D::Trend2D(const NRLib::Grid2D<double> & trend,
                 int                           reference1,
                 int                           reference2,
                 double                        dz1,
                 double                        dz2)
: trend_(trend)
{
  //Use this constructor when the increment needs to be specified

  ns1_ = static_cast<int>(trend_.GetNI());
  ns2_ = static_cast<int>(trend_.GetNJ());

  inv_s1_inc_ = 1/dz1;
  inv_s2_inc_ = 1/dz2;

  reference1_ = reference1;
  reference2_ = reference2;
}


Trend2D::Trend2D(const Trend2D & trend)
: Trend(),
  reference1_(trend.reference1_),
  reference2_(trend.reference2_),
  trend_(trend.trend_),
  ns1_(trend.ns1_),
  ns2_(trend.ns2_),
  inv_s1_inc_(trend.inv_s1_inc_),
  inv_s2_inc_(trend.inv_s2_inc_)
{
}


Trend2D::Trend2D(const std::string & filename,
                 int                 reference1,
                 int                 reference2)
{
  std::ifstream file;
  OpenRead(file, filename, std::ios::in | std::ios::binary);
  std::string line_string;
  std::getline(file, line_string);
  std::vector<std::string> tokens = NRLib::GetTokens(line_string);
  file.close();
  if (tokens.size() == 2) {//Test needed for hufsa project formats
    std::string err_txt;
    NRLib::ReadTrend2DPlainAscii(filename,
                          err_txt,
                          trend_);
    if (err_txt != "")
      throw NRLib::Exception(err_txt);
  }
  else{
    NRLib::SurfaceFileFormat file_format = FindSurfaceFileType(filename);
    if (file_format != NRLib::SURF_UNKNOWN) {
      NRLib::RegularSurface<double> flate(filename);
      trend_ = flate;
    }
  }

  reference1_ = reference1;
  reference2_ = reference2;
  ns1_ = static_cast<int>(trend_.GetNI());
  ns2_ = static_cast<int>(trend_.GetNJ());
  inv_s1_inc_ = ns1_ - 1;
  inv_s2_inc_ = ns2_ - 1;
}


Trend2D::Trend2D() {
}


Trend2D::~Trend2D() {
}


void
Trend2D::AddConstant(double c)
{
  for(size_t i = 0; i < trend_.GetNI(); i++)
    for(size_t j = 0; j < trend_.GetNJ(); j++)
      trend_(i, j) += c;

}


double
Trend2D::GetTrendElement(int i, int j, int k) const
{
  int ii = 0;
  int jj = 0;
  if(reference1_ == 1)
    ii = i;
  else if(reference1_ == 2)
    ii = j;
  else if(reference1_ == 3)
    ii = k;
  if(reference2_ == 1)
    jj = i;
  else if(reference2_==2)
    jj = j;
  else if(reference2_ == 3)
    jj = k;
  return trend_(ii, jj);

}


double
Trend2D::GetValue(double                s1,
                  double                s2,
                  double                s3,
                  const NRLib::Volume & volume) const {
  //Linear interpolation method
   double s1rel, s2rel;
   double val;
   if(reference1_ == 2 && reference2_ == 3) {     //Constant in x-direction
     s1rel = s2/volume.GetLY();
     s2rel = s3/volume.GetLZ();
   }
   else if(reference1_ == 1 && reference2_ == 3) { //Constant in y-direction
    s1rel = s1/volume.GetLX();
    s2rel = s3/volume.GetLZ();
   }
   else {                                         //Constant in z-direction
    s1rel = s1/volume.GetLX();
    s2rel = s2/volume.GetLY();
   }

   val = Trend::GetValue2D(s1rel, s2rel, trend_.GetStorage(), trend_.GetNI(), trend_.GetNJ(), 0);

   return val;
}

double
Trend2D::GetValue(double s1, double s2, double /*s3*/) const {
  //Bilinear interpolation method
  double i = s1*inv_s1_inc_;
  int i1   = static_cast<int>(floor(i));

  double j = s2*inv_s2_inc_;
  int j1   = static_cast<int>(floor(j));

  double t1 = i - i1;
  double t2 = j - j1;

  double val = 0;
  if (i1 < 0) {
    if (j1 < 0)
      val = trend_(0);
    else if (j1 + 1 > ns2_ - 1)
      val = trend_(0 + (ns2_ - 1)*ns1_);
    else
      val = t2*trend_(0 + (j1+1)*ns1_) + (1.0 - t2)*trend_(0 + j1*ns1_);
  }
  else if (i1 + 1 > ns1_ - 1) {
    if (j1 < 0)
      val = trend_(ns1_ - 1 + 0);
    else if (j1 + 1 > ns2_ - 1)
      val = trend_(ns1_ - 1 + (ns2_ - 1)*ns1_);
    else
      val = t2*trend_(ns1_ - 1 + (j1+1)*ns1_) + (1.0 - t2)*trend_(ns1_ - 1 + j1*ns1_);
  }
  else if (j1 < 0)// i1 is inside
    val = t1*trend_(i1+1 + 0) + (1.0 - t1)*trend_(i1 + 0);
  else if (j1 + 1 > ns2_ - 1)
    val = t1*trend_(i1+1 + (ns2_-1)*ns1_) + (1.0 - t1)*trend_(i1 + (ns2_-1)*ns1_);
  else { //i1 AND j1 are inside
    double val_s1_1 = t1*trend_(i1+1 + j1*ns1_) + (1.0 - t1)*trend_(i1 + j1*ns1_);
    double val_s1_2 = t1*trend_(i1+1 + (j1+1)*ns1_) + (1.0 - t1)*trend_(i1 + (j1+1)*ns1_);

    val = t2*val_s1_2 + (1-t2)*val_s1_1;
  }

  return val;
}


std::vector<double>
Trend2D::GetIncrement(void) const
{
  std::vector<double> increment(3);
  if(reference1_==1 && reference2_ ==2){ // constant in z
    increment[0] = 1/inv_s1_inc_;
    increment[1] = 1/inv_s2_inc_;
    increment[2] = 0;
  }
  else if(reference1_ == 2 && reference2_ == 3){ // constant in x
    increment[1] = 1/inv_s1_inc_;
    increment[2] = 1/inv_s2_inc_;
    increment[0] = 0;
  }
  else if(reference1_ == 1 && reference2_ ==3){ // constant in y
    increment[0] = 1/inv_s1_inc_;
    increment[2] = 1/inv_s2_inc_;
    increment[1] = 0;
  }
  return(increment);
}


std::vector<int>
Trend2D::GetTrendSize(void) const
{
  std::vector<int> size(3);
  if(reference1_==1 && reference2_ ==2){ // constant in z
    size[0] = ns1_;
    size[1] = ns2_;
    size[2] = 0;
  }
  else if(reference1_ == 2 && reference2_ == 3){ // constant in x
    size[1] = ns1_;
    size[2] = ns2_;
    size[0] = 0;
  }
  else if(reference1_ == 1 && reference2_ ==3){ // constant in y
    size[0] = ns1_;
    size[2] = ns2_;
    size[1] = 0;
  }
  return(size);
}


int
Trend2D::GetTrendDimension(void) const
{
  int dim = 2;
  return(dim);
}


double
Trend2D::GetMaxValue(void) const
{
  assert (trend_.GetN() > 0);
  double max = trend_(0);
  for (size_t i = 0; i < trend_.GetNI(); i++)
    for (size_t j = 0; j < trend_.GetNJ(); j++)
      if (trend_(i,j) > max)
        max = trend_(i,j);

  return(max);
}


double
Trend2D::GetMinValue(void) const
{
  assert (trend_.GetN() > 0);
  double min = trend_(0);
  for (size_t i = 0; i < trend_.GetNI(); i++)
    for (size_t j = 0; j < trend_.GetNJ(); j++)
      if (trend_(i,j) < min)
        min = trend_(i,j);

  return(min);
}


double
Trend2D::GetMeanValue(void) const
{
  assert (trend_.GetN() > 0);
  double mean = trend_(0);
  for (size_t i = 0; i < trend_.GetNI(); i++)
    for (size_t j = 0; j < trend_.GetNJ(); j++)
      mean += trend_(i,j);

  return(mean/(trend_.GetNI()*trend_.GetNJ()));
}

void
Trend2D::ScaleToAverage(double factor) {
  //normalize trend to have expectation = factor.
  std::vector<int> dim = GetTrendSize();
  double sum_trend = 0.0;
  for (int i = 0; i < dim[0]; i++)
    for (int j = 0; j < dim[1]; j++)
      for (int k = 0; k < dim[2]; k++)
        sum_trend += GetTrendElement(i, j, k);

  double avg_trend = sum_trend / (dim[0] * dim[1]);
  double scale     = factor/avg_trend;
  for (size_t i = 0; i < trend_.GetNI(); i++)
    for (size_t j = 0; j < trend_.GetNJ(); j++)
        trend_(i, j) *= scale;
}

bool
Trend2D::CheckForNegative() const
{
  bool ok = true;
  for (size_t i = 0; i < trend_.GetNI(); i++)
    for (size_t j = 0; j < trend_.GetNJ(); j++)
      if (trend_(i,j) < 0)
        ok = false;

  return ok;
}


//-------------------------------------------------------------------//

Trend3D::Trend3D(NRLib::Grid<double> & values)
  : Trend()
{

  trend_ = values;
}


Trend3D::Trend3D(const Trend3D & trend)
: Trend(),
  trend_(trend.trend_)
{
}


Trend3D::Trend3D(const std::string & file_name,
                 bool                binary,
                 Endianess           file_format)
  : Trend()
{
  std::string err_txt;
  if(binary == false)
    ReadTrend3DPlainAscii(file_name,
                          err_txt,
                          trend_);
  else
    ReadTrend3DBinary(file_name,
                      file_format,
                      err_txt,
                      trend_);

  if(err_txt !="")
    throw NRLib::Exception(err_txt);

}


Trend3D::~Trend3D()
{
}


void
Trend3D::AddConstant(double c)
{
  for(size_t i = 0; i < trend_.GetNI(); i++)
    for(size_t j = 0; j < trend_.GetNJ(); j++)
      for(size_t k = 0; k < trend_.GetNK(); k++)
        trend_(i, j, k) += c;

}


double
Trend3D::GetValue(double s1, double s2, double s3) const
{
  assert(s1 >= 0 && s2 >= 0 && s3 >= 0);
  assert(s1 <= 1 && s2 <= 1 && s3 <= 1);

  //Linear interpolation method between two bilinear interpolations (above, below)
  size_t ni         = trend_.GetNI();
  size_t nj         = trend_.GetNJ();
  size_t nk         = trend_.GetNK();
  double inv_s3_inc = nk - 1.0;
  double k_frac     = s3*inv_s3_inc;
  size_t k0         = static_cast<size_t>(floor(k_frac));

  size_t top_index  = k0*ni*nj;
  size_t bot_index  = (k0 + 1)*ni*nj;
  //bilinear
  double val_top = Trend::GetValue2D(s1, s2, trend_.GetStorage(), ni, nj, top_index);
  double val_bot;
  if (k0 + 1 == nk) {
    val_bot = val_top;
  }
  else {
    val_bot = Trend::GetValue2D(s1, s2, trend_.GetStorage(), ni, nj, bot_index);
  }

  double t          = k_frac - k0;
  double val        = t*val_top + (1.0 - t)*val_bot;

  return val;
}


double
Trend3D::GetValue(double                s1,
                  double                s2,
                  double                s3,
                  const NRLib::Volume & volume) const {
  //Linear interpolation method
  //Uses global coordinates
  double s1rel;
  double s2rel;
  double s3rel;
  double val;
   s1rel = s1/volume.GetLX();
   s2rel = s2/volume.GetLY();
   s3rel = s3/volume.GetLZ();

   if (s1rel < 0)    //This check is made because GetValue does not accept location points outside volume
     s1rel = 0.0;
   if (s2rel < 0)
     s2rel = 0.0;
   if (s3rel < 0)
     s3rel = 0.0;

   if (s1rel > 1)    //This check is made because GetValue does not accept location points outside volume
     s1rel = 1.0;
   if (s2rel > 1)
     s2rel = 1.0;
   if (s3rel > 1)
     s3rel = 1.0;

   val = GetValue(s1rel, s2rel, s3rel);

   return val;
}


std::vector<double>
Trend3D::GetIncrement(void) const
{
  // Dummy function. Implement if needed
  std::vector<double> increment(3);
  increment[0] = 0;
  increment[1] = 0;
  increment[2] = 0;
  return(increment);
}


std::vector<int>
Trend3D::GetTrendSize(void) const
{
  std::vector<int> size(3);
  size[0] = static_cast<int>(trend_.GetNI());
  size[1] = static_cast<int>(trend_.GetNJ());
  size[2] = static_cast<int>(trend_.GetNK());
  return(size);
}


int
Trend3D::GetTrendDimension(void) const {
  int dim = 3;
  return(dim);
}


double
Trend3D::GetMaxValue(void) const
{
  assert (trend_.GetN() > 0);
  double max = trend_(0);
  for (size_t i = 0; i < trend_.GetNI(); i++)
    for (size_t j = 0; j < trend_.GetNJ(); j++)
      for (size_t k = 0; k < trend_.GetNK(); k++)
        if (trend_(i,j,k) > max)
          max = trend_(i,j,k);

  return(max);
}


double
Trend3D::GetMinValue(void) const
{
  assert (trend_.GetN() > 0);
  double min = trend_(0);
  for (size_t i = 0; i < trend_.GetNI(); i++)
    for (size_t j = 0; j < trend_.GetNJ(); j++)
      for (size_t k = 0; k < trend_.GetNK(); k++)
        if (trend_(i,j,k) < min)
          min = trend_(i,j,k);

  return(min);
}


double
Trend3D::GetMeanValue(void) const
{
  assert (trend_.GetN() > 0);
  double mean = trend_(0);
  for (size_t i = 0; i < trend_.GetNI(); i++)
    for (size_t j = 0; j < trend_.GetNJ(); j++)
      for (size_t k = 0; k < trend_.GetNK(); k++)
        mean += trend_(i,j, k);

  return(mean/(trend_.GetNI()*trend_.GetNJ()*trend_.GetNK()));
}

void
Trend3D::ScaleToAverage(double factor) {
  //normalize trend to have expectation = factor.
  std::vector<int> dim = GetTrendSize();
  double sum_trend = 0.0;
  for (int i = 0; i < dim[0]; i++) {
    for (int j = 0; j < dim[1]; j++) {
      for (int k = 0; k < dim[2]; k++) {
        sum_trend += GetTrendElement(i, j, k);
      }
    }
  }
  double avg_trend = sum_trend / (dim[0] * dim[1] * dim[2]);
  double scale = factor / avg_trend;
  double debug_sum = 0.0;
  for (size_t i = 0; i < trend_.GetNI(); i++) {
    for (size_t j = 0; j < trend_.GetNJ(); j++) {
      for (size_t k = 0; k < trend_.GetNK(); k++) {
        trend_(i, j, k) *= scale;
        debug_sum += trend_(i, j, k);
      }
    }
  }
}

bool
Trend3D::CheckForNegative() const
{
  bool ok = true;
  for (size_t i = 0; i < trend_.GetNI(); i++)
    for (size_t j = 0; j < trend_.GetNJ(); j++)
      for (size_t k = 0; k < trend_.GetNK(); k++)
        if (trend_(i, j,k) < 0)
          ok = false;

  return ok;
}
