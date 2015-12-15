#include <seismic_parameters.hpp>
#include <ctime>

#include <nrlib/math/constants.hpp>
#include <physics/wavelet.hpp>
#include <nrlib/eclipsegrid/eclipsegrid.hpp>
#include <nrlib/surface/regularsurfacerotated.hpp>
#include <nrlib/segy/segygeometry.hpp>
#include <nrlib/stormgrid/stormcontgrid.hpp>
#include "nrlib/geometry/interpolation.hpp"
#include "seismic_geometry.hpp"
#include <physics/zoeppritz.hpp>
#include <physics/zoeppritz_ps.hpp>
#include <physics/zoeppritz_pp.hpp>


SeismicParameters::SeismicParameters(ModelSettings *model_settings) 
 : missing_(-999.0)
{
  this->model_settings_ = model_settings;

  seismic_geometry_ = new SeismicGeometry();

  if (model_settings->GetNMOCorr()) {
    CalculateOffsetSpan();
  }
  else {
    CalculateAngleSpan();
  }

  segy_geometry_ = NULL;

  SetupWavelet();
  //time_t t1 = time(0);   // get time now
  ReadEclipseGrid();
  //PrintElapsedTime(t1, "reading eclipsegrid");
  //t1 = time(0);
  FindGeometry();


  seismic_output_ = new SeismicOutput(model_settings_);
  FindSurfaceGeometry();
  CreateGrids();
  //PrintElapsedTime(t1, "finding geometry and creating grids");
}

void SeismicParameters::CalculateAngleSpan() {
  theta_0_ = model_settings_->GetTheta0();
  dtheta_ = model_settings_->GetDTheta();
  theta_max_ = model_settings_->GetThetaMax();

  if (dtheta_ == 0) {
    ntheta_ = 1;
  }
  else {
    ntheta_ = static_cast<size_t>((theta_max_ - theta_0_) / dtheta_ + 1.01);
  }
  theta_vec_.resize(ntheta_);
  for (size_t i = 0; i < ntheta_; ++i) {
    theta_vec_[i] = theta_0_ + i*dtheta_;
  }
}


void SeismicParameters::CalculateOffsetSpan() {
  offset_0_ = model_settings_->GetOffset0();
  doffset_ = model_settings_->GetDOffset();
  offset_max_ = model_settings_->GetOffsetMax();

  if (doffset_ == 0) {
    noffset_ = 1;
  }
  else {
    noffset_ = size_t((offset_max_ - offset_0_) / doffset_) + 1;
  }
  offset_vec_.resize(noffset_);
  for (size_t i = 0; i < noffset_; ++i) {
    offset_vec_[i] = offset_0_ + i*doffset_;
  }
}

void SeismicParameters::SetupWavelet() {
  if (model_settings_->GetRicker()) {
    double peakF = model_settings_->GetPeakFrequency();
    wavelet_ = new Wavelet(peakF);
  }
  else {
    std::string wavelet_file_format = model_settings_->GetWaveletFileFormat();
    std::string wavelet_file_name = model_settings_->GetWaveletFileName();
    wavelet_ = new Wavelet(wavelet_file_name, wavelet_file_format);
  }
  wavelet_scale_ = model_settings_->GetWaveletScale();
}

void SeismicParameters::SetSegyGeometry(const NRLib::SegyGeometry &geometry)
{
  if (segy_geometry_ != NULL) {
    delete segy_geometry_;
  }
  segy_geometry_ = new NRLib::SegyGeometry(geometry);
}

void SeismicParameters::FindLoopIndeces(int               &n_xl,
                                        int               &il_min,
                                        int               &il_max,
                                        int               &il_step,
                                        int               &xl_min,
                                        int               &xl_max,
                                        int               &xl_step,
                                        bool              &segy)
{
  if (segy_geometry_ == NULL) {
    il_min  = 0;
    il_max  = seismic_geometry_->nx();
    il_step = 1;
    xl_min  = 0;
    xl_max  = seismic_geometry_->ny()-1;
    xl_step = 1;
    n_xl = xl_max;
    segy = false;
  }
  else {
    segy_geometry_->FindILXLGeometry();
    il_min  = segy_geometry_->GetMinIL();
    il_max  = segy_geometry_->GetMaxIL();
    il_step = segy_geometry_->GetILStep();
    xl_min  = segy_geometry_->GetMinXL();
    xl_max  = segy_geometry_->GetMaxXL();
    xl_step = segy_geometry_->GetXLStep();
    n_xl = (xl_max - xl_min + 1) / xl_step;
    segy = true;
  }
}



void SeismicParameters::FindVrms(std::vector<double>       &vrms_vec,
                                 std::vector<double>       &vrms_vec_reg,
                                 const std::vector<double> &twt_vec,
                                 const std::vector<double> &twt_vec_reg,
                                 const std::vector<double> &v_vec,
                                 double                     const_v,
                                 size_t                     i,
                                 size_t                     j,
                                 bool                       include_regular) const
{
  double v_w = model_settings_->GetVw();
  double z_w = model_settings_->GetZw();

  double v_over;
  double twt_w = 2000 * z_w / v_w;
  double tmp, tmp0;
  size_t nk = twt_vec.size();

  //generate vrms_vec - only for each layer in reservoir
  if (twt_vec[0] == -999.0) {
    for (size_t k = 0; k < nk; ++k) {
      vrms_vec[k] = -999.0;
    }
  }
  else {
    v_over = 2000 * ((*zgrid_)(i, j, 0) - z_w) / (twt_vec[0] - 2000 * z_w / v_w);
    tmp0 = v_w*v_w*twt_w + v_over*v_over*(twt_vec[0] - twt_w);
    for (size_t k = 0; k < nk; ++k) {
      tmp = tmp0;
      for (size_t l = 1; l <= k; ++l) {
        tmp += v_vec[l]* v_vec[l]*(twt_vec[l] - twt_vec[l-1]);
      }
      vrms_vec[k] = std::sqrt(tmp / twt_vec[k]);
    }
  }

  if (include_regular) {
    //generate vrms_vec_reg - including overburden and underburden
    
    double twt_wavelet = wavelet_->GetTwtWavelet();
    double vrms_under  = vrms_vec[nk - 1]*vrms_vec[nk - 1]*twt_vec[nk - 1] + const_v * const_v * twt_wavelet;
    vrms_under *= (1 / (twt_vec[nk - 1] + twt_wavelet));
    vrms_under  = std::sqrt(vrms_under);

    //sample vrms regularly:
    std::vector<double> twt_vec_in(nk + 2), vrms_vec_in(nk + 2);
    twt_vec_in [0] = twt_w;
    vrms_vec_in[0] = v_w;
    twt_vec_in [1] = twt_vec[0];
    vrms_vec_in[1] = vrms_vec[0];
    size_t index = 2;
    for (size_t k = 0; k < nk; ++k) {
      if (twt_vec[k] != twt_vec_in[index - 1]) {
        twt_vec_in [index] = twt_vec[k];
        vrms_vec_in[index] = vrms_vec[k];
        ++index;
      }
    }
    twt_vec_in.resize(index + 1);
    vrms_vec_in.resize(index + 1);
    twt_vec_in [index] = twt_vec_in[index - 1] + twt_wavelet;
    vrms_vec_in[index] = vrms_under;

    vrms_vec_reg = NRLib::Interpolation::Interpolate1D(twt_vec_in, vrms_vec_in, twt_vec_reg, "linear");
  }  
}



void SeismicParameters::FindReflections(NRLib::Grid2D<double>       &r_vec,
                                        const std::vector<double>   &theta_vec,
                                        size_t                       i,
                                        size_t                       j){



  bool                ps_seis  = model_settings_->GetPSSeismic();

  Zoeppritz *zoeppritz = NULL;
  if (ps_seis) {
    zoeppritz = new ZoeppritzPS();
  } else {
    zoeppritz = new ZoeppritzPP();
  }

  double diffvp, meanvp, diffvs, meanvs, diffrho, meanrho;
  std::vector<double> vp_vec(bottom_k_-top_k_+3), vs_vec(bottom_k_-top_k_+3), rho_vec(bottom_k_-top_k_+3);

  for (size_t theta = 0; theta < theta_vec.size(); ++theta) {
    for (size_t k = top_k_; k <= bottom_k_ + 2; k++) {
      vp_vec[k - top_k_]  = (*vpgrid_)(i, j, (k - top_k_));
      vs_vec[k - top_k_]  = (*vsgrid_)(i, j, (k - top_k_));
      rho_vec[k - top_k_] = (*rhogrid_)(i, j, (k - top_k_));
    }
    zoeppritz->ComputeConstants(theta_vec[theta]);
    for (size_t k = top_k_; k <= bottom_k_ + 1; k++) {
        diffvp = vp_vec[k - top_k_ + 1] - vp_vec[k - top_k_];
        meanvp = 0.5 *  (vp_vec[k - top_k_ + 1] + vp_vec[k - top_k_]);
        diffvs = vs_vec[k - top_k_ + 1] - vs_vec[k - top_k_];
        meanvs = 0.5 *  (vs_vec[k - top_k_ + 1] + vs_vec[k - top_k_]);
        diffrho = rho_vec[k - top_k_ + 1] - rho_vec[k - top_k_];
        meanrho = 0.5 * (rho_vec[k - top_k_ + 1] + rho_vec[k - top_k_]);
        r_vec(k - top_k_, theta) = static_cast<float>(zoeppritz->GetReflection(diffvp, meanvp, diffrho, meanrho, diffvs, meanvs));
    }
  }
  delete zoeppritz;
}


void SeismicParameters::FindNMOReflections(NRLib::Grid2D<double>       &r_vec,
                                           const NRLib::Grid2D<double> &theta,
                                           size_t                       i,
                                           size_t                       j){

  bool                ps_seis  = model_settings_->GetPSSeismic();

  Zoeppritz *zoeppritz = NULL;
  if (ps_seis) {
    zoeppritz = new ZoeppritzPS();
  } else {
    zoeppritz = new ZoeppritzPP();
  }

  double diffvp, meanvp, diffvs, meanvs, diffrho, meanrho;
  std::vector<double> vp_vec(bottom_k_-top_k_+3), vs_vec(bottom_k_-top_k_+3), rho_vec(bottom_k_-top_k_+3);

  for (size_t off = 0; off < theta.GetNJ(); ++off) {
    for (size_t k = top_k_; k <= bottom_k_ + 2; k++) {
      vp_vec[k - top_k_]  = (*vpgrid_)(i, j, (k - top_k_));
      vs_vec[k - top_k_]  = (*vsgrid_)(i, j, (k - top_k_));
      rho_vec[k - top_k_] = (*rhogrid_)(i, j, (k - top_k_));
      }
    for (size_t k = top_k_; k <= bottom_k_ + 1; k++) {
        diffvp = vp_vec[k - top_k_ + 1] - vp_vec[k - top_k_];
        meanvp = 0.5 *  (vp_vec[k - top_k_ + 1] + vp_vec[k - top_k_]);
        diffvs = vs_vec[k - top_k_ + 1] - vs_vec[k - top_k_];
        meanvs = 0.5 *  (vs_vec[k - top_k_ + 1] + vs_vec[k - top_k_]);
        diffrho = rho_vec[k - top_k_ + 1] - rho_vec[k - top_k_];
        meanrho = 0.5 * (rho_vec[k - top_k_ + 1] + rho_vec[k - top_k_]);
        zoeppritz->ComputeConstants(theta(k, off));
        r_vec(k - top_k_, off) = static_cast<float>(zoeppritz->GetReflection(diffvp, meanvp, diffrho, meanrho, diffvs, meanvs));
      }
    }
  delete zoeppritz;
}

void SeismicParameters::FindMaxTwtIndex(size_t & i_max,
                                        size_t & j_max,
                                        double & max_value)
{
  max_value = 0;
  size_t k_max = (*twtgrid_).GetNK() - 1;
  for (size_t i = 0; i < (*twtgrid_).GetNI(); ++i) {
    for (size_t j = 0; j < (*twtgrid_).GetNJ(); ++j) {
      if ((*twtgrid_)(i, j, k_max) > max_value) {
        max_value = (*twtgrid_)(i, j, k_max);
        i_max = i;
        j_max = j;
      }
    }
  }
}

void SeismicParameters::GenerateTwt0AndZ0(std::vector<double> &twt_0,
                                          std::vector<double> &z_0,
                                          std::vector<double> &twts_0, 
                                          size_t              &time_samples_stretch,
                                          bool                 ps_seis)
{
  if (model_settings_->GetNMOCorr()){
    twt_0 = GenerateTwt0ForNMO(time_samples_stretch, ps_seis);
    z_0   = GenerateZ0ForNMO();
    if (model_settings_->GetTwtFileName() != "") {
      twts_0 = GenerateTWT0Shift(twt_0[0], time_samples_stretch);
    }
  }
  else {
    double tmin = seismic_geometry_->t0();
    double dt   = seismic_geometry_->dt();
    size_t nz   = seismic_geometry_->nz();
    size_t nt   = seismic_geometry_->nt();
    twt_0.resize(nt);
    for (size_t i = 0; i < nt; ++i){
      twt_0[i] = tmin + i*dt;
    }
    double zmin = seismic_geometry_->z0();
    double dz   = seismic_geometry_->dz();
    z_0.resize(nz);
    //find min z in a sample
    size_t nzmin = static_cast<size_t>(floor(zmin) / dz + 0.5);
    double zmin_sampl = nzmin *dz;

    for (size_t i = 0; i < nz; ++i){
      z_0[i] = zmin_sampl + i*dz;
    }

    if (model_settings_->GetTwtFileName() != "") {
      twts_0 = GenerateTWT0Shift(twt_0[0], twt_0.size());
    }
  }
}

std::vector<double> SeismicParameters::GenerateTwt0ForNMO(size_t & nt_stretch,
                                                          bool ps_seis)
{
  //Account for stretch by making twt0 sufficiently long. Stretch upwards is also taken into account
  //through "xtra_samples_top". 
  //Max twt value and location is found from twtgrid and twtx is calculated for the largest offset.
  //"Time samples stretch" is number of samples in nmo-corrected seismic trace.

  size_t i_max, j_max;
  double twt_max, twtx_max;
  size_t nt                      = seismic_geometry_->nt();
  double dt                      = seismic_geometry_->dt();
  double t0                      = seismic_geometry_->t0();
  std::vector<double> constvp    = model_settings_->GetConstVp();
  double twt_wavelet             = wavelet_->GetTwtWavelet();
  size_t nzrefl                  = seismic_geometry_->zreflectorcount();

  //find max from twgrid, and index of max twt value
  FindMaxTwtIndex(i_max, j_max, twt_max);

  //find max TWTX for highest offset in order to find the highest TWT value to sample seismic
  if (ps_seis) { //------------PS seismic------------
    std::vector<double> vrms_pp_vec(nzrefl), vrms_ss_vec(nzrefl), dummy;
    std::vector<double> vp_vec(nzrefl), vs_vec(nzrefl), twt_pp_vec(nzrefl), twt_ss_vec(nzrefl);

    for (size_t k = 0; k < nzrefl; ++k) {
      twt_pp_vec[k] = (*twtppgrid_)(i_max, j_max, k);
      twt_ss_vec[k] = (*twtssgrid_)(i_max, j_max, k);
      vp_vec[k]     = (*vpgrid_)(i_max, j_max, k);
      vs_vec[k]     = (*vsgrid_)(i_max, j_max, k);
    }

    FindVrms(vrms_pp_vec, dummy, twt_pp_vec, dummy, vp_vec, 1.0, i_max, j_max, false);
    FindVrms(vrms_ss_vec, dummy, twt_ss_vec, dummy, vs_vec, 1.0, i_max, j_max, false);

    double vrms_pp     = vrms_pp_vec[nzrefl - 1];
    double vrms_ss     = vrms_ss_vec[nzrefl - 1];
    double twt_pp_max  = twt_pp_vec[nzrefl - 1];
    double twt_ss_max  = twt_ss_vec[nzrefl - 1];
    double offset_max = offset_0_ + doffset_*(noffset_ - 1);
    double tmp         = offset_max / (vrms_pp_vec[nzrefl - 1]*twt_pp_max / 1000);
    double start_value = atan(tmp);
    if (start_value >= 1.0)
      start_value = 0.99;
    double dU    = vrms_ss * twt_ss_max / 2000;
    double dD    = vrms_pp* twt_pp_max / 2000;
    double vr    = vrms_ss / vrms_pp;
    size_t n_it  = 10;
    double y_out = FindSinThetaPSWithNewtonsMethod(start_value,
                                                   offset_max,
                                                   dU,
                                                   dD,
                                                   vr,
                                                   0.00001,
                                                   n_it);
    double theta_ss  = asin(vr*y_out);
    double theta_pp  = asin(y_out);
    double offset_pp = tan(theta_pp)*dD;
    double offset_ss = tan(theta_ss)*dU;

    double twtx_pp = std::sqrt(twt_pp_max * twt_pp_max / 4 + 1000 * 1000 * (offset_pp * offset_pp) / (vrms_pp * vrms_pp));
    double twtx_ss = std::sqrt(twt_ss_max * twt_ss_max / 4 + 1000 * 1000 * (offset_ss * offset_ss) / (vrms_ss * vrms_ss));
    twtx_max       = twtx_pp + twtx_ss;
  }
  else {  //------------PP seismic------------
    //find max vrms in index
    std::vector<double> vrms_vec(nzrefl), vp_vec(nzrefl), twt_vec(nzrefl), dummy;
    for (size_t k = 0; k < nzrefl; ++k) {
      twt_vec[k] = (*twtgrid_)(i_max, j_max, k);
      vp_vec[k] = (*vpgrid_) (i_max, j_max, k);
    }
    FindVrms(vrms_vec, dummy, twt_vec, dummy, vp_vec, 1.0, i_max, j_max, false);
    double vrms_max_t = vrms_vec[vrms_vec.size() - 1];

    //find max offset
    double offset_max = offset_0_ + doffset_*(noffset_ - 1);
    twtx_max = std::sqrt(twt_max*twt_max + 1000 * 1000 * offset_max*offset_max / (vrms_max_t*vrms_max_t));
  } //---------------------------

  double sf            = twtx_max / twt_max;                   //stretch factor. NO wavelet in twtx_max and twt_max
  double tmin          = t0;                                   //min twt sample. Includes wavlet
  size_t nt_top        = 0;                                    //samples on top due to stretch
  double tmax_stretch  = seismic_geometry_->tmax();            //max twt sample due to stretch. Includes wavelet
  nt_stretch           = nt;                                   //samples in nmo corrected seismic - include stretch top and bot
  size_t nt_seis       = nt;                                   //number of samples in prenmo seis (twt_0)
  twtx_max            += twt_wavelet;                          //add one wavelet (as no wavelet is included here)

  if (sf > 1) {
    tmin         -= sf * 2 * twt_wavelet;
    nt_top        = static_cast<size_t>((t0 - tmin) / dt);
    tmax_stretch += sf * 6 * twt_wavelet;
    nt_stretch    = static_cast<size_t>(floor((tmax_stretch - tmin) / dt + 0.5));
    nt_seis       = static_cast<size_t>(floor((twtx_max     - tmin) / dt + 0.5));
  }
 
  twt_0_.resize(nt_seis);
  for (size_t i = 0; i < nt_seis; ++i){
    twt_0_[i] = (t0 - nt_top * dt) + i*dt;
  }
  if (nt_stretch > twt_0_.size()){
    nt_stretch = twt_0_.size();
  }
  return twt_0_;
}

std::vector<double>  SeismicParameters::GenerateZ0ForNMO(){
  size_t nz                      = seismic_geometry_->nz();
  double zmin                    = seismic_geometry_->z0();
  double dz                      = seismic_geometry_->dz();

  std::vector<double> constvp    = model_settings_->GetConstVp();
  double twt_wavelet             = wavelet_->GetTwtWavelet();
  double z_top_wavelet           = twt_wavelet * constvp[0] / 2000;
  double z_bot_wavelet           = twt_wavelet * constvp[2] / 2000;

  double tmax                    = seismic_geometry_->tmax();
  double twt_0_max               = twt_0_[twt_0_.size()-1];
  double sf                      = (twt_0_max - twt_wavelet) / (tmax - twt_wavelet);
  double min_z                   = zmin;
  double max_z                   = zmin + (nz-1)*dz;

  //account for stretch above and below reservoir
  if (sf > 1) {
    min_z -= sf * 2 * z_top_wavelet;
    max_z += sf * 2 * z_bot_wavelet;
  }

  size_t nz_seis                 = static_cast<size_t>(std::ceil((max_z - min_z)/dz));

  //find min z in a sample
  size_t nzmin = static_cast<size_t>(floor(min_z) / dz + 0.5);
  double min_z_sampl = nzmin *dz;

  z_0_.resize(nz_seis);
  for (size_t i = 0; i < nz_seis; ++i){
    z_0_[i] = min_z_sampl + i*dz;
  }
  return z_0_;
}


std::vector<double>  SeismicParameters::GenerateTWT0Shift(double twt_0_min,
                                                          size_t n_samples)
{
  size_t i_max, j_max;
  double max_twt_value;

  //find max twt value
  FindMaxTwtIndex(i_max, j_max, max_twt_value);

  size_t k_max  = (*twt_timeshift_).GetNK() - 1;
  double ts_0   = (*twt_timeshift_)(i_max, j_max, 0);
  double ts_max = (*twt_timeshift_)(i_max, j_max, k_max);
  k_max         = (*twtgrid_).GetNK() - 1;
  double t_0    = (*twtgrid_)(i_max, j_max, 0);
  double t_max  = (*twtgrid_)(i_max, j_max, k_max);

  double dt                      = seismic_geometry_->dt();
  std::vector<double> twt_0_s;

  double delta_top = ts_0 - t_0;
  double delta_bot = ts_max - t_max;
  
  size_t n_samples_top = 0;
  size_t n_samples_bot = 0;
  if (delta_top < 0) { //shift upwards...include more samples in the top
    n_samples_top = static_cast<size_t>(std::ceil((-1*delta_top)/dt));
  }
  if (delta_bot > 0) { //shift downwards...include more samples below
    n_samples_bot = static_cast<size_t>(std::ceil(delta_bot/dt));
  }

  size_t n_samples_tot = n_samples_bot + n_samples + n_samples_top;
  double twts_min = twt_0_min - n_samples_top *dt;

  twt_0_s.resize(n_samples_tot);

  for (size_t k = 0; k < n_samples_tot; ++k){
    twt_0_s[k] = twts_min + k*dt;
  }
  return twt_0_s;
}


void SeismicParameters::FindPSNMOThetaAndOffset(NRLib::Grid2D<double>     &thetagrid,
                                                NRLib::Grid2D<double>     &offset_down_grid,
                                                NRLib::Grid2D<double>     &offset_up_grid,
                                                const std::vector<double> &twt_pp_vec,
                                                const std::vector<double> &twt_ss_vec,
                                                const std::vector<double> &vrms_pp_vec,
                                                const std::vector<double> &vrms_ss_vec,
                                                const std::vector<double> &offset,
                                                bool                       save_theta)
{
  double tol = 0.000001;
  size_t n_it = 10;
  size_t n_it_avg = 0;
  double theta_up, theta_down, y_out, dU, dD, vr;
  double start_value = 0.1;

  for (size_t off = 0; off < offset.size(); off++) {
    double tmp = offset[off] / (vrms_pp_vec[0]*twt_pp_vec[0] / 1000);
    start_value = atan(tmp);
    if (start_value >= 1.0)
      start_value = 0.99;
    for (size_t k = 0; k < twt_pp_vec.size(); k++) {
      dU = vrms_ss_vec[k] * twt_ss_vec[k] / 2000;
      dD = vrms_pp_vec[k] * twt_pp_vec[k] / 2000;
      vr = vrms_ss_vec[k] / vrms_pp_vec[k];
      n_it = 20;
      y_out = FindSinThetaPSWithNewtonsMethod(start_value,
                                              offset[off],
                                              dU,
                                              dD,
                                              vr,
                                              tol,
                                              n_it);
      n_it_avg  += n_it;
      theta_up   = asin(vr*y_out);
      theta_down = asin(y_out);
      if (save_theta) {
        thetagrid(k, off) = theta_down;
      }      
      offset_down_grid(k, off) = tan(theta_down)*dD;
      offset_up_grid(k, off)   = tan(theta_up)  *dU;
      start_value = y_out;
    }
  }
  size_t n_values = offset.size() * twt_pp_vec.size();
  double it_avg = static_cast<double>(n_it_avg) / static_cast<double>(n_values);
}


double SeismicParameters::FindSinThetaPSWithNewtonsMethod(double start_value,
                                                          double offset,
                                                          double dU,
                                                          double dD,
                                                          double vr,
                                                          double tol,
                                                          size_t &n_it)
{
  double y_old = start_value;
  double y_new, f_y, f_der_y;

  for (size_t i = 0; i < n_it; ++i) {
    f_y = -offset + dD*y_old / sqrt(1.0 - std::pow(y_old, 2)) + dU*vr*y_old / sqrt(1.0 - std::pow(vr, 2) * std::pow(y_old, 2));
    f_der_y = dD / (std::pow((1.0 - y_old), (3 / 2))) + dU*vr / (std::pow((1.0 - std::pow(vr, 2) * std::pow(y_old, 2)), (3 / 2)));

    if (f_der_y == 0) {
      std::cout << "failure in newtons method: zero derivative.";
      return 0;
    }
    y_new = y_old - f_y / f_der_y;

    if (std::abs(y_new) > 1.0) {
      std::cout << "failure in newtons method: Value > 1.0 : y_old = " << y_old << ", y_new = " << y_new << ". New value: y_new = 0.1 suggested.\n";
      y_new = 0.1;
    }

    if (std::abs(y_new - y_old) < tol) {
      n_it = i + 1;
      return y_new;
    }
    y_old = y_new;
  }
  return y_new;

}



void SeismicParameters::ReadEclipseGrid() {
    std::string filename = model_settings_->GetEclipseFileName();

    printf("Start reading Eclipsegrid from file.\n");
    eclipse_grid_ = new NRLib::EclipseGrid(filename);
    printf("Eclipsegrid read.\n");

    std::vector<std::string> names = model_settings_->GetParameterNames();
    if (!eclipse_grid_->HasParameter(names[0])) {
        std::cout << "Parameter " + names[0] + " is not found in Eclipse grid\n";
        exit(0);
    }
    if (!eclipse_grid_->HasParameter(names[1])) {
        std::cout << "Parameter " + names[1] + " is not found in Eclipse grid\n";
        exit(0);
    }
    if (!eclipse_grid_->HasParameter(names[2])) {
        std::cout << "Parameter " + names[2] + " is not found in Eclipse grid\n";
        exit(0);
    }
    std::vector<std::string> extra_parameter_names = model_settings_->GetExtraParameterNames();
    for (size_t i = 0; i < extra_parameter_names.size(); ++i) {
        if (!eclipse_grid_->HasParameter(extra_parameter_names[i])) {
            std::cout << "Parameter " + extra_parameter_names[i] + " is not found in Eclipse grid\n";
            exit(0);
        }
    }
}


void SeismicParameters::DeleteEclipseGrid() {
  delete eclipse_grid_;
}

void SeismicParameters::DeleteElasticParameterGrids() {
  delete vpgrid_;
  delete vsgrid_;
  delete rhogrid_;
}

void SeismicParameters::DeleteExtraParameterGrids() {
  delete extra_parameter_grid_;
}

void SeismicParameters::DeleteZandRandTWTGrids() {
  if (twtgrid_ != NULL)
    delete twtgrid_;
  if (twtssgrid_ != NULL)
    delete twtssgrid_;
  if (twtppgrid_ != NULL)
    delete twtppgrid_;
  delete zgrid_;
  delete rgridvec_;
  if (twt_timeshift_ != NULL)
    delete twt_timeshift_;
}

void SeismicParameters::DeleteVrmsGrid() {
  delete vrmsgrid_;
}

void SeismicParameters::DeleteWavelet() {
  delete wavelet_;
}

void SeismicParameters::DeleteGeometryAndOutput() {
  delete seismic_geometry_;
  delete segy_geometry_;
  delete seismic_output_;
  delete model_settings_;
}




void SeismicParameters::FindGeometry() {
  seismic_geometry_->setDxDy(model_settings_->GetDx(), model_settings_->GetDy());
  seismic_geometry_->setDz(model_settings_->GetDz());
  seismic_geometry_->setDt(model_settings_->GetDt());

  const NRLib::EclipseGeometry &geometry = eclipse_grid_->GetGeometry();


  if (model_settings_->GetAreaFromSegy() != "") {
    std::cout << "Area from <area-from-segy>.\n";
    int scalcoloc = 71;
    NRLib::TraceHeaderFormat::coordSys_t coord = NRLib::TraceHeaderFormat::UTM;
    NRLib::TraceHeaderFormat *thf = new NRLib::TraceHeaderFormat(scalcoloc, model_settings_->GetUtmxIn(), model_settings_->GetUtmyIn(), model_settings_->GetIL0In(), model_settings_->GetXL0In(), coord);
    double z0 = 0.0;
    NRLib::Volume *volume = NULL;
    std::vector<NRLib::TraceHeaderFormat *> thfvec;
    thfvec.push_back(thf);

    NRLib::SegY segy(model_settings_->GetAreaFromSegy(), static_cast<float>(z0), thfvec);
    segy.ReadAllTraces(volume, z0);
    segy.CreateRegularGrid();
    const NRLib::SegyGeometry *temp_segy_geometry = segy.GetGeometry();

    segy_geometry_ = new NRLib::SegyGeometry(temp_segy_geometry);
    segy_geometry_->WriteGeometry();
    segy_geometry_->WriteILXL();

    double x0    = temp_segy_geometry->GetX0();
    double y0    = temp_segy_geometry->GetY0();
    double lx    = temp_segy_geometry->Getlx();
    double ly    = temp_segy_geometry->Getly();
    double angle = temp_segy_geometry->GetAngle();
    double dx    = temp_segy_geometry->GetDx();
    double dy    = temp_segy_geometry->GetDy();

    seismic_geometry_->setGeometry(x0, y0, lx, ly, angle);
    seismic_geometry_->setDxDy(dx, dy);
  }
  else if (model_settings_->GetAreaGiven()) {
    std::cout << "Area from <area>.\n";
    double x0 = model_settings_->GetX0();
    double y0 = model_settings_->GetY0();
    double lx = model_settings_->GetLx();
    double ly = model_settings_->GetLy();
    double angle = model_settings_->GetAngle();
    seismic_geometry_->setGeometry(x0, y0, lx, ly, angle);
  }
  else if (model_settings_->GetAreaFromSurface() != "") {
    std::cout << "Area from <area-from-surface>.\n";
    NRLib::RegularSurfaceRotated<double> toptime_rot = NRLib::RegularSurfaceRotated<double>(model_settings_->GetAreaFromSurface());
    double x0 = toptime_rot.GetXRef();
    double y0 = toptime_rot.GetYRef();
    double lx = toptime_rot.GetLengthX();
    double ly = toptime_rot.GetLengthY();
    double angle = toptime_rot.GetAngle();
    seismic_geometry_->setGeometry(x0, y0, lx, ly, angle);
  }
  else {
    std::cout << "Area from Eclipsegrid.\n";
    double x0, y0, lx, ly, angle;
    geometry.FindEnclosingVolume(x0, y0, lx, ly, angle);
    seismic_geometry_->setGeometry(x0, y0, lx, ly, angle);
  }
}

void SeismicParameters::FindSurfaceGeometry() {
  const NRLib::EclipseGeometry &geometry = eclipse_grid_->GetGeometry();

  double dx = seismic_geometry_->dx();
  double dy = seismic_geometry_->dy();

  double lxsurf = seismic_geometry_->xsurfacelength();
  double lysurf = seismic_geometry_->ysurfacelength();

  double xmin = seismic_geometry_->xmin();
  double ymin = seismic_geometry_->ymin();

  size_t nxsurfec = seismic_geometry_->nxsurfaceeclipse();
  size_t nysurfec = seismic_geometry_->nysurfaceeclipse();

  bool const_top_given = true;
  if (model_settings_->GetTopTimeSurfaceFile() != "") {
    NRLib::RegularSurfaceRotated<double> top_time_rotated = NRLib::RegularSurfaceRotated<double>(model_settings_->GetTopTimeSurfaceFile());
    double topmin = top_time_rotated.Min();
    top_time_ = NRLib::RegularSurface<double>(xmin - dx, ymin - dy, lxsurf + 2 * dx, lysurf + 2 * dy, nxsurfec + 2, nysurfec + 2, topmin);
    top_time_.SetMissingValue(top_time_rotated.GetMissingValue());
    for (size_t i = 0; i < top_time_.GetNI(); i++) {
      for (size_t j = 0; j < top_time_.GetNJ(); j++) {
        double x, y;
        top_time_.GetXY(i, j, x, y);
        double value = top_time_rotated.GetZ(x, y);
        top_time_(i, j) = value;
      }
    }

    bot_time_ = NRLib::RegularSurface<double>(xmin - dx, ymin - dy, lxsurf + 2 * dx, lysurf + 2 * dy, nxsurfec + 2, nysurfec + 2, top_time_.Max());
    const_top_given = false;
  } 
  else {
    double t1 = model_settings_->GetTopTimeConstant();
    top_time_ = NRLib::RegularSurface<double>(xmin - dx, ymin - dy, lxsurf + 2 * dx, lysurf + 2 * dy, nxsurfec + 2, nysurfec + 2, t1);
    bot_time_ = NRLib::RegularSurface<double>(xmin - dx, ymin - dy, lxsurf + 2 * dx, lysurf + 2 * dy, nxsurfec + 2, nysurfec + 2, t1);
  }

  topeclipse_ = NRLib::RegularSurface<double>(xmin - dx, ymin - dy, lxsurf + 2 * dx, lysurf + 2 * dy, nxsurfec + 2, nysurfec + 2, -999.0);
  boteclipse_ = NRLib::RegularSurface<double>(xmin - dx, ymin - dy, lxsurf + 2 * dx, lysurf + 2 * dy, nxsurfec + 2, nysurfec + 2, -999.0);

  top_k_    = geometry.FindTopLayer();
  bottom_k_ = geometry.FindBottomLayer();

  seismic_geometry_->setZReflectorCount(static_cast<size_t>(bottom_k_ + 2 - top_k_));

  NRLib::Grid2D<double> values(nxsurfec + 2, nysurfec + 2, 0.0);
  if (model_settings_->GetUseCornerpointInterpol()) {
    geometry.FindLayerSurfaceCornerpoint(values, top_k_, 0, topeclipse_.GetDX(), topeclipse_.GetDY(), xmin - dx, ymin - dy, 0.0, 0);
  }
  else {
    geometry.FindLayerSurface(values, top_k_, 0, topeclipse_.GetDX(), topeclipse_.GetDY(), xmin - dx, ymin - dy, 0.0, 0);
  }

  for (size_t i = 0; i < topeclipse_.GetNI(); i++) {
    for (size_t j = 0; j < topeclipse_.GetNJ(); j++) {
      topeclipse_(i, j) = values(i, j);
    }
  }

  if (model_settings_->GetUseCornerpointInterpol()) {
    geometry.FindLayerSurfaceCornerpoint(values, bottom_k_, 1, boteclipse_.GetDX(), boteclipse_.GetDY(), xmin - dx, ymin - dy, 0.0, 0);
  }
  else {
    geometry.FindLayerSurface(values, bottom_k_, 1, boteclipse_.GetDX(), boteclipse_.GetDY(), xmin - dx, ymin - dy, 0.0, 0);
  }

  for (size_t i = 0; i < boteclipse_.GetNI(); i++) {
    for (size_t j = 0; j < boteclipse_.GetNJ(); j++) {
      boteclipse_(i, j) = values(i, j);
    }
  }

  if (model_settings_->GetOutputDepthSurfaces()) {
    seismic_output_->WriteDepthSurfaces(topeclipse_, boteclipse_);
  }


  double min, max, d1 = 1e10, d2 = 0;
  bool found;
  for (size_t i = 0; i < geometry.GetNI(); ++i) {
    for (size_t j = 0; j < geometry.GetNJ(); ++j) {
      found = false;
      min = geometry.FindZTopInCellActiveColumn(i, j, top_k_, found);
      if (found && min < d1) {
        d1 = min;
      }
      found = false;
      max = geometry.FindZBotInCellActiveColumn(i, j, bottom_k_, found);
      if (found && max > d2) {
        d2 = max;
      }
    }
  }

  if (const_top_given) {
    double const_v = model_settings_->GetConstVp()[0];
    if (model_settings_->GetPSSeismic()) {
      double const_vs = model_settings_->GetConstVs()[0];
      const_v = 2 / (1 / const_v + 1 / const_vs);
    }
    double t1 = model_settings_->GetTopTimeConstant();
    for (size_t i = 0; i < top_time_.GetNI(); i++)
      for (size_t j = 0; j < top_time_.GetNJ(); j++) {
        top_time_(i, j) = t1 + 2000.0 * (topeclipse_(i, j) - d1) / const_v;
        bot_time_(i, j) = top_time_(i, j);
      }
  }
  double twt_wavelet   = wavelet_->GetTwtWavelet();
  std::vector<double> constvp = model_settings_->GetConstVp();
  double z_top_wavelet = twt_wavelet * constvp[0] / 2000;
  double z_bot_wavelet = twt_wavelet * constvp[2] / 2000;
  topeclipse_.Add(-1 * z_top_wavelet); // add one wavelet length to bot and subtract from top
  boteclipse_.Add(     z_bot_wavelet);

  d1 -= twt_wavelet;
  d2 += twt_wavelet;

  seismic_geometry_->setZRange(d1, d2);
}

void SeismicParameters::CreateGrids() {
  size_t nx = seismic_geometry_->nx();
  size_t ny = seismic_geometry_->ny();
  size_t nzrefl = seismic_geometry_->zreflectorcount();

  NRLib::Volume volume = seismic_geometry_->createDepthVolume();

  zgrid_    = new NRLib::StormContGrid(volume, nx, ny, nzrefl,     0.0);
  vpgrid_   = new NRLib::StormContGrid(volume, nx, ny, nzrefl + 1, missing_);
  vsgrid_   = new NRLib::StormContGrid(volume, nx, ny, nzrefl + 1, missing_);
  rhogrid_  = new NRLib::StormContGrid(volume, nx, ny, nzrefl + 1, missing_);
  twtgrid_  = new NRLib::StormContGrid(volume, nx, ny, nzrefl,     0.0);
  
  if (model_settings_->GetNMOCorr() && model_settings_->GetPSSeismic()) {
    twtssgrid_ = new NRLib::StormContGrid(volume, nx, ny, nzrefl, 0.0);
    twtppgrid_ = new NRLib::StormContGrid(volume, nx, ny, nzrefl, 0.0);
  }
  else {
    twtssgrid_ = NULL;
    twtppgrid_ = NULL;
  }

  if (model_settings_->GetNMOCorr() && model_settings_->GetOutputVrms()) {
    vrmsgrid_ = new NRLib::StormContGrid(volume, nx, ny, nzrefl, 0.0); //dimensions??
  }
  else { 
    vrmsgrid_ = NULL;
  }

  if (model_settings_->GetOutputReflections()){
    if (model_settings_->GetWhiteNoise())
      rgridvec_ = new std::vector<NRLib::StormContGrid>(2);
    else
      rgridvec_ = new std::vector<NRLib::StormContGrid>(1);
  }
  else {
    rgridvec_ = NULL;
  }

  NRLib::StormContGrid rgrid(volume, nx, ny, nzrefl, 0.0);

  if (model_settings_->GetOutputReflections()) {
    (*rgridvec_)[0] = rgrid;
    if (model_settings_->GetWhiteNoise()) {
      (*rgridvec_)[1] = rgrid;
    }
  }

  std::vector<std::string> extra_parameter_names = model_settings_->GetExtraParameterNames();
  std::vector<double> extra_parameter_default_values = model_settings_->GetExtraParameterDefaultValues();
  extra_parameter_grid_ = new std::vector<NRLib::StormContGrid>(extra_parameter_names.size());
  for (size_t i = 0; i < extra_parameter_names.size(); ++i) {
    (*extra_parameter_grid_)[i] = NRLib::StormContGrid(volume, nx, ny, nzrefl + 1, extra_parameter_default_values[i]);
  }

  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      for (size_t epi = 0; epi < extra_parameter_names.size(); ++epi) {
        (*extra_parameter_grid_)[epi](i, j, nzrefl) = 0.0;
      }
    }
  }
  
  if (model_settings_->GetTwtFileName() != "") {
    twt_timeshift_ = new NRLib::StormContGrid(model_settings_->GetTwtFileName());
    if ((*twt_timeshift_).GetNI() != nx || (*twt_timeshift_).GetNJ() != ny || (*twt_timeshift_).GetNK() != nzrefl) {
      printf("TWT timeshift from file has wrong dimension. Aborting. \n");
      exit(1);
    }
  }
  else {
      twt_timeshift_ = NULL;
  }
}

void SeismicParameters::PrintElapsedTime(time_t start_time, std::string work)
{
  time_t end_time = time(0);   // get time now
  size_t seconds = static_cast<size_t>(difftime(end_time, start_time));

  size_t hours = static_cast<int>(seconds / 3600);
  seconds = seconds % 3600;
  size_t minutes = static_cast<int>(seconds / 60);
  seconds = seconds % 60;

  std::string hours_s = NRLib::ToString(hours);
  std::string zeros = std::string(2 - hours_s.length(), '0');
  hours_s = zeros + hours_s;
  std::string minutes_s = NRLib::ToString(minutes);
  zeros = std::string(2 - minutes_s.length(), '0');
  minutes_s = zeros + minutes_s;
  std::string seconds_s = NRLib::ToString(seconds);
  zeros = std::string(2 - seconds_s.length(), '0');
  seconds_s = zeros + seconds_s;

  std::cout << "Time " << work << ": "
    << hours_s << ':'
    << minutes_s << ':'
    << seconds_s
    << "\n";
}


