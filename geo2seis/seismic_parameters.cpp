#include <seismic_parameters.hpp>

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


SeismicParameters::SeismicParameters(ModelSettings *model_settings) {
    this->model_settings_ = model_settings;

    seismic_geometry_ = new SeismicGeometry();

    if (model_settings->GetNMOCorr()) {
      calculateOffsetSpan();
    }
    else {
      calculateAngleSpan();
    }

    segy_geometry_ = NULL;

    setupWavelet();
    readEclipseGrid();
    findGeometry();

    seismic_output_ = new SeismicOutput(model_settings_);
    findSurfaceGeometry();

    createGrids();
}

void SeismicParameters::calculateAngleSpan() {
    theta_0_ = model_settings_->GetTheta0();
    dtheta_ = model_settings_->GetDTheta();
    theta_max_ = model_settings_->GetThetaMax();

    if (dtheta_ == 0) {
        ntheta_ = 1;
    } else {
        ntheta_ = size_t((theta_max_ - theta_0_) / dtheta_) + 1;
        dtheta_ = (theta_max_ - theta_0_) / (ntheta_ - 1);
    }
    theta_vec_.resize(ntheta_);
    for (size_t i = 0; i < ntheta_; ++i) {
      theta_vec_[i] = theta_0_ + i*dtheta_;
    }
}


void SeismicParameters::calculateOffsetSpan() {
    offset_0_ = model_settings_->GetOffset0();
    doffset_ = model_settings_->GetDOffset();
    offset_max_ = model_settings_->GetOffsetMax();

    if (doffset_ == 0) {
        noffset_ = 1;
    } else {
        noffset_ = size_t((offset_max_ - offset_0_) / doffset_) + 1;
        doffset_ = (offset_max_ - offset_0_) / (noffset_ - 1);
    }
    offset_vec_.resize(noffset_);
    for (size_t i = 0; i < noffset_; ++i) {
      offset_vec_[i] = offset_0_ + i*doffset_;
    }
}

void SeismicParameters::setupWavelet() {
    if (model_settings_->GetRicker()) {
        double peakF = model_settings_->GetPeakFrequency();
        wavelet_ = new Wavelet(peakF);
    } else {
        std::string wavelet_file_format = model_settings_->GetWaveletFileFormat();
        std::string wavelet_file_name = model_settings_->GetWaveletFileName();
        wavelet_ = new Wavelet(wavelet_file_name, wavelet_file_format);
    }
    wavelet_scale_ = model_settings_->GetWaveletScale();
}

void SeismicParameters::setSegyGeometry(const NRLib::SegyGeometry &geometry)
{
  if (segy_geometry_ != NULL) {
    delete segy_geometry_;
  }
  segy_geometry_ = new NRLib::SegyGeometry(geometry);
}

void SeismicParameters::getSeisLimits(size_t               n_twt_0,
                                      std::vector<double>  vrms_vec,
                                      std::vector<double>  offset_vec,
                                      std::vector<size_t> &n_min,
                                      std::vector<size_t> &n_max)
{
  size_t nt                      = seismic_geometry_->nt();
  double dt                      = seismic_geometry_->dt();
  double tmin                    = seismic_geometry_->t0();
  double tmax                    = seismic_geometry_->tmax();
  double vrms_nk = vrms_vec[vrms_vec.size() - 1];
  double vrms_0  = vrms_vec[0];

  for (size_t i = 0; i < offset_vec.size(); ++i) {
    double offset = offset_vec[i];
    double twtx_max                = std::sqrt(tmax*tmax + 1000*1000*offset*offset/(vrms_nk*vrms_nk));
    double twtx_min                = std::sqrt(tmin*tmin + 1000*1000*offset*offset/(vrms_0*vrms_0));
    if (twtx_min > tmin){
      double test = (twtx_min - tmin)/dt;
      n_min[i] = static_cast<size_t>((twtx_min - tmin)/dt); //static_cast cuts decimal. Ok, because target is (x - 0.5).
    }
    double test = (twtx_max - tmin)/dt;
    n_max[i] = std::min(static_cast<size_t>((twtx_max - tmin)/dt),n_twt_0-1);
  }
}


void SeismicParameters::findVrmsPos(std::vector<double>       &vrms_vec,
                                    std::vector<double>       &vrms_vec_reg,
                                    const std::vector<double> &twt_0,
                                    size_t                    i,
                                    size_t                    j,
                                    bool                      include_regular)
{
  double v_w = model_settings_->GetVw();
  double z_w = model_settings_->GetZw();

  double v_over;
  double twt_w = 2000*z_w/v_w;
  double tmp, tmp0;
  size_t nk = (*twtgrid_).GetNK();

  //generate vrms_vec - only for each layer in reservoir
  if ((*twtgrid_)(i,j,0) == -999.0) {
    for (size_t k = 0; k < (*twtgrid_).GetNK(); ++k) {
      vrms_vec[k] = -999.0;
    }
  }
  else {
    v_over = 2000*((*zgrid_)(i,j,0) - z_w)/((*twtgrid_)(i,j,0) - 2000*z_w/v_w);
    tmp0 = v_w*v_w*twt_w + v_over*v_over*((*twtgrid_)(i,j,0) - twt_w);
    for (size_t k = 0; k < nk; ++k) {
      tmp = tmp0;
      for (size_t l = 1; l <= k; ++l) {
        tmp += (*vpgrid_)(i,j,l)* (*vpgrid_)(i,j,l)*((*twtgrid_)(i,j,l) - (*twtgrid_)(i,j,l-1));
      }
      vrms_vec[k] = std::sqrt(tmp / (*twtgrid_)(i,j,k));
    }
  }
  
  if (include_regular) {
    //generate vrms_vec_reg - including overburden and underburden
    std::vector<double> constvp    = model_settings_->GetConstVp();
    double twt_wavelet = 2000 / constvp[2] * wavelet_->GetDepthAdjustmentFactor();  
    double vrms_under  = vrms_vec[nk-1]*vrms_vec[nk-1]*(*twtgrid_)(i,j,(nk-1)) + constvp[2]*constvp[2]*twt_wavelet;
    vrms_under        *= (1/((*twtgrid_)(i,j,(nk-1))+twt_wavelet));
    vrms_under        = std::sqrt(vrms_under);

    //sample vrms regularly:
    std::vector<double> twt_vec_in(nk+2), vrms_vec_in(nk+2);
    twt_vec_in[0]  = twt_w;
    vrms_vec_in[0] = v_w;
    twt_vec_in[1]  = (*twtgrid_)(i,j,0);
    vrms_vec_in[1] = vrms_vec[0];
    size_t index   = 2;
    for (size_t k = 0; k < nk; ++k) {
      if ((*twtgrid_)(i,j,k) != twt_vec_in[index-1]) {
        twt_vec_in[index]  = (*twtgrid_)(i,j,k);
        vrms_vec_in[index] = vrms_vec[k];
        ++index;
      }
    }
    twt_vec_in.resize(index+1);
    vrms_vec_in.resize(index+1);
    twt_vec_in[index]  = twt_vec_in[index-1]+twt_wavelet;
    vrms_vec_in[index] = vrms_under;

    vrms_vec_reg = NRLib::Interpolation::Interpolate1D(twt_vec_in, vrms_vec_in, twt_0, "linear");
  }
}



void SeismicParameters::findReflections(NRLib::Grid2D<double>       &r_vec,
                                        const std::vector<double>   &theta_vec,
                                        size_t                       i,
                                        size_t                       j){



  bool ps_seis   = model_settings_->GetPSSeismic();

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
      diffvp  =         vp_vec[k-top_k_ + 1] - vp_vec[k-top_k_];
      meanvp  = 0.5 *  (vp_vec[k-top_k_ + 1] + vp_vec[k-top_k_]);
      diffvs  =         vs_vec[k-top_k_ + 1] - vs_vec[k-top_k_];
      meanvs  = 0.5 *  (vs_vec[k-top_k_ + 1] + vs_vec[k-top_k_]);
      diffrho =        rho_vec[k-top_k_ + 1] - rho_vec[k-top_k_];
      meanrho = 0.5 * (rho_vec[k-top_k_ + 1] + rho_vec[k-top_k_]);
      r_vec(k - top_k_, theta) = static_cast<float>(zoeppritz->GetReflection(diffvp, meanvp, diffrho, meanrho, diffvs, meanvs));
    }
  }
  delete zoeppritz;
}


void SeismicParameters::findNMOReflections(NRLib::Grid2D<double>       &r_vec,
                                           const NRLib::Grid2D<double> &theta_vec,
                                           const std::vector<double>   &offset_vec,
                                           size_t                       i,
                                           size_t                       j){

  bool ps_seis   = model_settings_->GetPSSeismic();

  Zoeppritz *zoeppritz = NULL;
  if (ps_seis) {
    zoeppritz = new ZoeppritzPS();
  } else {
    zoeppritz = new ZoeppritzPP();
  }

  double diffvp, meanvp, diffvs, meanvs, diffrho, meanrho;
  std::vector<double> vp_vec(bottom_k_-top_k_+3), vs_vec(bottom_k_-top_k_+3), rho_vec(bottom_k_-top_k_+3);

  for (size_t off = 0; off < offset_vec.size(); ++off) {
    for (size_t k = top_k_; k <= bottom_k_ + 2; k++) {
      vp_vec[k - top_k_]  = (*vpgrid_)(i, j, (k - top_k_));
      vs_vec[k - top_k_]  = (*vsgrid_)(i, j, (k - top_k_));
      rho_vec[k - top_k_] = (*rhogrid_)(i, j, (k - top_k_));
    }
    for (size_t k = top_k_; k <= bottom_k_ + 1; k++) {
      diffvp  =         vp_vec[k-top_k_ + 1] - vp_vec[k-top_k_];
      meanvp  = 0.5 *  (vp_vec[k-top_k_ + 1] + vp_vec[k-top_k_]);
      diffvs  =         vs_vec[k-top_k_ + 1] - vs_vec[k-top_k_];
      meanvs  = 0.5 *  (vs_vec[k-top_k_ + 1] + vs_vec[k-top_k_]);
      diffrho =        rho_vec[k-top_k_ + 1] - rho_vec[k-top_k_];
      meanrho = 0.5 * (rho_vec[k-top_k_ + 1] + rho_vec[k-top_k_]);
      zoeppritz->ComputeConstants(theta_vec(k, off));
      r_vec(k - top_k_, off) = static_cast<float>(zoeppritz->GetReflection(diffvp, meanvp, diffrho, meanrho, diffvs, meanvs));
    }
  }
  delete zoeppritz;
}


std::vector<double> SeismicParameters::generateTWT_0(){
  size_t i_max, j_max, k_max;
  double x, y;
  size_t nt                      = seismic_geometry_->nt();
  double tmin                    = seismic_geometry_->t0();
  double dt                      = seismic_geometry_->dt();
  std::vector<double> constvp    = model_settings_->GetConstVp();

  //find max TWT for highest offset in order to find the highest TWT value to sample seismic
  double max_twt_value           = bot_time_.MaxNode(i_max, j_max);
  bot_time_.GetXY(i_max, j_max, x, y);
  double bot_y_value_twt = bot_time_.GetZ(x,y) -  2000 / constvp[2] * wavelet_->GetDepthAdjustmentFactor();
  (*twtgrid_).FindIndex(x, y, bot_y_value_twt, i_max, j_max, k_max);
  std::vector<double> vrms_vec((*twtgrid_).GetNK()), vrms_dummy, twt_0_dummy;
  findVrmsPos(vrms_vec, vrms_dummy, twt_0_dummy, i_max, j_max, false);
  double vrms_max_t              = vrms_vec[vrms_vec.size() - 1];
  double offset_max              = offset_0_+doffset_*(noffset_ - 1);
  double twtx_max                = std::sqrt(max_twt_value*max_twt_value + 1000*1000*offset_max*offset_max/(vrms_max_t*vrms_max_t));
  size_t nt_seis                 = nt;
  if (twtx_max > tmin + nt*dt) {
    nt_seis = static_cast<size_t>(std::ceil((twtx_max - tmin)/dt));
  }
  twt_0_.resize(nt_seis);
  for (size_t i = 0; i < nt_seis; ++i){
    twt_0_[i] = tmin + (0.5 + i)*dt;
  }
  return twt_0_;
}

std::vector<double>  SeismicParameters::generateZ_0(){
  size_t nz                      = seismic_geometry_->nz();
  double zmin                    = seismic_geometry_->z0();
  double dz                      = seismic_geometry_->dz();

  double tmax                    = seismic_geometry_->tmax();
  double twt_0_max               = twt_0_[twt_0_.size()-1];
  double factor                  = 2 * twt_0_max / tmax;
  double max_z                   = zmin + (nz-1)*dz + factor * wavelet_->GetDepthAdjustmentFactor();

  size_t nz_seis                 = static_cast<size_t>(std::ceil((max_z - zmin)/dz));

  z_0_.resize(nz_seis);
  for (size_t i = 0; i < nz_seis; ++i){
    z_0_[i] = zmin + (0.5 + i)*dz;
  }
  return z_0_;
}

void SeismicParameters::readEclipseGrid() {
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


void SeismicParameters::deleteEclipseGrid() {
  delete eclipse_grid_;
}

void SeismicParameters::deleteElasticParameterGrids() {
  delete vpgrid_;
  delete vsgrid_;
  delete rhogrid_;  
}

void SeismicParameters::deleteExtraParameterGrids() {
  delete extra_parameter_grid_;
}

void SeismicParameters::deleteZandRandTWTGrids() {
  delete twtgrid_;
  delete zgrid_;
  delete rgridvec_;
}

void SeismicParameters::deleteVrmsGrid() {
  delete vrmsgrid_;
}

void SeismicParameters::deleteWavelet() {
  delete wavelet_;
}

void SeismicParameters::deleteGeometryAndOutput() {
  delete seismic_geometry_;
  delete segy_geometry_;
  delete seismic_output_;
  delete model_settings_;
}

bool SeismicParameters::GetTimeOutput() {
  bool output = (model_settings_->GetOutputSeismicTime()
              || model_settings_->GetOutputTimeSegy()
              || model_settings_->GetOutputSeismicStackTimeStorm()      
              || model_settings_->GetOutputSeismicStackTimeSegy()
              || model_settings_->GetOutputPrenmoTimeSegy());
  return output;
}

bool SeismicParameters::GetDepthOutput() {
  return (model_settings_->GetOutputSeismicDepth()               
       || model_settings_->GetOutputDepthSegy()
       || model_settings_->GetOutputSeismicStackDepthStorm()     
       || model_settings_->GetOutputSeismicStackDepthSegy());
}

bool SeismicParameters::GetTimeshiftOutput() {
  return (model_settings_->GetOutputSeismicTimeshift()           
       || model_settings_->GetOutputTimeshiftSegy()
       || model_settings_->GetOutputSeismicStackTimeShiftStorm() 
       || model_settings_->GetOutputSeismicStackTimeShiftSegy());
}
bool SeismicParameters::GetStackOutput() {
  return (model_settings_->GetOutputSeismicStackTimeStorm()      || model_settings_->GetOutputSeismicStackTimeSegy()
       || model_settings_->GetOutputSeismicStackTimeShiftStorm() || model_settings_->GetOutputSeismicStackTimeShiftSegy()
       || model_settings_->GetOutputSeismicStackDepthStorm()     || model_settings_->GetOutputSeismicStackDepthSegy());
}
bool SeismicParameters::GetSegyOutput() {
  return (model_settings_->GetOutputTimeSegy()
       || model_settings_->GetOutputSeismicStackTimeSegy()
       || model_settings_->GetOutputDepthSegy()
       || model_settings_->GetOutputSeismicStackDepthSegy()
       || model_settings_->GetOutputTimeshiftSegy()
       || model_settings_->GetOutputSeismicStackTimeShiftSegy()
       || model_settings_->GetOutputPrenmoTimeSegy());
}
bool SeismicParameters::GetTimeStormOutput() {
  return (model_settings_->GetOutputSeismicTime()      || model_settings_->GetOutputSeismicStackTimeStorm());
}
bool SeismicParameters::GetDepthStormOutput() {
  return (model_settings_->GetOutputSeismicDepth()     || model_settings_->GetOutputSeismicStackDepthStorm());
}
bool SeismicParameters::GetTimeshiftStormOutput() {
  return (model_settings_->GetOutputSeismicTimeshift() || model_settings_->GetOutputSeismicStackTimeShiftStorm());
}

bool SeismicParameters::GetStormOutput() {
  return (GetTimeStormOutput() || GetTimeshiftStormOutput() || GetDepthStormOutput());
}



void SeismicParameters::findGeometry() {
  seismic_geometry_->setDxDy(model_settings_->GetDx(), model_settings_->GetDy());
  seismic_geometry_->setDz(model_settings_->GetDz());
  seismic_geometry_->setDt(model_settings_->GetDt());

  const NRLib::EclipseGeometry &geometry = eclipse_grid_->GetGeometry();

  if (model_settings_->GetAreaGiven()) {
    double x0 = model_settings_->GetX0();
    double y0 = model_settings_->GetY0();
    double lx = model_settings_->GetLx();
    double ly = model_settings_->GetLy();
    double angle = model_settings_->GetAngle();
    seismic_geometry_->setGeometry(x0, y0, lx, ly, angle);

  }
  else if (model_settings_->GetAreaFromSurface() != "") {
    NRLib::RegularSurfaceRotated<double> toptime_rot = NRLib::RegularSurfaceRotated<double>(model_settings_->GetAreaFromSurface());
    double x0 = toptime_rot.GetXRef();
    double y0 = toptime_rot.GetYRef();
    double lx = toptime_rot.GetLengthX();
    double ly = toptime_rot.GetLengthY();
    double angle = toptime_rot.GetAngle();
    seismic_geometry_->setGeometry(x0, y0, lx, ly, angle);

  }
  else if (model_settings_->GetAreaFromSegy() != "") {
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
  else {
    double x0, y0, lx, ly, angle;
    geometry.FindEnclosingVolume(x0, y0, lx, ly, angle);
    seismic_geometry_->setGeometry(x0, y0, lx, ly, angle);
  }

  if (model_settings_->GetAreaGiven() && model_settings_->GetAreaFromSurface() != "") {
    printf("WARNING! Area defined in two different ways. The area specified by the area command is used.\n");
  }
}

void SeismicParameters::findSurfaceGeometry() {
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
  } else {
    double t1 = model_settings_->GetTopTimeConstant();
    top_time_ = NRLib::RegularSurface<double>(xmin - dx, ymin - dy, lxsurf + 2 * dx, lysurf + 2 * dy, nxsurfec + 2, nysurfec + 2, t1);
    bot_time_ = NRLib::RegularSurface<double>(xmin - dx, ymin - dy, lxsurf + 2 * dx, lysurf + 2 * dy, nxsurfec + 2, nysurfec + 2, t1);

  }

  topeclipse_ = NRLib::RegularSurface<double>(xmin - dx, ymin - dy, lxsurf + 2 * dx, lysurf + 2 * dy, nxsurfec + 2, nysurfec + 2, -999.0);
  boteclipse_ = NRLib::RegularSurface<double>(xmin - dx, ymin - dy, lxsurf + 2 * dx, lysurf + 2 * dy, nxsurfec + 2, nysurfec + 2, -999.0);

  top_k_ = geometry.FindTopLayer();
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
  } else {
    geometry.FindLayerSurface(values, bottom_k_, 1, boteclipse_.GetDX(), boteclipse_.GetDY(), xmin - dx, ymin - dy, 0.0, 0);
  }

  for (size_t i = 0; i < boteclipse_.GetNI(); i++) {
    for (size_t j = 0; j < boteclipse_.GetNJ(); j++) {
      boteclipse_(i, j) = values(i, j);
    }
  }

  if (model_settings_->GetOutputDepthSurfaces()) {
    seismic_output_->writeDepthSurfaces(topeclipse_, boteclipse_);
  }

  double d1 = topeclipse_.Min();
  double d2 = boteclipse_.Max();

  if (const_top_given) {
    double t1 = model_settings_->GetTopTimeConstant();
    std::vector<double> const_vp = model_settings_->GetConstVp();
    for (size_t i = 0; i < top_time_.GetNI(); i++)
      for (size_t j = 0; j < top_time_.GetNJ(); j++) {
        top_time_(i, j) = t1 + 2000.0 * (topeclipse_(i, j) - d1) / const_vp[0];
        bot_time_(i, j) = top_time_(i, j);
      }
  }

  topeclipse_.Add(-1 * wavelet_->GetDepthAdjustmentFactor()); // add one wavelet length to bot and subtract from top
  boteclipse_.Add(wavelet_->GetDepthAdjustmentFactor());
  d1 = topeclipse_.Min();
  d2 = boteclipse_.Max();

  seismic_geometry_->setZRange(d1, d2);
}

void SeismicParameters::createGrids() {
  size_t nx = seismic_geometry_->nx();
  size_t ny = seismic_geometry_->ny();
  size_t nzrefl = seismic_geometry_->zreflectorcount();

  NRLib::Volume volume = seismic_geometry_->createDepthVolume();

  zgrid_    = new NRLib::StormContGrid(volume, nx, ny, nzrefl);
  vpgrid_   = new NRLib::StormContGrid(volume, nx, ny, nzrefl + 1);
  vsgrid_   = new NRLib::StormContGrid(volume, nx, ny, nzrefl + 1);
  rhogrid_  = new NRLib::StormContGrid(volume, nx, ny, nzrefl + 1);
  twtgrid_  = new NRLib::StormContGrid(volume, nx, ny, nzrefl);
  if (model_settings_->GetNMOCorr()){
    vrmsgrid_ = new NRLib::StormContGrid(volume, nx, ny, nzrefl); //dimensions??
    if (model_settings_->GetWhiteNoise())
      rgridvec_ = new std::vector<NRLib::StormContGrid>(2);
    else
      rgridvec_ = new std::vector<NRLib::StormContGrid>(1);
  }
  else {
    rgridvec_ = new std::vector<NRLib::StormContGrid>(ntheta_);
    vrmsgrid_ = NULL;
  }
  NRLib::StormContGrid rgrid(volume, nx, ny, nzrefl);

  std::vector<std::string> extra_parameter_names = model_settings_->GetExtraParameterNames();

  std::vector<double> extra_parameter_default_values = model_settings_->GetExtraParameterDefaultValues();
  extra_parameter_grid_ = new std::vector<NRLib::StormContGrid>(extra_parameter_names.size());
  for (size_t i = 0; i < extra_parameter_names.size(); ++i) {
    (*extra_parameter_grid_)[i] = NRLib::StormContGrid(volume, nx, ny, nzrefl + 1);
  }

  std::vector<double> const_vp  = model_settings_->GetConstVp();
  std::vector<double> const_vs  = model_settings_->GetConstVs();
  std::vector<double> const_rho = model_settings_->GetConstRho();

  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      for (size_t k = 0; k < nzrefl; k++) {
        (*zgrid_)(i, j, k) = 0.0;
        (*vpgrid_)(i, j, k) = static_cast<float>(const_vp[1]);
        (*vsgrid_)(i, j, k) = static_cast<float>(const_vs[1]);
        (*rhogrid_)(i, j, k) = static_cast<float>(const_rho[1]);
        (*twtgrid_)(i, j, k) = 0.0;
        rgrid(i, j, k) = 0.0;
        if (model_settings_->GetNMOCorr()){
          (*vrmsgrid_)(i, j, k) = 0.0;
        }
        for (size_t epi = 0; epi < extra_parameter_names.size(); ++epi) {
          (*extra_parameter_grid_)[epi](i, j, k) = static_cast<float>(extra_parameter_default_values[epi]);
        }
      }

      (*vpgrid_)(i, j, nzrefl) = static_cast<float>(const_vp[2]);
      (*vsgrid_)(i, j, nzrefl) = static_cast<float>(const_vs[2]);
      (*rhogrid_)(i, j, nzrefl) = static_cast<float>(const_rho[2]);
      for (size_t epi = 0; epi < extra_parameter_names.size(); ++epi) {
        (*extra_parameter_grid_)[epi](i, j, nzrefl) = 0.0;
      }
    }
  }
  if (model_settings_->GetNMOCorr()){
    if (model_settings_->GetWhiteNoise()){
      (*rgridvec_)[0] = rgrid;
      (*rgridvec_)[1] = rgrid;
    }
    else {
      (*rgridvec_)[0] = rgrid;
    }
  }
  else {
    for (size_t i = 0; i < ntheta_; i++) {
      (*rgridvec_)[i] = rgrid;
    }
  }
}
