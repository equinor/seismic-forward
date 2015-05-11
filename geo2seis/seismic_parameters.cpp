#include <seismic_parameters.hpp>

#include <nrlib/math/constants.hpp>
#include <physics/wavelet.hpp>
#include <nrlib/eclipsegrid/eclipsegrid.hpp>
#include <nrlib/surface/regularsurfacerotated.hpp>
#include <nrlib/segy/segygeometry.hpp>
#include <nrlib/stormgrid/stormcontgrid.hpp>
#include "seismic_geometry.hpp"


SeismicParameters::SeismicParameters(ModelSettings *model_settings) {
    this->model_settings_ = model_settings;

    seismic_geometry_ = new SeismicGeometry();
    seismic_output_ = new SeismicOutput(model_settings_);

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

std::vector<double> SeismicParameters::twt_0(){
  size_t i_max, j_max, k_max;
  double x, y;
  size_t nt                      = seismic_geometry_->nt();
  double tmin                    = seismic_geometry_->t0();
  double dt                      = seismic_geometry_->dt();
  std::vector<double> constvp    = model_settings_->GetConstVp();

  
  double max_twt_value           = bot_time_.MaxNode(i_max, j_max);  
  bot_time_.GetXY(i_max, j_max, x, y);
  (*twtgrid_).FindIndex(x, y, bot_time_.GetZ(x,y), i_max, j_max, k_max);
  //double max_twt_value           = (*twtgrid_)( i_max, j_max, (*twtgrid_).GetNK()-1); //need half wavelet for NMO correction
  double vrms_max_t              = (*vrmsgrid_)(i_max, j_max, (*vrmsgrid_).GetNK()-1);  
  double offset_max              = offset_0_+doffset_*noffset_;

  double twtx_max                = std::sqrt(max_twt_value*max_twt_value + 2000*2000*offset_max*offset_max/vrms_max_t);      
  twtx_max                      += 2000 / constvp[2] * wavelet_->GetDepthAdjustmentFactor();
  
  size_t nt_seis                 = nt;
  if (twtx_max > tmin + nt*dt) {
    nt_seis = std::ceil((twtx_max - tmin)/dt);
  }
  twt_0_.resize(nt_seis);
  for (size_t i = 0; i < nt_seis; ++i){
    twt_0_[i] = tmin + (0.5 + i)*dt;
  }
  return twt_0_;
}

std::vector<double>  SeismicParameters::z_0(){
  size_t nz                      = seismic_geometry_->nz();
  double zmin                    = seismic_geometry_->z0();
  double dz                      = seismic_geometry_->dz();

  z_0_.resize(nz); 
  for (size_t i = 0; i < nz; ++i){
    z_0_[i] = zmin + (0.5 + i)*dz;
  }  
  return z_0_;
}

std::vector<double> SeismicParameters::offset_vec(){
  offset_vec_.resize(noffset_);
  for (size_t i = 0; i < noffset_; ++i) {
    offset_vec_[i] = offset_0_ + i*doffset_;
  }
  return offset_vec_;
}

void SeismicParameters::readEclipseGrid() {
    std::string filename = model_settings_->GetEclipseFileName();

    printf("Start reading Eclipsegrid from file\n");
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

void SeismicParameters::deleteParameterGrids() {
  delete vpgrid_;
  delete vsgrid_;
  delete rhogrid_;
  delete extra_parameter_grid_;
}

void SeismicParameters::deleteZandRandTWTGrids() {
  delete twtgrid_;
  delete zgrid_;
  delete rgridvec_;
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

    double x0 = temp_segy_geometry->GetX0();
    double y0 = temp_segy_geometry->GetY0();
    double lx = temp_segy_geometry->Getlx();
    double ly = temp_segy_geometry->Getly();
    double angle = temp_segy_geometry->GetAngle();
    double dx = temp_segy_geometry->GetDx();
    double dy = temp_segy_geometry->GetDy();

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
    rgridvec_ = new std::vector<NRLib::StormContGrid>(1);
    twtxgrid_ = new NRLib::StormContGrid(volume, nx, ny, nzrefl);
    thetagrid_ = new NRLib::StormContGrid(volume, nx, ny, nzrefl);
  }
  else {
    rgridvec_ = new std::vector<NRLib::StormContGrid>(ntheta_);
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
          (*twtxgrid_)(i, j, k) = 0.0;
          (*thetagrid_)(i, j, k) = 0.0;
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
    (*rgridvec_)[0] = rgrid;
  }
  else {
    for (size_t i = 0; i < ntheta_; i++) {
      (*rgridvec_)[i] = rgrid;
    }
  }
}
