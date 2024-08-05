#include "nrlib/extrapolation/extrapolategrid2d.hpp"
#include "nrlib/eclipsegrid/eclipsegrid.hpp"
#include "nrlib/geometry/interpolation.hpp"

#include <tbb/concurrent_queue.h>
//#include "tbb/compat/thread"

#include "seismic_parameters.hpp"
#include "wavelet.hpp"


//------------------------------------------------------------------
SeismicParameters::SeismicParameters(ModelSettings * model_settings)
//------------------------------------------------------------------
 : missing_(-999.0f)
{
  model_settings_   = model_settings;

  seismic_output_   = new SeismicOutput(model_settings);
  segy_geometry_    = NULL;

  theta_vec_        = model_settings->GetThetaVec();
  offset_vec_       = model_settings->GetOffsetVec();
  wavelet_scale_    = model_settings->GetWaveletScale();

  SetupWavelet(wavelet_,
               model_settings->GetOutputWavelet(),
               model_settings->GetRicker(),
               model_settings->GetPeakFrequency(),
               model_settings->GetWaveletLength(),
               model_settings->GetWaveletLengthFactor(),
               model_settings->GetDt(),
               model_settings->GetWaveletFileName(),
               model_settings->GetWaveletFileFormat(),
               model_settings->GetPrefix());

  ReadEclipseGrid(eclipse_grid_,
                  model_settings->GetEclipseFileName(),
                  model_settings->GetParameterNames(),
                  model_settings->GetExtraParameterNames());

  FindGeometry(seismic_geometry_,
               segy_geometry_,
               eclipse_grid_->GetGeometry(),
               model_settings);

  FindTopAndBaseSurfaces(top_time_,
                         bot_time_,
                         topeclipse_,
                         boteclipse_,
                         top_k_,
                         bottom_k_,
                         seismic_geometry_,
                         eclipse_grid_->GetGeometry(),
                         wavelet_,
                         model_settings);

  CreateGrids(seismic_geometry_,
              model_settings);
}


//------------------------------------------------------------------------
void SeismicParameters::SetupWavelet(Wavelet           *& wavelet,
                                     const bool           write_wavelet,
                                     const bool           use_ricker,
                                     const double         peakF,
                                     const double         length,
                                     const double         length_factor,
                                     const double         dt,
                                     const std::string  & file_name,
                                     const std::string  & file_format,
                                     const std::string  & prefix)
//------------------------------------------------------------------------
{
  if (use_ricker) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nMaking Ricker wavelet with peak frequency %.1f Hz\n", peakF);
    wavelet = new Wavelet(peakF,
                          dt,
                          length,
                          length_factor,
                          write_wavelet,
                          prefix);
  }
  else {
    bool error = false;
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nReading wavelet from file \'%s\'\n", file_name.c_str());
    wavelet = new Wavelet(file_name,
                          file_format,
                          dt,
                          length,
                          length_factor,
                          write_wavelet,
                          prefix,
                          error);
    if (error) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\nWavelet format \'%s\' is not recognized", file_format.c_str());
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\nAborting ...\n");
      exit(1);
    }
  }
}

//----------------------------------------------------------------------------------------------
void SeismicParameters::ReadEclipseGrid(NRLib::EclipseGrid             *& eclipse_grid,
                                        const std::string               & filename,
                                        const std::vector<std::string>  & names,
                                        const std::vector<std::string>  & extra_parameter_names)
//----------------------------------------------------------------------------------------------
{
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nReading Eclipse grid from file \'%s'\n", filename.c_str());

  eclipse_grid = new NRLib::EclipseGrid(filename);

  bool error = false;

  for (size_t i = 0; i < names.size(); ++i) {
    if (!eclipse_grid->HasParameter(names[i])) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nParameter \'%s\' is not found!\n", names[i].c_str());
      error = true;
    }
  }

  for (size_t i = 0; i < extra_parameter_names.size(); ++i) {
    if (!eclipse_grid->HasParameter(extra_parameter_names[i])) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nParameter \'%s\' is not found!\n", extra_parameter_names[i].c_str());
      error = true;
    }
  }

  if (error) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nAborting ...\n");
    exit(1);
  }
}


//------------------------------------------------------------------------------------
void SeismicParameters::FindGeometry(SeismicGeometry              *& seismic_geometry,
                                     NRLib::SegyGeometry          *& segy_geometry,
                                     const NRLib::EclipseGeometry  & eclipse_geometry,
                                     ModelSettings                 * model_settings)
//------------------------------------------------------------------------------------
{
  seismic_geometry = new SeismicGeometry();
  seismic_geometry->setDxDy(model_settings->GetDx(), model_settings->GetDy());
  seismic_geometry->setDz(model_settings->GetDz());
  seismic_geometry->setDt(model_settings->GetDt());

  double x0, y0, lx, ly, angle;
  std::string text;

  if (model_settings->GetAreaFromSegy() != "") {
    text = "SegY";
    NRLib::TraceHeaderFormat::coordSys_t coord = NRLib::TraceHeaderFormat::UTM;
    NRLib::TraceHeaderFormat *thf = new NRLib::TraceHeaderFormat(model_settings->GetScalcoLoc(),
                                                                 model_settings->GetUtmxLoc(),
                                                                 model_settings->GetUtmyLoc(),
                                                                 model_settings->GetIL0Loc(),
                                                                 model_settings->GetXL0Loc(),
                                                                 model_settings->GetStartTimeLoc(),
                                                                 coord);
    double z0 = 0.0;
    NRLib::Volume * volume = NULL;
    std::vector<NRLib::TraceHeaderFormat *> thfvec;
    thfvec.push_back(thf);

    NRLib::SegY segy(model_settings->GetAreaFromSegy(),
                     static_cast<float>(z0),
                     thfvec);
    segy.ReadAllTraces(volume, z0);
    segy.CreateRegularGrid();
    const NRLib::SegyGeometry * tmp_segy_geometry = segy.GetGeometry();

    segy_geometry = new NRLib::SegyGeometry(tmp_segy_geometry);
    segy_geometry->WriteGeometry();
    segy_geometry->WriteILXL();

    x0        = tmp_segy_geometry->GetX0();
    y0        = tmp_segy_geometry->GetY0();
    lx        = tmp_segy_geometry->Getlx();
    ly        = tmp_segy_geometry->Getly();
    angle     = tmp_segy_geometry->GetAngle();
    double dx = tmp_segy_geometry->GetDx();
    double dy = tmp_segy_geometry->GetDy();

    seismic_geometry->setDxDy(dx, dy);
  }
  else if (model_settings->GetAreaGiven()) {
    text = "model file";
    x0    = model_settings->GetX0();
    y0    = model_settings->GetY0();
    lx    = model_settings->GetLx();
    ly    = model_settings->GetLy();
    angle = model_settings->GetAngle();
  }
  else if (model_settings->GetAreaFromSurface() != "") {
    text = "surface";
    NRLib::RegularSurfaceRotated<double> toptime_rot = NRLib::RegularSurfaceRotated<double>(model_settings->GetAreaFromSurface());
    x0    = toptime_rot.GetXRef();
    y0    = toptime_rot.GetYRef();
    lx    = toptime_rot.GetLengthX();
    ly    = toptime_rot.GetLengthY();
    angle = toptime_rot.GetAngle();
  }
  else {
    text = "Eclipse grid";
    eclipse_geometry.FindEnclosingVolume(x0, y0, lx, ly, angle);
  }
  seismic_geometry->setGeometry(x0, y0, lx, ly, angle);

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\n                                     x            y           lx         ly     angle");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\n-------------------------------------------------------------------------------------");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nArea from %-15s %12.2f %12.2f   %10.2f %10.2f    %6.2f\n",
                              text.c_str(), x0, y0, lx, ly, (180.0/NRLib::Pi)*angle);
}

//----------------------------------------------------------------------------------------------
void SeismicParameters::FindTopAndBaseSurfaces(NRLib::RegularSurface<double> & top_time,
                                               NRLib::RegularSurface<double> & bot_time,
                                               NRLib::RegularSurface<double> & topeclipse,
                                               NRLib::RegularSurface<double> & boteclipse,
                                               size_t                        & top_k,
                                               size_t                        & bot_k,
                                               SeismicGeometry               * seismic_geometry,
                                               const NRLib::EclipseGeometry  & eclipse_geometry,
                                               Wavelet                       * wavelet,
                                               ModelSettings                 * model_settings)
//----------------------------------------------------------------------------------------------
{
  //
  // This setup below defines a rectangular non-rotated surface that
  // totally encloses the Eclipse grid. In addition, there is a small
  // safety margin of dx and dy. Strictly speaking dx and dy should
  // have been calculated the same way as lxsurf and lysurf
  //
  double xmin     = seismic_geometry->xmin();
  double ymin     = seismic_geometry->ymin();
  double dx       = seismic_geometry->dx();
  double dy       = seismic_geometry->dy();
  double lxsurf   = seismic_geometry->xsurfacelength();
  double lysurf   = seismic_geometry->ysurfacelength();
  size_t nxsurfec = seismic_geometry->nxsurfaceeclipse();
  size_t nysurfec = seismic_geometry->nysurfaceeclipse();

  double x0       = xmin - dx;
  double y0       = ymin - dy;
  double lx       = lxsurf + 2*dx;
  double ly       = lysurf + 2*dy;
  size_t nx       = nxsurfec + 2;
  size_t ny       = nysurfec + 2;

  //
  // Finding top and base time surfaces
  //
  bool const_top_given = true;
  if (model_settings->GetTopTimeSurfaceFile() != "") {
    const_top_given = false;
    NRLib::RegularSurfaceRotated<double> top_time_rotated = NRLib::RegularSurfaceRotated<double>(model_settings->GetTopTimeSurfaceFile());
    double topmin = top_time_rotated.Min();

    top_time = NRLib::RegularSurface<double>(x0, y0, lx, ly, nx, ny, topmin);
    top_time.SetMissingValue(top_time_rotated.GetMissingValue());
    for (size_t i = 0; i < top_time.GetNI(); i++) {
      for (size_t j = 0; j < top_time.GetNJ(); j++) {
        double x, y;
        top_time.GetXY(i, j, x, y);
        top_time(i, j) = top_time_rotated.GetZ(x, y);
      }
    }
    bot_time = NRLib::RegularSurface<double>(x0, y0, lx, ly, nx, ny, top_time.Max());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nTaking top time surface from file \'%s\'", model_settings->GetTopTimeSurfaceFile().c_str());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nSetting base time surface to maximum of top time : %.2f\n", top_time.Max());
  }
  else {
    double t1 = model_settings->GetTopTimeConstant();
    top_time = NRLib::RegularSurface<double>(x0, y0, lx, ly, nx, ny, t1);
    bot_time = NRLib::RegularSurface<double>(x0, y0, lx, ly, nx, ny, t1);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nSetting top and base time surfaces to     : %8.2f\n", t1);
  }

  //
  // Layers for which at least one cell is active and non-collapsed
  //
  top_k = eclipse_geometry.FindTopLayer();
  bot_k = eclipse_geometry.FindBottomLayer();

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nTop layer of Eclipse grid                 : %4d"  , top_k);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nBase layer of Eclipse grid                : %4d\n", bot_k);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nNumber of layers needed in z-grid         : %4d\n",bot_k + 2 - top_k);

  seismic_geometry->setZReflectorCount(static_cast<size_t>(bot_k + 2 - top_k));

  //
  // Finding top and base Eclipse surfaces in depth (regular and non-rotated)
  //
  topeclipse = NRLib::RegularSurface<double>(x0, y0, lx, ly, nx, ny, missing_);
  boteclipse = NRLib::RegularSurface<double>(x0, y0, lx, ly, nx, ny, missing_);

  double etdx       = topeclipse.GetDX();
  double etdy       = topeclipse.GetDY();
  double ebdx       = boteclipse.GetDX(); // equals etdx
  double ebdy       = boteclipse.GetDY(); // equals etdy
  double angle      = 0.0;
  bool   cornerpt   = model_settings->GetUseCornerpointInterpol();
  bool   fixed      = model_settings->GetUseFixedTriangularization();
  bool   bilinear   = model_settings->GetUseBilinearInterpolation();
  bool   faults     = model_settings->GetCornerpointInterpolationAtFaults();
  bool   horizontal = true; // Always use horizontal interpolation for top and base surface

  if (cornerpt)
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nFinding Eclipse top and base surfaces using cornerpoint interpolation.\n");
  else
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nFinding Eclipse top and base surfaces (not corner point interpolation).\n");

  NRLib::Grid2D<double> tvalues(nx, ny, 0.0);
  NRLib::Grid2D<double> bvalues(nx, ny, 0.0);

  eclipse_geometry.FindLayer(tvalues, top_k, 0, etdx, etdy, x0, y0, 0.0, cornerpt, faults, bilinear, fixed, horizontal, missing_);
  eclipse_geometry.FindLayer(bvalues, bot_k, 1, ebdx, ebdy, x0, y0, 0.0, cornerpt, faults, bilinear, fixed, horizontal, missing_);

  for (size_t i = 0; i < topeclipse.GetNI(); i++) {    // This fixes small, insignificant defects
    for (size_t j = 0; j < topeclipse.GetNJ(); j++) {
      topeclipse(i, j) = tvalues(i, j);
      boteclipse(i, j) = bvalues(i, j);
      if (topeclipse(i, j) > boteclipse(i, j)) {
        topeclipse(i, j) = boteclipse(i, j);
      }
    }
  }

  if (model_settings->GetOutputDepthSurfaces()) {
    seismic_output_->WriteDepthSurfaces(topeclipse, boteclipse);
  }

  double min;
  double max;
  double d1 = 1e10;
  double d2 =  0.0;

  bool found;
  for (size_t i = 0; i < eclipse_geometry.GetNI(); ++i) {
    for (size_t j = 0; j < eclipse_geometry.GetNJ(); ++j) {
      found = false;
      min = eclipse_geometry.FindZTopInCellActiveColumn(i, j, top_k, found);
      if (found && min < d1) {
        d1 = min;
      }
      found = false;
      max = eclipse_geometry.FindZBotInCellActiveColumn(i, j, bot_k, found);
      if (found && max > d2) {
        d2 = max;
      }
    }
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nEclipse grid minimum value                : %8.2f", d1);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nEclipse grid maximum value                : %8.2f\n", d2);

  //
  // Finding top and base Eclipse surfaces in time
  //

  // NBNB-PAL: Testen nedenfor er ikke god ettersom en flate spesifisert fra fil kan være konstant
  if (const_top_given) {
    std::string text = "PP";
    double const_v   = model_settings->GetConstVp()[0];
    double const_vs  = model_settings->GetConstVs()[0];
    double t1        = model_settings->GetTopTimeConstant();
    if (model_settings->GetPSSeismic()) {
      const_v = 2/(1/const_v + 1/const_vs);
      text = "PS";
    }
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nUsing %s velocity                         : %8.2f\n", text.c_str(), const_v);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nCalculating Eclipse top time surface\n");

    for (size_t i = 0; i < top_time.GetNI(); i++)
      for (size_t j = 0; j < top_time.GetNJ(); j++) {
        top_time(i, j) = t1 + 2000.0*(topeclipse(i, j) - d1)/const_v;
        bot_time(i, j) = top_time(i, j);
      }
  }

  //
  // Finding additional grid size due to wavelet length
  //
  double twt_wavelet = wavelet->GetTwtLength();
  std::vector<double> constvp = model_settings->GetConstVp();

  double z_top_wavelet = twt_wavelet*constvp[0]/2000;
  double z_bot_wavelet = twt_wavelet*constvp[2]/2000;

  std::string text = "PP";
  if (model_settings->GetPSSeismic()) {
    std::vector<double> constvs = model_settings->GetConstVs();
    double vel_top = 2/(1/constvp[0] + 1/constvs[0]);
    double vel_bot = 2/(1/constvp[2] + 1/constvs[2]);
    z_top_wavelet = twt_wavelet*vel_top/2000;
    z_bot_wavelet = twt_wavelet*vel_bot/2000;
    text = "PS";
  }
  model_settings->SetZWaveletTop(z_top_wavelet);
  model_settings->SetZWaveletBot(z_bot_wavelet);

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nWavelet time length is                    : %8.2f\n", twt_wavelet*2.0);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nEclipse top surface lift due to wavelet   : %8.2f"  , z_top_wavelet);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nEclipse base surface drop due to wavelet  : %8.2f\n", z_bot_wavelet);

  topeclipse.Add(-1 * z_top_wavelet); // add one wavelet length to bot and subtract from top
  boteclipse.Add(     z_bot_wavelet);

  double factor = model_settings_->GetPaddingFactorSeismicModelling();

  d1 -= z_top_wavelet*factor;
  d2 += z_bot_wavelet*factor;

  seismic_geometry->setZRange(d1, d2);

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nGrid minimum value                        : %8.2f (Eclipse grid minimum - %.1f wavelet)"  , d1, factor/2.0);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\nGrid maximum value                        : %8.2f (Eclipse grid maximum + %.1f wavelet)\n", d2, factor/2.0);
}

//---------------------------------------------------------------------
void SeismicParameters::CreateGrids(SeismicGeometry * seismic_geometry,
                                    ModelSettings   * model_settings)
//---------------------------------------------------------------------
{
  size_t                           nx                 = seismic_geometry->nx();
  size_t                           ny                 = seismic_geometry->ny();
  size_t                           nzrefl             = seismic_geometry->zreflectorcount();
  NRLib::Volume                    volume             = seismic_geometry->createDepthVolume();

  bool                             nmo_corr           = model_settings->GetNMOCorr();
  bool                             ps_seismic         = model_settings->GetPSSeismic();
  bool                             noise              = model_settings->GetAddNoiseToReflCoef();
  bool                             output_vrms        = model_settings->GetOutputVrms();
  bool                             output_refl        = model_settings->GetOutputReflections();
  const std::vector<std::string> & xtr_par_names      = model_settings->GetExtraParameterNames();
  const std::vector<double>      & xtr_par_def_values = model_settings->GetExtraParameterDefaultValues();
  const std::string              & twt_filename       = model_settings->GetTwtFileName();

  zgrid_   = new NRLib::StormContGrid(volume, nx, ny, nzrefl, 0.0);
  twtgrid_ = new NRLib::StormContGrid(volume, nx, ny, nzrefl, 0.0);
  vpgrid_  = new NRLib::StormContGrid(volume, nx, ny, nzrefl + 1, missing_);
  vsgrid_  = new NRLib::StormContGrid(volume, nx, ny, nzrefl + 1, missing_);
  rhogrid_ = new NRLib::StormContGrid(volume, nx, ny, nzrefl + 1, missing_);

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nMaking regular grids:");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n  z, TWT               %4d x %4d x %4d : %10d", nx, ny, nzrefl    , nx * ny * nzrefl);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n  Vp, Vs, Rho          %4d x %4d x %4d : %10d", nx, ny, nzrefl + 1, nx * ny * (nzrefl + 1));

  if (nmo_corr && ps_seismic) {
    twtssgrid_ = new NRLib::StormContGrid(volume, nx, ny, nzrefl, 0.0);
    twtppgrid_ = new NRLib::StormContGrid(volume, nx, ny, nzrefl, 0.0);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n  TWTss, TWTpp       %4d x %4d x %4d : %10d", nx, ny, nzrefl, nx * ny * nzrefl);
  } else {
    twtssgrid_ = NULL;
    twtppgrid_ = NULL;
  }

  if (nmo_corr && output_vrms) {
    vrmsgrid_ = new NRLib::StormContGrid(volume, nx, ny, nzrefl, 0.0);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\n  Vrms                 %4d x %4d x %4d : %10d", nx, ny, nzrefl, nx*ny*nzrefl);
  }
  else
    vrmsgrid_ = NULL;

  if (output_refl) {
    if (noise) {
      rgridvec_ = new std::vector<NRLib::StormContGrid>(2);
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\n  Refl.coef/Noise 2 x %4d x %4d x %4d : %10d", nx, ny, nzrefl, nx*ny*nzrefl*2);
    }
    else {
      rgridvec_ = new std::vector<NRLib::StormContGrid>(1);
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\n  Refl. coef.          %4d x %4d x %4d : %10d", nx, ny, nzrefl, nx*ny*nzrefl);
    }
    NRLib::StormContGrid rgrid(volume, nx, ny, nzrefl, 0.0);

    (*rgridvec_)[0] = rgrid;
    if (noise) {
      (*rgridvec_)[1] = rgrid;
    }
  }
  else {
    rgridvec_ = NULL;
  }

  if (twt_filename != "") {
    twt_timeshift_ = new NRLib::StormContGrid(model_settings->GetTwtFileName());
    if ((*twt_timeshift_).GetNI() != nx || (*twt_timeshift_).GetNJ() != ny || (*twt_timeshift_).GetNK() != nzrefl) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\nTWT timeshift from file has wrong dimension.\nAborting...\n");
      exit(1);
    }
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\n  TWT time shift     %4d x %4d x %4d : %10d", nx, ny, nzrefl, nx*ny*nzrefl);
  }
  else {
    twt_timeshift_ = NULL;
  }

  size_t n = xtr_par_names.size();
  if (n > 1) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n  Extra par. (%d)       %4d x %4d x %4d : %10d\n", n, nx, ny, nzrefl, nx * ny * nzrefl);
    extra_parameter_grid_.resize(n);
    for (size_t p = 0; p < n; ++p) {
      extra_parameter_grid_[p] = new NRLib::StormContGrid(volume, nx, ny, nzrefl + 1, static_cast<float>(xtr_par_def_values[p]));
      for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
          (*extra_parameter_grid_[p])(i, j, nzrefl) = 0.0;
        }
      }
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nFilling in default value %.2f for extra parameter \'%s\'",
                                  xtr_par_def_values[p], xtr_par_names[p].c_str());
    }
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n\n");
}


//========================================== BELOW IS UNPROCESSED ==========================================

//---------------------------------------------------------------------------
void SeismicParameters::SetSegyGeometry(const NRLib::SegyGeometry * geometry)
//---------------------------------------------------------------------------
{
  if (segy_geometry_ != NULL) {
    delete segy_geometry_;
  }
  segy_geometry_ = new NRLib::SegyGeometry(geometry);
}

void SeismicParameters::FindLoopIndeces(int  & n_xl,
                                        int  & il_min,
                                        int  & il_max,
                                        int  & il_step,
                                        int  & xl_min,
                                        int  & xl_max,
                                        int  & xl_step,
                                        bool & segy)
{
  if (segy_geometry_ == NULL) {
    il_min  = 0;
    il_max  = static_cast<int>(seismic_geometry_->nx());
    il_step = 1;
    xl_min  = 0;
    xl_max  = static_cast<int>(seismic_geometry_->ny()-1);
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

//-------------------------------------------------------------------------
void SeismicParameters::FindVrms(std::vector<double>       & vrms_vec,
                                 const std::vector<double> & twt_vec,
                                 const std::vector<double> & v_vec,
                                 double                      z_res) const
//-------------------------------------------------------------------------
{
  double undef  = -999.0;
  double v_w    = model_settings_->GetVw(); // Velocity of water
  double z_w    = model_settings_->GetZw(); // Depth to seabed
  double twt_w  = 2000 * z_w / v_w;
  size_t nk     = twt_vec.size();

  //generate vrms_vec - only for each layer in reservoir
  if (twt_vec[0] == undef) {
    for (size_t k = 0; k < nk; ++k) {
      vrms_vec[k] = undef;
    }
  }
  else {
    double v_over = 2000 * (z_res - z_w) / (twt_vec[0] - twt_w);
    double tmp0   = v_w*v_w*twt_w + v_over*v_over*(twt_vec[0] - twt_w);
    for (size_t k = 0; k < nk; ++k) {
      double tmp = tmp0;
      for (size_t l = 1; l <= k; ++l) {
        tmp += v_vec[l]* v_vec[l]*(twt_vec[l] - twt_vec[l-1]);
      }
      vrms_vec[k] = std::sqrt(tmp / twt_vec[k]);
    }
  }
}

//---------------------------------------------------------------------------------------------
void SeismicParameters::FindVrmsRegular(const std::vector<double> & vrms_vec,
                                        std::vector<double>       & vrms_vec_reg,
                                        const std::vector<double> & twt_vec,
                                        const std::vector<double> & twt_vec_reg,
                                        const std::vector<double> & v_vec,
                                        double                      const_v,
                                        double                      twt_wavelet_extrapol) const
//---------------------------------------------------------------------------------------------
{
  //generate vrms_vec_reg - including overburden and underburden

  size_t nk          = twt_vec.size();
  double v_w         = model_settings_->GetVw(); // Velocity of water
  double z_w         = model_settings_->GetZw(); // Depth to seabed
  double twt_w       = 2000 * z_w / v_w;
  double vrms_under  = vrms_vec[nk - 1]*vrms_vec[nk - 1]*twt_vec[nk - 1] + const_v * const_v * twt_wavelet_extrapol;
  vrms_under        *= (1 / (twt_vec[nk - 1] + twt_wavelet_extrapol));
  vrms_under         = std::sqrt(vrms_under);

  //sample vrms regularly:
  std::vector<double> twt_vec_in(nk + 2);
  std::vector<double> vrms_vec_in(nk + 2);
  twt_vec_in[0]  = twt_w;
  twt_vec_in[1]  = twt_vec[0];
  vrms_vec_in[0] = v_w;
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
  twt_vec_in [index] = twt_vec_in[index - 1] + twt_wavelet_extrapol;
  vrms_vec_in[index] = vrms_under;

  vrms_vec_reg = NRLib::Interpolation::Interpolate1D(twt_vec_in, vrms_vec_in, twt_vec_reg, "linear");
}

//---------------------------------------------------------------------------------
void SeismicParameters::FindReflections(NRLib::Grid2D<double>          & r_vec,
                                           const NRLib::Grid2D<double> & theta_vec,
                                           size_t                        i,
                                           size_t                        j)
//---------------------------------------------------------------------------------
{
  bool ps_seis   = model_settings_->GetPSSeismic();

  std::vector<double> vp_vec (bottom_k_ - top_k_ + 3);
  std::vector<double> vs_vec (bottom_k_ - top_k_ + 3);
  std::vector<double> rho_vec(bottom_k_ - top_k_ + 3);

  for (size_t off = 0; off < theta_vec.GetNJ(); ++off) {
    for (size_t k = top_k_; k <= bottom_k_ + 2; k++) {
      size_t kk = k - top_k_;
      vp_vec [kk] = (*vpgrid_ )(i, j, kk);
      vs_vec [kk] = (*vsgrid_ )(i, j, kk);
      rho_vec[kk] = (*rhogrid_)(i, j, kk);
    }
    for (size_t k = top_k_; k <= bottom_k_ + 1; k++) {
      size_t kk        = k - top_k_;
      double diffvp    =      vp_vec [kk + 1] - vp_vec [kk];
      double meanvp    = 0.5*(vp_vec [kk + 1] + vp_vec [kk]);
      double diffvs    =      vs_vec [kk + 1] - vs_vec [kk];
      double meanvs    = 0.5*(vs_vec [kk + 1] + vs_vec [kk]);
      double diffrho   =      rho_vec[kk + 1] - rho_vec[kk];
      double meanrho   = 0.5*(rho_vec[kk + 1] + rho_vec[kk]);
      double refl;
      double theta     = theta_vec(kk, off);
      double vpvs      = meanvs/meanvp;
      double sin2theta = sin(theta)*sin(theta);

      if (ps_seis) {
        double sin_theta = sin(theta);
        double cos_theta = cos(theta);
        double cos_phi   = cos(asin(sin_theta*vpvs));
        double a2        = 2*sin_theta*vpvs*(vpvs*sin2theta/cos_phi - cos_theta);
        double a3        = sin_theta*(-0.5 + sin2theta*vpvs*vpvs - cos_theta*cos_phi*vpvs)/cos_phi;
        refl             = a2*diffvs/meanvs + a3*diffrho/meanrho;

      }
      else {
        double tan2theta = tan(theta)*tan(theta);
        double a1        =  0.5*(1.0 + tan2theta);
        double a2        = -4.0*sin2theta*vpvs*vpvs;
        double a3        =  0.5*(1.0 + a2);
        refl             = a1*diffvp/meanvp + a2*diffvs/meanvs + a3*diffrho/meanrho;
      }
      r_vec(kk, off) = static_cast<float>(refl);
    }
  }
}

//----------------------------------------------------------
void SeismicParameters::FindMaxTwtIndex(size_t & i_max,
                                        size_t & j_max,
                                        double & max_value)
//----------------------------------------------------------
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

//-----------------------------------------------------------------------------------
void SeismicParameters::GenerateTwt0AndZ0(ModelSettings       * model_settings,
                                          std::vector<double> & twt_0,
                                          std::vector<double> & z_0,
                                          std::vector<double> & twts_0,
                                          size_t              & time_samples_stretch,
                                          bool                  ps_seis)
//-----------------------------------------------------------------------------------
{
  size_t nt             = seismic_geometry_->nt();
  double dt             = seismic_geometry_->dt();
  double t0             = seismic_geometry_->t0();
  size_t nzrefl         = seismic_geometry_->zreflectorcount();
  size_t nz             = seismic_geometry_->nz();
  double twt_wavelet    = wavelet_->GetTwtLength();              // Half a wavelet
  double stretch_factor = -999.0;

  if (model_settings->GetNMOCorr() && !model_settings->GetOffsetWithoutStretch()){
    GenerateTwt0ForNMO(twt_0,
                       time_samples_stretch,
                       stretch_factor,
                       ps_seis,
                       twt_wavelet,
                       nt,
                       dt,
                       t0,
                       nzrefl);

    GenerateZ0ForNMO(z_0,
                     stretch_factor);

    if (model_settings_->GetTwtFileName() != "") {
      twts_0 = GenerateTWT0Shift(twt_0[0], time_samples_stretch);
    }
  }
  else {
    twt_0.resize(nt);
    for (size_t i = 0; i < nt; ++i){
      twt_0[i] = t0 + i*dt;
    }
    double zmin       = seismic_geometry_->z0();
    double dz         = seismic_geometry_->dz();
    size_t nzmin      = static_cast<size_t>(floor(zmin) / dz + 0.5); // Find min z in a sample
    double zmin_sampl = nzmin *dz;
    z_0.resize(nz);

    for (size_t i = 0; i < nz; ++i){
      z_0[i] = zmin_sampl + i*dz;
    }

    if (model_settings_->GetTwtFileName() != "") {
      twts_0 = GenerateTWT0Shift(twt_0[0], twt_0.size());
    }
  }
}

//------------------------------------------------------------------------------
void SeismicParameters::GenerateTwt0ForNMO(std::vector<double> & twt_0,
                                           size_t              & nt_stretch,
                                           double              & stretch_factor,
                                           const bool            ps_seis,
                                           const double          twt_wavelet,
                                           const size_t          nt,
                                           const double          dt,
                                           const double          t0,
                                           const size_t          nzrefl)
//------------------------------------------------------------------------------
{
  // Account for stretch by making twt0 sufficiently long. Stretch upwards is also taken into account
  // through "xtra_samples_top".
  // Max twt value and location is found from twtgrid and twtx is calculated for the largest offset.
  // "Time samples stretch" is number of samples in nmo-corrected seismic trace.

  size_t i_max;
  size_t j_max;
  double twt_max;
  double twtx_max;
  size_t noffset    = offset_vec_.size();
  double offset_max = offset_vec_[noffset - 1];

  FindMaxTwtIndex(i_max, j_max, twt_max);   // Find max from twgrid, and index of max twt value

  // Find max TWTX for highest offset in order to find the highest TWT value to sample seismic

  //
  //------------PS seismic------------
  //
  if (ps_seis) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGenerate Twt0 for NMO assuming PS seismic\n");
    std::vector<double> vp_vec(nzrefl);
    std::vector<double> vs_vec(nzrefl);
    std::vector<double> twt_pp_vec(nzrefl);
    std::vector<double> twt_ss_vec(nzrefl);
    std::vector<double> vrms_pp_vec(nzrefl);
    std::vector<double> vrms_ss_vec(nzrefl);

    for (size_t k = 0; k < nzrefl; ++k) {
      twt_pp_vec[k] = (*twtppgrid_)(i_max, j_max, k);
      twt_ss_vec[k] = (*twtssgrid_)(i_max, j_max, k);
      vp_vec[k]     = (*vpgrid_)(i_max, j_max, k);
      vs_vec[k]     = (*vsgrid_)(i_max, j_max, k);
    }

    FindVrms(vrms_pp_vec, twt_pp_vec, vp_vec, (*zgrid_)(i_max, j_max, 0));
    FindVrms(vrms_ss_vec, twt_ss_vec, vs_vec, (*zgrid_)(i_max, j_max, 0));

    double vrms_pp     = vrms_pp_vec[nzrefl - 1];
    double vrms_ss     = vrms_ss_vec[nzrefl - 1];
    double twt_pp_max  = twt_pp_vec [nzrefl - 1];
    double twt_ss_max  = twt_ss_vec [nzrefl - 1];
    double tmp         = offset_max / (vrms_pp_vec[nzrefl - 1]*twt_pp_max / 1000);
    double start_value = atan(tmp);
    if (start_value >= 1.0)
      start_value = 0.99;
    double dU        = vrms_ss * twt_ss_max / 2000;
    double dD        = vrms_pp * twt_pp_max / 2000;
    double vr        = vrms_ss / vrms_pp;
    size_t n_it      = 10;
    double y_out     = FindSinThetaPSWithNewtonsMethod(start_value,
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

    double twtx_pp   = std::sqrt(twt_pp_max * twt_pp_max / 4 + 1000 * 1000 * (offset_pp * offset_pp) / (vrms_pp * vrms_pp));
    double twtx_ss   = std::sqrt(twt_ss_max * twt_ss_max / 4 + 1000 * 1000 * (offset_ss * offset_ss) / (vrms_ss * vrms_ss));
    twtx_max         = twtx_pp + twtx_ss;
  }
  //
  //------------PP seismic------------
  //
  else {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGenerate Twt0 for NMO assuming PP seismic.\n");
    std::vector<double> vrms_vec(nzrefl);
    std::vector<double> vp_vec(nzrefl);
    std::vector<double> twt_vec(nzrefl);
    for (size_t k = 0; k < nzrefl; ++k) {
      twt_vec[k] = (*twtgrid_)(i_max, j_max, k);
      vp_vec[k]  = (*vpgrid_) (i_max, j_max, k);
    }
    FindVrms(vrms_vec, twt_vec, vp_vec, (*zgrid_)(i_max, j_max, 0));  // Find max vrms in index
    double vrms_max_t = vrms_vec[vrms_vec.size() - 1];
    twtx_max = std::sqrt(twt_max*twt_max + 1000 * 1000 * offset_max*offset_max / (vrms_max_t*vrms_max_t));
  }

  stretch_factor       = twtx_max / twt_max;                          // stretch factor. NO wavelet in twtx_max and twt_max
  double tmin          = t0;                                          // min twt sample. Includes wavlet
  size_t ishift        =  0;                                          // samples on top due to stretch
  double tmax_stretch  = seismic_geometry_->tmax();                   // max twt sample due to stretch. Includes wavelet
  size_t nt_seis       = nt;                                          // number of samples in prenmo seis (twt_0)
  nt_stretch           = nt;                                          // samples in nmo corrected seismic - include stretch top and bot
  twtx_max            += twt_wavelet;                                 // add one wavelet (as no wavelet is included here)

  if (stretch_factor > 1) {
    tmin         -= stretch_factor*2*twt_wavelet;                     // Subtracting a stetched wavlet
    tmax_stretch += stretch_factor*6*twt_wavelet;                     // Add three stetched wavlets
    ishift        = static_cast<size_t>(floor((t0 - tmin)/dt + 0.5)); // Difference in start time between non NMO and NMO
    nt_stretch    = static_cast<size_t>(ceil((tmax_stretch - tmin) / dt + 0.5));
    nt_seis       = static_cast<size_t>(ceil((twtx_max     - tmin) / dt + 0.5));
  }

  twt_0.resize(nt_seis);
  for (size_t i = 0; i < nt_seis; ++i){
    twt_0[i] = (t0 - ishift * dt) + i*dt;
  }
  if (nt_stretch > nt_seis){
    nt_stretch = nt_seis;
  }

  if (stretch_factor > 1) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "   stretch factor :  %7.5f  (twtx_max / twt_max)\n", stretch_factor);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "   wavelet length :  %7.2f\n", twt_wavelet*2);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "   nt_stretch     :  %7d\n"  , nt_stretch);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "   nt_seis        :  %7d\n"  , nt_seis);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "   tmax_stretch   :  %7.2f\n", tmax_stretch);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "   twtx_max       :  %7.2f\n", twtx_max);
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "   twt0_min       :  %7.2f\n",twt_0[0]);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "   twt0_max       :  %7.2f\n",twt_0[nt_seis - 1]);
}

//----------------------------------------------------------------------------
void SeismicParameters::GenerateZ0ForNMO(std::vector<double> & z_0,
                                         const double          stretch_factor)
//----------------------------------------------------------------------------
{
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGenerate Z0 for NMO.\n");

  double z_top_wavelet           = model_settings_->GetZWaveletTop();
  double z_bot_wavelet           = model_settings_->GetZWaveletBot();

  double zmin                    = seismic_geometry_->z0();
  size_t nz                      = seismic_geometry_->nz();
  double dz                      = seismic_geometry_->dz();

  double min_z                   = zmin;
  double max_z                   = zmin + (nz - 1)*dz;

  // Account for stretch above and below reservoir

  if (stretch_factor > 1) {
    min_z -= stretch_factor * 2 * z_top_wavelet;
    max_z += stretch_factor * 2 * z_bot_wavelet;
  }

  size_t nz_seis     = static_cast<size_t>(std::floor((max_z - min_z)/dz + 0.5)) + 1;
  size_t nzmin       = static_cast<size_t>(floor(min_z) / dz + 0.5);                  // Find min z in a sample
  double min_z_sampl = nzmin *dz;

  z_0.resize(nz_seis);
  for (size_t i = 0; i < nz_seis; ++i){
    z_0[i] = min_z_sampl + i*dz;
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "   z0_min        :  %7.2f\n", min_z_sampl);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "   z0_max        :  %7.2f\n", z_0[nz_seis - 1]);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "   dz            :  %7.2f\n", dz);
}

//-------------------------------------------------------------------------
std::vector<double> SeismicParameters::GenerateTWT0Shift(double twt_0_min,
                                                         size_t n_samples)
//-------------------------------------------------------------------------
{
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGenerate TWT0 shift.\n");

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

  double dt     = seismic_geometry_->dt();
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

//-------------------------------------------------------------------------------------------
void SeismicParameters::FindPSNMOThetaAndOffset(NRLib::Grid2D<double>     & thetagrid,
                                                NRLib::Grid2D<double>     & offset_down_grid,
                                                NRLib::Grid2D<double>     & offset_up_grid,
                                                const std::vector<double> & twt_pp_vec,
                                                const std::vector<double> & twt_ss_vec,
                                                const std::vector<double> & vrms_pp_vec,
                                                const std::vector<double> & vrms_ss_vec,
                                                const std::vector<double> & offset,
                                                bool                        save_theta)
//-------------------------------------------------------------------------------------------
{
  NRLib::LogKit::LogFormatted(NRLib::LogKit::High, "\nFinding PS NMO angles and offsets\n");

  double tol         = 0.000001;
  size_t n_it        = 10;
  size_t n_it_avg    = 0;
  double start_value = 0.1;

  for (size_t off = 0; off < offset.size(); off++) {
    double tmp = offset[off] / (vrms_pp_vec[0]*twt_pp_vec[0] / 1000);
    start_value = atan(tmp);
    if (start_value >= 1.0)
      start_value = 0.99;
    for (size_t k = 0; k < twt_pp_vec.size(); k++) {
      double dU = vrms_ss_vec[k] * twt_ss_vec[k] / 2000;
      double dD = vrms_pp_vec[k] * twt_pp_vec[k] / 2000;
      double vr = vrms_ss_vec[k] / vrms_pp_vec[k];
      n_it = 20;
      double y_out = FindSinThetaPSWithNewtonsMethod(start_value,
                                                     offset[off],
                                                     dU,
                                                     dD,
                                                     vr,
                                                     tol,
                                                     n_it);
      n_it_avg  += n_it;
      double theta_up   = asin(vr*y_out);
      double theta_down = asin(y_out);
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

//-----------------------------------------------------------------------------
double SeismicParameters::FindSinThetaPSWithNewtonsMethod(double   start_value,
                                                          double   offset,
                                                          double   dU,
                                                          double   dD,
                                                          double   vr,
                                                          double   tol,
                                                          size_t & n_it)
//-----------------------------------------------------------------------------
{
  NRLib::LogKit::LogFormatted(NRLib::LogKit::High, "\nFind sin(theta) PS with Newtons method.\n");
  double y_old = start_value;
  double y_new, f_y, f_der_y;

  for (size_t i = 0; i < n_it; ++i) {
    f_y     = -offset + dD*y_old / sqrt(1.0 - std::pow(y_old, 2)) + dU*vr*y_old / sqrt(1.0 - std::pow(vr, 2) * std::pow(y_old, 2));
    f_der_y = dD / (std::pow((1.0 - y_old), (3 / 2))) + dU*vr / (std::pow((1.0 - std::pow(vr, 2) * std::pow(y_old, 2)), (3 / 2)));

    if (f_der_y == 0) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\nFailure in newtons method: zero derivative.\n");
      return 0;
    }
    y_new = y_old - f_y / f_der_y;

    if (std::abs(y_new) > 1.0) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\nFailure in newtons method: Value > 1.0 : y_old = %.2f , y_new = %.2f. New value: y_new = 0.1 suggested.\n",y_old,y_new);
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

//-----------------------------------------
void SeismicParameters::DeleteEclipseGrid()
//-----------------------------------------
{
  delete eclipse_grid_;
}

//---------------------------------------------------
void SeismicParameters::DeleteElasticParameterGrids()
//---------------------------------------------------
{
  delete vpgrid_;
  delete vsgrid_;
  delete rhogrid_;
}

//---------------------------------------------------
void SeismicParameters::DeleteExtraParameterGrids()
//---------------------------------------------------
{
  for (size_t i = 0; i < extra_parameter_grid_.size(); ++i)
    delete extra_parameter_grid_[i];
}

//----------------------------------------------
void SeismicParameters::DeleteZandRandTWTGrids()
//----------------------------------------------
{
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

//----------------------------------------------
void SeismicParameters::DeleteVrmsGrid()
//----------------------------------------------
{
  delete vrmsgrid_;
}

//----------------------------------------------
void SeismicParameters::DeleteWavelet()
//----------------------------------------------
{
  delete wavelet_;
}

//----------------------------------------------
void SeismicParameters::DeleteGeometryAndOutput()
//----------------------------------------------
{
  delete seismic_geometry_;
  delete segy_geometry_;
  delete seismic_output_;
  delete model_settings_;
}

//--------------------------------------------------------------
void SeismicParameters::PrintElapsedTime(time_t      start_time,
                                         std::string work)
//--------------------------------------------------------------
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

//-------------------------------------------------------------------------------------
tbb::concurrent_queue<Trace*> SeismicParameters::FindTracesInForward(size_t & n_traces)
//-------------------------------------------------------------------------------------
{
  tbb::concurrent_queue<Trace*> traces;
  int n_xl, il_min, il_max, il_step, xl_min, xl_max, xl_step;
  bool ilxl_loop = false;
  //find index min and max and whether loop over i,j or il,xl:
  FindLoopIndeces(n_xl, il_min, il_max, il_step, xl_min, xl_max, xl_step, ilxl_loop);
  //NRLib::SegyGeometry *geometry = GetSegyGeometry();
  int il_steps = 0;
  int xl_steps = 0;
  //----------------------LOOP OVER I,J OR IL,XL---------------------------------
  size_t job_number = 0;
  for (int il = il_min; il <= il_max; il += il_step) {
    ++il_steps;
    xl_steps = 0;
    for (int xl = xl_min; xl <= xl_max; xl += xl_step) {
      ++xl_steps;
      size_t i, j;
      double x, y;
      if (ilxl_loop) { //loop over il,xl, find corresponding x,y,i,j
        segy_geometry_->FindXYFromILXL(il, xl, x, y);
        segy_geometry_->FindIndex(x, y, i, j);
      }
      else { //loop over i,j, no segy output
        i = il;
        j = xl;
        x = 0;
        y = 0;
      }
      Trace * trace = new Trace(job_number, x, y, i, j);
      traces.push(trace);
      ++job_number;
    }
  }
  n_traces = job_number;
  return traces;
}

//--------------------------------------------------------------
void SeismicParameters::MonitorInitialize(size_t   n_traces,
                                          float  & monitor_size,
                                          float  & next_monitor)
//--------------------------------------------------------------
{
  monitor_size = static_cast<float>(n_traces) * 0.02f;
  if (monitor_size < 1.0f)
    monitor_size = 1.0f;
  next_monitor = monitor_size;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";
}

//----------------------------------------------------
void SeismicParameters::Monitor(size_t   trace,
                                float    monitor_size,
                                float  & next_monitor)
//----------------------------------------------------
{
  if (trace + 1 >= static_cast<size_t>(next_monitor)) {
    next_monitor += monitor_size;
    std::cout << "^";
    fflush(stdout);
    if (next_monitor > monitor_size * 51) {
      std::cout << "\n";
    }
  }
}

//-----------------------------------------------------------------------------------------------
std::vector<double> SeismicParameters::SplineInterp1D(const std::vector<double> & x_in,
                                                      const std::vector<double> & y_in,
                                                      const std::vector<double> & x_out,
                                                      double                      extrap_value)
//-----------------------------------------------------------------------------------------------
{
  std::vector<double> x_in_copy(x_in.size());
  std::vector<double> y_in_copy(y_in.size());
  x_in_copy[0] = x_in[0];
  y_in_copy[0] = y_in[0];
  size_t index = 1;
  for (size_t i = 1; i < x_in.size(); ++i) {
    if (x_in[i] != x_in[i - 1]) {
      x_in_copy[index] = x_in[i];
      y_in_copy[index] = y_in[i];
      ++index;
    }
  }
  x_in_copy.resize(index);
  y_in_copy.resize(index);
  if (index == 1) {
    std::vector<double> y_out(x_out.size(), extrap_value);
    return y_out;
  }
  return NRLib::Interpolation::Interpolate1D(x_in_copy, y_in_copy, x_out, "spline", extrap_value);
}

//------------------------------------------------------------------------------------
std::vector<double> SeismicParameters::LinInterp1D(const std::vector<double> & x_in,
                                                   const std::vector<double> & y_in,
                                                   const std::vector<double> & x_out)
//------------------------------------------------------------------------------------

{
  std::vector<double> x_in_copy(x_in.size());
  std::vector<double> y_in_copy(y_in.size());
  x_in_copy[0] = x_in[0];
  y_in_copy[0] = y_in[0];
  size_t index = 1; //sjekk om denne skal være 0???
  for (size_t i = 1; i < x_in.size(); ++i) {
    if (x_in[i] != x_in[i - 1]) {
      x_in_copy[index] = x_in[i];
      y_in_copy[index] = y_in[i];
      ++index;
    }
  }
  x_in_copy.resize(index);
  y_in_copy.resize(index);
  if (index == 1) {
    std::vector<double> y_out(x_out.size(), 0);
    return y_out;
  }
  return NRLib::Interpolation::Interpolate1D(x_in_copy, y_in_copy, x_out, "linear");
}
