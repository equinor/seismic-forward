#include "output.hpp"


#include <seismic_geometry.hpp>
#include "modelsettings.hpp"

Output::Output(SeismicParameters &seismic_parameters,
                     std::vector<double> twt_0,
                     std::vector<double> z_0,
                     std::vector<double> twts_0,
                     std::vector<double> offset_vec,
                     size_t time_samples_stretch)
  : segy_ok_(false),
    time_segy_ok_(false),
    prenmo_time_segy_ok_(false),
    time_stack_segy_ok_(false),
    depth_segy_ok_(false),
    depth_stack_segy_ok_(false),
    timeshift_segy_ok_(false),
    timeshift_stack_segy_ok_(false),
    twtx_segy_ok_(false),
    twt_0_(twt_0),
    z_0_(z_0),
    twts_0_(twts_0),
    offset_vec_(offset_vec)
{
  size_t nx = seismic_parameters.seismicGeometry()->nx();
  size_t ny = seismic_parameters.seismicGeometry()->ny();
  size_t nz = seismic_parameters.seismicGeometry()->nz();
  size_t nt = seismic_parameters.seismicGeometry()->nt();

  NRLib::Volume volume   = seismic_parameters.seismicGeometry()->createDepthVolume();
  NRLib::Volume volume_t = seismic_parameters.seismicGeometry()->createTimeVolume();
  bool nmo = seismic_parameters.modelSettings()->GetNMOCorr();

  if (seismic_parameters.GetSegyOutput()) {
    seismic_parameters.seismicOutput()->SetSegyGeometry(seismic_parameters, volume_t, nx, ny);
    segy_ok_ = seismic_parameters.seismicOutput()->CheckUTMPrecision(seismic_parameters, volume_t, nx, ny);
  }
  if (segy_ok_) {
    if (seismic_parameters.modelSettings()->GetOutputTimeSegy()) {
      std::string filename        = "seismic_time";
      time_segy_ok_            = seismic_parameters.seismicOutput()->PrepareSegy(time_segy_, twt_0_, time_samples_stretch, filename, seismic_parameters, offset_vec_, offset_vec_.size(), true, nmo);
    }
    if (seismic_parameters.modelSettings()->GetOutputPrenmoTimeSegy()) {
      std::string filename        = "seismic_prenmo_time";
      prenmo_time_segy_ok_         = seismic_parameters.seismicOutput()->PrepareSegy(prenmo_time_segy_, twt_0_, twt_0_.size(), filename, seismic_parameters, offset_vec_, offset_vec_.size(), true, nmo);
    }
    if (seismic_parameters.modelSettings()->GetOutputSeismicStackTimeSegy()) {
      std::string filename        = "seismic_time_stack";
      time_stack_segy_ok_      = seismic_parameters.seismicOutput()->PrepareSegy(time_stack_segy_, twt_0_, time_samples_stretch, filename, seismic_parameters, offset_vec_, 1, true, nmo);
    }
    if (seismic_parameters.modelSettings()->GetOutputDepthSegy()) {
      std::string filename        = "seismic_depth";
      depth_segy_ok_           = seismic_parameters.seismicOutput()->PrepareSegy(depth_segy_, z_0_, z_0_.size(), filename, seismic_parameters, offset_vec_, offset_vec_.size(), false, nmo);
    }
    if (seismic_parameters.modelSettings()->GetOutputSeismicStackDepthSegy()) {
      std::string filename        = "seismic_depth_stack";
      depth_stack_segy_ok_     = seismic_parameters.seismicOutput()->PrepareSegy(depth_stack_segy_, z_0_, z_0_.size(), filename, seismic_parameters, offset_vec_, 1, false, nmo);
    }
    if (seismic_parameters.modelSettings()->GetOutputTimeshiftSegy()) {
      std::string filename        = "seismic_timeshift";
      timeshift_segy_ok_       = seismic_parameters.seismicOutput()->PrepareSegy(timeshift_segy_, twts_0_, twts_0_.size(), filename, seismic_parameters, offset_vec_, offset_vec_.size(), true, nmo);
    }
    if (seismic_parameters.modelSettings()->GetOutputSeismicStackTimeShiftSegy()) {
      std::string filename        = "seismic_timeshift_stack";
      timeshift_stack_segy_ok_ = seismic_parameters.seismicOutput()->PrepareSegy(timeshift_stack_segy_, twts_0_, twts_0_.size(), filename, seismic_parameters, offset_vec_, 1, true, nmo);
    }
    if (seismic_parameters.modelSettings()->GetOutputTwtOffset()) {
      std::string filename        = "twt_offset";
      twtx_segy_ok_                = seismic_parameters.seismicOutput()->PrepareSegy(twtx_segy_, twt_0_, twt_0_.size(), filename, seismic_parameters, offset_vec_, offset_vec_.size(), true, nmo);
    }
  }
  
  //prepare grid if output of seismic in storm is requested
  if (seismic_parameters.GetTimeStormOutput()) {
    if (nmo) {
      NRLib::Volume volume_t_nmo = NRLib::Volume(seismic_parameters.seismicGeometry()->x0(),
                                                 seismic_parameters.seismicGeometry()->y0(),
                                                 twt_0_[0] - seismic_parameters.seismicGeometry()->dt()/2,
                                                 seismic_parameters.seismicGeometry()->xlength(),
                                                 seismic_parameters.seismicGeometry()->ylength(),
                                                 (time_samples_stretch * seismic_parameters.seismicGeometry()->dt()),
                                                 seismic_parameters.seismicGeometry()->angle());

      timegrid_      = new NRLib::StormContGrid(volume_t_nmo, nx, ny, time_samples_stretch);
    }
    else {
      timegrid_      = new NRLib::StormContGrid(volume_t, nx, ny, nt);
    }
  }
  if (seismic_parameters.GetTimeshiftStormOutput()) {
    NRLib::Volume volume_ts_nmo = NRLib::Volume(seismic_parameters.seismicGeometry()->x0(),
                                                seismic_parameters.seismicGeometry()->y0(),
                                                twts_0_[0] - seismic_parameters.seismicGeometry()->dt()/2,
                                                seismic_parameters.seismicGeometry()->xlength(),
                                                seismic_parameters.seismicGeometry()->ylength(),
                                                (twts_0_.size() * seismic_parameters.seismicGeometry()->dt()),
                                                seismic_parameters.seismicGeometry()->angle());
    timeshiftgrid_ = new NRLib::StormContGrid(volume_ts_nmo, nx, ny, twts_0_.size());
  }
  if (seismic_parameters.GetDepthStormOutput()) {
    if (nmo) {
      NRLib::Volume volume_nmo = NRLib::Volume(seismic_parameters.seismicGeometry()->x0(),
                                               seismic_parameters.seismicGeometry()->y0(),
                                               z_0_[0] - seismic_parameters.seismicGeometry()->dz()/2,
                                               seismic_parameters.seismicGeometry()->xlength(),
                                               seismic_parameters.seismicGeometry()->ylength(),
                                               (z_0_.size() * seismic_parameters.seismicGeometry()->dz()),
                                               seismic_parameters.seismicGeometry()->angle());
      depthgrid_     = new NRLib::StormContGrid(volume_nmo, nx, ny, z_0_.size());
    }
    else {
      depthgrid_     = new NRLib::StormContGrid(volume, nx, ny, nz);
    }
  }
}

void Output::AddTrace(SeismicParameters     &seismic_parameters,
                         NRLib::Grid2D<double> &timegrid_pos,
                         NRLib::Grid2D<double> &prenmo_timegrid_pos,
                         NRLib::Grid2D<double> &timegrid_stack_pos,
                         NRLib::Grid2D<double> &depthgrid_pos,
                         NRLib::Grid2D<double> &depthgrid_stack_pos,
                         NRLib::Grid2D<double> &timeshiftgrid_pos,
                         NRLib::Grid2D<double> &timeshiftgrid_stack_pos,
                         NRLib::Grid2D<double> &twtx_reg,
                         double                 x,
                         double                 y,
                         size_t                 i,
                         size_t                 j)
{
  std::vector<double>               zero_vec(1);
  zero_vec[0] = 0;
  size_t nz = seismic_parameters.seismicGeometry()->nz();
  size_t nt = seismic_parameters.seismicGeometry()->nt();
  bool nmo = seismic_parameters.modelSettings()->GetNMOCorr();

  //write seismic time
  if (time_segy_ok_){
    seismic_parameters.seismicOutput()->WriteSegyGather(timegrid_pos, time_segy_, twt_0_, offset_vec_, true, x, y, nmo);
  }
  if (prenmo_time_segy_ok_){
    seismic_parameters.seismicOutput()->WriteSegyGather(prenmo_timegrid_pos, prenmo_time_segy_, twt_0_, offset_vec_, true, x, y, nmo);
  }
  if (time_stack_segy_ok_) {
    seismic_parameters.seismicOutput()->WriteSegyGather(timegrid_stack_pos, time_stack_segy_, twt_0_, zero_vec, true, x, y, nmo);
  }
  //write twtx
  if (twtx_segy_ok_) {
    seismic_parameters.seismicOutput()->WriteSegyGather(twtx_reg, twtx_segy_, twt_0_, offset_vec_, true, x, y, nmo);
  }
  //write seismic depth
  if (depth_segy_ok_) {
    seismic_parameters.seismicOutput()->WriteSegyGather(depthgrid_pos, depth_segy_, z_0_, offset_vec_, false, x, y, nmo);
  }
  if (depth_stack_segy_ok_) {
    seismic_parameters.seismicOutput()->WriteSegyGather(depthgrid_stack_pos, depth_stack_segy_, z_0_, zero_vec, false, x, y, nmo);
  }
  //write seismic timeshift
  if (timeshift_segy_ok_) {
    seismic_parameters.seismicOutput()->WriteSegyGather(timeshiftgrid_pos, timeshift_segy_, twts_0_, offset_vec_, true, x, y, nmo);
  }

  if (timeshift_stack_segy_ok_){
    seismic_parameters.seismicOutput()->WriteSegyGather(timeshiftgrid_stack_pos, timeshift_stack_segy_, twts_0_, zero_vec, true, x, y, nmo);
  }
  //save to storm grid for output, print storm when finish loop
  if (seismic_parameters.GetTimeStormOutput()) {
    for (size_t k = 0; k < timegrid_->GetNK(); ++k){
      (*timegrid_)(i, j, k) = float(timegrid_stack_pos(k,0));
    }
  }
  if (seismic_parameters.GetDepthStormOutput()) {
    for (size_t k = 0; k < depthgrid_->GetNK(); ++k){
      (*depthgrid_)(i, j, k) = float(depthgrid_stack_pos(k,0));
    }
  }
  if (seismic_parameters.GetTimeshiftStormOutput()) {
    for (size_t k = 0; k < timeshiftgrid_->GetNK(); ++k){
      (*timeshiftgrid_)(i, j, k) = float(timeshiftgrid_stack_pos(k,0));
    }
  }
}


void Output::AddZeroTrace(SeismicParameters     &seismic_parameters,
                           double                 x, 
                           double                 y,
                           size_t                 i,
                           size_t                 j)
{
  std::vector<double>               zero_vec(1);
  zero_vec[0] = 0;
  size_t nz = seismic_parameters.seismicGeometry()->nz();
  size_t nt = seismic_parameters.seismicGeometry()->nt();
  bool nmo = seismic_parameters.modelSettings()->GetNMOCorr();

  if (time_segy_ok_){
    seismic_parameters.seismicOutput()->WriteZeroSegyGather(time_segy_, offset_vec_, x, y, nmo);
  }
  if (prenmo_time_segy_ok_){
    seismic_parameters.seismicOutput()->WriteZeroSegyGather(prenmo_time_segy_, offset_vec_, x, y, nmo);
  }
  if (time_stack_segy_ok_) {
    seismic_parameters.seismicOutput()->WriteZeroSegyGather(time_stack_segy_, zero_vec, x, y, nmo);
  }
  if (depth_segy_ok_) {
    seismic_parameters.seismicOutput()->WriteZeroSegyGather(depth_segy_, offset_vec_, x, y, nmo);
  }
  if (depth_stack_segy_ok_) {
    seismic_parameters.seismicOutput()->WriteZeroSegyGather(depth_stack_segy_, zero_vec, x, y, nmo);
  }
  if (timeshift_segy_ok_) {
    seismic_parameters.seismicOutput()->WriteZeroSegyGather(timeshift_segy_, offset_vec_, x, y, nmo);
  }
  if (timeshift_stack_segy_ok_){
    seismic_parameters.seismicOutput()->WriteZeroSegyGather(timeshift_stack_segy_, zero_vec, x, y, nmo);
  }
  if (twtx_segy_ok_) {
    seismic_parameters.seismicOutput()->WriteZeroSegyGather(twtx_segy_, offset_vec_, x, y, nmo);
  }
  if (seismic_parameters.GetTimeStormOutput()) {
    for (size_t k = 0; k < timegrid_->GetNK(); ++k){
      (*timegrid_)(i, j, k) = 0.0;
    }
  }
  if (seismic_parameters.GetDepthStormOutput()) {
    for (size_t k = 0; k < depthgrid_->GetNK(); ++k){
      (*depthgrid_)(i, j, k) = 0.0;
    }
  }
  if (seismic_parameters.GetTimeshiftStormOutput()) {
    for (size_t k = 0; k < timeshiftgrid_->GetNK(); ++k){
      (*timeshiftgrid_)(i, j, k) = 0.0;
    }
  }
}

void Output::WriteSeismicStorm(SeismicParameters     &seismic_parameters)
{
  if (seismic_parameters.GetTimeStormOutput()) {
    seismic_parameters.seismicOutput()->WriteSeismicTimeStorm(seismic_parameters, (*timegrid_), 0, true);
    timegrid_ = NULL;
  }
  if (seismic_parameters.GetDepthStormOutput()) {
    seismic_parameters.seismicOutput()->WriteSeismicDepthStorm(seismic_parameters, (*depthgrid_), 0, true);
    depthgrid_ = NULL;
  }
  if (seismic_parameters.GetTimeshiftStormOutput()) {
    seismic_parameters.seismicOutput()->WriteSeismicTimeshiftStorm(seismic_parameters, (*timeshiftgrid_), 0, true);
    timeshiftgrid_ = NULL;
  }
  //write reflections
  if (seismic_parameters.modelSettings()->GetOutputReflections()) {
    seismic_parameters.seismicOutput()->WriteReflections(seismic_parameters, offset_vec_[0]);
  }
}
