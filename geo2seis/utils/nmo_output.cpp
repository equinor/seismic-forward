#include "nmo_output.hpp"


#include <seismic_geometry.hpp>
#include "modelsettings.hpp"


NMOOutput::NMOOutput(SeismicParameters &seismic_parameters,
                     std::vector<double> twt_0,
                     std::vector<double> z_0,
                     std::vector<double> twts_0,
                     size_t time_samples_stretch)
  : segy_ok_(false),
    nmo_time_segy_ok_(false),
    prenmo_time_segy_ok_(false),
    nmo_time_stack_segy_ok_(false),
    nmo_depth_segy_ok_(false),
    nmo_depth_stack_segy_ok_(false),
    nmo_timeshift_segy_ok_(false),
    nmo_timeshift_stack_segy_ok_(false),
    twtx_segy_ok_(false),
    twt_0_(twt_0),
    z_0_(z_0),
    twts_0_(twts_0)
{
  size_t nx = seismic_parameters.seismicGeometry()->nx();
  size_t ny = seismic_parameters.seismicGeometry()->ny();
  size_t nz = seismic_parameters.seismicGeometry()->nz();
  size_t nt = seismic_parameters.seismicGeometry()->nt();
  std::vector<double> & offset_vec = seismic_parameters.GetOffsetVec();
  NRLib::Volume volume   = seismic_parameters.seismicGeometry()->createDepthVolume();
  NRLib::Volume volume_t = seismic_parameters.seismicGeometry()->createTimeVolume();

  if (seismic_parameters.GetSegyOutput()) {
    seismic_parameters.seismicOutput()->setSegyGeometry(seismic_parameters, volume_t, nx, ny);
    segy_ok_ = seismic_parameters.seismicOutput()->checkUTMPrecision(seismic_parameters, volume_t, nx, ny);
  }
  if (segy_ok_) {
    if (seismic_parameters.modelSettings()->GetOutputTimeSegy()) {
      std::string filename        = "seismic_time";
      nmo_time_segy_ok_            = seismic_parameters.seismicOutput()->prepareSegy(nmo_time_segy_, twt_0_, filename, seismic_parameters, offset_vec, offset_vec.size(), true, true);
    }
    if (seismic_parameters.modelSettings()->GetOutputPrenmoTimeSegy()) {
      std::string filename        = "prenmo_seismic_time";
      prenmo_time_segy_ok_         = seismic_parameters.seismicOutput()->prepareSegy(prenmo_time_segy_, twt_0_, filename, seismic_parameters, offset_vec, offset_vec.size(), true, true);
    }
    if (seismic_parameters.modelSettings()->GetOutputSeismicStackTimeSegy()) {
      std::string filename        = "seismic_time_stack";
      nmo_time_stack_segy_ok_      = seismic_parameters.seismicOutput()->prepareSegy(nmo_time_stack_segy_, twt_0_, filename, seismic_parameters, offset_vec, 1, true, true);
    }
    if (seismic_parameters.modelSettings()->GetOutputDepthSegy()) {
      std::string filename        = "seismic_depth";
      nmo_depth_segy_ok_           = seismic_parameters.seismicOutput()->prepareSegy(nmo_depth_segy_, z_0_, filename, seismic_parameters, offset_vec, offset_vec.size(), false, true);
    }
    if (seismic_parameters.modelSettings()->GetOutputSeismicStackDepthSegy()) {
      std::string filename        = "seismic_depth_stack";
      nmo_depth_stack_segy_ok_     = seismic_parameters.seismicOutput()->prepareSegy(nmo_depth_stack_segy_, z_0_, filename, seismic_parameters, offset_vec, 1, false, true);
    }
    if (seismic_parameters.modelSettings()->GetOutputTimeshiftSegy()) {
      std::string filename        = "seismic_timeshift";
      nmo_timeshift_segy_ok_       = seismic_parameters.seismicOutput()->prepareSegy(nmo_timeshift_segy_, twts_0_, filename, seismic_parameters, offset_vec, offset_vec.size(), true, true);
    }
    if (seismic_parameters.modelSettings()->GetOutputSeismicStackTimeShiftSegy()) {
      std::string filename        = "seismic_timeshift_stack";
      nmo_timeshift_stack_segy_ok_ = seismic_parameters.seismicOutput()->prepareSegy(nmo_timeshift_stack_segy_, twts_0_, filename, seismic_parameters, offset_vec, 1, true, true);
    }
    if (seismic_parameters.modelSettings()->GetOutputTwtOffset()) {
      std::string filename        = "twt_offset";
      twtx_segy_ok_                = seismic_parameters.seismicOutput()->prepareSegy(twtx_segy_, twt_0_, filename, seismic_parameters, offset_vec, offset_vec.size(), true, true);
    }
  }
  
  //prepare grid if output of seismic in storm is requested
  if (seismic_parameters.GetTimeStormOutput()) {
    NRLib::Volume volume_t_nmo = NRLib::Volume(seismic_parameters.seismicGeometry()->x0(),
                                               seismic_parameters.seismicGeometry()->y0(),
                                               twt_0_[0] - seismic_parameters.seismicGeometry()->dt()/2,
                                               seismic_parameters.seismicGeometry()->xlength(),
                                               seismic_parameters.seismicGeometry()->ylength(),
                                               (time_samples_stretch * seismic_parameters.seismicGeometry()->dt()),
                                               seismic_parameters.seismicGeometry()->angle());

    timegrid_      = new NRLib::StormContGrid(volume_t_nmo, nx, ny, time_samples_stretch);
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
    NRLib::Volume volume_nmo = NRLib::Volume(seismic_parameters.seismicGeometry()->x0(),
                                             seismic_parameters.seismicGeometry()->y0(),
                                             z_0_[0] - seismic_parameters.seismicGeometry()->dz()/2,
                                             seismic_parameters.seismicGeometry()->xlength(),
                                             seismic_parameters.seismicGeometry()->ylength(),
                                             (z_0_.size() * seismic_parameters.seismicGeometry()->dz()),
                                             seismic_parameters.seismicGeometry()->angle());
    depthgrid_     = new NRLib::StormContGrid(volume_nmo, nx, ny, z_0_.size());
  }
}

void NMOOutput::AddTrace(SeismicParameters     &seismic_parameters,
                         NRLib::Grid2D<double> &timegrid_pos,
                         NRLib::Grid2D<double> &nmo_timegrid_pos,
                         NRLib::Grid2D<double> &nmo_timegrid_stack_pos,
                         NRLib::Grid2D<double> &nmo_depthgrid_pos,
                         NRLib::Grid2D<double> &nmo_depthgrid_stack_pos,
                         NRLib::Grid2D<double> &nmo_timeshiftgrid_pos,
                         NRLib::Grid2D<double> &nmo_timeshiftgrid_stack_pos,
                         NRLib::Grid2D<double> &twtx_reg,
                         double                 x, 
                         double                 y,
                         size_t                 i,
                         size_t                 j)
{
  const std::vector<double> & offset_vec = seismic_parameters.GetOffsetVec();
  std::vector<double>               zero_vec(1);
  zero_vec[0] = 0;
  size_t nz = seismic_parameters.seismicGeometry()->nz();
  size_t nt = seismic_parameters.seismicGeometry()->nt();

  //write seismic time
  if (nmo_time_segy_ok_){
    seismic_parameters.seismicOutput()->writeSegyGather(nmo_timegrid_pos, nmo_time_segy_, twt_0_, offset_vec, true, x,y);
  }
  if (prenmo_time_segy_ok_){
    seismic_parameters.seismicOutput()->writeSegyGather(timegrid_pos, prenmo_time_segy_, twt_0_, offset_vec, true, x,y);
  }
  if (nmo_time_stack_segy_ok_) {
    seismic_parameters.seismicOutput()->writeSegyGather(nmo_timegrid_stack_pos, nmo_time_stack_segy_, twt_0_, zero_vec, true, x,y);
  }
  //write twtx
  if (twtx_segy_ok_) {
    seismic_parameters.seismicOutput()->writeSegyGather(twtx_reg, twtx_segy_, twt_0_, offset_vec, true, x,y);
  }
  //write seismic depth
  if (nmo_depth_segy_ok_) {
    seismic_parameters.seismicOutput()->writeSegyGather(nmo_depthgrid_pos, nmo_depth_segy_, z_0_, offset_vec, false, x,y);
  }
  if (nmo_depth_stack_segy_ok_) {
    seismic_parameters.seismicOutput()->writeSegyGather(nmo_depthgrid_stack_pos, nmo_depth_stack_segy_, z_0_, zero_vec, false, x,y);
  }
  //write seismic timeshift
  if (nmo_timeshift_segy_ok_) {
    seismic_parameters.seismicOutput()->writeSegyGather(nmo_timeshiftgrid_pos, nmo_timeshift_segy_, twts_0_, offset_vec, true, x,y);
  }

  if (nmo_timeshift_stack_segy_ok_){
    seismic_parameters.seismicOutput()->writeSegyGather(nmo_timeshiftgrid_stack_pos, nmo_timeshift_stack_segy_, twts_0_, zero_vec, true, x,y);
  }
  //save to storm grid for output, print storm when finish loop
  if (seismic_parameters.GetTimeStormOutput()) {
    for (size_t k = 0; k < timegrid_->GetNK(); ++k){
      (*timegrid_)(i, j, k) = float(nmo_timegrid_stack_pos(k,0));
    }
  }
  if (seismic_parameters.GetDepthStormOutput()) {
    for (size_t k = 0; k < depthgrid_->GetNK(); ++k){
      (*depthgrid_)(i, j, k) = float(nmo_depthgrid_stack_pos(k,0));
    }
  }
  if (seismic_parameters.GetTimeshiftStormOutput()) {
    for (size_t k = 0; k < timeshiftgrid_->GetNK(); ++k){
      (*timeshiftgrid_)(i, j, k) = float(nmo_timeshiftgrid_stack_pos(k,0));
    }
  }
}


void NMOOutput::AddZeroTrace(SeismicParameters     &seismic_parameters,
                           double                 x, 
                           double                 y,
                           size_t                 i,
                           size_t                 j)
{
  std::vector<double> & offset_vec = seismic_parameters.GetOffsetVec();
  std::vector<double>               zero_vec(1);
  zero_vec[0] = 0;
  size_t nz = seismic_parameters.seismicGeometry()->nz();
  size_t nt = seismic_parameters.seismicGeometry()->nt();

  if (nmo_time_segy_ok_){
    seismic_parameters.seismicOutput()->writeZeroSegyGather(nmo_time_segy_, offset_vec, x,y);
  }
  if (prenmo_time_segy_ok_){
    seismic_parameters.seismicOutput()->writeZeroSegyGather(prenmo_time_segy_, offset_vec, x,y);
  }
  if (nmo_time_stack_segy_ok_) {
    seismic_parameters.seismicOutput()->writeZeroSegyGather(nmo_time_stack_segy_, zero_vec, x,y);
  }
  if (nmo_depth_segy_ok_) {
    seismic_parameters.seismicOutput()->writeZeroSegyGather(nmo_depth_segy_, offset_vec, x,y);
  }
  if (nmo_depth_stack_segy_ok_) {
    seismic_parameters.seismicOutput()->writeZeroSegyGather(nmo_depth_stack_segy_, zero_vec, x,y);
  }
  if (nmo_timeshift_segy_ok_) {
    seismic_parameters.seismicOutput()->writeZeroSegyGather(nmo_timeshift_segy_, offset_vec, x,y);
  }
  if (nmo_timeshift_stack_segy_ok_){
    seismic_parameters.seismicOutput()->writeZeroSegyGather(nmo_timeshift_stack_segy_, zero_vec, x,y);
  }
  if (twtx_segy_ok_) {
    seismic_parameters.seismicOutput()->writeZeroSegyGather(twtx_segy_, offset_vec, x,y);
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

void NMOOutput::WriteSeismicStorm(SeismicParameters     &seismic_parameters)
{
  if (seismic_parameters.GetTimeStormOutput()) {
    seismic_parameters.seismicOutput()->writeNMOSeismicTimeStorm(seismic_parameters, (*timegrid_), 0, true);
    timegrid_ = NULL;
  }
  if (seismic_parameters.GetDepthStormOutput()) {
    seismic_parameters.seismicOutput()->writeNMOSeismicDepthStorm(seismic_parameters, (*depthgrid_), 0, true);
    depthgrid_ = NULL;
  }
  if (seismic_parameters.GetTimeshiftStormOutput()) {
    seismic_parameters.seismicOutput()->writeNMOSeismicTimeshiftStorm(seismic_parameters, (*timeshiftgrid_), 0, true);
    timeshiftgrid_ = NULL;
  }
  //write reflections
  if (seismic_parameters.modelSettings()->GetOutputReflections()) {
    seismic_parameters.seismicOutput()->writeNMOReflections(seismic_parameters, 0.0);
  }
}