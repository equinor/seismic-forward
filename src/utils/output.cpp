#include "utils/result_trace.hpp"

#include "seismic_geometry.hpp"
#include "modelsettings.hpp"
#include "output.hpp"

Output::Output(SeismicParameters   & seismic_parameters,
               ModelSettings       * model_settings,
               std::vector<double>   twt_0,
               std::vector<double>   z_0,
               std::vector<double>   twts_0,
               std::vector<double>   offset_vec,
               size_t                time_samples_stretch)
  : segy_ok_                ( false      ),
    time_segy_ok_           ( false      ),
    prenmo_time_segy_ok_    ( false      ),
    time_stack_segy_ok_     ( false      ),
    depth_segy_ok_          ( false      ),
    depth_stack_segy_ok_    ( false      ),
    timeshift_segy_ok_      ( false      ),
    timeshift_stack_segy_ok_( false      ),
    twtx_segy_ok_           ( false      ),
    twt_0_                  ( twt_0      ),
    z_0_                    ( z_0        ),
    twts_0_                 ( twts_0     ),
    offset_vec_             ( offset_vec ),
    timegrid_               ( NULL       ),
    timeshiftgrid_          ( NULL       ),
    depthgrid_              ( NULL       )
{
  SeismicGeometry * seismic_geometry = seismic_parameters.GetSeismicGeometry();
  SeismicOutput   * seismic_output   = seismic_parameters.GetSeismicOutput();

  size_t            nx               = seismic_geometry->nx();
  size_t            ny               = seismic_geometry->ny();
  size_t            nz               = seismic_geometry->nz();
  size_t            nt               = seismic_geometry->nt();
  NRLib::Volume     volume           = seismic_geometry->createDepthVolume();
  NRLib::Volume     volume_t         = seismic_geometry->createTimeVolume();

  bool              nmo              = model_settings->GetNMOCorr();
  size_t            n_off            = offset_vec_.size();
  double            x0               = seismic_geometry->x0();
  double            y0               = seismic_geometry->y0();
  double            dt               = seismic_geometry->dt();
  double            dz               = seismic_geometry->dz();
  double            xlength          = seismic_geometry->xlength();
  double            ylength          = seismic_geometry->ylength();
  double            angle            = seismic_geometry->angle();

  if (model_settings->GetSegyOutput()) {
    seismic_output->SetSegyGeometry(seismic_parameters, volume_t, nx, ny);
    segy_ok_ = seismic_output->CheckUTMPrecision(seismic_parameters, volume_t, nx, ny);
  }

  if (segy_ok_) {
    if (model_settings->GetOutputTimeSegy()                 ) time_segy_ok_            = seismic_output->PrepareSegy(time_segy_           , twt_0_ , time_samples_stretch, "seismic_time"           , seismic_parameters, offset_vec_, n_off, true , nmo);
    if (model_settings->GetOutputSeismicStackTimeSegy()     ) time_stack_segy_ok_      = seismic_output->PrepareSegy(time_stack_segy_     , twt_0_ , time_samples_stretch, "seismic_time_stack"     , seismic_parameters, offset_vec_, 1    , true , nmo);
    if (model_settings->GetOutputPrenmoTimeSegy()           ) prenmo_time_segy_ok_     = seismic_output->PrepareSegy(prenmo_time_segy_    , twt_0_ , twt_0_.size()       , "seismic_time_prenmo"    , seismic_parameters, offset_vec_, n_off, true , nmo);
    if (model_settings->GetOutputTwtOffset()                ) twtx_segy_ok_            = seismic_output->PrepareSegy(twtx_segy_           , twt_0_ , twt_0_.size()       , "twt_offset"             , seismic_parameters, offset_vec_, n_off, true , nmo);
    if (model_settings->GetOutputDepthSegy()                ) depth_segy_ok_           = seismic_output->PrepareSegy(depth_segy_          , z_0_   , z_0_.size()         , "seismic_depth"          , seismic_parameters, offset_vec_, n_off, false, nmo);
    if (model_settings->GetOutputSeismicStackDepthSegy()    ) depth_stack_segy_ok_     = seismic_output->PrepareSegy(depth_stack_segy_    , z_0_   , z_0_.size()         , "seismic_depth_stack"    , seismic_parameters, offset_vec_, 1    , false, nmo);
    if (model_settings->GetOutputTimeshiftSegy()            ) timeshift_segy_ok_       = seismic_output->PrepareSegy(timeshift_segy_      , twts_0_, twts_0_.size()      , "seismic_timeshift"      , seismic_parameters, offset_vec_, n_off, true , nmo);
    if (model_settings->GetOutputSeismicStackTimeShiftSegy()) timeshift_stack_segy_ok_ = seismic_output->PrepareSegy(timeshift_stack_segy_, twts_0_, twts_0_.size()      , "seismic_timeshift_stack", seismic_parameters, offset_vec_, 1    , true , nmo);
  }

  //prepare grid if output of seismic in storm is requested

  if (model_settings->GetTimeOutput()) {
    if (nmo) {
      NRLib::Volume volume_t_nmo = NRLib::Volume(x0, y0, twt_0_[0] - dt/2, xlength, ylength, time_samples_stretch * dt, angle);
      timegrid_ = new NRLib::StormContGrid(volume_t_nmo, nx, ny, time_samples_stretch);
    }
    else {
      timegrid_ = new NRLib::StormContGrid(volume_t, nx, ny, nt);
    }
  }

  if (model_settings->GetTimeshiftOutput()) {
    NRLib::Volume volume_ts_nmo = NRLib::Volume(x0, y0, twts_0_[0] - dt/2, xlength, ylength, (twts_0_.size() * dt), angle);
    timeshiftgrid_ = new NRLib::StormContGrid(volume_ts_nmo, nx, ny, twts_0_.size());
  }

  if (model_settings->GetDepthOutput()) {
    if (nmo) {
      NRLib::Volume volume_nmo = NRLib::Volume(x0, y0, (z_0_[0] - dz/2), xlength, ylength, z_0_.size()*dz, angle);
      depthgrid_ = new NRLib::StormContGrid(volume_nmo, nx, ny, z_0_.size());
    }
    else {
      depthgrid_ = new NRLib::StormContGrid(volume    , nx, ny, nz);
    }
  }
}

//--------------------
Output::~Output(void)
//--------------------
{
  if (timegrid_ != NULL)
    delete timegrid_;

  if (timeshiftgrid_ != NULL)
    delete timeshiftgrid_;

  if (depthgrid_ != NULL)
    delete depthgrid_;
}

//---------------------------------------------------
void Output::AddTrace(ResultTrace   * result_trace,
                      ModelSettings * model_settings,
                      SeismicOutput * seismic_output)
//---------------------------------------------------
{
  std::vector<double> zero_vec(1, 0.0);
  bool                nmo            = model_settings->GetNMOCorr();
  double              x              = result_trace->GetX();
  double              y              = result_trace->GetY();
  size_t              i              = result_trace->GetI();
  size_t              j              = result_trace->GetJ();

  if (result_trace->GetIsEmpty()) {
    if (time_segy_ok_           ) seismic_output->WriteZeroSegyGather(time_segy_             , offset_vec_, x, y, nmo);
    if (prenmo_time_segy_ok_    ) seismic_output->WriteZeroSegyGather(prenmo_time_segy_      , offset_vec_, x, y, nmo);
    if (time_stack_segy_ok_     ) seismic_output->WriteZeroSegyGather(time_stack_segy_       , zero_vec   , x, y, nmo);
    if (depth_segy_ok_          ) seismic_output->WriteZeroSegyGather(depth_segy_            , offset_vec_, x, y, nmo);
    if (depth_stack_segy_ok_    ) seismic_output->WriteZeroSegyGather(depth_stack_segy_      , zero_vec   , x, y, nmo);
    if (timeshift_segy_ok_      ) seismic_output->WriteZeroSegyGather(timeshift_segy_        , offset_vec_, x, y, nmo);
    if (timeshift_stack_segy_ok_) seismic_output->WriteZeroSegyGather(timeshift_stack_segy_  , zero_vec   , x, y, nmo);
    if (twtx_segy_ok_           ) seismic_output->WriteZeroSegyGather(twtx_segy_             , offset_vec_, x, y, nmo);
  }
  else {
    if (time_segy_ok_           ) seismic_output->WriteSegyGather(result_trace->GetTimeTrace()          , time_segy_           , twt_0_ , offset_vec_, true , x, y, nmo);
    if (prenmo_time_segy_ok_    ) seismic_output->WriteSegyGather(result_trace->GetPreNMOTimeTrace()    , prenmo_time_segy_    , twt_0_ , offset_vec_, true , x, y, nmo);
    if (time_stack_segy_ok_     ) seismic_output->WriteSegyGather(result_trace->GetTimeStackTrace()     , time_stack_segy_     , twt_0_ , zero_vec   , true , x, y, nmo);
    if (twtx_segy_ok_           ) seismic_output->WriteSegyGather(result_trace->GetTWTxReg()            , twtx_segy_           , twt_0_ , offset_vec_, true , x, y, nmo);
    if (depth_segy_ok_          ) seismic_output->WriteSegyGather(result_trace->GetDepthTrace()         , depth_segy_          , z_0_   , offset_vec_, false, x, y, nmo);
    if (depth_stack_segy_ok_    ) seismic_output->WriteSegyGather(result_trace->GetDepthStackTrace()    , depth_stack_segy_    , z_0_   , zero_vec   , false, x, y, nmo);
    if (timeshift_segy_ok_      ) seismic_output->WriteSegyGather(result_trace->GetTimeShiftTrace()     , timeshift_segy_      , twts_0_, offset_vec_, true , x, y, nmo);
    if (timeshift_stack_segy_ok_) seismic_output->WriteSegyGather(result_trace->GetTimeShiftStackTrace(), timeshift_stack_segy_, twts_0_, zero_vec   , true , x, y, nmo);
  }

  //
  // Save to storm grid for output, print storm when finish loop
  //
  if (model_settings->GetTimeOutput()) {
    for (size_t k = 0 ; k < timegrid_->GetNK() ; ++k) {
      if (result_trace->GetIsEmpty() || !time_stack_segy_ok_)
        (*timegrid_)(i, j, k) = 0.0f;
      else
        (*timegrid_)(i, j, k) = float(result_trace->GetTimeStackTrace()(k,0));
    }
  }
  if (model_settings->GetTimeshiftOutput()) {
    for (size_t k = 0 ; k < timeshiftgrid_->GetNK() ; ++k) {
      if (result_trace->GetIsEmpty() || !timeshift_stack_segy_ok_)
        (*timeshiftgrid_)(i, j, k) = 0.0f;
      else
        (*timeshiftgrid_)(i, j, k) = float(result_trace->GetTimeShiftStackTrace()(k,0));
    }
  }
  if (model_settings->GetDepthOutput()) {
    for (size_t k = 0 ; k < depthgrid_->GetNK() ; ++k) {
      if (result_trace->GetIsEmpty() || !depth_stack_segy_ok_)
        (*depthgrid_)(i, j, k) = 0.0f;
      else
        (*depthgrid_)(i, j, k) = float(result_trace->GetDepthStackTrace()(k,0));
    }
  }
}

void Output::WriteStatisticsForSeismic(ModelSettings * model_settings)
{
  if (model_settings->GetTimeOutput() || model_settings->GetTimeshiftOutput() || model_settings->GetDepthOutput()) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nStatistics for generated seismic data.\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nType               Avg         Min         Max");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n----------------------------------------------\n");
  }
  float min, max, avg;
  if (model_settings->GetTimeOutput()) {
    timegrid_->GetAvgMinMaxWithMissing(avg, min , max, timegrid_->GetMissingCode());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Time       %11.4f %11.4f %11.4f\n", avg, min, max);
  }
  if (model_settings->GetTimeshiftOutput()) {
    timeshiftgrid_->GetAvgMinMaxWithMissing(avg, min , max, timeshiftgrid_->GetMissingCode());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Time shift %11.4f %11.4f %11.4f\n", avg, min, max);
  }
  if (model_settings->GetDepthOutput()) {
    depthgrid_->GetAvgMinMaxWithMissing(avg, min , max, depthgrid_->GetMissingCode());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Depth      %11.4f %11.4f %11.4f\n", avg, min, max);
  }
}

void Output::WriteSeismicStorm(ModelSettings                     * model_settings,
                               SeismicOutput                     * seismic_output,
                               std::vector<NRLib::StormContGrid> & rgrids)
{
  if (model_settings->GetTimeStormOutput())      { seismic_output->WriteSeismicTimeStorm     (*timegrid_     , 0, true); timegrid_      = NULL; }
  if (model_settings->GetDepthStormOutput())     { seismic_output->WriteSeismicDepthStorm    (*depthgrid_    , 0, true); depthgrid_     = NULL; }
  if (model_settings->GetTimeshiftStormOutput()) { seismic_output->WriteSeismicTimeshiftStorm(*timeshiftgrid_, 0, true); timeshiftgrid_ = NULL; }
  if (model_settings->GetOutputReflections())    { seismic_output->WriteReflections(rgrids, offset_vec_[0]); }
}
