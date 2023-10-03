#include "result_trace.hpp"
#include <seismic_geometry.hpp>

//------------------------------------------------------
ResultTrace::ResultTrace(ModelSettings * model_settings,
                         const size_t    nzrefl,
                         const size_t    nt,
                         const size_t    nz,
                         const size_t    ntwts0,
                         const size_t    nt_stretch,
                         const size_t    noff)
//------------------------------------------------------
  : empty_(false)
{
  if (model_settings->GetNMOCorr()) {
    twtx_reg_.           Resize(nt    , noff);
    twtx_.               Resize(nzrefl, noff);
    theta_.              Resize(nzrefl, noff);
    refl_.               Resize(nzrefl, noff);
    prenmo_timegrid_pos_.Resize(nt    , noff);

    if (model_settings->GetPSSeismic()) {
      offset_pp_.    Resize(nzrefl, noff);
      offset_ss_.    Resize(nzrefl, noff);
      offset_pp_reg_.Resize(nt    , noff);
      offset_ss_reg_.Resize(nt    , noff);
    }
  }
  timegrid_pos_.Resize(nt_stretch, noff);

  if (model_settings->GetStackOutput() || model_settings->GetStormOutput()) {
    timegrid_stack_pos_.Resize(nt_stretch, 1);
  }
  if (model_settings->GetTimeshiftOutput()) {
    timeshiftgrid_pos_.      Resize(ntwts0, noff);
    timeshiftgrid_stack_pos_.Resize(ntwts0, 1);
  }
  if (model_settings->GetDepthOutput()){
    depthgrid_pos_.      Resize(nz, noff);
    depthgrid_stack_pos_.Resize(nz, 1);
  }
}

void ResultTrace::SetJobID(Trace *trace)
{
  x_          = trace->GetX();
  y_          = trace->GetY();
  i_          = trace->GetI();
  j_          = trace->GetJ();
  job_number_ = trace->GetJobNumber();
}
