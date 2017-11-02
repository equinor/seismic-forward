#include "result_trace.hpp"
#include <seismic_geometry.hpp>


ResultTrace::ResultTrace(SeismicParameters         & seismic_parameters,
                         std::vector<double>         twt_0,
                         std::vector<double>         z_0,
                         std::vector<double>         twts_0,
                         size_t                      time_samples_stretch,
                         const std::vector<double> & offset_vec)
  : empty_(false)
{
  size_t nzrefl = seismic_parameters.GetSeismicGeometry()->zreflectorcount();
  if (seismic_parameters.GetModelSettings()->GetNMOCorr()) {

    twtx_reg_.           Resize(twt_0.size(), offset_vec.size());
    twtx_.               Resize(nzrefl,       offset_vec.size());
    theta_.              Resize(nzrefl,       offset_vec.size());
    refl_.               Resize(nzrefl,       offset_vec.size());
    prenmo_timegrid_pos_.Resize(twt_0.size(), offset_vec.size());

    if (seismic_parameters.GetModelSettings()->GetPSSeismic()) {
      offset_pp_.    Resize(nzrefl, offset_vec.size());
      offset_ss_.    Resize(nzrefl, offset_vec.size());
      offset_pp_reg_.Resize(twt_0.size(), offset_vec.size());
      offset_ss_reg_.Resize(twt_0.size(), offset_vec.size());
    }
  }
  timegrid_pos_.Resize(time_samples_stretch, offset_vec.size());

  if (seismic_parameters.GetStackOutput() || seismic_parameters.GetStormOutput()) {
    timegrid_stack_pos_.Resize(time_samples_stretch, 1);
  }
  if (seismic_parameters.GetTimeshiftOutput()) {
    timeshiftgrid_pos_.      Resize(twts_0.size(), offset_vec.size());
    timeshiftgrid_stack_pos_.Resize(twts_0.size(), 1);
  }
  if (seismic_parameters.GetDepthOutput()){
    depthgrid_pos_.      Resize(z_0.size(), offset_vec.size());
    depthgrid_stack_pos_.Resize(z_0.size(), 1);
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
