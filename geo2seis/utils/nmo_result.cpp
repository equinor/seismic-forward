#include "nmo_result.hpp"



NMOResult::NMOResult(SeismicParameters         &seismic_parameters,
                     std::vector<double>        twt_0,
                     std::vector<double>        z_0,
                     std::vector<double>        twts_0,
                     size_t                     time_samples_stretch,
                     const std::vector<double> &offset_vec)
  : empty_(false)
{
  twtx_reg_.        Resize(twt_0.size(),         offset_vec.size());
  timegrid_pos_.    Resize(twt_0.size(),         offset_vec.size());
  nmo_timegrid_pos_.Resize(time_samples_stretch, offset_vec.size());

  if (seismic_parameters.GetStackOutput() || seismic_parameters.GetStormOutput()) {
    nmo_timegrid_stack_pos_.Resize(time_samples_stretch, 1);
  }
  if (seismic_parameters.GetTimeshiftOutput()) {
    nmo_timeshiftgrid_pos_.      Resize(twts_0.size(), offset_vec.size());
    nmo_timeshiftgrid_stack_pos_.Resize(twts_0.size(), 1);
  }
  if (seismic_parameters.GetDepthOutput()){
    nmo_depthgrid_pos_.      Resize(z_0.size(), offset_vec.size());
    nmo_depthgrid_stack_pos_.Resize(z_0.size(), 1);
  }
}


void NMOResult::SetJobID(Trace *trace)
{
  x_          = trace->GetX();
  y_          = trace->GetY();
  i_          = trace->GetI();
  j_          = trace->GetJ();
  job_number_ = trace->GetJobNumber();
}
