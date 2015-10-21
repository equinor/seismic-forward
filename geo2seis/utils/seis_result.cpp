#include "seis_result.hpp"


//#include "modelsettings.hpp"

SeisResult::SeisResult(SeismicParameters         &seismic_parameters,
                       std::vector<double>        twt_0,
                       std::vector<double>        z_0,
                       std::vector<double>        twts_0,
                       const std::vector<double> &theta_vec)
  : empty_(false)
{
  timegrid_pos_.Resize(twt_0.size(), theta_vec.size());

  if (seismic_parameters.GetStackOutput() || seismic_parameters.GetStormOutput()) {
    timegrid_stack_pos_.Resize(twt_0.size(), 1);
  }
  if (seismic_parameters.GetTimeshiftOutput()) {
    timeshiftgrid_pos_.      Resize(twts_0.size(), theta_vec.size());
    timeshiftgrid_stack_pos_.Resize(twts_0.size(), 1);
  }
  if (seismic_parameters.GetDepthOutput()){
    depthgrid_pos_.      Resize(z_0.size(), theta_vec.size());
    depthgrid_stack_pos_.Resize(z_0.size(), 1);
  }
}


void SeisResult::SetJobID(Trace *trace)
{
  x_          = trace->GetX();
  y_          = trace->GetY();
  i_          = trace->GetI();
  j_          = trace->GetJ();
  job_number_ = trace->GetJobNumber();
}
