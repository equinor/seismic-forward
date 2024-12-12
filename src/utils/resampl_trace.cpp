
#include "resampl_trace.hpp"

ResamplTrace::ResamplTrace(std::vector<NRLib::Grid2D<double>> & traces)
  : traces_(traces),
    empty_(false)
{
}

void ResamplTrace::SetJobID(Trace * trace)
{
  i_          = trace->GetI();
  j_          = trace->GetJ();
  job_number_ = trace->GetJobNumber();
}
