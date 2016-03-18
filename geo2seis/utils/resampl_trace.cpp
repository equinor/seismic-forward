#include "resampl_trace.hpp"

ResamplTrace::ResamplTrace(std::vector<NRLib::Grid2D<double> > &traces)
  : empty_(false),
    traces_(traces)
{
}

void ResamplTrace::SetJobID(Trace *trace)
{
  x_ = trace->GetX();
  y_ = trace->GetY();
  i_ = trace->GetI();
  j_ = trace->GetJ();
  job_number_ = trace->GetJobNumber();
}

