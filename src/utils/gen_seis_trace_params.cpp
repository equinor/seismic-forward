#include "gen_seis_trace_params.hpp"

GenSeisTraceParams::GenSeisTraceParams(SeismicParameters                           &seismic_parameters,
                                       const std::vector<double>                   &twt_0,
                                       const std::vector<double>                   &z_0,
                                       const std::vector<double>                   &twts_0,
                                       const std::vector<double>                   &theta_vec,
                                       const std::vector<double>                   &offset_vec,
                                       tbb::concurrent_queue<ResultTrace*>         &empty_queue,
                                       tbb::concurrent_bounded_queue<ResultTrace*> &result_queue,
                                       tbb::concurrent_queue<Trace*>               &seismic_traces,
                                       size_t                                       n_traces,
                                       size_t                                       time_samples_stretch)
: seismic_parameters(seismic_parameters),
  twt_0(twt_0),
  z_0(z_0),
  twts_0(twts_0),
  theta_vec(theta_vec),
  offset_vec(offset_vec),
  empty_queue(empty_queue),
  result_queue(result_queue),
  seismic_traces(seismic_traces),
  n_traces(n_traces),
  time_samples_stretch(time_samples_stretch)
{
}

