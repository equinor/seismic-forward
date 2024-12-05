#include "gen_seis_trace_params.hpp"

GenSeisTraceParams::GenSeisTraceParams(tbb::concurrent_queue<ResultTrace*>         & empty_queue,
                                       tbb::concurrent_bounded_queue<ResultTrace*> & result_queue,
                                       tbb::concurrent_queue<Trace*>               & seismic_traces,
                                       size_t                                        n_traces)
: empty_queue(empty_queue),
  result_queue(result_queue),
  seismic_traces(seismic_traces),
  n_traces(n_traces)
{
}
