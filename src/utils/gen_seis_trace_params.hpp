#ifndef GEN_SEIS_TRACE_PARAMS_HPP
#define GEN_SEIS_TRACE_PARAMS_HPP

#include "utils/result_trace.hpp"
#include "utils/trace.hpp"

#include "seismic_parameters.hpp"
#include "modelsettings.hpp"

#include "tbb/concurrent_queue.h"

#include <vector>

class GenSeisTraceParams
{
  public:
  GenSeisTraceParams(tbb::concurrent_queue<ResultTrace*>         & empty_queue,
                     tbb::concurrent_bounded_queue<ResultTrace*> & result_queue,
                     tbb::concurrent_queue<Trace*>               & seismic_traces,
                     size_t                                        n_traces);

    ~GenSeisTraceParams() {};
    tbb::concurrent_queue<ResultTrace*>           empty_queue;
    tbb::concurrent_bounded_queue<ResultTrace*>   result_queue;
    tbb::concurrent_queue<Trace*>                 seismic_traces;
    size_t                                        n_traces;
};

#endif
