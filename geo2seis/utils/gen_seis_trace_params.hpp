#ifndef GEN_SEIS_TRACE_PARAMS_HPP
#define GEN_SEIS_TRACE_PARAMS_HPP

#include <vector>
#include "modelsettings.hpp"
#include <seismic_parameters.hpp>
#include "tbb/concurrent_queue.h"
#include <utils/trace.hpp>
#include "utils/seis_output.hpp"
#include "utils/nmo_output.hpp"
#include "utils/result_trace.hpp"

class GenSeisTraceParams {
  public:
    GenSeisTraceParams(SeismicParameters                           &seismic_parameters,
                       const std::vector<double>                   &twt_0,
                       const std::vector<double>                   &z_0,
                       const std::vector<double>                   &twts_0,
                       const std::vector<double>                   &theta_vec,
                       const std::vector<double>                   &offset_vec,
                       tbb::concurrent_queue<ResultTrace*>         &empty_queue,
                       tbb::concurrent_bounded_queue<ResultTrace*> &result_queue,
                       tbb::concurrent_queue<Trace*>               &seismic_traces,
                       size_t                                       n_traces,
                       size_t                                       time_samples_stretch);

    ~GenSeisTraceParams() {};
    void SetNTraces(size_t n_trace) { n_traces = n_trace; };


    SeismicParameters                           seismic_parameters;
    const std::vector<double>                   twt_0;
    const std::vector<double>                   z_0;
    const std::vector<double>                   twts_0;
    const std::vector<double>                   theta_vec;
    const std::vector<double>                   offset_vec;
    tbb::concurrent_queue<ResultTrace*>         empty_queue;
    tbb::concurrent_bounded_queue<ResultTrace*> result_queue;
    tbb::concurrent_queue<Trace*>               seismic_traces;
    size_t                                      n_traces;
    size_t                                      time_samples_stretch;

};

#endif
