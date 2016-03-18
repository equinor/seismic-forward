#ifndef GEN_RESAMPL_PARAM_HPP
#define GEN_RESAMPL_PARAM_HPP

//#include <seismic_output.hpp>
#include <seismic_parameters.hpp>
//#include "modelsettings.hpp"
//#include "nrlib/segy/segygeometry.hpp"
//#include <seismic_geometry.hpp>
//#include <vector>
#include <nrlib/stormgrid/stormcontgrid.hpp>
#include <nrlib/surface/regularsurface.hpp>
#include "tbb/concurrent_queue.h"
#include <utils/resampl_trace.hpp>
#include <utils/trace.hpp>

//class SeismicParameters;
//class ResampleTrace;

class GenResamplParam {
  public:
    GenResamplParam(SeismicParameters                   &seismic_parameters,
                    const std::vector<double>           &time_or_depth_vec_reg,
                    const NRLib::StormContGrid          &time_or_depth_grid,
                    const NRLib::RegularSurface<double> &toptime,
                    size_t                               n_samples,
                    size_t                               n_traces,
                    bool                                 time,
                    tbb::concurrent_queue<ResamplTrace*>         &empty_queue,
                    tbb::concurrent_bounded_queue<ResamplTrace*> &result_queue,
                    tbb::concurrent_queue<Trace*>                &traces);

    ~GenResamplParam() {};

    SeismicParameters             seismic_parameters;
    std::vector<double>           time_or_depth_vec_reg;
    NRLib::StormContGrid          time_or_depth_grid;
    NRLib::RegularSurface<double> toptime;
    size_t                        n_samples;
    size_t                        n_traces;
    bool                          time;
    tbb::concurrent_queue<ResamplTrace*>         empty_queue;
    tbb::concurrent_bounded_queue<ResamplTrace*> result_queue;
    tbb::concurrent_queue<Trace*>                traces;

  };

#endif
