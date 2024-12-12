#ifndef GEN_RESAMPL_PARAM_HPP
#define GEN_RESAMPL_PARAM_HPP

#include <seismic_parameters.hpp>
#include <nrlib/stormgrid/stormcontgrid.hpp>
#include <nrlib/surface/regularsurface.hpp>
#include "tbb/concurrent_queue.h"
#include <utils/resampl_trace.hpp>
#include <utils/trace.hpp>

class GenResamplParam {
  public:
  GenResamplParam(tbb::concurrent_queue<ResamplTrace*>         & empty_queue,
                  tbb::concurrent_bounded_queue<ResamplTrace*> & result_queue);

  ~GenResamplParam() {};

  tbb::concurrent_queue<ResamplTrace*>         empty_queue;
  tbb::concurrent_bounded_queue<ResamplTrace*> result_queue;

};

#endif
