#include "gen_resampl_param.hpp"
#include <utils/trace.hpp>
#include <utils/resampl_trace.hpp>
#include "tbb/concurrent_queue.h"


GenResamplParam::GenResamplParam(tbb::concurrent_queue<ResamplTrace*>         & empty_queue,
                                 tbb::concurrent_bounded_queue<ResamplTrace*> & result_queue,
                                 tbb::concurrent_queue<Trace*>                & traces)
  :
    empty_queue(empty_queue),
    result_queue(result_queue),
    traces(traces)
{
}
