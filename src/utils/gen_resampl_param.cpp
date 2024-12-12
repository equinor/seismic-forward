#include "gen_resampl_param.hpp"
#include <utils/trace.hpp>
#include <utils/resampl_trace.hpp>
#include "tbb/concurrent_queue.h"


GenResamplParam::GenResamplParam(SeismicParameters                            & seismic_parameters,
                                 // const std::vector<double>                    & time_or_depth_vec_reg,
                                 //NRLib::StormContGrid                         & time_or_depth_grid,
                                 const NRLib::RegularSurface<double>          & toptime,
                                 size_t                                         n_samples,
                                 size_t                                         n_traces,
                                 bool                                           time,
                                 tbb::concurrent_queue<ResamplTrace*>         & empty_queue,
                                 tbb::concurrent_bounded_queue<ResamplTrace*> & result_queue,
                                 tbb::concurrent_queue<Trace*>                & traces)
  : seismic_parameters(seismic_parameters),
//    time_or_depth_vec_reg(time_or_depth_vec_reg),
//    time_or_depth_grid(&time_or_depth_grid),
    toptime(toptime),
    n_samples(n_samples),
    n_traces(n_traces),
    time(time),
    empty_queue(empty_queue),
    result_queue(result_queue),
    traces(traces)
{
}
