#ifndef RESAMPL_TRACE_HPP
#define RESAMPL_TRACE_HPP

#include "nrlib/grid/grid2d.hpp"
#include "seismic_output.hpp"

#include "utils/trace.hpp"

#include <vector>

class ResamplTrace
{
public:

  ResamplTrace(std::vector<NRLib::Grid2D<double>> & traces);

  const std::vector<NRLib::Grid2D<double>> & GetTraces() const { return traces_     ;}
  size_t GetI()                                          const { return i_          ;}
  size_t GetJ()                                          const { return j_          ;}
  bool   GetIsEmpty()                                    const { return empty_      ;}
  size_t GetJobNumber()                                  const { return job_number_ ;}

  void SetIsEmpty(bool empty)                                  { empty_ = empty     ;}
  void SetJobID(Trace * trace);

private:

  const std::vector<NRLib::Grid2D<double>> & traces_;

  size_t                                     i_;
  size_t                                     j_;
  size_t                                     job_number_;
  bool                                       empty_;

};

#endif
