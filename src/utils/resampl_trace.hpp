#ifndef RESAMPL_TRACE_HPP
#define RESAMPL_TRACE_HPP

#include <utils/trace.hpp>
#include <vector>
#include "nrlib/grid/grid2d.hpp"
#include "seismic_output.hpp"

class ResamplTrace {

public:

  ResamplTrace(std::vector<NRLib::Grid2D<double> > &traces);

  void SetJobID(Trace *trace);
  void SetIsEmpty(bool empty) { empty_ = empty; };
  std::vector<NRLib::Grid2D<double> > &GetTraces() { return traces_; }
  size_t GetI() { return i_; };
  size_t GetJ() { return j_; };
  double GetX() { return x_; };
  double GetY() { return y_; };
  bool   GetIsEmpty() { return empty_; };
  size_t GetJobNumber() { return job_number_; };

private:

  std::vector<NRLib::Grid2D<double> > traces_;

  double                x_;
  double                y_;
  size_t                i_;
  size_t                j_;
  size_t                job_number_;
  bool                  empty_;

};

#endif

