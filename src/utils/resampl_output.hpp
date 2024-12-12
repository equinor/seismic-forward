#ifndef RESAMPL_OUTPUT_HPP
#define RESAMPL_OUTPUT_HPP

#include "nrlib/stormgrid/stormcontgrid.hpp"
#include "nrlib/segy/segygeometry.hpp"

#include "seismic_parameters.hpp"
#include "seismic_geometry.hpp"

class ResamplOutput
{
public:
  ResamplOutput(const bool   segy_ok,
                const bool   time,
                const size_t n_samples);

  ~ResamplOutput(void);

  void AddResampleCase(std::string            filename,
                       NRLib::StormContGrid & input_grid,
                       bool                   time,
                       std::vector<double>  & time_or_depth_vec_reg,
                       SeismicParameters    & seismic_parameters);

  void AddTrace(SeismicOutput                            * seismic_output,
                std::vector<double>                      & time_or_depth_vec_reg,
                const std::vector<NRLib::Grid2D<double>> & traces,
                const double                               x,
                const double                               y,
                const bool                                 time);

  std::vector<NRLib::Grid2D<double>> & GetTraces()     { return traces_     ;}
  std::vector<NRLib::StormContGrid*>   GetInputGrid()  { return input_grid_ ;}

private:
  bool                                segy_ok_;
  double                              n_samples_;
  bool                                time_;

  std::vector<bool>                   segy_files_ok_;
  std::vector<NRLib::SegY*>           segy_files_;

  std::vector<NRLib::Grid2D<double>>  traces_;
  std::vector<NRLib::StormContGrid*>  input_grid_;

};

#endif
