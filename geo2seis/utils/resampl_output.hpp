#ifndef RESAMPL_OUTPUT_HPP
#define RESAMPL_OUTPUT_HPP

#include <seismic_parameters.hpp>
#include "nrlib/segy/segygeometry.hpp"
#include <seismic_geometry.hpp>
#include <nrlib/stormgrid/stormcontgrid.hpp>

class ResamplOutput {
public:
  ResamplOutput(SeismicParameters &seismic_parameters, bool time, size_t n_samples);


  void AddResampleCase(std::string            filename,
                        NRLib::StormContGrid &input_grid,
                        bool                  time,
                        std::vector<double>  &time_or_depth_vec_reg,
                        SeismicParameters    &seismic_parameters);

  void AddTrace(SeismicParameters                   &seismic_parameters,
                std::vector<double>                 &time_or_depth_vec_reg,
                std::vector<NRLib::Grid2D<double> > &traces,
                double                               x,
                double                               y);

  std::vector<NRLib::Grid2D<double> > &GetTraces()     { return traces_;      }
  std::vector<NRLib::StormContGrid*>  GetInputGrid()  { return input_grid_;  }

private:
  bool                      segy_ok_;
  std::vector<bool>         segy_files_ok_;
  std::vector<NRLib::SegY*> segy_files_;
  //Til Paal: Hvis du finner en bedre løsning på dette er det fint:
  NRLib::SegY segy_1_, segy_2_, segy_3_, segy_4_, segy_5_, segy_6_, segy_7_, segy_8_, segy_9_, segy_10_;
  std::vector<NRLib::Grid2D<double> > traces_;
  double                              n_samples_;
  bool                                time_;
  std::vector<NRLib::StormContGrid*>   input_grid_;

};

#endif
