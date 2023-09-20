#ifndef SEISMIC_OUTPUT_HPP
#define SEISMIC_OUTPUT_HPP

#include "nrlib/stormgrid/stormcontgrid.hpp"
#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/segy/segygeometry.hpp"
#include "nrlib/segy/traceheader.hpp"

#include "seismic_parameters.hpp"
#include "seismic_geometry.hpp"
#include "modelsettings.hpp"

#include <stdio.h>
#include <string>
#include <vector>

class SeismicParameters;

class SeismicOutput
{
public:
  SeismicOutput(ModelSettings *model_settings);

  void SetSegyGeometry(SeismicParameters   & seismic_parameters,
                       const NRLib::Volume & vol,
                       size_t                nx,
                       size_t                ny);

  bool CheckUTMPrecision(SeismicParameters   & seismic_parameters,
                         const NRLib::Volume & vol,
                         size_t                nx,
                         size_t                ny);

  bool PrepareSegy(NRLib::SegY               & segyout,
                   const std::vector<double> & twt_0,
                   size_t                      n_samples,
                   std::string                 fileName,
                   SeismicParameters         & seismic_parameters,
                   const std::vector<double> & offset_vec,
                   size_t                      n_traces_per_ensamble,
                   bool                        time,
                   bool                        nmo);

  void WriteSegyGather(const NRLib::Grid2D<double> & data_gather,
                       NRLib::SegY                 & segyout,
                       const std::vector<double>   & twt_0,
                       const std::vector<double>   & offset_vec,
                       bool                          time,
                       double                        x,
                       double                        y,
                       bool                          nmo);

  void WriteZeroSegyGather(NRLib::SegY               & segyout,
                           const std::vector<double>   offset_vec,
                           double                      x,
                           double                      y,
                           bool                        nmo);

  void WriteDepthSurfaces(const NRLib::RegularSurface<double> & top_eclipse,
                          const NRLib::RegularSurface<double> & bottom_eclipse);

  void WriteTimeSurfaces(SeismicParameters &seismic_parameters);

  void WriteVpVsRho(SeismicParameters &seismic_parameters);
  void WriteZValues(SeismicParameters &seismic_parameters);
  void WriteTwt(SeismicParameters &seismic_parameters);
  void WriteVrms(SeismicParameters &seismic_parameters, std::string name_pp_or_ps = "");

  void WriteSeismicTimeStorm     (NRLib::StormContGrid & timegrid     , double offset, bool is_stack = false);
  void WriteSeismicDepthStorm    (NRLib::StormContGrid & depthgrid    , double offset, bool is_stack = false);
  void WriteSeismicTimeshiftStorm(NRLib::StormContGrid & timeshiftgrid, double offset, bool is_stack = false);
  void WriteReflections(std::vector<NRLib::StormContGrid> & rgridvec, double angle_or_offset);

  void PrintVector(std::vector<double> vec, std::string filename);
  void PrintVectorSizeT(std::vector<size_t> vec, std::string filename);
  void PrintMatrix(NRLib::Grid2D<double> matrix, std::string filename);

private:

  double                   top_time_window_;
  double                   bot_time_window_;
  bool                     time_window_;
  bool                     depth_window_;
  double                   top_depth_window_;
  double                   bot_depth_window_;

  std::string              prefix_;
  std::string              suffix_;

  std::vector<std::string> extra_parameter_names_;

  int                      inline_start_;
  int                      xline_start_;
  std::string              inline_direction_;
  short                    scalco_;
  bool                     xline_x_axis_;
  int                      inline_step_;
  int                      xline_step_;

  NRLib::TraceHeaderFormat thf_;
};

#endif
