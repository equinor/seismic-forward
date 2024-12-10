#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include "seismic_parameters.hpp"

class ResultTrace;

class Output
{

public:
  Output(SeismicParameters   & seismic_parameters,
         ModelSettings       * model_settings,
         std::vector<double>   twt_0,
         std::vector<double>   z_0,
         std::vector<double>   twts_0,
         std::vector<double>   offset_vec,
         size_t                time_samples_stretch);

  ~Output(void);

  void AddTrace(const ResultTrace   & result_trace,
                const ModelSettings & model_settings,
                SeismicOutput       * seismic_output);

  void AddZeroTraceToStormGrid(NRLib::StormContGrid & grid,
                               const size_t           i,
                               const size_t           j);

  void AddTraceToStormGrid(NRLib::StormContGrid        & grid,
                           const NRLib::Grid2D<double> & trace,
                           const size_t                  i,
                           const size_t                  j);

  void WriteStatisticsForSeismic(ModelSettings * model_settings);

  void WriteSeismicStorm(ModelSettings                     * model_settings,
                         SeismicOutput                     * seismic_output,
                         std::vector<NRLib::StormContGrid> & rgrids);

  bool GetDepthSegyOk(void)          const { return depth_segy_ok_          ;}
  bool GetDepthStackSegyOk(void)     const { return depth_stack_segy_ok_    ;}
  bool GetTimeshiftSegyOk(void)      const { return timeshift_segy_ok_      ;}
  bool GetTimeshiftStackSegyOk(void) const { return timeshift_stack_segy_ok_;}

private:

  bool                   segy_ok_;
  bool                   time_segy_ok_;
  bool                   prenmo_time_segy_ok_;
  bool                   time_stack_segy_ok_;
  bool                   depth_segy_ok_;
  bool                   depth_stack_segy_ok_;
  bool                   timeshift_segy_ok_;
  bool                   timeshift_stack_segy_ok_;
  bool                   twtx_segy_ok_;

  NRLib::SegY            time_segy_;
  NRLib::SegY            prenmo_time_segy_;
  NRLib::SegY            time_stack_segy_;
  NRLib::SegY            depth_segy_;
  NRLib::SegY            depth_stack_segy_;
  NRLib::SegY            timeshift_segy_;
  NRLib::SegY            timeshift_stack_segy_;
  NRLib::SegY            twtx_segy_;

  NRLib::StormContGrid * timegrid_;
  NRLib::StormContGrid * timeshiftgrid_;
  NRLib::StormContGrid * depthgrid_;
  std::vector<double>    twt_0_;
  std::vector<double>    z_0_;
  std::vector<double>    twts_0_;
  std::vector<double>    offset_vec_;
};

#endif
