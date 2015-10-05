#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include <seismic_parameters.hpp>

class Output {
  public:
    Output(SeismicParameters &seismic_parameters,
              std::vector<double> twt_0,
              std::vector<double> z_0,
              std::vector<double> twts_0,
              std::vector<double> offset_vec,
              size_t time_samples_stretch);

    void AddTrace(SeismicParameters     &seismic_parameters,
                  NRLib::Grid2D<double> &timegrid_pos,
                  NRLib::Grid2D<double> &pre_nmo_timegrid_pos,
                  NRLib::Grid2D<double> &timegrid_stack_pos,
                  NRLib::Grid2D<double> &depthgrid_pos,
                  NRLib::Grid2D<double> &depthgrid_stack_pos,
                  NRLib::Grid2D<double> &timeshiftgrid_pos,
                  NRLib::Grid2D<double> &timeshiftgrid_stack_pos,
                  NRLib::Grid2D<double> &twtx_reg,
                  double                 x, 
                  double                 y,
                  size_t                 i,
                  size_t                 j);
    void AddZeroTrace(SeismicParameters     &seismic_parameters,
                     double                 x, 
                     double                 y,
                     size_t                 i,
                     size_t                 j);
    void WriteSeismicStorm(SeismicParameters     &seismic_parameters);

    const bool GetDepthSegyOk(void)          { return depth_segy_ok_;};
    const bool GetDepthStackSegyOk(void)     { return depth_stack_segy_ok_;};
    const bool GetTimeshiftSegyOk(void)      { return timeshift_segy_ok_;};
    const bool GetTimeshiftStackSegyOk(void) { return timeshift_stack_segy_ok_;};



  private:
      bool segy_ok_;                    
      bool time_segy_ok_;           
      bool prenmo_time_segy_ok_;        
      bool time_stack_segy_ok_;     
      bool depth_segy_ok_;          
      bool depth_stack_segy_ok_;    
      bool timeshift_segy_ok_;      
      bool timeshift_stack_segy_ok_;
      bool twtx_segy_ok_;               
      NRLib::SegY time_segy_, prenmo_time_segy_, time_stack_segy_, depth_segy_, depth_stack_segy_;
      NRLib::SegY timeshift_segy_, timeshift_stack_segy_, twtx_segy_;

      NRLib::StormContGrid *timegrid_;
      NRLib::StormContGrid *timeshiftgrid_;
      NRLib::StormContGrid *depthgrid_;
      std::vector<double> twt_0_;
      std::vector<double> z_0_;
      std::vector<double> twts_0_;
      std::vector<double> offset_vec_;




};

#endif
