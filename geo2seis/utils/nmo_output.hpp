#ifndef NMO_OUTPUT_HPP
#define NMO_OUTPUT_HPP

#include <seismic_parameters.hpp>

class NMOOutput {
  public:
    NMOOutput(SeismicParameters &seismic_parameters,
            std::vector<double> twt_0,
            std::vector<double> z_0);

    void AddTrace(SeismicParameters     &seismic_parameters,
                  NRLib::Grid2D<double> &timegrid_pos,
                  NRLib::Grid2D<double> &nmo_timegrid_pos,
                  NRLib::Grid2D<double> &nmo_timegrid_stack_pos,
                  NRLib::Grid2D<double> &nmo_depthgrid_pos,
                  NRLib::Grid2D<double> &nmo_depthgrid_stack_pos,
                  NRLib::Grid2D<double> &nmo_timeshiftgrid_pos,
                  NRLib::Grid2D<double> &nmo_timeshiftgrid_stack_pos,
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

    bool GetNMODepthSegyOk(void)          { return nmo_depth_segy_ok_;};
    bool GetNMODepthStackSegyOk(void)     { return nmo_depth_stack_segy_ok_;};
    bool GetNMOTimeshiftSegyOk(void)      { return nmo_timeshift_segy_ok_;};
    bool GetNMOTimeshiftStackSegyOk(void) { return nmo_timeshift_stack_segy_ok_;};



  private:
      bool segy_ok_;                    
      bool nmo_time_segy_ok_;           
      bool prenmo_time_segy_ok_;        
      bool nmo_time_stack_segy_ok_;     
      bool nmo_depth_segy_ok_;          
      bool nmo_depth_stack_segy_ok_;    
      bool nmo_timeshift_segy_ok_;      
      bool nmo_timeshift_stack_segy_ok_;
      bool twtx_segy_ok_;               
      NRLib::SegY nmo_time_segy_, prenmo_time_segy_, nmo_time_stack_segy_, nmo_depth_segy_, nmo_depth_stack_segy_;
      NRLib::SegY nmo_timeshift_segy_, nmo_timeshift_stack_segy_, twtx_segy_;

      NRLib::StormContGrid *timegrid_;
      NRLib::StormContGrid *timeshiftgrid_;
      NRLib::StormContGrid *depthgrid_;
      std::vector<double> twt_0_;
      std::vector<double> z_0_;

      //NRLib::Grid2D<double> nmo_timegrid_stack_pos_;
      //NRLib::Grid2D<double> nmo_timeshiftgrid_pos_;
      //NRLib::Grid2D<double> nmo_timeshiftgrid_stack_pos_;
      //NRLib::Grid2D<double> nmo_depthgrid_pos_;
      //NRLib::Grid2D<double> nmo_depthgrid_stack_pos_;


};

#endif
