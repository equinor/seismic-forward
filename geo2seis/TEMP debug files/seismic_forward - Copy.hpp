#ifndef SEISMIC_FORWARD_HPP
#define SEISMIC_FORWARD_HPP

#include <stdio.h>
#include <string>
#include <vector>
#include "modelsettings.hpp"
#include <seismic_parameters.hpp>

#include <nrlib/stormgrid/stormcontgrid.hpp>
#include <nrlib/surface/regularsurface.hpp>
//#include "nrlib/geometry/interpolation.hpp"


class SeismicForward {
  public:
    static void seismicForward(SeismicParameters &seismic_parameters);

  private:
    static void generateSeismic(std::vector<NRLib::StormContGrid> &rgridvec,
                                NRLib::StormContGrid &twtgrid,
                                NRLib::StormContGrid &zgrid,
                                NRLib::StormContGrid &twt_timeshift,
                                std::vector<NRLib::StormContGrid> &timegridvec,
                                std::vector<NRLib::StormContGrid> &depthgridvec,
                                std::vector<NRLib::StormContGrid> &timeshiftgridvec,
                                Wavelet *wavelet,
                                double dt,
                                NRLib::RegularSurface<double> &bot,
                                NRLib::RegularSurface<double> &toptime,
                                double t0, double dz, double z0,
                                std::vector<double> &constvp,
                                double waveletScale,
                                bool time_output,
                                bool depth_output,
                                bool timeshift_output,
                                SeismicParameters &seismic_parameters);


    static void generateSeismicPos(std::vector<std::vector<double> > &timegrid_pos,
                                   std::vector<std::vector<double> > refl_pos,
                                   std::vector<std::vector<double> > twtx_pos,
                                   NRLib::StormContGrid              &zgrid,
                                   NRLib::RegularSurface<double>     &toptime,
                                   Wavelet                           *wavelet,
                                   double                            waveletScale,
                                   std::vector<double>               offset,
                                   double                            t0,
                                   double                            dt,
                                   size_t                            i,
                                   size_t                            j, 
                                   std::vector<size_t>               n_min,
                                   std::vector<size_t>               n_max,
                                     SeismicParameters &seismic_parameters);


    static void generateSeismicOnFile(std::vector<NRLib::StormContGrid> &rgridvec,
                                      NRLib::StormContGrid &twtgrid,
                                      NRLib::StormContGrid &zgrid,
                                      NRLib::StormContGrid &twt_timeshift,
                                      Wavelet *wavelet,
                                      double dt,
                                      int nt, int nz, int nx, int ny,
                                      NRLib::RegularSurface<double> &bot,
                                      NRLib::RegularSurface<double> &toptime,
                                      double t0, double dz, double z0,
                                      std::vector<double> &constvp,
                                      double waveletScale,
                                      bool time_output,
                                      bool depth_output,
                                      bool timeshift_output);

    static double findTFromZ(double z, 
                            std::vector<double> &zvec, 
                            std::vector<double> &tvec);

    static void   convertSeis(std::vector<double>               twt_vec,
                              std::vector<double>               twt_0, 
                              std::vector<double>               zgrid, 
                              std::vector<double>               z_0, 
                              std::vector<std::vector<double> > seismic,
                              std::vector<std::vector<double> > &conv_seismic);

    static void   nmoCorrInterpol1Pos(std::vector<double>                t_in, 
                                      std::vector<std::vector<double> >  data_in, 
                                      std::vector<std::vector<double> >  t_out, 
                                      std::vector<std::vector<double> > &data_out);


    static void   findReflectionsPos(SeismicParameters                &seismic_parameters,  
                                    std::vector<std::vector<double> > &r_vec, 
                                    std::vector<std::vector<double> >  theta_vec, 
                                    std::vector<double>                offset, 
                                    size_t                             i, 
                                    size_t                             j);



    static void   findThetaPos(std::vector<std::vector<double> > &thetagrid, 
                               std::vector<double>                twt_vec, 
                               std::vector<double>                vrms_vec, 
                               std::vector<double>                offset);


    static void   findTWTxPos(std::vector<std::vector<double> > &twtx_grid,  
                              std::vector<double>                twt_vec, 
                              std::vector<double>                vrms_vec, 
                              std::vector<double>                offset);


    static void   findLoopIndeces(SeismicParameters &seismic_parameters,
                                  int               &n_xl,
                                  int               &il_min,
                                  int               &il_max,
                                  int               &il_step,
                                  int               &xl_min,
                                  int               &xl_max,
                                  int               &xl_step,
                                  bool              &segy);
    
};

#endif
