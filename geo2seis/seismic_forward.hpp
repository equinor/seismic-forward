#ifndef SEISMIC_FORWARD_HPP
#define SEISMIC_FORWARD_HPP

#include <stdio.h>
#include <string>
#include <vector>
#include "modelsettings.hpp"
#include <seismic_parameters.hpp>

#include <nrlib/stormgrid/stormcontgrid.hpp>
#include <nrlib/surface/regularsurface.hpp>


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
            bool timeshift_output);

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

    static double findTFromZ(double z, std::vector<double> &zvec, std::vector<double> &tvec);
};

#endif
