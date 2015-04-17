#ifndef SEISMIC_PARAMETERS_HPP
#define SEISMIC_PARAMETERS_HPP

#include <stdio.h>
#include <string>
#include <vector>
#include <nrlib/surface/regularsurface.hpp>
#include <nrlib/volume/volume.hpp>
#include "modelsettings.hpp"
#include "seismic_output.hpp"

class Wavelet;
class ModelSettings;
class SeismicGeometry;
class SeismicOutput;

namespace NRLib {
    class EclipseGrid;
    class EclipseGeometry;
    class SegyGeometry;
    class StormContGrid;
}


class SeismicParameters {
  public:
    SeismicParameters(ModelSettings *model_settings);

    ~SeismicParameters() {
    };

    NRLib::StormContGrid &zGrid()  { return *zgrid; };
    NRLib::StormContGrid &vpGrid()  { return *vpgrid; };
    NRLib::StormContGrid &vsGrid()  { return *vsgrid; };
    NRLib::StormContGrid &rhoGrid()  { return *rhogrid; };
    NRLib::StormContGrid &twtGrid()  { return *twtgrid; };
    std::vector<NRLib::StormContGrid> &rGrids()  { return *rgridvec; };
    std::vector<NRLib::StormContGrid> &extraParametersGrids()  { return *extra_parameter_grid; };
    NRLib::EclipseGrid &eclipseGrid()  { return *eclipse_grid; };

    size_t topK() { return top_k; }
    size_t bottomK() { return bottom_k; }
    double theta0() { return theta_0; }
    double dTheta() { return dtheta; }
    size_t nTheta() { return ntheta; }

    NRLib::RegularSurface<double> &topTime() { return top_time; };
    NRLib::RegularSurface<double> &bottomTime() { return bot_time; };

    NRLib::RegularSurface<double> &topEclipse() { return topeclipse; };
    NRLib::RegularSurface<double> &bottomEclipse() { return boteclipse; };


    ModelSettings *model_settings;

    SeismicOutput* seismicOutput() { return seismic_output; };
    SeismicGeometry* seismicGeometry() { return seismic_geometry; };

    NRLib::SegyGeometry *segyGeometry() { return segy_geometry; };
    Wavelet* wavelet() { return _wavelet; };
    double waveletScale() { return  _wavelet_scale; };

  private:
    SeismicGeometry *seismic_geometry;
    SeismicOutput *seismic_output;

    size_t ntheta;
    double theta_0;
    double dtheta;
    double theta_max;

    Wavelet *_wavelet;
    double _wavelet_scale;

    NRLib::EclipseGrid *eclipse_grid;

    size_t top_k;
    size_t bottom_k;

    NRLib::RegularSurface<double> top_time;
    NRLib::RegularSurface<double> bot_time;

    NRLib::RegularSurface<double> topeclipse;
    NRLib::RegularSurface<double> boteclipse;

    NRLib::SegyGeometry *segy_geometry;

    NRLib::StormContGrid *zgrid;
    NRLib::StormContGrid *vpgrid;
    NRLib::StormContGrid *vsgrid;
    NRLib::StormContGrid *rhogrid;
    NRLib::StormContGrid *twtgrid;
    std::vector<NRLib::StormContGrid> *rgridvec;
    std::vector<NRLib::StormContGrid> *extra_parameter_grid;

    void setupWavelet();

    void readEclipseGrid();

    void findGeometry();

    void findSurfaceGeometry();

    void calculateAngleSpan();

    void createGrids();
};

#endif
