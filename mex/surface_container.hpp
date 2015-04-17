#ifndef __SURFACE_CONTAINER_HPP__
#define __SURFACE_CONTAINER_HPP__

#include <string>
#include <vector>
#include <mex.h>
#include "modelsettings.hpp"
#include <seismic_parameters.hpp>
#include <nrlib/stormgrid/stormcontgrid.hpp>


class SurfaceContainer {

  public:
    SurfaceContainer(std::string path);

    ~SurfaceContainer() {
        mxDestroyArray(_zgrid);
        mxDestroyArray(_twtgrid);
        mxDestroyArray(_vpgrid);
        mxDestroyArray(_vsgrid);
        mxDestroyArray(_rhogrid);

        mxDestroyArray(_toptime);
        mxDestroyArray(_bottime);
        mxDestroyArray(_topeclipse);
        mxDestroyArray(_boteclipse);

        delete _seismic_parameters;
        delete _model_settings;
    };

    mxArray * vpgrid();
    mxArray * vsgrid();
    mxArray * rhogrid();
    mxArray * zgrid();
    mxArray * twtgrid();

    mxArray * toptime();
    mxArray * bottime();
    mxArray * topeclipse();
    mxArray * boteclipse();

  private:
    mxArray * createArrayFromGrid(NRLib::StormContGrid & grid);

    mxArray * createArrayFromGrid(NRLib::RegularSurface<double> & grid);

    void copyValuesFromGridToArray(NRLib::StormContGrid & grid, mxArray * array);

    void copyValuesFromGridToArray(NRLib::RegularSurface<double> & grid, mxArray * array);


    ModelSettings * _model_settings;
    SeismicParameters * _seismic_parameters;


    mxArray * _vpgrid;
    mxArray * _vsgrid;
    mxArray * _rhogrid;
    mxArray * _zgrid;
    mxArray * _twtgrid;
    mxArray * _toptime;
    mxArray * _bottime;
    mxArray * _topeclipse;
    mxArray * _boteclipse;
};

#endif