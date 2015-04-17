#ifndef SEISMIC_GEOMETRY_HPP
#define SEISMIC_GEOMETRY_HPP

#include <stddef.h>
#include <nrlib/volume/volume.hpp>

class SeismicGeometry {
  public:
    SeismicGeometry();

    NRLib::Volume createDepthVolume();
    NRLib::Volume createTimeVolume();

    double dx();
    double dy();
    double dz();
    double dt();

    double x0();
    double y0();
    double z0();
    double t0();

    double xlength();
    double ylength();
    double zlength();
    double tlength();

    double angle();

    double xmin();
    double ymin();

    double xsurfacelength();
    double ysurfacelength();

    size_t nx();
    size_t ny();
    size_t nz();
    size_t nt();

    size_t nxsurfaceeclipse();
    size_t nysurfaceeclipse();

    size_t zreflectorcount();

    void setDxDy(double dx, double dy);
    void setDz(double dz);
    void setDt(double dt);
    void setGeometry(double x0, double y0, double lx, double ly, double angle);
    void setZRange(double z_min, double z_max);
    void setZReflectorCount(size_t nzrefl);

    void setTRange(double t_min, double t_max);
    void setNt(size_t nt);

  private:

    void printValues();

    double _dx;
    double _dy;
    double _dz;
    double _dt;

    double _x0;
    double _y0;
    double _z0;
    double _t0;

    double _z_max;
    double _t_max;

    double _lx;
    double _ly;

    double _angle;

    double _lx_surf;
    double _ly_surf;

    size_t _z_reflector_count;
    size_t _nt;
};

#endif

