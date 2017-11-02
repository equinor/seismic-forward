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
    double tmax();

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

    double dx_;
    double dy_;
    double dz_;
    double dt_;

    double x0_;
    double y0_;
    double z0_;
    double t0_;

    double z_max_;
    double t_max_;

    double lx_;
    double ly_;

    double angle_;

    double lx_surf_;
    double ly_surf_;

    size_t z_reflector_count_;
    size_t nt_;
};

#endif

