#include <seismic_geometry.hpp>
#include <nrlib/math/constants.hpp>
#include <math.h>


SeismicGeometry::SeismicGeometry() {
  x0_ = 0.0;
  y0_ = 0.0;
  z0_ = 0.0;
  t0_ = 0.0;

  lx_ = 0.0;
  ly_ = 0.0;

  z_max_ = 0.0;
  t_max_ = 0.0;

  angle_ = 0.0;
  lx_surf_ = 0.0;
  ly_surf_ = 0.0;

  dx_ = 0.0;
  dy_ = 0.0;
  dz_ = 0.0;
  dt_ = 0.0;

  z_reflector_count_ = 0;
  nt_ = 0;
}

double SeismicGeometry::x0() {
  return x0_;
}

double SeismicGeometry::y0() {
  return y0_;
}

double SeismicGeometry::xlength() {
  return lx_;
}

double SeismicGeometry::ylength() {
  return ly_;
}

double SeismicGeometry::angle() {
  return angle_;
}

double SeismicGeometry::dx() {
  return dx_;
}

double SeismicGeometry::dy() {
  return dy_;
}

double SeismicGeometry::dz() {
  return dz_;
}

double SeismicGeometry::dt() {
  return dt_;
}

double SeismicGeometry::xmin() {
  double xmin;
  if (angle_ < -NRLib::PiHalf) {
    xmin = x0_ - lx_ * fabs(cos(angle_));
  } else if (angle_ < 0) {
    xmin = x0_;
  } else if (angle_ < NRLib::PiHalf) {
    xmin = x0_ - ly_ * cos(NRLib::PiHalf - angle_);
  } else {
    xmin = x0_ - ly_ * cos(angle_ - NRLib::PiHalf) - lx_ * cos(NRLib::Pi - angle_);
  }
  return xmin;
}

double SeismicGeometry::ymin() {
  double ymin;
  if (angle_ < -NRLib::PiHalf) {
    ymin = y0_ - ly_surf_;
  } else if (angle_ < 0) {
    ymin = y0_ - lx_ * sin(-angle_);
  } else if (angle_ < NRLib::PiHalf) {
    ymin = y0_;
  } else {
    ymin = y0_ - ly_ * sin(angle_ - NRLib::PiHalf);
  }
  return ymin;
}

double SeismicGeometry::xsurfacelength() {
  return lx_surf_;
}

double SeismicGeometry::ysurfacelength() {
  return ly_surf_;
}

size_t SeismicGeometry::nx() {
  return static_cast<size_t>(floor(xlength() / dx()));
}

size_t SeismicGeometry::ny() {
  return static_cast<size_t>(floor(ylength() / dy()));
}

size_t SeismicGeometry::nxsurfaceeclipse() {
  return static_cast<size_t>(xsurfacelength() / dx());
}

size_t SeismicGeometry::nysurfaceeclipse() {
  return static_cast<size_t>(ysurfacelength() / dy());
}

double SeismicGeometry::zlength() {
  return nz() * dz();
}

double SeismicGeometry::z0() {
  return z0_;
}

size_t SeismicGeometry::nz() {
  return static_cast<size_t>(floor((z_max_ - z0_) / dz_));
}

size_t SeismicGeometry::zreflectorcount() {
  return z_reflector_count_;
}

size_t SeismicGeometry::nt() {
  return nt_;
}

void SeismicGeometry::setGeometry(double x0, double y0, double lx, double ly, double angle) {
  x0_ = x0;
  y0_ = y0;
  lx_ = lx;
  ly_ = ly;

  if (angle > NRLib::Pi) {
    angle = angle - 2 * NRLib::Pi;
  }
  if (angle < -NRLib::Pi) {
    angle = angle + 2 * NRLib::Pi;
  }

  angle_ = angle;

  lx_surf_ = lx * fabs(cos(angle)) + ly * sin(fabs(angle));
  ly_surf_ = lx * sin(fabs(angle)) + ly * fabs(cos(angle));
}

void SeismicGeometry::setDxDy(double dx, double dy) {
  dx_ = dx;
  dy_ = dy;
}

void SeismicGeometry::setDz(double dz) {
  dz_ = dz;
}

void SeismicGeometry::setDt(double dt) {
  dt_ = dt;
}

void SeismicGeometry::setNt(size_t nt) {
  nt_ = nt;
}

void SeismicGeometry::setZReflectorCount(size_t nzrefl) {
  z_reflector_count_ = nzrefl;
}

void SeismicGeometry::setZRange(double z_min, double z_max) {
  z0_ = z_min;
  z_max_ = z_max;
}

double SeismicGeometry::t0() {
  return t0_;
}

double SeismicGeometry::tmax() {
  return t_max_;
}

void SeismicGeometry::setTRange(double t_min, double t_max) {
  t0_ = t_min;
  t_max_ = t_max;
}

double SeismicGeometry::tlength() {
  return nt() * dt();
}

NRLib::Volume SeismicGeometry::createDepthVolume() {
  return NRLib::Volume(x0_, y0_, z0_, lx_, ly_, zlength(), angle_);
}

NRLib::Volume SeismicGeometry::createTimeVolume() {
  return NRLib::Volume(x0_, y0_, t0_, lx_, ly_, tlength(), angle_);
}

void SeismicGeometry::printValues() {
  printf("x0: %f lx: %f nx: %lu dx: %f\n", x0(), xlength(), nx(), dx());
  printf("y0: %f ly: %f ny: %lu dy: %f\n", y0(), ylength(), ny(), dy());
  printf("z0: %f lz: %f nz: %lu dz: %f\n", z0(), zlength(), nz(), dz());
  printf("t0: %f lt: %f nt: %lu dt: %f\n", t0(), tlength(), nt(), dt());
  printf("zrefl: %lu lx_surf: %f ly_surf: %f angle: %f\n", zreflectorcount(), xsurfacelength(), ysurfacelength(), angle());
  printf("nxse: %lu nyse: %lu\n", nxsurfaceeclipse(), nysurfaceeclipse());
}
