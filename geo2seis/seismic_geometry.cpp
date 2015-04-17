#include <seismic_geometry.hpp>
#include <nrlib/math/constants.hpp>
#include <math.h>


SeismicGeometry::SeismicGeometry() {
    _x0 = 0.0;
    _y0 = 0.0;
    _z0 = 0.0;
    _t0 = 0.0;

    _lx = 0.0;
    _ly = 0.0;

    _z_max = 0.0;
    _t_max = 0.0;

    _angle = 0.0;
    _lx_surf = 0.0;
    _ly_surf = 0.0;

    _dx = 0.0;
    _dy = 0.0;
    _dz = 0.0;
    _dt = 0.0;

    _z_reflector_count = 0;
    _nt = 0;
}

double SeismicGeometry::x0() {
    return _x0;
}

double SeismicGeometry::y0() {
    return _y0;
}

double SeismicGeometry::xlength() {
    return _lx;
}

double SeismicGeometry::ylength() {
    return _ly;
}

double SeismicGeometry::angle() {
    return _angle;
}

double SeismicGeometry::dx() {
    return _dx;
}

double SeismicGeometry::dy() {
    return _dy;
}

double SeismicGeometry::dz() {
    return _dz;
}

double SeismicGeometry::dt() {
    return _dt;
}

double SeismicGeometry::xmin() {
    double xmin;
    if (_angle < -NRLib::PiHalf) {
        xmin = _x0 - _lx * fabs(cos(_angle));
    } else if (_angle < 0) {
        xmin = _x0;
    } else if (_angle < NRLib::PiHalf) {
        xmin = _x0 - _ly * cos(NRLib::PiHalf - _angle);
    } else {
        xmin = _x0 - _ly * cos(_angle - NRLib::PiHalf) - _lx * cos(NRLib::Pi - _angle);
    }
    return xmin;
}

double SeismicGeometry::ymin() {
    double ymin;
    if (_angle < -NRLib::PiHalf) {
        ymin = _y0 - _ly_surf;
    } else if (_angle < 0) {
        ymin = _y0 - _lx * sin(-_angle);
    } else if (_angle < NRLib::PiHalf) {
        ymin = _y0;
    } else {
        ymin = _y0 - _ly * sin(_angle - NRLib::PiHalf);
    }
    return ymin;
}

double SeismicGeometry::xsurfacelength() {
    return _lx_surf;
}

double SeismicGeometry::ysurfacelength() {
    return _ly_surf;
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
    return _z0;
}

size_t SeismicGeometry::nz() {
    return static_cast<size_t>(floor((_z_max - _z0) / _dz));
}

size_t SeismicGeometry::zreflectorcount() {
    return _z_reflector_count;
}

size_t SeismicGeometry::nt() {
    return _nt;
}

void SeismicGeometry::setGeometry(double x0, double y0, double lx, double ly, double angle) {
    _x0 = x0;
    _y0 = y0;
    _lx = lx;
    _ly = ly;

    if (angle > NRLib::Pi) {
        angle = angle - 2 * NRLib::Pi;
    }
    if (angle < -NRLib::Pi) {
        angle = angle + 2 * NRLib::Pi;
    }

    _angle = angle;

    _lx_surf = lx * fabs(cos(angle)) + ly * sin(fabs(angle));
    _ly_surf = lx * sin(fabs(angle)) + ly * fabs(cos(angle));
}

void SeismicGeometry::setDxDy(double dx, double dy) {
    _dx = dx;
    _dy = dy;
}

void SeismicGeometry::setDz(double dz) {
    _dz = dz;
}

void SeismicGeometry::setDt(double dt) {
    _dt = dt;
}

void SeismicGeometry::setNt(size_t nt) {
    _nt = nt;
}

void SeismicGeometry::setZReflectorCount(size_t nzrefl) {
    _z_reflector_count = nzrefl;
}

void SeismicGeometry::setZRange(double z_min, double z_max) {
    _z0 = z_min;
    _z_max = z_max;
}

double SeismicGeometry::t0() {
    return _t0;
}

void SeismicGeometry::setTRange(double t_min, double t_max) {
    _t0 = t_min;
    _t_max = t_max;
}

double SeismicGeometry::tlength() {
    return nt() * dt();
}

NRLib::Volume SeismicGeometry::createDepthVolume() {
    return NRLib::Volume(_x0, _y0, _z0, _lx, _ly, zlength(), _angle);
}

NRLib::Volume SeismicGeometry::createTimeVolume() {
    return NRLib::Volume(_x0, _y0, _t0, _lx, _ly, tlength(), _angle);
}

void SeismicGeometry::printValues() {
    printf("x0: %f lx: %f nx: %u dx: %f\n", x0(), xlength(), nx(), dx());
    printf("y0: %f ly: %f ny: %u dy: %f\n", y0(), ylength(), ny(), dy());
    printf("z0: %f lz: %f nz: %u dz: %f\n", z0(), zlength(), nz(), dz());
    printf("t0: %f lt: %f nt: %u dt: %f\n", t0(), tlength(), nt(), dt());
    printf("zrefl: %u lx_surf: %f ly_surf: %f angle: %f\n", zreflectorcount(), xsurfacelength(), ysurfacelength(), angle());
    printf("nxse: %u nyse: %u\n", nxsurfaceeclipse(), nysurfaceeclipse());
}
