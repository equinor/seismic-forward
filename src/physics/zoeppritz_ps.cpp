
#include <stdio.h>
#include "zoeppritz_ps.hpp"
#include <math.h>

ZoeppritzPS::ZoeppritzPS(void) {

}

ZoeppritzPS::~ZoeppritzPS() {

}

double ZoeppritzPS::GetReflection(double diffvp,
                                  double meanvp,
                                  double diffrho,
                                  double meanrho,
                                  double diffvs,
                                  double meanvs,
                                  double theta)
{
  double sin2theta = sin(theta) * sin(theta);
  double sin_theta = sin(theta);
  double cos_theta = cos(theta);

  double cos_phi = cos(asin(meanvs * sin_theta / meanvp));

  double a2 = 2*sin_theta*(meanvs*meanvs*sin2theta/ (cos_phi*meanvp*meanvp) - meanvs*cos_theta/meanvp);
  double a3 = sin_theta*(-0.5 + meanvs*meanvs*sin2theta/(meanvp*meanvp) - meanvs*cos_theta*cos_phi/meanvp)/cos_phi;
  return a2*diffvs/meanvs + a3*diffrho/meanrho;

}
