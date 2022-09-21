
#include "zoeppritz_pp.hpp"
#include <math.h>

ZoeppritzPP::ZoeppritzPP(void) {

}

ZoeppritzPP::~ZoeppritzPP() {

}

double ZoeppritzPP::GetReflection(double diffvp,
                                  double meanvp,
                                  double diffrho,
                                  double meanrho,
                                  double diffvs,
                                  double meanvs,
                                  double theta)
{
  double sin2theta = sin(theta) * sin(theta);
  double tan2theta = tan(theta) * tan(theta);
  double vpvs      = meanvs/meanvp;

  double a1 =  0.5*(1.0 + tan2theta);
  double a2 = -4.0*sin2theta*vpvs*vpvs;
  double a3 =  0.5*(1.0 + a2);
  return  (a1*diffvp/meanvp + a2*diffvs/meanvs + a3*diffrho/meanrho);
}
