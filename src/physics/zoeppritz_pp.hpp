#ifndef ZOEPPRITZPP_HPP
#define ZOEPPRITZPP_HPP

#include "zoeppritz.hpp"

class ZoeppritzPP : public Zoeppritz {
  public:
    ZoeppritzPP(void);

    ~ZoeppritzPP();

  virtual double GetReflection(double diffvp,
                               double meanvp,
                               double diffrho,
                               double meanrho,
                               double diffvs = 0.0,
                               double meanvs = 1.0,
                               double theta  = 0.0);
};

#endif
