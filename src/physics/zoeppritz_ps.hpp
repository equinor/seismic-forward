#ifndef ZOEPPRITZPS_HPP
#define ZOEPPRITZPS_HPP

#include "zoeppritz.hpp"

class ZoeppritzPS : public Zoeppritz {

  public:
    ZoeppritzPS(void);

    ~ZoeppritzPS();

    virtual double GetReflection(double diffvp,
                                 double meanvp,
                                 double diffrho,
                                 double meanrho,
                                 double diffvs = 0.0,
                                 double meanvs = 1.0,
                                 double theta  = 0.0);
};

#endif
