
#ifndef ZOEPPRITZ_HPP
#define ZOEPPRITZ_HPP

class Zoeppritz {
  public:
    Zoeppritz(void);

    virtual ~Zoeppritz();

    virtual double GetReflection(double diffvp,
                                 double meanvp,
                                 double diffrho,
                                 double meanrho,
                                 double diffvs = 0.0,
                                 double meanvs = 1.0,
                                 double theta  = 0.0) = 0;

  protected :
};

#endif
