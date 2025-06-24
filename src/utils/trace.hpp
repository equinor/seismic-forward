#ifndef TRACE_HPP
#define TRACE_HPP

#include <string>

class Trace {
  public:
    Trace(size_t job_number,
          double x,
          double y,
          size_t i,
          size_t j);

    size_t GetI()         const { return i_          ;}
    size_t GetJ()         const { return j_          ;}
    double GetX()         const { return x_          ;}
    double GetY()         const { return y_          ;}
    size_t GetJobNumber() const { return job_number_ ;}

  private:
    size_t job_number_;
    double x_;
    double y_;
    size_t i_;
    size_t j_;
};
#endif
