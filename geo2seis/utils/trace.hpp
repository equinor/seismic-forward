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

    size_t GetI() { return i_; };
    size_t GetJ() { return j_; };
    double GetX() { return x_; };
    double GetY() { return y_; };
    size_t GetJobNumber() { return job_number_; };

  private:
    size_t job_number_;
    double x_;
    double y_;
    size_t i_;
    size_t j_;
};
#endif
