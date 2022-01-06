#include "trace.hpp"

Trace::Trace(size_t job_number,
             double x,
             double y,
             size_t i,
             size_t j)
: job_number_(job_number),
  x_(x),
  y_(y),
  i_(i),
  j_(j)
{
}
