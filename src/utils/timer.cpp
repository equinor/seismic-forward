#include "nrlib/iotools/logkit.hpp"

#include "utils/timer.hpp"

#include <cstdlib>
#if defined (unix) || defined (__unix) || defined (__unix__)
#include <sys/times.h>
#endif

//----------------
Timer::Timer(void)
//----------------
{
  clock_         = clock();
  time_          = time(0);
  initial_clock_ = clock_;
  initial_time_  = time_;
}

//-----------------
Timer::~Timer(void)
//-----------------
{
}

//---------------------
void Timer::reset(void)
//---------------------
{
  clock_ = clock();
  time_  = time(0);
}

//-------------------------------------
double Timer::CPU(ClockType type) const
//-------------------------------------
{
  double cpu = -999.0;
  switch (type) {
  case TOTAL: {
    cpu = static_cast<double>(clock() - initial_clock_)/CLOCKS_PER_SEC;
    break;
  }
  case STOPWATCH: {
    cpu = static_cast<double>(clock() - clock_)/CLOCKS_PER_SEC;
    break;
  }
  default: {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\nBug in Timer::CPU. Unknown clock type.\n");
    std::exit(1);
  }
  }
  return cpu;
}

//---------------------------------------
double Timer::Clock(ClockType type) const
//---------------------------------------
{
  switch (type) {
  case TOTAL:
    return static_cast<double>(time(nullptr) - initial_time_);

  case STOPWATCH:
    return static_cast<double>(time(nullptr) - time_);

  default:
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\nBug in Timer::Clock. Unknown clock type.\n");
    std::exit(EXIT_FAILURE);
  }
}
