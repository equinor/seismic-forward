#ifndef TIMER_H
#define TIMER_H

#include <string>
#include <ctime>

class Timer
{
public:
  Timer(void);
  ~Timer(void);

  enum ClockType {TOTAL, STOPWATCH};               ///< Timer type

  double  Clock(ClockType type=STOPWATCH) const;  ///< Time since start (type=TOTAL) or time since last reset (type=STOPWATCH)
  double  CPU(ClockType   type=STOPWATCH) const;  ///< Time since start (type=TOTAL) or time since last reset (type=STOPWATCH)
  void    reset(void);                            ///< Reset the timer to zero.

private:

  clock_t clock_;                                 ///< Reset clock tick. Ticks at last reset.
  clock_t initial_clock_;                         ///< Initial clock ticks. Ticks at start of process.
  time_t  time_;                                  ///< Reset time. Time at last reset.
  time_t  initial_time_;                          ///< Initial time. Time at start of process.
};

#ifndef DOXYGEN_SKIP
#define TIMER_END
#elif !defined TIMER_END
#error timer.h is part of a cyclic dependency structure
#endif
#endif
