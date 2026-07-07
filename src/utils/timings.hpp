#ifndef TIMINGS_H
#define TIMINGS_H

#include <string>

#include "nrlib/iotools/logkit.hpp"

#include "utils/timer.hpp"
class Timings
{
public:
  static void    reportAll(double threshold);
  static void    reportTotal(void);

  static void    setTimeTotal(Timer & timer);
  static void    setTimeLoadCornerPointGrid(Timer & timer);
  static void    setTimeFindZValues(Timer & timer);
  static void    setTimeFindElasticParameters(Timer & timer);
  static void    setTimeForwardModelling(Timer & timer);
  static void    setTimeWriteSegy(Timer & timer);
  static void    addTimeDummy(Timer & timer);

private:
  static void    report(const std::string & text,
                        double              cpuThis,
                        double              wallThis,
                        double              threshold);

  static void    calculateRest(void);

  static double  c_total_;
  static double  w_total_;

  static double  c_rest_;
  static double  w_rest_;

  static double  c_load_cornerpoint_grid_;
  static double  w_load_cornerpoint_grid_;

  static double  c_find_zvalues_;
  static double  w_find_zvalues_;

  static double  c_find_elastic_parameters_;
  static double  w_find_elastic_parameters_;

  static double  c_forward_modelling_;
  static double  w_forward_modelling_;

  static double  c_write_segy_;
  static double  w_write_segy_;

  static double  c_dummy_;
  static double  w_dummy_;
};
#endif
