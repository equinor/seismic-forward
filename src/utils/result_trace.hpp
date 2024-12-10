#ifndef RESULT_TRACE_HPP
#define RESULT_TRACE_HPP

#include "nrlib/grid/grid2d.hpp"

#include <cstddef>

class SeismicParameters;
class ModelSettings;
class Trace;

class ResultTrace
{

public:
  ResultTrace(const SeismicParameters & seismic_parameters,
              const ModelSettings     & model_settings,
              const Trace             & trace,
              const size_t              nzrefl,
              const size_t              ntwt0,
              const size_t              nz0,
              const size_t              ntwts0,
              const size_t              nt_stretch,
              const size_t              noff);

  size_t                        GetI()                   const { return i_                       ;}
  size_t                        GetJ()                   const { return j_                       ;}
  double                        GetX()                   const { return x_                       ;}
  double                        GetY()                   const { return y_                       ;}
  bool                          GetIsEmpty()             const { return empty_                   ;}
  size_t                        GetJobNumber()           const { return job_number_              ;}

  const NRLib::Grid2D<double> & GetTimeTrace()           const { return timegrid_pos_            ;}
  const NRLib::Grid2D<double> & GetPreNMOTimeTrace()     const { return prenmo_timegrid_pos_     ;}
  const NRLib::Grid2D<double> & GetTimeStackTrace()      const { return timegrid_stack_pos_      ;}
  const NRLib::Grid2D<double> & GetDepthTrace()          const { return depthgrid_pos_           ;}
  const NRLib::Grid2D<double> & GetDepthStackTrace()     const { return depthgrid_stack_pos_     ;}
  const NRLib::Grid2D<double> & GetTimeShiftTrace()      const { return timeshiftgrid_pos_       ;}
  const NRLib::Grid2D<double> & GetTimeShiftStackTrace() const { return timeshiftgrid_stack_pos_ ;}
  const NRLib::Grid2D<double> & GetTWTxReg()             const { return twtx_reg_                ;}
  const NRLib::Grid2D<double> & GetTWTx()                const { return twtx_                    ;}
  const NRLib::Grid2D<double> & GetTheta()               const { return theta_                   ;}
  const NRLib::Grid2D<double> & GetRefl()                const { return refl_                    ;}
  const NRLib::Grid2D<double> & GetOffsetPP()            const { return offset_pp_               ;}
  const NRLib::Grid2D<double> & GetOffsetSS()            const { return offset_ss_               ;}
  const NRLib::Grid2D<double> & GetOffsetPPReg()         const { return offset_pp_reg_           ;}
  const NRLib::Grid2D<double> & GetOffsetSSReg()         const { return offset_ss_reg_           ;}

  NRLib::Grid2D<double>       & GetTimeTrace()                 { return timegrid_pos_            ;}
  NRLib::Grid2D<double>       & GetPreNMOTimeTrace()           { return prenmo_timegrid_pos_     ;}
  NRLib::Grid2D<double>       & GetTimeStackTrace()            { return timegrid_stack_pos_      ;}
  NRLib::Grid2D<double>       & GetDepthTrace()                { return depthgrid_pos_           ;}
  NRLib::Grid2D<double>       & GetDepthStackTrace()           { return depthgrid_stack_pos_     ;}
  NRLib::Grid2D<double>       & GetTimeShiftTrace()            { return timeshiftgrid_pos_       ;}
  NRLib::Grid2D<double>       & GetTimeShiftStackTrace()       { return timeshiftgrid_stack_pos_ ;}
  NRLib::Grid2D<double>       & GetTWTxReg()                   { return twtx_reg_                ;}
  NRLib::Grid2D<double>       & GetTWTx()                      { return twtx_                    ;}
  NRLib::Grid2D<double>       & GetTheta()                     { return theta_                   ;}
  NRLib::Grid2D<double>       & GetRefl()                      { return refl_                    ;}
  NRLib::Grid2D<double>       & GetOffsetPP()                  { return offset_pp_               ;}
  NRLib::Grid2D<double>       & GetOffsetSS()                  { return offset_ss_               ;}
  NRLib::Grid2D<double>       & GetOffsetPPReg()               { return offset_pp_reg_           ;}
  NRLib::Grid2D<double>       & GetOffsetSSReg()               { return offset_ss_reg_           ;}

private:

  bool GenerateTraceOk(const SeismicParameters & seismic_parameters,
                       const ModelSettings     & model_settings,
                       const size_t              i,
                       const size_t              j);

  NRLib::Grid2D<double> timegrid_pos_;
  NRLib::Grid2D<double> prenmo_timegrid_pos_;
  NRLib::Grid2D<double> timegrid_stack_pos_;
  NRLib::Grid2D<double> depthgrid_pos_;
  NRLib::Grid2D<double> depthgrid_stack_pos_;
  NRLib::Grid2D<double> timeshiftgrid_pos_;
  NRLib::Grid2D<double> timeshiftgrid_stack_pos_;
  NRLib::Grid2D<double> twtx_reg_;
  NRLib::Grid2D<double> twtx_;
  NRLib::Grid2D<double> theta_;
  NRLib::Grid2D<double> refl_;

  NRLib::Grid2D<double> offset_pp_;
  NRLib::Grid2D<double> offset_ss_;
  NRLib::Grid2D<double> offset_pp_reg_;
  NRLib::Grid2D<double> offset_ss_reg_;

  const double          x_;
  const double          y_;
  const size_t          i_;
  const size_t          j_;
  const size_t          job_number_;

  bool                  empty_;
};
#endif
