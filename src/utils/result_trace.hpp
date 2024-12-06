#ifndef RESULT_TRACE_HPP
#define RESULT_TRACE_HPP

#include "seismic_parameters.hpp"

class Trace;

class ResultTrace
{
  public:

    ResultTrace(ModelSettings * model_settings,
                Trace         * trace,
                const size_t    nzrefl,
                const size_t    ntwt0,
                const size_t    nz0,
                const size_t    ntwts0,
                const size_t    nt_stretch,
                const size_t    noff);

    void SetIsEmpty(bool empty) { empty_ = empty; };

    NRLib::Grid2D<double> & GetTimeTrace()           { return timegrid_pos_            ;}
    NRLib::Grid2D<double> & GetPreNMOTimeTrace()     { return prenmo_timegrid_pos_     ;}
    NRLib::Grid2D<double> & GetTimeStackTrace()      { return timegrid_stack_pos_      ;}
    NRLib::Grid2D<double> & GetDepthTrace()          { return depthgrid_pos_           ;}
    NRLib::Grid2D<double> & GetDepthStackTrace()     { return depthgrid_stack_pos_     ;}
    NRLib::Grid2D<double> & GetTimeShiftTrace()      { return timeshiftgrid_pos_       ;}
    NRLib::Grid2D<double> & GetTimeShiftStackTrace() { return timeshiftgrid_stack_pos_ ;}
    NRLib::Grid2D<double> & GetTWTxReg()             { return twtx_reg_                ;}
    NRLib::Grid2D<double> & GetTWTx()                { return twtx_                    ;}
    NRLib::Grid2D<double> & GetTheta()               { return theta_                   ;}
    NRLib::Grid2D<double> & GetRefl()                { return refl_                    ;}
    NRLib::Grid2D<double> & GetOffsetPP()            { return offset_pp_               ;}
    NRLib::Grid2D<double> & GetOffsetSS()            { return offset_ss_               ;}
    NRLib::Grid2D<double> & GetOffsetPPReg()         { return offset_pp_reg_           ;}
    NRLib::Grid2D<double> & GetOffsetSSReg()         { return offset_ss_reg_           ;}

    size_t                  GetI()                   { return i_                       ;}
    size_t                  GetJ()                   { return j_                       ;}
    double                  GetX()                   { return x_                       ;}
    double                  GetY()                   { return y_                       ;}
    bool                    GetIsEmpty()             { return empty_                   ;}
    size_t                  GetJobNumber()           { return job_number_              ;}

  private:

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

    double                x_;
    double                y_;
    size_t                i_;
    size_t                j_;
    size_t                job_number_;
    bool                  empty_;

};
#endif
