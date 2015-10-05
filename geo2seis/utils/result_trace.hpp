#ifndef RESULT_TRACE_HPP
#define RESULT_TRACE_HPP

#include <seismic_parameters.hpp>
#include <utils/trace.hpp>

class ResultTrace {

  public:

    ResultTrace(SeismicParameters         &seismic_parameters,
                std::vector<double>        twt_0,
                std::vector<double>        z_0,
                std::vector<double>        twts_0,
                size_t                     time_samples_stretch,
                const std::vector<double> &offset_vec);

    void SetJobID(Trace *trace);
    void SetIsEmpty(bool empty) { empty_ = empty; };

    NRLib::Grid2D<double> &GetTimeTrace()              { return timegrid_pos_; };
    NRLib::Grid2D<double> &GetPreNMOTimeTrace()        { return prenmo_timegrid_pos_; };
    NRLib::Grid2D<double> &GetTimeStackTrace()         { return timegrid_stack_pos_; };
    NRLib::Grid2D<double> &GetDepthTrace()             { return depthgrid_pos_; };
    NRLib::Grid2D<double> &GetDepthStackTrace()        { return depthgrid_stack_pos_; };
    NRLib::Grid2D<double> &GetTimeShiftTrace()         { return timeshiftgrid_pos_; };
    NRLib::Grid2D<double> &GetTimeShiftStackTrace()    { return timeshiftgrid_stack_pos_; };
    NRLib::Grid2D<double> &GetTWTx()                   { return twtx_reg_; };
    size_t GetI() { return i_; };
    size_t GetJ() { return j_; };
    double GetX() { return x_; };
    double GetY() { return y_; };
    bool   GetIsEmpty() { return empty_; };
    size_t GetJobNumber() { return job_number_; };

  private:

    NRLib::Grid2D<double> timegrid_pos_;
    NRLib::Grid2D<double> prenmo_timegrid_pos_;
    NRLib::Grid2D<double> timegrid_stack_pos_;
    NRLib::Grid2D<double> depthgrid_pos_;
    NRLib::Grid2D<double> depthgrid_stack_pos_;
    NRLib::Grid2D<double> timeshiftgrid_pos_;
    NRLib::Grid2D<double> timeshiftgrid_stack_pos_;
    NRLib::Grid2D<double> twtx_reg_;
    double                x_;
    double                y_;
    size_t                i_;
    size_t                j_;
    size_t                job_number_;
    bool                  empty_;

};
#endif
