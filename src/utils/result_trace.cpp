
#include "seismic_parameters.hpp"
#include "result_trace.hpp"
#include "trace.hpp"

//--------------------------------------------------------------------
ResultTrace::ResultTrace(const SeismicParameters & seismic_parameters,
                         const ModelSettings     & model_settings,
                         const Trace             & trace,
                         const size_t              nzrefl,
                         const size_t              nt,
                         const size_t              nz,
                         const size_t              ntwts0,
                         const size_t              nt_stretch,
                         const size_t              noff)
//--------------------------------------------------------------------
  : x_         ( trace.GetX()         ),
    y_         ( trace.GetY()         ),
    i_         ( trace.GetI()         ),
    j_         ( trace.GetJ()         ),
    job_number_( trace.GetJobNumber() ),
    empty_     ( false                )
{
  if (model_settings.GetNMOCorr()) {
    twtx_reg_.           Resize(nt    , noff);
    twtx_.               Resize(nzrefl, noff);
    theta_.              Resize(nzrefl, noff);
    refl_.               Resize(nzrefl, noff);
    prenmo_timegrid_pos_.Resize(nt    , noff);

    if (model_settings.GetPSSeismic()) {
      offset_pp_.    Resize(nzrefl, noff);
      offset_ss_.    Resize(nzrefl, noff);
      offset_pp_reg_.Resize(nt    , noff);
      offset_ss_reg_.Resize(nt    , noff);
    }
  }
  timegrid_pos_.Resize(nt_stretch, noff);

  if (model_settings.GetStackOutput() || model_settings.GetStormOutput()) {
    timegrid_stack_pos_.Resize(nt_stretch, 1);
  }
  if (model_settings.GetTimeshiftOutput()) {
    timeshiftgrid_pos_.      Resize(ntwts0, noff);
    timeshiftgrid_stack_pos_.Resize(ntwts0, 1);
  }
  if (model_settings.GetDepthOutput()){
    depthgrid_pos_.      Resize(nz, noff);
    depthgrid_stack_pos_.Resize(nz, 1);
  }

  empty_ = GenerateTraceOk(seismic_parameters,
                           model_settings,
                           trace.GetI(),
                           trace.GetJ());
}

//-----------------------------------------------------------------------------
bool ResultTrace::GenerateTraceOk(const SeismicParameters & seismic_parameters,
                                  const ModelSettings     & model_settings,
                                  const size_t              i,
                                  const size_t              j)
//-----------------------------------------------------------------------------
{
  const double                 const_vp    = model_settings.GetConstVp ()[1];
  const double                 const_vs    = model_settings.GetConstVs ()[1];
  const double                 const_rho   = model_settings.GetConstRho()[1];

  const NRLib::StormContGrid & vpgrid      = seismic_parameters.GetVpGrid();
  const NRLib::StormContGrid & vsgrid      = seismic_parameters.GetVsGrid();
  const NRLib::StormContGrid & rhogrid     = seismic_parameters.GetRhoGrid();
  const NRLib::StormContGrid & twtgrid     = seismic_parameters.GetTwtGrid();

  bool generate = false;

  if (twtgrid(i, j, 0) != -999) {
    size_t nk = vpgrid.GetNK();

    for (size_t k = 1; k < nk - 1 ; k++) {
      if (generate == false) {
        if (vpgrid(i, j, k) != const_vp){
          generate = true;
        }
        if (vsgrid(i, j, k) != const_vs){
          generate = true;
        }
        if (rhogrid(i, j, k) != const_rho){
          generate = true;
        }
      }
      else{
        break;
      }
    }
  }
  return !generate;
}
