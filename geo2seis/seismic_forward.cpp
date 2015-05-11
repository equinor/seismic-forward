#include <seismic_forward.hpp>
#include <seismic_regridding.hpp>

#include <physics/wavelet.hpp>

#include <seismic_geometry.hpp>
#include <seismic_output.hpp>
#include <physics/zoeppritz.hpp>
#include <physics/zoeppritz_ps.hpp>
#include <physics/zoeppritz_pp.hpp>



void SeismicForward::seismicForward(SeismicParameters &seismic_parameters) {
  ModelSettings * model_settings = seismic_parameters.modelSettings();
  size_t nx = seismic_parameters.seismicGeometry()->nx();
  size_t ny = seismic_parameters.seismicGeometry()->ny();
  size_t nz = seismic_parameters.seismicGeometry()->nz();
  size_t nt = seismic_parameters.seismicGeometry()->nt();

  size_t nzrefl = seismic_parameters.seismicGeometry()->zreflectorcount();
  size_t ntheta = seismic_parameters.nTheta();


  double dt = seismic_parameters.seismicGeometry()->dt();
  
  NRLib::StormContGrid              &zgrid    = seismic_parameters.zGrid();
  NRLib::StormContGrid              &twtgrid  = seismic_parameters.twtGrid();
  std::vector<NRLib::StormContGrid> &rgridvec = seismic_parameters.rGrids(); //hold


  NRLib::RegularSurface<double> &toptime = seismic_parameters.topTime();
  NRLib::RegularSurface<double> &bottime = seismic_parameters.bottomTime(); //hold

  NRLib::Volume volume   = seismic_parameters.seismicGeometry()->createDepthVolume();
  NRLib::Volume volume_t = seismic_parameters.seismicGeometry()->createTimeVolume();

  double tmin            = seismic_parameters.seismicGeometry()->t0();
  double dz              = seismic_parameters.seismicGeometry()->dz();
  double d1              = seismic_parameters.seismicGeometry()->z0();
  double wavelet_scale   = seismic_parameters.waveletScale();
  Wavelet * wavelet      = seismic_parameters.wavelet();

  std::vector<double> constvp = seismic_parameters.modelSettings()->GetConstVp();

  NRLib::RegularSurface<double> &bottom_eclipse = seismic_parameters.bottomEclipse();


  bool time_output      = (model_settings->GetOutputSeismicTime()                || model_settings->GetOutputTimeSegy()  
                       || model_settings->GetOutputSeismicStackTimeStorm()       || model_settings->GetOutputSeismicStackTimeSegy());
  bool depth_output     = (model_settings->GetOutputSeismicDepth()               || model_settings->GetOutputDepthSegy() 
                       || model_settings->GetOutputSeismicStackDepthStorm()      || model_settings->GetOutputSeismicStackDepthSegy());
  bool timeshift_output = (model_settings->GetOutputSeismicTimeshift()           || model_settings->GetOutputTimeshiftSegy()
                       || model_settings->GetOutputSeismicStackTimeShiftStorm()  || model_settings->GetOutputSeismicStackTimeShiftSegy());  
  bool stack_output     = (model_settings->GetOutputSeismicStackTimeStorm()      || model_settings->GetOutputSeismicStackTimeSegy()
                       || model_settings->GetOutputSeismicStackTimeShiftStorm()  || model_settings->GetOutputSeismicStackTimeShiftSegy()
                       || model_settings->GetOutputSeismicStackDepthStorm()      || model_settings->GetOutputSeismicStackDepthSegy());  



  NRLib::StormContGrid *twt_timeshift;
  if (model_settings->GetTwtFileName() != "") {
      twt_timeshift = new NRLib::StormContGrid(model_settings->GetTwtFileName());
      if (twtgrid.GetNI() != nx || twtgrid.GetNJ() != ny || twtgrid.GetNK() != nzrefl) {
          printf("TWT from file has wrong dimension. Aborting. \n");
          exit(1);
      }
  } else {
      twt_timeshift = NULL;
  }



    //-------------------------------NMO CORR SEIS-----------------------------------------
  if (seismic_parameters.modelSettings()->GetNMOCorr()){

    if (time_output || depth_output || timeshift_output) { //OK?
      
      
      unsigned long seed             = seismic_parameters.modelSettings()->GetSeed();
      double        deviation        = seismic_parameters.modelSettings()->GetStandardDeviation();

      NRLib::StormContGrid &vrmsgrid = seismic_parameters.vrmsGrid();
      NRLib::StormContGrid &twtxgrid = seismic_parameters.twtxGrid();


      //find max twt for seismic grid - should handle the largest offset.
      std::vector<double> twt_0      = seismic_parameters.twt_0();
      std::vector<double> z_0        = seismic_parameters.z_0();
      std::vector<double> offset_vec = seismic_parameters.offset_vec();

 
      if (true) {

        double nzrefl = twtgrid.GetNK();
        std::vector<double>               twt_vec(nzrefl),   vrms_vec(nzrefl);
        std::vector<std::vector<double> > theta_pos(nzrefl), refl_pos(nzrefl), twtx_pos(nzrefl);
        std::vector<double>               n_offset(offset_vec.size());
        for (size_t k = 0; k < nzrefl; ++k) {
          theta_pos[k] = n_offset;
          refl_pos[k]  = n_offset;
          twtx_pos[k]  = n_offset;
        }
        std::vector<std::vector<double> > twtx_pos_reg(twt_0.size()),     timegrid_pos(twt_0.size());
        std::vector<std::vector<double> > nmo_timegrid_pos(twt_0.size()), nmo_timegrid_stack_pos(twt_0.size());
        for (size_t k = 0; k < twt_0.size(); ++k) {
          twtx_pos_reg[k]           = n_offset; 
          timegrid_pos[k]           = n_offset;
          nmo_timegrid_pos[k]       = n_offset;
          nmo_timegrid_stack_pos[k] = n_offset;
        }

        size_t i = 260; //for loop 782
        size_t j = 550; //for loop 1408
        //----------------------BEGIN GEN SEIS WITH NMO FOR I,J---------------------------------
        for (size_t k = 0; k < nzrefl; ++k) {
          twt_vec[k]   = twtgrid(i,j,k);
          vrms_vec[k]  = vrmsgrid(i,j,k);
        }
        
        //sample vrms to regular grid:
        printf("sample vrms\n");
        std::vector<double> vrms_vec_reg = interpol1(twt_vec, vrms_vec, twt_0);
     
        //find theta ------------ for each reflection for each offset:
        printf("find theta \n");
        findThetaPos(theta_pos, twt_vec, vrms_vec, offset_vec);

        //find reflection coeff - for each reflection for each offset:
        printf("find reflections \n");
        findReflectionsPos(seismic_parameters, refl_pos, theta_pos, offset_vec, i, j);
        
        //add noise to reflections
        if (model_settings->GetWhiteNoise()) {
          printf("add noise to reflections \n");
          SeismicRegridding::addNoiseToReflectionsPos(seed+i, deviation, refl_pos); //nb, make unique seed when i and j loop is made
        }

        //find twtx ------------- for each reflection for each offset:
        printf("find twtx \n");
        findTWTxPos(twtx_pos, twt_vec, vrms_vec, offset_vec);

        //generate seismic
        printf("generate seismic \n");
        generateSeismicPos(timegrid_pos,
                           refl_pos,
                           twtx_pos,
                           zgrid,
                           toptime,
                           wavelet,
                           wavelet_scale,
                           offset_vec,
                           tmin,
                           dt,
                           i,
                           j);
        
        //sample twtx to regular grid - for each reflection for each offset:
        printf("find twtx regular \n");
        findTWTxPos(twtx_pos_reg, twt_0, vrms_vec_reg, offset_vec);
        
        //NMO correction:
        printf("NMO correction \n");
        nmoCorrInterpol1Pos(twt_0, timegrid_pos, twtx_pos_reg, nmo_timegrid_pos);
        
        if (model_settings->GetOutputSeismicStackTimeSegy() || model_settings->GetOutputSeismicStackTimeStorm()) {
          float noffset_inv = static_cast<float>(1.0 / offset_vec.size());
          for (size_t off = 0; off < offset_vec.size(); ++off) {
            for (size_t k = 0; k < twt_0.size(); ++k) {
              nmo_timegrid_stack_pos[k][off] += noffset_inv * nmo_timegrid_pos[k][off];
            }
          }
        }

        //debug print
        seismic_parameters.seismicOutput()->printVector(twt_0, "twt_0.txt");  
        seismic_parameters.seismicOutput()->printVector(twt_vec, "twt_vec.txt");  
        seismic_parameters.seismicOutput()->printVector(vrms_vec, "vrms_vec.txt");  
        seismic_parameters.seismicOutput()->printVector(vrms_vec_reg, "vrms_vec_reg.txt");   
        seismic_parameters.seismicOutput()->printMatrix(theta_pos, "theta_pos.txt");
        seismic_parameters.seismicOutput()->printMatrix(refl_pos, "refl_pos.txt");
        seismic_parameters.seismicOutput()->printMatrix(twtx_pos, "twtx_pos.txt");
        seismic_parameters.seismicOutput()->printMatrix(timegrid_pos, "timegrid_pos.txt");
        seismic_parameters.seismicOutput()->printMatrix(twtx_pos_reg, "twtx_pos_reg.txt");
        seismic_parameters.seismicOutput()->printMatrix(nmo_timegrid_pos, "nmo_timegrid_pos.txt");



      }
      if (false) {

        double offset_0                = seismic_parameters.offset0();
        double doffset                 = seismic_parameters.dOffset();
        double noffset                 = seismic_parameters.nOffset();


        NRLib::StormContGrid              vrms_reg     = NRLib::StormContGrid(volume_t, nx, ny, twt_0.size()); //hold
        NRLib::StormContGrid              twtxgrid_reg = NRLib::StormContGrid(volume_t, nx, ny, twt_0.size()); //hold
        std::vector<NRLib::StormContGrid> timegridvec(1);                                                      //hold
        timegridvec[0]                                 = NRLib::StormContGrid(volume_t, nx, ny, twt_0.size()); //hold
        NRLib::StormContGrid              nmo_timegrid = NRLib::StormContGrid(volume_t, nx, ny, twt_0.size()); 
        NRLib::StormContGrid              nmo_timegrid_stack(0,0,0); 
        if (stack_output == true) {
          nmo_timegrid_stack = NRLib::StormContGrid(volume_t, nx, ny, twt_0.size());
        }
        //grid under kan opprettes og slettes før og etter der de trengs - må vurderes mtp minne, men usikkert hva som lønner seg..
        NRLib::StormContGrid zt_reg(0,0,0);
        NRLib::StormContGrid nmo_depthgrid(0,0,0);
        NRLib::StormContGrid nmo_depthgrid_stack(0,0,0);
        if (depth_output) {
          zt_reg              = NRLib::StormContGrid(volume, nx, ny, twt_0.size()); //volume_t?
          nmo_depthgrid       = NRLib::StormContGrid(volume, nx, ny, nz);
        }
        NRLib::StormContGrid tst_reg(0,0,0); 
        NRLib::StormContGrid nmo_timeshiftgrid(0,0,0);
        NRLib::StormContGrid nmo_timeshiftgrid_stack(0,0,0);
        if (timeshift_output){
          tst_reg                 = NRLib::StormContGrid(volume_t, nx, ny, twt_0.size()); //volume_t?
          nmo_timeshiftgrid       = NRLib::StormContGrid(volume_t, nx, ny, nt);
        }

        printf("sample vrms to a regular grid in dt: \n");
        regSamplInterpol1(twtgrid, vrmsgrid, twt_0, vrms_reg);

        double offset   = offset_0;
        //------------------LOOP OVER OFFSET------------------------------
        for (size_t i = 0; i < noffset; ++i) {
          
          printf("find theta \n");
          findTheta(seismic_parameters, offset);
          seismic_parameters.seismicOutput()->writeTheta(seismic_parameters, offset);
          printf("find reflections \n");
          //seismic_parameters.seismicOutput()->writeNMOReflections(seismic_parameters, offset, true);
          findReflections(seismic_parameters);
          printf("Reflections found \n");

          std::vector<NRLib::StormContGrid> &rgridvec = seismic_parameters.rGrids();
          unsigned long seed                          = seismic_parameters.modelSettings()->GetSeed();
          double        deviation                     = seismic_parameters.modelSettings()->GetStandardDeviation();
          if (model_settings->GetOutputReflections()) {
            seismic_parameters.seismicOutput()->writeNMOReflections(seismic_parameters, offset, false);
          }
          if (model_settings->GetWhiteNoise()) {
            SeismicRegridding::addNoiseToReflections(seed+i, deviation, rgridvec);
            if (model_settings->GetOutputReflections()) {
              seismic_parameters.seismicOutput()->writeNMOReflections(seismic_parameters, offset, true);
            }
          }
          printf("find twtx \n");
          findTWTx(twtxgrid, vrmsgrid, twtgrid, offset);
          printf("generate seismic \n");
          generateSeismic(timegridvec,
                          rgridvec,
                          twtxgrid,
                          zgrid,
                          toptime,
                          wavelet,
                          wavelet_scale,
                          tmin,
                          dt);
          printf("find twtx on a regular grid in dt \n");
          //twtx on a regular grid in dt
          findRegTWTx(twtxgrid_reg, vrms_reg, twt_0, offset); //should this be calculated directly from twtx? 
          printf("NMO correction \n");
          //NMO correction:
          nmoCorrInterpol1(twt_0, timegridvec[0], twtxgrid_reg, nmo_timegrid);

          if (model_settings->GetOutputTimeSegy()) {
            seismic_parameters.seismicOutput()->writeNMOSeismicTimeSegy(seismic_parameters, nmo_timegrid, offset);
          }
          if (model_settings->GetOutputSeismicTime()) {
            seismic_parameters.seismicOutput()->writeNMOSeismicTimeStorm(seismic_parameters, nmo_timegrid, offset);
          }
          if (model_settings->GetOutputSeismicStackTimeSegy() || model_settings->GetOutputSeismicStackTimeStorm()) {
            float noffset_inv = static_cast<float>(1.0 / noffset);
            for (size_t i = 0; i < nx; ++i) {
              for (size_t j = 0; j < ny; ++j) {
                for (size_t k = 0; k < nt; ++k) {
                  nmo_timegrid_stack(i, j, k) += noffset_inv * nmo_timegrid(i, j, k);
                }
              }
            }
          }

          if (depth_output) {
            regSamplInterpol1(twtgrid, zgrid, twt_0, zt_reg);
            regSamplInterpol1(zt_reg, nmo_timegrid, z_0, nmo_depthgrid);
            if (model_settings->GetOutputDepthSegy()) {

            }
            if (model_settings->GetOutputSeismicDepth()) {

            }
          }

          if (timeshift_output) {
            regSamplInterpol1(twtgrid, *twt_timeshift, twt_0, tst_reg);
            regSamplInterpol1(tst_reg, nmo_timegrid, twt_0, nmo_timeshiftgrid);
            if (model_settings->GetOutputTimeshiftSegy()) {

            }
            if (model_settings->GetOutputSeismicTimeshift()) {

            }
          }
        

          //print vrms,vrms_reg, twtx, twtx_reg? other?
        
          offset += doffset;
        }
        //------------------END LOOP OVER OFFSET------------------------------
        //write stack time
        if (model_settings->GetOutputSeismicStackTimeSegy()) {
          //seismic_parameters.seismicOutput()->writeNMOSeismicTimeSegy(seismic_parameters, nmo_timegrid_stack, offset);
        }
        if (model_settings->GetOutputSeismicStackTimeStorm()) {
          //seismic_parameters.seismicOutput()->writeNMOSeismicTimeStorm(seismic_parameters, nmo_timegrid_stack, offset, true);
        }

        if (model_settings->GetOutputSeismicStackDepthSegy()     || model_settings->GetOutputSeismicStackDepthStorm()) {
          nmo_depthgrid_stack = NRLib::StormContGrid(volume, nx, ny, nz);
          regSamplInterpol1(twtgrid, zgrid, twt_0, zt_reg);
          regSamplInterpol1(zt_reg, nmo_timegrid_stack, z_0, nmo_depthgrid_stack);
          if (model_settings->GetOutputSeismicStackDepthSegy()) {

          }
          if (model_settings->GetOutputSeismicStackDepthStorm()) {

          }
          //print and delete
        }
        if (model_settings->GetOutputSeismicStackTimeShiftSegy() || model_settings->GetOutputSeismicStackTimeShiftStorm()) {
          nmo_timeshiftgrid_stack = NRLib::StormContGrid(volume_t, nx, ny, nt);
          regSamplInterpol1(twtgrid, *twt_timeshift, twt_0, tst_reg);
          regSamplInterpol1(tst_reg, nmo_timegrid_stack, twt_0, nmo_timeshiftgrid_stack);
          if (model_settings->GetOutputSeismicStackTimeShiftSegy()) {

          }
          if (model_settings->GetOutputSeismicStackTimeShiftStorm()) {

          }
          //print and delete
        }
      }
    }
  }
  else {

    //  cout << "Number of GB ram needed without writing to file:" << memory_use/1000000000.0 << "\n";
    //  double minimum_mem_use = static_cast<double>(4.0*(nx*ny*nzrefl*(2+ntheta)));
    //  cout << "Number of GB ram needed with writing to file:" << minimum_mem_use/1000000000.0 << "\n";
    double memory_use     = static_cast<double>(4.0 * (nx * ny * nzrefl * (2 + ntheta) + nx * ny * nz * ntheta * depth_output + nx * ny * nt * ntheta * time_output + (0.5 * nx * ny * nz))); // last term is a buffer
    double memory_limit   = model_settings->GetMemoryLimit();

    ///-----------------------------------------------------------------------------------------------
    ///----------------------------MEMORY USE < MEMORY LIMIT------------------------------------------
    ///-----------------------------------------------------------------------------------------------
    if (memory_use < memory_limit) {
      // ---- Create grids for seismic data
      std::vector<NRLib::StormContGrid> timegridvec(ntheta);
      std::vector<NRLib::StormContGrid> depthgridvec(ntheta);
      std::vector<NRLib::StormContGrid> timeshiftgridvec(ntheta);

      NRLib::StormContGrid timegrid(0, 0, 0);
      if (time_output) {
        timegrid = NRLib::StormContGrid(volume_t, nx, ny, nt);
      }
      NRLib::StormContGrid timeshiftgrid(0, 0, 0);
      if (timeshift_output) {
        timeshiftgrid = NRLib::StormContGrid(volume_t, nx, ny, nt);
      }
      NRLib::StormContGrid depthgrid(0, 0, 0);
      if (depth_output) {
        depthgrid = NRLib::StormContGrid(volume, nx, ny, nz);
      }

      for (size_t i = 0; i < ntheta; i++) {
        timegridvec[i] = timegrid;
        depthgridvec[i] = depthgrid;
        timeshiftgridvec[i] = timeshiftgrid;
      }

      //-----------------Generate seismic----------------------------
      if (model_settings->GetPSSeismic()) {
        printf("\nGenerating PS seismic:\n");
      }

      generateSeismicOLD(rgridvec,
                      twtgrid,
                      zgrid,
                      *twt_timeshift,
                      timegridvec,
                      depthgridvec,
                      timeshiftgridvec,
                      wavelet,
                      dt,
                      bottom_eclipse,
                      toptime,
                      tmin, dz, d1,
                      constvp,
                      wavelet_scale,
                      time_output,
                      depth_output,
                      timeshift_output);


      seismic_parameters.deleteZandRandTWTGrids();
      if (twt_timeshift != NULL) {
        delete twt_timeshift;
      }


      if (model_settings->GetOutputTimeSegy()) {
        seismic_parameters.seismicOutput()->writeSeismicTimeSegy(seismic_parameters, timegridvec);
      }

      if (model_settings->GetOutputSeismicTime()) {
        seismic_parameters.seismicOutput()->writeSeismicTimeStorm(seismic_parameters, timegridvec);
      }

      if (model_settings->GetOutputTimeshiftSegy()) {
        seismic_parameters.seismicOutput()->writeSeismicTimeshiftSegy(seismic_parameters, timeshiftgridvec);
      }

      if (model_settings->GetOutputSeismicTimeshift()) {
        seismic_parameters.seismicOutput()->writeSeismicTimeshiftStorm(seismic_parameters, timeshiftgridvec);
      }

      if (model_settings->GetOutputDepthSegy()) {
        seismic_parameters.seismicOutput()->writeSeismicDepthSegy(seismic_parameters, depthgridvec);
      }

      if (model_settings->GetOutputSeismicDepth()) {
        seismic_parameters.seismicOutput()->writeSeismicDepthStorm(seismic_parameters, depthgridvec);
      }

      if (model_settings->GetOutputSeismicStackTimeStorm() || model_settings->GetOutputSeismicStackTimeSegy()) {
        seismic_parameters.seismicOutput()->writeSeismicStackTime(seismic_parameters, timegridvec);
      }

      if (model_settings->GetOutputSeismicStackTimeShiftStorm() || model_settings->GetOutputSeismicStackTimeShiftSegy()) {
        seismic_parameters.seismicOutput()->writeSeismicStackTimeshift(seismic_parameters, timeshiftgridvec);
      }


      if (model_settings->GetOutputSeismicStackDepthStorm() || model_settings->GetOutputSeismicStackDepthSegy()) {
        seismic_parameters.seismicOutput()->writeSeismicStackDepth(seismic_parameters, depthgridvec);
      }
    }
  
    ///-----------------------------------------------------------------------------------------------
    ///----------------------------END(MEMORY USE < MEMORY LIMIT)-------------------------------------
    ///-----------------------------------------------------------------------------------------------

    ///-----------------------------------------------------------------------------------------------
    ///-------------------------------MEMORY USE TOO BIG----------------------------------------------
    ///-----------------------------------------------------------------------------------------------
    else {

      generateSeismicOnFile(rgridvec,
                            twtgrid,
                            zgrid,
                            *twt_timeshift,
                            wavelet,
                            dt,
                            nt, nz, nx, ny,
                            bottom_eclipse,
                            toptime,
                            tmin, dz, d1,
                            constvp,
                            wavelet_scale,
                            time_output,
                            depth_output,
                            timeshift_output);



      seismic_parameters.deleteZandRandTWTGrids();
      if (twt_timeshift != NULL) {
        delete twt_timeshift;
      }
      seismic_parameters.seismicOutput()->writeSeismicTimeSeismicOnFile(seismic_parameters, time_output);
      seismic_parameters.seismicOutput()->writeSeismicDepthSeismicOnFile(seismic_parameters, depth_output);
    }
    //-----------------------------------------------------------------------------------------------
    //-------------------------------END( MEMORY USE TOO BIG)----------------------------------------
    //-----------------------------------------------------------------------------------------------
    seismic_parameters.deleteWavelet();
    seismic_parameters.deleteGeometryAndOutput();
  }
}


void SeismicForward::nmoCorrInterpol1(std::vector<double> &t_in, NRLib::StormContGrid &data_in, NRLib::StormContGrid t_out, NRLib::StormContGrid &data_out){
  
  std::vector<double> data_vec_in(data_in.GetNK()), t_vec_out(data_in.GetNK()), data_vec_out(data_in.GetNK());
  
  for (size_t i = 0; i < data_in.GetNI(); i++) {
    for (size_t j = 0; j < data_in.GetNJ(); j++) {
      for (size_t k = 0; k < data_in.GetNK(); k++) {
        data_vec_in[k] = data_in(i,j,k);
        t_vec_out[k]   = t_out(i,j,k);
      }
      data_vec_out = interpol1(t_in, data_vec_in, t_vec_out);
      for (size_t k = 0; k < data_in.GetNK(); k++) {
        data_out(i,j,k) = data_vec_out[k];
      }
    }
  }
}

void SeismicForward::nmoCorrInterpol1Pos(std::vector<double>                t_in, 
                                         std::vector<std::vector<double> >  data_in, 
                                         std::vector<std::vector<double> >  t_out, 
                                         std::vector<std::vector<double> > &data_out){
  
  std::vector<double> data_vec_in(data_in.size()), t_vec_out(data_in.size()), data_vec_out(data_in.size());

  for (size_t off = 0; off < data_in[0].size(); off++) {
    for (size_t k = 0; k < data_in.size(); k++) {
      data_vec_in[k] = data_in[k][off];
      t_vec_out[k]   = t_out[k][off];
    }
    data_vec_out = interpol1(t_in, data_vec_in, t_vec_out);
    for (size_t k = 0; k < data_in.size(); k++) {
      data_out[k][off] = data_vec_out[k];
    }
  }
}




void SeismicForward::regSamplInterpol1(NRLib::StormContGrid &t_in, NRLib::StormContGrid &data_in, std::vector<double> t_out, NRLib::StormContGrid &data_out){
  
  std::vector<double> data_vec_in(data_in.GetNK()-1), t_vec_in(data_in.GetNK()-1), data_vec_out(data_in.GetNK()-1);
  
  for (size_t i = 0; i < data_in.GetNI(); i++) {
    for (size_t j = 0; j < data_in.GetNJ(); j++) {
      for (size_t k = 0; k < data_in.GetNK()-1; k++) { //dimension? 
        data_vec_in[k] = data_in(i,j,k);
        t_vec_in[k]    = t_in(i,j,k);
      }
      data_vec_out = interpol1(t_vec_in, data_vec_in, t_out);
      for (size_t k = 0; k < data_in.GetNK(); k++) {
        data_out(i,j,k) = data_vec_out[k];
      }
    }
  }
}


std::vector<double> SeismicForward::interpol1(std::vector<double> x_in, std::vector<double> y_in, std::vector<double> x_out){
  
  std::vector<double> dx(x_in.size()), dy(x_in.size()), slope(x_in.size()), intercept(x_in.size());
  std::vector<double> y_out(x_out.size());

  for (size_t i = 0; i < x_in.size() - 1; i++) {
    if (x_in[i+1] == x_in[i] && y_in[i+1] == y_in[i] && i > 0) {
      dx[i] = dx[i-1];
      dy[i] = dy[i-1];
    }
    else {
      dx[i] = x_in[i+1] - x_in[i];
      dy[i] = y_in[i+1] - y_in[i];
    }
    slope[i] = dy[i] / dx[i];
    intercept[i] = y_in[i] - x_in[i] * slope[i];
  }
  dx[x_in.size()-1]        = dx[x_in.size()-2];
  dy[x_in.size()-1]        = dy[x_in.size()-2];
  slope[x_in.size()-1]     = slope[x_in.size()-2];
  intercept[x_in.size()-1] = intercept[x_in.size()-2];
  for (size_t i = 0; i < y_out.size(); i++) {
    size_t idx = findNearestNeighbourIndex(x_out[i], x_in);
    if (idx == -1)
      idx = 0;
    y_out[i] = slope[idx] * x_out[i] + intercept[idx];
  }

  return y_out;
}

size_t SeismicForward::findNearestNeighbourIndex(double x, std::vector<double> x_in){
  double dist = 1e10;
  int idx = -1;
  for (size_t i = 0; i < x_in.size(); ++i) {
    float newDist = x - x_in[i];
    if ( newDist > 0 && newDist < dist ) {
      dist = newDist;
      idx = i;
    }
  }
  return idx;
}

//move to seismic_parameters.cpp
void SeismicForward::findReflections(SeismicParameters &seismic_parameters){

  NRLib::StormContGrid &vpgrid                = seismic_parameters.vpGrid();
  NRLib::StormContGrid &vsgrid                = seismic_parameters.vsGrid();
  NRLib::StormContGrid &rhogrid               = seismic_parameters.rhoGrid();
  std::vector<NRLib::StormContGrid> &rgridvec = seismic_parameters.rGrids();
  NRLib::StormContGrid &thetagrid             = seismic_parameters.thetaGrid();
  NRLib::StormContGrid &twtgrid               = seismic_parameters.twtGrid();


  size_t topk    = seismic_parameters.topK();
  size_t botk    = seismic_parameters.bottomK();
  bool ps_seis   = seismic_parameters.modelSettings()->GetPSSeismic();

  Zoeppritz *zoeppritz = NULL;
  if (ps_seis) {
    zoeppritz = new ZoeppritzPS();
  } else {
    zoeppritz = new ZoeppritzPP();
  }

  double diffvp, meanvp, diffvs, meanvs, diffrho, meanrho;
  
  std::vector<double> vp_vec(botk-topk+3), vs_vec(botk-topk+3), rho_vec(botk-topk+3);
  for (size_t i = 0; i < rgridvec[0].GetNI(); i++) {
    for (size_t j = 0; j <  rgridvec[0].GetNJ(); j++) {
      if (twtgrid(i, j, 0) != -999.0) {
        for (size_t k = topk; k <= botk + 2; k++) {         
          vp_vec[k - topk]  = vpgrid(i, j, (k - topk));
          vs_vec[k - topk]  = vsgrid(i, j, (k - topk));
          rho_vec[k - topk] = rhogrid(i, j, (k - topk));
        }
        for (size_t k = topk; k <= botk + 1; k++) {  
          zoeppritz->ComputeConstants(thetagrid(i,j,k));
          diffvp  =         vp_vec[k-topk + 1] - vp_vec[k-topk];
          meanvp  = 0.5 *  (vp_vec[k-topk + 1] + vp_vec[k-topk]);
          diffvs  =         vs_vec[k-topk + 1] - vs_vec[k-topk];
          meanvs  = 0.5 *  (vs_vec[k-topk + 1] + vs_vec[k-topk]);
          diffrho =        rho_vec[k-topk + 1] - rho_vec[k-topk];
          meanrho = 0.5 * (rho_vec[k-topk + 1] + rho_vec[k-topk]);
          //diffvp =        vpgrid(i, j, (k - topk) + 1) - vpgrid(i, j, (k - topk));
          //meanvp = 0.5 * (vpgrid(i, j, (k - topk) + 1) + vpgrid(i, j, (k - topk)));
          //diffvs = vsgrid(i, j, (k - topk) + 1) - vsgrid(i, j, (k - topk));
          //meanvs = 0.5 * (vsgrid(i, j, (k - topk) + 1) + vsgrid(i, j, (k - topk)));
          //diffrho = rhogrid(i, j, (k - topk) + 1) - rhogrid(i, j, (k - topk));
          //meanrho = 0.5 * (rhogrid(i, j, (k - topk) + 1) + rhogrid(i, j, (k - topk)));
          rgridvec[0](i, j, k - topk) = static_cast<float>(zoeppritz->GetReflection(diffvp, meanvp, diffrho, meanrho, diffvs, meanvs));
          
        } 
        std::cout << i << " " << j << " " << "\n";
        //printf("%d %d\n",i, k);
      }            
    }
  }
}


//move to seismic_parameters.cpp
void SeismicForward::findReflectionsPos(SeismicParameters &seismic_parameters,  
                                        std::vector<std::vector<double> > &r_vec, 
                                        std::vector<std::vector<double> > theta_vec, 
                                        std::vector<double> offset, 
                                        size_t i, 
                                        size_t j){

  NRLib::StormContGrid &vpgrid             = seismic_parameters.vpGrid();
  NRLib::StormContGrid &vsgrid             = seismic_parameters.vsGrid();
  NRLib::StormContGrid &rhogrid            = seismic_parameters.rhoGrid();

  size_t topk    = seismic_parameters.topK();
  size_t botk    = seismic_parameters.bottomK();
  bool ps_seis   = seismic_parameters.modelSettings()->GetPSSeismic();

  Zoeppritz *zoeppritz = NULL;
  if (ps_seis) {
    zoeppritz = new ZoeppritzPS();
  } else {
    zoeppritz = new ZoeppritzPP();
  }

  double diffvp, meanvp, diffvs, meanvs, diffrho, meanrho;
  std::vector<double> vp_vec(botk-topk+3), vs_vec(botk-topk+3), rho_vec(botk-topk+3);

  for (size_t off = 0; off < offset.size(); ++off) {
    for (size_t k = topk; k <= botk + 2; k++) {         
      vp_vec[k - topk]  = vpgrid(i, j, (k - topk));
      vs_vec[k - topk]  = vsgrid(i, j, (k - topk));
      rho_vec[k - topk] = rhogrid(i, j, (k - topk));
    }
    for (size_t k = topk; k <= botk + 1; k++) { 
      diffvp  =         vp_vec[k-topk + 1] - vp_vec[k-topk];
      meanvp  = 0.5 *  (vp_vec[k-topk + 1] + vp_vec[k-topk]);
      diffvs  =         vs_vec[k-topk + 1] - vs_vec[k-topk];
      meanvs  = 0.5 *  (vs_vec[k-topk + 1] + vs_vec[k-topk]);
      diffrho =        rho_vec[k-topk + 1] - rho_vec[k-topk];
      meanrho = 0.5 * (rho_vec[k-topk + 1] + rho_vec[k-topk]);
      zoeppritz->ComputeConstants(theta_vec[k][off]);
      //diffvp  =        vpgrid(i, j, (k - topk) + 1) - vpgrid(i, j, (k - topk));
      //meanvp  = 0.5 * (vpgrid(i, j, (k - topk) + 1) + vpgrid(i, j, (k - topk)));
      //diffvs  =        vsgrid(i, j, (k - topk) + 1) - vsgrid(i, j, (k - topk));
      //meanvs  = 0.5 * (vsgrid(i, j, (k - topk) + 1) + vsgrid(i, j, (k - topk)));
      //diffrho =        rhogrid(i, j, (k - topk) + 1) - rhogrid(i, j, (k - topk));
      //meanrho = 0.5 * (rhogrid(i, j, (k - topk) + 1) + rhogrid(i, j, (k - topk)));
      r_vec[k - topk][off] = static_cast<float>(zoeppritz->GetReflection(diffvp, meanvp, diffrho, meanrho, diffvs, meanvs));
    }
  }
}

//move to seismic_parameters.cpp
void SeismicForward::findTheta(SeismicParameters &seismic_parameters, double offset){
  NRLib::StormContGrid &twtgrid             = seismic_parameters.twtGrid();
  NRLib::StormContGrid &vrmsgrid            = seismic_parameters.vrmsGrid();
  NRLib::StormContGrid &thetagrid           = seismic_parameters.thetaGrid();

  for (size_t i = 0; i < thetagrid.GetNI(); i++) {
    for (size_t j = 0; j < thetagrid.GetNJ(); j++) {
      for (size_t k = 0; k < thetagrid.GetNK(); k++) {
        double tmp = offset * (vrmsgrid(i,j,k)*twtgrid(i,j,k));
        thetagrid(i,j,k) = atan(tmp);
      }
    }
  }
}

void SeismicForward::findThetaPos(std::vector<std::vector<double> > & thetagrid, std::vector<double> twt_vec, std::vector<double> vrms_vec, std::vector<double> offset){

  for (size_t off = 0; off < offset.size(); off++) {
    for (size_t k = 0; k < twt_vec.size(); k++) {
      double tmp = offset[off] / (std::sqrt(vrms_vec[k])*twt_vec[k] / 2000); //sqrt av vrms???
      thetagrid[k][off] = atan(tmp); //HOLD, sjekkk dimensjon
    }
  }
}

//move to seismic_parameters.cpp?
void SeismicForward::findRegTWTx(NRLib::StormContGrid &twtxgrid,  
                                 NRLib::StormContGrid &vrmsgrid, 
                                 std::vector<double>  twt, 
                                 double               offset){
  double twtx;
  for (size_t i = 0; i < twtxgrid.GetNI(); ++i) {
    for (size_t j = 0; j < twtxgrid.GetNJ(); ++j) {
      for (size_t k = 0; k < vrmsgrid.GetNK(); ++k) {
        twtx = twt[k]*twt[k] + 2000*2000*(offset*offset/vrmsgrid(i,j,k));
        twtxgrid(i,j,k) = std::sqrt(twtx);
      }
    }
  }
}

void SeismicForward::findTWTx(NRLib::StormContGrid &twtxgrid,  
                              NRLib::StormContGrid &vrmsgrid, 
                              NRLib::StormContGrid &twtgrid,  
                              double                offset){
  double twtx;
  for (size_t i = 0; i < twtxgrid.GetNI(); ++i) {
    for (size_t j = 0; j < twtxgrid.GetNJ(); ++j) {
      for (size_t k = 0; k < vrmsgrid.GetNK(); ++k) {
        twtx = twtgrid(i,j,k)*twtgrid(i,j,k) + 2000*2000*(offset*offset/vrmsgrid(i,j,k));
        twtxgrid(i,j,k) = std::sqrt(twtx);
      }
    }
  }
}

//move to seismic_parameters.cpp?
void SeismicForward::findRegTWTxPos(std::vector<std::vector<double> > &twtx_pos, 
                                    std::vector<double>                twt_vec, 
                                    std::vector<double>                vrms_vec, 
                                    std::vector<double>                offset){
  double twtx;
  for (size_t off = 0; off < offset.size(); ++off) {
    for (size_t k = 0; k < vrms_vec.size(); ++k) {
      twtx = twt_vec[k]*twt_vec[k] + 2000*2000*offset[off]*offset[off]/vrms_vec[k];
      twtx_pos[k][off] = std::sqrt(twtx); 
    }
  }
}



void SeismicForward::findTWTxPos(std::vector<std::vector<double> > &twtx_grid,  
                                 std::vector<double>                twt_vec, 
                                 std::vector<double>                vrms_vec, 
                                 std::vector<double>                offset){
  double twtx;
  for (size_t off = 0; off < offset.size(); ++off) {
    for (size_t k = 0; k < twt_vec.size(); ++k) {
      twtx = twt_vec[k]*twt_vec[k] + 2000*2000*(offset[off]*offset[off]/vrms_vec[k]);
      twtx_grid[k][off] = std::sqrt(twtx); 
      //twtx_grid[k][off] = twt_vec[k]; //NBNB test for not NMO corr
      //std::cout << twtx_grid[k][off] << " " << k << " " << off << "\n";
      //std::cout << twt_vec[k] << " " << vrms_vec[k] << " " << offset[off] << "\n";
      //std::cout << off << " " << k << " " << twtx << "\n";
    }
    
  }
}



void SeismicForward::generateSeismic(std::vector<NRLib::StormContGrid> &timegridvec,
                                     std::vector<NRLib::StormContGrid> &rgridvec,
                                     NRLib::StormContGrid              &twtgrid,
                                     NRLib::StormContGrid              &zgrid,
                                     NRLib::RegularSurface<double>     &toptime,
                                     Wavelet *wavelet,
                                     double waveletScale,
                                     double t0,
                                     double dt)
{
  size_t nx = timegridvec[0].GetNI();
  size_t ny = timegridvec[0].GetNJ();
  size_t nt = timegridvec[0].GetNK();
  size_t nc = rgridvec[0].GetNK();
  std::vector<double> seis(rgridvec.size());
  double x, y, z;
  double topt        = 0.0; //does this need to be set to zero?
  double rickerLimit = wavelet->GetDepthAdjustmentFactor();

  printf("\nComputing synthetic seismic:");
  float monitorSize = std::max(1.0f, static_cast<float>(nx * ny) * 0.02f);
  float nextMonitor = monitorSize;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";
  
  //  double start, stop;
  //  start=omp_get_wtime();
  //  #ifdef WITH_OMP
  // int n_threads = omp_get_num_procs();
  // int n_threads = 1;
  //  printf("Antall tr�der: %d \n", n_threads);
  //#pragma omp parallel for private(x,y,z,topt)  num_threads(n_threads)
  //#endif

  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      zgrid.FindCenterOfCell(i, j, 0, x, y, z);
      topt = toptime.GetZ(x, y);
      if (toptime.IsMissing(topt) == false) {
        double t = t0 + 0.5 * dt;

        for (size_t k = 0; k < nt; k++) {
          for (size_t l = 0; l < rgridvec.size(); l++) {
            seis[l] = 0.0;
          }
          for (size_t kk = 0; kk < nc; kk++) {
            if (fabs(twtgrid(i, j, kk) - t) < rickerLimit) {
              double ricker = waveletScale * wavelet->FindWaveletPoint(twtgrid(i, j, kk) - t);
              for (size_t l = 0; l < rgridvec.size(); l++) {
                seis[l] += rgridvec[l](i, j, kk) * ricker;
              }
            }
          }
          for (size_t l = 0; l < rgridvec.size(); l++) {
            timegridvec[l](i, j, k) = static_cast<float>(seis[l]);
          }
          t = t + dt;
        }
      } 
      else {
        for (size_t k = 0; k < nt; k++){
          for (size_t l = 0; l < rgridvec.size(); l++) {
            timegridvec[l](i, j, k) = 0.0;
          }
        }
      }
      
      if (ny * i + j + 1 >= static_cast<size_t>(nextMonitor)) {
        nextMonitor += monitorSize;
        std::cout << "^";
        fflush(stdout);
      }
    }
  }
  std::cout << "\n";
  // stop=omp_get_wtime();
  // double par_time=stop-start;

  // printf("Time: %f \n", par_time);
}


void SeismicForward::generateSeismicPos(std::vector<std::vector<double> > &timegrid_pos,
                                        std::vector<std::vector<double> > refl_pos,
                                        std::vector<std::vector<double> > twtx_pos,
                                        NRLib::StormContGrid              &zgrid,
                                        NRLib::RegularSurface<double>     &toptime,
                                        Wavelet                           *wavelet,
                                        double                            waveletScale,
                                        std::vector<double>               offset,
                                        double                            t0,
                                        double                            dt,
                                        size_t                            i,
                                        size_t                            j)
{

  size_t nt = timegrid_pos.size(); //hold, sjekk
  size_t nc = refl_pos.size();     // hold, sjekk
  double seis;
  double x, y, z;
  double topt        = 0.0; //does this need to be set to zero?
  double rickerLimit = wavelet->GetDepthAdjustmentFactor();

  //std::cout << "nc = " << nc << " nt = " << nt << "\n";

  zgrid.FindCenterOfCell(i, j, 0, x, y, z);
  topt = toptime.GetZ(x, y);
  if (toptime.IsMissing(topt) == false) {
    for (size_t off = 0; off < offset.size(); off++) {
     //std::cout << "offset= " << offset[off] << "\n";
      double t = t0 + 0.5 * dt;

      for (size_t k = 0; k < nt; k++) {
        seis = 0.0;
        for (size_t kk = 0; kk < nc; kk++) {
          //std::cout << twtx_pos[kk][off] << " " << kk << " " << off << " " << t << "\n";
          if (fabs(twtx_pos[kk][off] - t) < rickerLimit) {            
            double ricker = waveletScale * wavelet->FindWaveletPoint(twtx_pos[kk][off] - t);
            seis += refl_pos[kk][off] * ricker;    
            //std::cout << "ricker = " << ricker << "\n";
          }
        }
        //std::cout << "seis = " << seis <<  "\n";
        timegrid_pos[k][off] = static_cast<float>(seis);
        //std::cout << " k= " << k << " nt= " << nt << "\n";
        t = t + dt;
      }
    } 
  }
  else {
    for (size_t k = 0; k < nt; k++){
      for (size_t off = 0; off < offset.size(); off++) {
        timegrid_pos[k][off] = 0.0;
      }
    }
  }
}


void SeismicForward::generateSeismicOLD(std::vector<NRLib::StormContGrid> &rgridvec,
                                     NRLib::StormContGrid &twtgrid,
                                     NRLib::StormContGrid &zgrid,
                                     NRLib::StormContGrid &twt_timeshift,
                                     std::vector<NRLib::StormContGrid> &timegridvec,
                                     std::vector<NRLib::StormContGrid> &depthgridvec,
                                     std::vector<NRLib::StormContGrid> &timeshiftgridvec,
                                     Wavelet *wavelet,
                                     double dt,
                                     NRLib::RegularSurface<double> &bot,
                                     NRLib::RegularSurface<double> &toptime,
                                     double t0, double dz, double z0,
                                     std::vector<double> &constvp,
                                     double waveletScale,
                                     bool time_output,
                                     bool depth_output,
                                     bool timeshift_output) 
{
  size_t nt = timegridvec[0].GetNK();
  size_t nc = rgridvec[0].GetNK();
  size_t nz = depthgridvec[0].GetNK();
  std::vector<double> seis(rgridvec.size());
  double x, y, z;
  size_t nx, ny;
  double rickerLimit = wavelet->GetDepthAdjustmentFactor();

  if (time_output == true) {
    nx = timegridvec[0].GetNI();
    ny = timegridvec[0].GetNJ();
  } else {
    nx = depthgridvec[0].GetNI();
    ny = depthgridvec[0].GetNJ();
  }

  printf("\nComputing synthetic seismic:");
  float monitorSize = std::max(1.0f, static_cast<float>(nx * ny) * 0.02f);
  float nextMonitor = monitorSize;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";
  double topt = 0.0;
  //  double start, stop;
  //  start=omp_get_wtime();
  //  #ifdef WITH_OMP
  // int n_threads = omp_get_num_procs();
  // int n_threads = 1;
  //  printf("Antall tr�der: %d \n", n_threads);
  //#pragma omp parallel for private(x,y,z,topt)  num_threads(n_threads)
  //#endif
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      zgrid.FindCenterOfCell(i, j, 0, x, y, z);
      topt = toptime.GetZ(x, y);
      if (time_output == true) {
        if (toptime.IsMissing(topt) == false) {
          double t = t0 + 0.5 * dt;

          for (size_t k = 0; k < nt; k++) {
            for (size_t l = 0; l < rgridvec.size(); l++) {
              seis[l] = 0.0;
            }
            for (size_t kk = 0; kk < nc; kk++) {
              if (fabs(twtgrid(i, j, kk) - t) < rickerLimit) {
                double ricker = waveletScale * wavelet->FindWaveletPoint(twtgrid(i, j, kk) - t);
                for (size_t l = 0; l < rgridvec.size(); l++) {
                  seis[l] += rgridvec[l](i, j, kk) * ricker;
                }
              }
            }
            for (size_t l = 0; l < rgridvec.size(); l++) {
              timegridvec[l](i, j, k) = static_cast<float>(seis[l]);
            }
            t = t + dt;
          }
        } 
        else {
          for (size_t k = 0; k < nt; k++){
            for (size_t l = 0; l < rgridvec.size(); l++) {
              timegridvec[l](i, j, k) = 0.0;
            }
          }
        }
      }

      if (depth_output == true) {
        if (toptime.IsMissing(topt) == false) {
          double z = z0 + 0.5 * dz;

          std::vector<double> zvec(zgrid.GetNK() + 2);
          std::vector<double> twt(zgrid.GetNK() + 2);

          zvec[0] = 0.0;
          twt[0] = 0.0;
          for (size_t k = 1; k <= zgrid.GetNK(); k++) {
            zvec[k] = zgrid(i, j, k - 1);
            twt[k] = twtgrid(i, j, k - 1);
          }
          zvec[zvec.size() - 1] = bot.GetZ(x, y);
          twt[twt.size() - 1] = twt[twt.size() - 2] + 2000.0 * (zvec[zvec.size() - 1] - zvec[zvec.size() - 2]) / constvp[2];


          for (size_t k = 0; k < nz; k++) {
            for (size_t l = 0; l < rgridvec.size(); l++) {
              seis[l] = 0.0;
            }
            double t = findTFromZ(z, zvec, twt);
            for (size_t kk = 0; kk < nc; kk++) {
              if (fabs(twtgrid(i, j, kk) - t) < rickerLimit) {
                double ricker = waveletScale * wavelet->FindWaveletPoint(twtgrid(i, j, kk) - t);
                for (size_t l = 0; l < rgridvec.size(); l++) {
                  seis[l] += rgridvec[l](i, j, kk) * ricker;
                }
              }
            }
            for (size_t l = 0; l < rgridvec.size(); l++) {
              depthgridvec[l](i, j, k) = static_cast<float>(seis[l]);
            }
            z = z + dz;

          }
        } else {
          for (size_t k = 0; k < nz; k++) {
            for (size_t l = 0; l < rgridvec.size(); l++) {
              depthgridvec[l](i, j, k) = 0.0;
            }
          }
        }
      }

      if (timeshift_output == true) {
        if (toptime.IsMissing(topt) == false) {
          double tt = t0 + 0.5 * dt;

          std::vector<double> twt_shift(twt_timeshift.GetNK() + 1);
          std::vector<double> twt(twtgrid.GetNK() + 1);

          twt_shift[0] = 0.0;
          twt[0] = 0.0;
          for (size_t k = 1; k <= twt_timeshift.GetNK(); k++) {
            twt_shift[k] = twt_timeshift(i, j, k - 1);
            twt[k] = twtgrid(i, j, k - 1);
          }

          for (size_t k = 0; k < nt; k++) {
            for (size_t l = 0; l < rgridvec.size(); l++) {
              seis[l] = 0.0;
            }
            double t = findTFromZ(tt, twt_shift, twt);
            for (size_t kk = 0; kk < nc; kk++) {
              if (fabs(twtgrid(i, j, kk) - t) < rickerLimit) {
                double ricker = waveletScale * wavelet->FindWaveletPoint(twtgrid(i, j, kk) - t);
                for (size_t l = 0; l < rgridvec.size(); l++) {
                  seis[l] += rgridvec[l](i, j, kk) * ricker;
                }
              }
            }
            for (size_t l = 0; l < rgridvec.size(); l++) {
              timeshiftgridvec[l](i, j, k) = static_cast<float>(seis[l]);
            }
            tt = tt + dt;

          }
        } 
        else {
          for (size_t k = 0; k < nt; k++) {
            for (size_t l = 0; l < rgridvec.size(); l++) {
              timeshiftgridvec[l](i, j, k) = 0.0;
            }
          }
        }
      }

      if (ny * i + j + 1 >= static_cast<size_t>(nextMonitor)) {
        nextMonitor += monitorSize;
        std::cout << "^";
        fflush(stdout);
      }
    }
  }
  std::cout << "\n";
  // stop=omp_get_wtime();
  // double par_time=stop-start;

  // printf("Time: %f \n", par_time);
}


void SeismicForward::generateSeismicOnFile(std::vector<NRLib::StormContGrid> &rgridvec,
                                           NRLib::StormContGrid &twtgrid,
                                           NRLib::StormContGrid &zgrid,
                                           NRLib::StormContGrid &twt_timeshift,
                                           Wavelet *wavelet,
                                           double dt,
                                           int nt, int nz, int nx, int ny,
                                           NRLib::RegularSurface<double> &bot,
                                           NRLib::RegularSurface<double> &toptime,
                                           double t0, double dz, double z0,
                                           std::vector<double> &constvp,
                                           double waveletScale,
                                           bool time_output,
                                           bool depth_output,
                                           bool timeshift_output) 
{

  //size_t nt = timegridvec[0].GetNK();
  size_t nc = rgridvec[0].GetNK();
  // size_t nz = depthgridvec[0].GetNK();
  std::vector<double> seis(rgridvec.size());
  //double dz;
  double x, y, z;
  //size_t nx, ny;
  //double rickerLimit = 1150.0/wavelet->GetPeakFrequency();
  double rickerLimit = wavelet->GetDepthAdjustmentFactor();
  //  if(time_output == true) {
  //    nx = timegridvec[0].GetNI();
  //    ny = timegridvec[0].GetNJ();
  //  }
  //  else {
  //    nx = depthgridvec[0].GetNI();
  //    ny = depthgridvec[0].GetNJ();
  //  }
  int n_angles = rgridvec.size();
  std::ofstream *time_file = new std::ofstream[n_angles];
  std::ofstream *depth_file = new std::ofstream[n_angles];
  std::ofstream *timeshift_file = new std::ofstream[n_angles];
  if (time_output == true) {

    for (int i = 0; i < n_angles; i++) {
      std::string filename = "time_" + NRLib::ToString(i);
      NRLib::OpenWrite(time_file[i], filename, std::ios::out | std::ios::binary);
    }
  }
  if (depth_output == true) {

    for (int i = 0; i < n_angles; i++) {
      std::string filename = "depth_" + NRLib::ToString(i);
      // printf("%s\n", filename);
      NRLib::OpenWrite(depth_file[i], filename, std::ios::out | std::ios::binary);
    }
  }
  if (timeshift_output == true) {

    for (int i = 0; i < n_angles; i++) {
      std::string filename = "timeshift_" + NRLib::ToString(i);
      NRLib::OpenWrite(timeshift_file[i], filename, std::ios::out | std::ios::binary);
    }
  }
  printf("\nComputing synthetic seismic:");
  float monitorSize = std::max(1.0f, static_cast<float>(nx * ny) * 0.02f);
  float nextMonitor = monitorSize;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";
  double topt = 0.0;

  for (size_t i = 0; i < static_cast<size_t>(nx); i++) {
    for (size_t j = 0; j < static_cast<size_t>(ny); j++) {
      zgrid.FindCenterOfCell(i, j, 0, x, y, z);
      topt = toptime.GetZ(x, y);
      if (time_output == true) {
        if (toptime.IsMissing(topt) == false) {
          double t = t0 + 0.5 * dt;

          for (size_t k = 0; k < static_cast<size_t>(nt); k++) {
            for (size_t l = 0; l < rgridvec.size(); l++) {
              seis[l] = 0.0;
            }
            for (size_t kk = 0; kk < nc; kk++) {
              if (fabs(twtgrid(i, j, kk) - t) < rickerLimit) {
                double ricker = waveletScale * wavelet->FindWaveletPoint(twtgrid(i, j, kk) - t);
                for (size_t l = 0; l < rgridvec.size(); l++) {
                  seis[l] += rgridvec[l](i, j, kk) * ricker;
                }
              }
            }
            for (size_t l = 0; l < rgridvec.size(); l++) {
              NRLib::WriteBinaryFloat(time_file[l], static_cast<float>(seis[l]));
            }
            t = t + dt;
          }
        } else
          for (size_t k = 0; k < static_cast<size_t>(nt); k++)
            for (size_t l = 0; l < rgridvec.size(); l++) {
              NRLib::WriteBinaryFloat(time_file[l], 0.0);
            }

      }

      if (depth_output == true) {
        if (toptime.IsMissing(topt) == false) {
          double z = z0 + 0.5 * dz;

          std::vector<double> zvec(zgrid.GetNK() + 2);
          std::vector<double> twt(zgrid.GetNK() + 2);


          zvec[0] = 0.0;
          twt[0] = 0.0;
          for (size_t k = 1; k <= zgrid.GetNK(); k++) {
            zvec[k] = zgrid(i, j, k - 1);
            twt[k] = twtgrid(i, j, k - 1);
          }
          zvec[zvec.size() - 1] = bot.GetZ(x, y);
          twt[twt.size() - 1] = twt[twt.size() - 2] + 2000.0 * (zvec[zvec.size() - 1] - zvec[zvec.size() - 2]) / constvp[2];


          for (size_t k = 0; k < static_cast<size_t>(nz); k++) {
            for (size_t l = 0; l < rgridvec.size(); l++) {
              seis[l] = 0.0;
            }
            double t = findTFromZ(z, zvec, twt);
            for (size_t kk = 0; kk < nc; kk++) {
              if (fabs(twtgrid(i, j, kk) - t) < rickerLimit) {
                double ricker = waveletScale * wavelet->FindWaveletPoint(twtgrid(i, j, kk) - t);
                for (size_t l = 0; l < rgridvec.size(); l++) {
                  seis[l] += rgridvec[l](i, j, kk) * ricker;
                }
              }
            }
            for (size_t l = 0; l < rgridvec.size(); l++) {
              NRLib::WriteBinaryFloat(depth_file[l], static_cast<float>(seis[l]));
            }
            z = z + dz;

          }
        } else
          for (size_t k = 0; k < static_cast<size_t>(nz); k++)
            for (size_t l = 0; l < rgridvec.size(); l++) {
              NRLib::WriteBinaryFloat(depth_file[l], 0.0);
            }
      }

      if (timeshift_output == true) {
        if (toptime.IsMissing(topt) == false) {
          double tt = t0 + 0.5 * dt;

          std::vector<double> twt_shift(twt_timeshift.GetNK() + 1);
          std::vector<double> twt(twtgrid.GetNK() + 1);

          twt_shift[0] = 0.0;
          twt[0] = 0.0;
          for (size_t k = 1; k <= twt_timeshift.GetNK(); k++) {
            twt_shift[k] = twt_timeshift(i, j, k - 1);
            twt[k] = twtgrid(i, j, k - 1);
          }

          for (size_t k = 0; k < static_cast<size_t>(nt); k++) {
            for (size_t l = 0; l < rgridvec.size(); l++) {
              seis[l] = 0.0;
            }
            double t = findTFromZ(tt, twt_shift, twt);
            for (size_t kk = 0; kk < nc; kk++) {
              if (fabs(twtgrid(i, j, kk) - t) < rickerLimit) {
                double ricker = waveletScale * wavelet->FindWaveletPoint(twtgrid(i, j, kk) - t);
                for (size_t l = 0; l < rgridvec.size(); l++) {
                  seis[l] += rgridvec[l](i, j, kk) * ricker;
                }
              }
            }
            for (size_t l = 0; l < rgridvec.size(); l++) {
              NRLib::WriteBinaryFloat(timeshift_file[l], static_cast<float>(seis[l]));
            }
            tt = tt + dt;

          }
        } else
          for (size_t k = 0; k < static_cast<size_t>(nt); k++)
            for (size_t l = 0; l < rgridvec.size(); l++) {
              NRLib::WriteBinaryFloat(timeshift_file[l], 0.0);
            }

      }


      if (ny * i + j + 1 >= static_cast<size_t>(nextMonitor)) {
        nextMonitor += monitorSize;
        std::cout << "^";
        fflush(stdout);
      }

    }
  }
  std::cout << "\n";
  // printf("seismic done\n");
  delete[] time_file;
  delete[] depth_file;
  delete[] timeshift_file;
  // stop=omp_get_wtime();
  // double par_time=stop-start;

  // printf("Time: %f \n", par_time);
}

double SeismicForward::findTFromZ(double z, std::vector<double> &zvec, std::vector<double> &tvec) {
  double t;
  size_t i = 0;
  while (i < zvec.size() - 1 && z > zvec[i]) {
    i++;
  }

  if (i > 0) {
    double a = (zvec[i] - z) / (zvec[i] - zvec[i - 1]);
    t = a * tvec[i - 1] + (1 - a) * tvec[i];
  } else {
    t = tvec[0];
  }

  return t;
}