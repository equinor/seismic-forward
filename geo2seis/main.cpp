// $Id: main.cpp 37 2014-01-20 14:45:30Z anner $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// o  Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// o  Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
// SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <iostream>
#include <fstream>
#include <ctime>

#include "nrlib/eclipsegrid/eclipsegrid.hpp"
#include "nrlib/iotools/logkit.hpp"

#include "seismic_parameters.hpp"
#include "seismic_regridding.hpp"
#include "seismic_forward.hpp"
#include "xmlmodelfile.hpp"

#ifdef WITH_OMP
#include <omp.h>
#endif

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("A modelfile must be provided.\n");
    printf("Usage: %s modelfile\n", argv[0]);
    exit(1);
  }

  NRLib::LogKit::SetScreenLog(NRLib::LogKit::L_Low);
  NRLib::LogKit::StartBuffering();

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\n***************************************************************************************************");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\n*****                                                                                         *****");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\n*****                     Seismic Forward Modeling / Geo2Seis                                 *****");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\n*****                                version 4.2 beta                                         *****");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\n*****                 Copyright (c) 2017 by Norsk Regnesentral / Statoil                      *****");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\n*****                                                                                         *****");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\n***************************************************************************************************\n\n");

  std::string inputfile(argv[1]);
  XmlModelFile    modelFile(inputfile);
  ModelSettings * model_settings = modelFile.getModelSettings();

  if (!modelFile.getParsingFailed()) {
    NRLib::LogKit::SetFileLog("logfile.txt", model_settings->GetLogLevel());
    NRLib::LogKit::EndBuffering();

    //---find number of threads available and specified------------
    int n_threads       = static_cast<int>(model_settings->GetMaxThreads());
    int n_threads_avail = 1;
#ifdef WITH_OMP
    n_threads_avail = omp_get_num_procs();
#endif
    if (n_threads > n_threads_avail)
      n_threads = n_threads_avail;
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"Threads in use                            :   %3d / %3d\n",n_threads, n_threads_avail);

    NRLib::LogKit::WriteHeader("Model settings");
    model_settings->PrintSettings();

    time_t t1 = time(0);   // get time now
    NRLib::LogKit::WriteHeader("Setting up seismic parameters");
    SeismicParameters seismic_parameters = SeismicParameters(model_settings);
    NRLib::LogKit::WriteHeader("Load earth model");
    SeismicRegridding::MakeSeismicRegridding(seismic_parameters, n_threads);
    seismic_parameters.PrintElapsedTime(t1, "for preprocesses");
    NRLib::LogKit::WriteHeader("Forward modelling");
    SeismicForward::DoSeismicForward(seismic_parameters);
    //seismic_parameters.PrintElapsedTime(t1, "for total program");
  }
  else {
    printf("Press a key and then enter to continue.\n");
    int x;
    cin >> x;
  }
  NRLib::LogKit::EndLog();
}
