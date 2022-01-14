// $Id: main.cpp 37 2014-01-20 14:45:30Z anner $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// �    Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// �    Redistributions in binary form must reproduce the above copyright notice, this list of
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

#include "nrlib/iotools/fileio.hpp"
#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/eclipsegrid/eclipsegrid.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "modelsettings.hpp"
#include "xmlmodelfile.hpp"
#include "seismic_parameters.hpp"
#include "seismic_regridding.hpp"
#include "seismic_forward.hpp"




int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("A modelfile must be provided.\n");
        printf("Usage: %s modelfile\n", argv[0]);
        exit(1);
    }

    NRLib::LogKit::SetScreenLog(NRLib::LogKit::L_Low);
    NRLib::LogKit::StartBuffering();

    std::string inputfile(argv[1]);
    bool failedModelFile = false;
    std::cout << "************************************************************************ \n";
    std::cout << "                Seismic Forward Modeling / Geo2Seis                      \n";
    std::cout << "                           ver 4.1  2016                                 \n";
    std::cout << "                   Norsk Regnesentral & Statoil                          \n";
    std::cout << "************************************************************************ \n\n";

    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed()) {
        failedModelFile = true;
    }

    if (!failedModelFile) {
      time_t t1 = time(0);   // get time now
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters);
      seismic_parameters.PrintElapsedTime(t1, "for preprocesses");
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
