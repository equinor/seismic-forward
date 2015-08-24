// $Id: xmlmodelfile.cpp 41 2014-03-28 09:42:21Z vigsnes $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification
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

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "nrlib/exception/exception.hpp"
#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/iotools/fileio.hpp"
#include "xmlmodelfile.hpp"
#include "modelsettings.hpp"


XmlModelFile::XmlModelFile(const std::string &fileName) {
    modelSettings_ = new ModelSettings();
    failed_ = false;

    std::ifstream file;
    try {
        NRLib::OpenRead(file, fileName);
    } catch (NRLib::IOError e) {
        NRLib::LogKit::LogMessage(NRLib::LogKit::Error, "\nERROR: " + std::string(e.what()));
        exit(1);
    }

    //Remove all comments, since this convention is outside xml.
    std::string line;
    std::string active;
    std::string clean;
    while (std::getline(file, line)) {
        active = line.substr(0, line.find_first_of("#"));
        clean = clean + active + "\n";
    }
    file.close();

    TiXmlDocument doc;
    doc.Parse(clean.c_str());

    if (doc.Error() == true) {
        NRLib::LogKit::WriteHeader("Invalid XML file");
        NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\n%s is not a valid XML file. %s In line %d, column %d.",
                fileName.c_str(), doc.ErrorDesc(), doc.ErrorRow(), doc.ErrorCol());
        if (doc.ErrorId() == 9) {  // Not very robust check, but a start
            NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\nPossible cause: Mis-spelled or forgotten end tag.");
        }
        NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\nAborting\n");
        failed_ = true;
    } else {
        std::string errTxt = "";
        if (ParseSeismicForward(&doc, errTxt) == false) {
            errTxt = "'" + std::string(fileName) + "' is not a SeismicForward model file (lacks the <seismic-forward> keyword.)\n";
        }


        if (errTxt != "") {
            NRLib::LogKit::WriteHeader("Invalid model file");
            NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\n%s is not a valid model file:\n", fileName.c_str());
            NRLib::LogKit::LogMessage(NRLib::LogKit::Error, errTxt);
            NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\nAborting\n");
            failed_ = true;
        }

    }

}

XmlModelFile::~XmlModelFile() {

}

bool XmlModelFile::ParseSeismicForward(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("seismic-forward");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("elastic-param");
    legalCommands.push_back("angle");
    legalCommands.push_back("nmo-stretch");
    legalCommands.push_back("output-grid");
    legalCommands.push_back("wavelet");
    legalCommands.push_back("white-noise");
    legalCommands.push_back("output-parameters");
    legalCommands.push_back("timeshift-twt");
    legalCommands.push_back("ps-seismic");


    if (ParseNMOStretch(root, errTxt)){
      modelSettings_->SetNMOCorr(true);
    }
    ParseAngle(root, errTxt);
    ParseOutputGrid(root, errTxt);
    ParseElasticParam(root, errTxt);
    ParseWavelet(root, errTxt);

    if (ParseWhiteNoise(root, errTxt)) {
        modelSettings_->SetWhiteNoise();
    }

    std::string file_name;
    if (ParseValue(root, "timeshift-twt", file_name, errTxt) == true) {
        modelSettings_->SetTwtFileName(file_name);
    }

    std::string ps;
    modelSettings_->SetPSSeismic(false);
    if (ParseValue(root, "ps-seismic", ps, errTxt) == true) if (ps == "yes") {
        modelSettings_->SetPSSeismic(true);
    }

    ParseOutputParameters(root, errTxt);

    CheckForJunk(root, errTxt, legalCommands);
    return (true);
}

bool XmlModelFile::ParseElasticParam(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("elastic-param");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("eclipse-file");
    legalCommands.push_back("default-values");
    legalCommands.push_back("parameter-names");
    legalCommands.push_back("zero-thickness-limit");
    legalCommands.push_back("cornerpt-interpolation-in-depth");
    legalCommands.push_back("extra-parameters");

    std::string value;
    if (ParseValue(root, "eclipse-file", value, errTxt) == true) {
        modelSettings_->SetEclipseGrid(value);
    } else {
        errTxt += " An Eclipse file name must be given\n";
    }

    ParseDefaultValues(root, errTxt);
    ParseParameterNames(root, errTxt);

    double val;
    if (ParseValue(root, "zero-thickness-limit", val, errTxt) == true) {
        modelSettings_->SetZeroThicknessLimit(val);
    }

    bool bolval;
    if (ParseBool(root, "cornerpt-interpolation-in-depth", bolval, errTxt) == true) {
        modelSettings_->SetUseCornerpointInterpol(bolval);
    }

    size_t i = 0;
    while (ParseExtraParameters(root, errTxt)) {
        ++i;
    }

    CheckForJunk(root, errTxt, legalCommands);
    return true;
}

bool XmlModelFile::ParseDefaultValues(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("default-values");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("vp-top");
    legalCommands.push_back("vp-mid");
    legalCommands.push_back("vp-bot");
    legalCommands.push_back("vs-top");
    legalCommands.push_back("vs-mid");
    legalCommands.push_back("vs-bot");
    legalCommands.push_back("rho-top");
    legalCommands.push_back("rho-mid");
    legalCommands.push_back("rho-bot");

    double value;
    if (ParseValue(root, "vp-top", value, errTxt) == true) {
        modelSettings_->SetVpTop(value);
    } else {
        errTxt += "Value for Vp above reservoir is not given.\n";
    }

    if (ParseValue(root, "vp-mid", value, errTxt) == true) {
        modelSettings_->SetVpMid(value);
    } else {
        errTxt += "Value for Vp in missing cells is not given.\n";
    }

    if (ParseValue(root, "vp-bot", value, errTxt) == true) {
        modelSettings_->SetVpBot(value);
    } else {
        errTxt += "Value for Vp below reservoir is not given\n";
    }

    if (ParseValue(root, "vs-top", value, errTxt) == true) {
        modelSettings_->SetVsTop(value);
    } else {
        errTxt += "Value for Vs above reservoir is not given.\n";
    }

    if (ParseValue(root, "vs-mid", value, errTxt) == true) {
        modelSettings_->SetVsMid(value);
    } else {
        errTxt += "Value for Vs in missing cells is not given.\n";
    }

    if (ParseValue(root, "vs-bot", value, errTxt) == true) {
        modelSettings_->SetVsBot(value);
    } else {
        errTxt += "Value for Vs below reservoir is not given\n";
    }

    if (ParseValue(root, "rho-top", value, errTxt) == true) {
        modelSettings_->SetRhoTop(value);
    } else {
        errTxt += "Value for Rho above reservoir is not given.\n";
    }

    if (ParseValue(root, "rho-mid", value, errTxt) == true) {
        modelSettings_->SetRhoMid(value);
    } else {
        errTxt += "Value for rho in missing cells is not given.\n";
    }

    if (ParseValue(root, "rho-bot", value, errTxt) == true) {
        modelSettings_->SetRhoBot(value);
    } else {
        errTxt += "Value for rho below reservoir is not given\n";
    }

    CheckForJunk(root, errTxt, legalCommands);
    return true;

}

bool XmlModelFile::ParseParameterNames(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("parameter-names");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("vp");
    legalCommands.push_back("vs");
    legalCommands.push_back("rho");

    std::string value;
    if (ParseValue(root, "vp", value, errTxt) == true) {
        modelSettings_->SetVpName(value);
    } else {
        errTxt += "Name for Vp parameter in Eclipse grid is not given.\n";
    }

    if (ParseValue(root, "vs", value, errTxt) == true) {
        modelSettings_->SetVsName(value);
    } else {
        errTxt += "Name for Vs parameter in Eclipse grid is not given.\n";
    }
    if (ParseValue(root, "rho", value, errTxt) == true) {
        modelSettings_->SetRhoName(value);
    } else {
        errTxt += "Name for rho parameter in Eclipse grid is not given.\n";
    }

    CheckForJunk(root, errTxt, legalCommands);
    return true;

}

bool XmlModelFile::ParseExtraParameters(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("extra-parameters");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("name");
    legalCommands.push_back("default-value");

    std::string name;
    double value;
    if (ParseValue(root, "name", name, errTxt) && ParseValue(root, "default-value", value, errTxt)) {
        modelSettings_->AddExtraParameterName(name);
        modelSettings_->AddExtraParameterDefaultValue(value);
    } else
        errTxt += "One or more keyword under command <" + root->ValueStr() + ">  on line " + NRLib::ToString(root->Row()) + ", column "
                + NRLib::ToString(root->Column()) + " is not legal or missing. <name> and <default-value> are required.\n";


    CheckForJunk(root, errTxt, legalCommands, true);
    return true;

}

bool XmlModelFile::ParseNMOStretch(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("nmo-stretch");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("seafloor-depth");
    legalCommands.push_back("velocity-water");
    legalCommands.push_back("offset");

    double value;
    if (ParseValue(root, "seafloor-depth", value, errTxt) == true) {
      modelSettings_->SetZw(value);
    } else {
      errTxt += "Value for seafloor depth is not given.\n";
    }
    if (ParseValue(root, "velocity-water", value, errTxt) == true) {
      modelSettings_->SetVw(value);
    } else {
      errTxt += "Value for velocity in water is not given.\n";
    }
    ParseOffset(root, errTxt);

    CheckForJunk(root, errTxt, legalCommands, true);
    return true;
}

bool XmlModelFile::ParseOffset(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("offset");
    if (root == 0) {
        return (false);
    }
    std::vector<std::string> legalCommands;
    legalCommands.push_back("offset-0");
    legalCommands.push_back("doffset");
    legalCommands.push_back("offset-max");

    double value;
    if (ParseValue(root, "offset-0", value, errTxt) == true) {
      modelSettings_->SetOffset0(value);
    } else {
      errTxt += "Value for minimum offset is not given.\n";
    }
    if (ParseValue(root, "doffset", value, errTxt) == true) {
      modelSettings_->SetDOffset(value);
    } else {
      errTxt += "Value for offset increment is not given.\n";
    }
    if (ParseValue(root, "offset-max", value, errTxt) == true) {
      modelSettings_->SetOffsetMax(value);
    } else {
      errTxt += "Value for maximum offset is not given.\n";
    }

    CheckForJunk(root, errTxt, legalCommands, true);
    return true;
}

bool XmlModelFile::ParseAngle(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("angle");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("theta-0");
    legalCommands.push_back("dtheta");
    legalCommands.push_back("theta-max");

    double value;
    if (ParseValue(root, "theta-0", value, errTxt) == true) {
        modelSettings_->SetTheta0(value);
    } else {
        errTxt += "Value for minimum angle is not given.\n";
    }

    if (ParseValue(root, "dtheta", value, errTxt) == true) {
        modelSettings_->SetDTheta(value);
    } else {
        errTxt += "Value for angle increment is not given.\n";
    }

    if (ParseValue(root, "theta-max", value, errTxt) == true) {
        modelSettings_->SetThetaMax(value);
    } else {
        errTxt += "Value for maximum angle is not given.\n";
    }

    CheckForJunk(root, errTxt, legalCommands);
    return true;

}

bool XmlModelFile::ParseWavelet(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("wavelet");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("ricker");
    legalCommands.push_back("scale");
    legalCommands.push_back("from-file");

    if (ParseRicker(root, errTxt) == true) {
        modelSettings_->SetRicker(true);
    } else {
        modelSettings_->SetRicker(false);
    }
    double value;
    if (ParseValue(root, "scale", value, errTxt) == true) {
        modelSettings_->SetWaveletScale(value);
    }

    if (ParseWaveletFromFile(root, errTxt) == false && modelSettings_->GetRicker() == false) {
        errTxt += "No wavelet is given. Should be given as either <ricker> or <from-file>.\n";
    }

    CheckForJunk(root, errTxt, legalCommands);
    return true;
}

bool XmlModelFile::ParseRicker(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("ricker");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("peak-frequency");
    double value;
    if (ParseValue(root, "peak-frequency", value, errTxt) == true) {
        modelSettings_->SetPeakF(value);
    } else {
        errTxt += "Value for peak frequency in Ricker wavelet is not given.\n";
    }

    CheckForJunk(root, errTxt, legalCommands);
    return true;

}

bool XmlModelFile::ParseWaveletFromFile(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("from-file");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("format");
    legalCommands.push_back("file-name");

    std::string format, file_name;
    if (ParseValue(root, "format", format, errTxt) && ParseValue(root, "file-name", file_name, errTxt)) {
        modelSettings_->SetWaveletFileFormat(format);
        modelSettings_->SetWaveletFileName(file_name);
    } else
        errTxt += "One or more keyword under command <" + root->ValueStr() + ">  on line " + NRLib::ToString(root->Row()) + ", column "
                + NRLib::ToString(root->Column()) + " is not legal or missing. <format> and <file-name> are required.\n";


    CheckForJunk(root, errTxt, legalCommands);
    return true;

}

bool XmlModelFile::ParseOutputGrid(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("output-grid");
    if (root == 0) {
        return (false);
    }
    std::vector<std::string> legalCommands;
    legalCommands.push_back("depth-window");
    legalCommands.push_back("area");
    legalCommands.push_back("depth");
    legalCommands.push_back("cell-size");
    legalCommands.push_back("area-from-surface");
    legalCommands.push_back("area-from-segy");
    legalCommands.push_back("segy-indexes");
    legalCommands.push_back("test");
    legalCommands.push_back("utm-precision");
    legalCommands.push_back("time-window");


    if (ParseArea(root, errTxt) == true) {
        modelSettings_->SetAreaGiven(true);
    } else {
        modelSettings_->SetAreaGiven(false);
    }

    bool area_from_segy = false;
    std::string value;
    if (ParseValue(root, "area-from-surface", value, errTxt) == true) {
        modelSettings_->SetAreaFromSurface(value);
    } else if (ParseAreaFromSegy(root, errTxt) == true) {
        area_from_segy = true;
    }


    ParseDepth(root, errTxt);
    ParseCellSize(root, errTxt, area_from_segy);
    ParseSegyIndexes(root, errTxt, area_from_segy);

    double utm_prec;
    if (ParseValue(root, "utm-precision", utm_prec, errTxt)) {
        if (utm_prec == 0.0001) {
            modelSettings_->SetUtmPrecision(-10000);
        } else if (utm_prec == 0.001) {
            modelSettings_->SetUtmPrecision(-1000);
        } else if (utm_prec == 0.01) {
            modelSettings_->SetUtmPrecision(-100);
        } else if (utm_prec == 0.1) {
            modelSettings_->SetUtmPrecision(-10);
        } else if (utm_prec == 1 || utm_prec == 10 || utm_prec == 100 || utm_prec == 1000 || utm_prec == 10000) {
            modelSettings_->SetUtmPrecision(static_cast<short>(utm_prec));
        } else {
            errTxt += "Value given in <utm-precision> is not valid. Optional values are 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000 or 10000.\n";
        }
    }

    ParseTimeWindow(root, errTxt);
    ParseDepthWindow(root, errTxt);

    CheckForJunk(root, errTxt, legalCommands);
    return true;
}

bool XmlModelFile::ParseAreaFromSegy(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("area-from-segy");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("il0");
    legalCommands.push_back("xl0");
    legalCommands.push_back("utmxLoc");
    legalCommands.push_back("utmyLoc");
    legalCommands.push_back("filename");

    int value;
    if (ParseValue(root, "il0", value, errTxt) == true) {
        modelSettings_->SetIL0In(value);
    }
    if (ParseValue(root, "xl0", value, errTxt) == true) {
        modelSettings_->SetXL0In(value);
    }
    if (ParseValue(root, "utmxLoc", value, errTxt) == true) {
        modelSettings_->SetUtmxIn(value);
    }
    if (ParseValue(root, "utmyLoc", value, errTxt) == true) {
        modelSettings_->SetUtmyIn(value);
    }

    std::string filename;
    if (ParseValue(root, "filename", filename, errTxt) == true) {
        modelSettings_->SetAreaFromSegy(filename);
    }
    CheckForJunk(root, errTxt, legalCommands);
    return true;

}

bool XmlModelFile::ParseArea(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("area");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("x0");
    legalCommands.push_back("y0");
    legalCommands.push_back("lx");
    legalCommands.push_back("ly");
    legalCommands.push_back("angle");

    double value;
    if (ParseValue(root, "x0", value, errTxt) == true) {
        modelSettings_->SetX0(value);
    } else {
        errTxt += "Value for x0 in area is not given.\n";
    }

    if (ParseValue(root, "y0", value, errTxt) == true) {
        modelSettings_->SetY0(value);
    } else {
        errTxt += "Value for y0 in area is not given.\n";
    }

    if (ParseValue(root, "lx", value, errTxt) == true) {
        modelSettings_->SetLx(value);
    } else {
        errTxt += "Value for lx in area is not given.\n";
    }

    if (ParseValue(root, "ly", value, errTxt) == true) {
        modelSettings_->SetLy(value);
    } else {
        errTxt += "Value for ly in area is not given.\n";
    }

    if (ParseValue(root, "angle", value, errTxt) == true) {
        modelSettings_->SetAreaAngle(value);
    } else {
        errTxt += "Value for angle in area is not given.\n";
    }

    CheckForJunk(root, errTxt, legalCommands);
    return true;


}

bool XmlModelFile::ParseDepth(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("depth");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    //legalCommands.push_back("top-surface-file"); commands not used. probably never used.
    //legalCommands.push_back("bot-surface-file");
    //legalCommands.push_back("constant-top");
    //legalCommands.push_back("constant-bot");
    legalCommands.push_back("top-time-surface");
    legalCommands.push_back("top-time-constant");
    legalCommands.push_back("nlayers-from-file");

    std::string value;
    //bool topfile, botfile;
    //topfile = false;
    //botfile = false;
    //if(ParseValue(root, "top-surface-file", value, errTxt) == true){
    //  modelSettings_->SetTopSurfaceFile(value);
    //  topfile = true;
    //}
    //if(ParseValue(root, "bot-surface-file", value, errTxt) == true){
    //  modelSettings_->SetBotSurfaceFile(value);
    //  botfile = true;
    //}
    if (ParseValue(root, "top-time-surface", value, errTxt) == true) {
        modelSettings_->SetTopTimeSurface(value);
    }
    double val;
    if (ParseValue(root, "top-time-constant", val, errTxt) == true) {
        modelSettings_->SetTopTimeConstant(val);
    }


    //if(ParseValue(root, "constant-top", val, errTxt) == true){
    //  modelSettings_->SetConstantTop(val);
    //  if(topfile == true)
    //    errTxt+= "Top depth is given both as surface and constant.\n";
    //}
    //if(ParseValue(root, "constant-bot", val, errTxt) == true){
    //  modelSettings_->SetConstantBot(val);
    //  if(botfile == true)
    //    errTxt+= "Bottom depth is given both as surface and constant.\n";
    //}
    if (ParseValue(root, "nlayers-from-file", value, errTxt) == true) {
        modelSettings_->SetNLayersFile(value);


    }
    CheckForJunk(root, errTxt, legalCommands);
    return true;

}

bool XmlModelFile::ParseCellSize(TiXmlNode *node, std::string &errTxt, bool &area_from_segy) {
    TiXmlNode *root = node->FirstChildElement("cell-size");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("dx");
    legalCommands.push_back("dy");
    legalCommands.push_back("dz");
    legalCommands.push_back("dt");

    double value;
    if (ParseValue(root, "dx", value, errTxt) == true) {
        if (area_from_segy == true) {
          printf("WARNING:: Both <dx> and <area-from-segy> is specified, and dx will be taken from segy file. \n\n");
        }
        else {
          modelSettings_->SetDx(value);
        }
    }
    // else
    //   errTxt += "Value for dx is not given.\n";

    if (ParseValue(root, "dy", value, errTxt) == true) {
        if (area_from_segy == true) {
            printf("WARNING:: Both <dy> and <area-from-segy> is specified, and dy will be taken from segy file. \n\n");
        }
        else {
          modelSettings_->SetDy(value);
        }
    }
    // else
    //   errTxt += "Value for dy is not given.\n";

    if (ParseValue(root, "dz", value, errTxt) == true) {
        modelSettings_->SetDz(value);
    }
    // else
    //   errTxt += "Value for dz is not given.\n";

    if (ParseValue(root, "dt", value, errTxt) == true) {
        modelSettings_->SetDt(value);
    }
    // else
    //   errTxt += "Value for dt is not given.\n";


    CheckForJunk(root, errTxt, legalCommands);

    return true;
}

bool XmlModelFile::ParseTimeWindow(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("time-window");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("top");
    legalCommands.push_back("bot");

    double value;
    bool top_value = false;
    if (ParseValue(root, "top", value, errTxt) == true) {
        modelSettings_->SetTopTimeWindow(value);
        modelSettings_->SetTimeWindowSpecified(true);
        top_value = true;
    }
    if (ParseValue(root, "bot", value, errTxt) == true) {
        modelSettings_->SetBotTimeWindow(value);
        modelSettings_->SetTimeWindowSpecified(true);
        if (top_value == false) {
            printf("WARNING:: Value for <top> under <time-window> is not given. Is set to top reservoir.\n\n");
        }
    } else if (top_value == true) {
        printf("WARNING:: Value for <bot> under <time-window> is not given. Is set to bottom reservoir.\n\n");
    }

    CheckForJunk(root, errTxt, legalCommands);
    return true;
}

bool XmlModelFile::ParseDepthWindow(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("depth-window");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("top");
    legalCommands.push_back("bot");

    double value;
    bool top_value = false;
    if (ParseValue(root, "top", value, errTxt) == true) {
        modelSettings_->SetTopDepthWindow(value);
        modelSettings_->SetDepthWindowSpecified(true);
        top_value = true;
    }
    if (ParseValue(root, "bot", value, errTxt) == true) {
        modelSettings_->SetBotDepthWindow(value);
        modelSettings_->SetDepthWindowSpecified(true);
        if (top_value == false) {
            printf("WARNING:: Value for <top> under <depth-window> is not given. Is set to top reservoir.\n\n");
        }
    } else if (top_value == true) {
        printf("WARNING:: Value for <bot> under <depth-window> is not given. Is set to bottom reservoir.\n\n");
    }
    CheckForJunk(root, errTxt, legalCommands);
    return true;
}

bool XmlModelFile::ParseOutputParameters(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("output-parameters");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("elastic-parameters");
    legalCommands.push_back("zvalues");
    legalCommands.push_back("reflections");
    legalCommands.push_back("seismic-time");
    legalCommands.push_back("seismic-timeshift");
    legalCommands.push_back("seismic-depth");
    legalCommands.push_back("time-surfaces");
    legalCommands.push_back("depth-surfaces");
    legalCommands.push_back("twt");
    legalCommands.push_back("vrms");
    legalCommands.push_back("twt-offset");
    legalCommands.push_back("prefix");
    legalCommands.push_back("suffix");
    legalCommands.push_back("seismic-time-segy");
    legalCommands.push_back("seismic-timeshift-segy");
    legalCommands.push_back("seismic-depth-segy");
    legalCommands.push_back("elastic-parameters-time-segy");
    legalCommands.push_back("elastic-parameters-depth-segy");
    legalCommands.push_back("extra-parameters-time-segy");
    legalCommands.push_back("extra-parameters-depth-segy");
    legalCommands.push_back("seismic-stack");
    legalCommands.push_back("seismic-prenmo-time-segy");

    bool value;
    if (ParseBool(root, "elastic-parameters", value, errTxt) == true) {
        modelSettings_->SetOutputVp(value);
    }

    if (ParseBool(root, "zvalues", value, errTxt) == true) {
        modelSettings_->SetOutputZvalues(value);
    }

    if (ParseBool(root, "reflections", value, errTxt) == true) {
        modelSettings_->SetOutputReflections(value);
    }

    if (ParseBool(root, "seismic-time", value, errTxt) == true) {
        modelSettings_->SetOutputSeismicTime(value);
    }

    if (ParseBool(root, "seismic-timeshift", value, errTxt) == true) {
        if (value == true && modelSettings_->GetTwtFileName() == "") {
            errTxt += "The command <twt-timeshift> must be set if output parameter <seismic-timeshift> is given.\n";
        } else {
            modelSettings_->SetOutputSeismicTimeshift(value);
        }
    }

    if (ParseBool(root, "seismic-depth", value, errTxt) == true) {
        modelSettings_->SetOutputSeismicDepth(value);
    }

    if (ParseBool(root, "time-surfaces", value, errTxt) == true) {
        modelSettings_->SetOutputTimeSurfaces(value);
    }

    if (ParseBool(root, "depth-surfaces", value, errTxt) == true) {
        modelSettings_->SetOutputDepthSurfaces(value);
    }

    if (ParseBool(root, "twt", value, errTxt) == true) {
        modelSettings_->SetOutputTwt(value);
    }

    if (ParseBool(root, "vrms", value, errTxt) == true) {
        modelSettings_->SetOutputVrms(value);
    }

    if (ParseBool(root, "twt-offset", value, errTxt) == true) {
        modelSettings_->SetOutputTwtOffset(value);
    }

    std::string val;
    if (ParseValue(root, "prefix", val, errTxt) == true) {
        modelSettings_->SetPrefix(val);
    }

    if (ParseValue(root, "suffix", val, errTxt) == true) {
        modelSettings_->SetSuffix(val);
    }

    if (ParseBool(root, "seismic-time-segy", value, errTxt) == true) {
        modelSettings_->SetOutputTimeSegy(value);
    }

    if (ParseBool(root, "seismic-timeshift-segy", value, errTxt) == true) {
        if (value == true && modelSettings_->GetTwtFileName() == "") {
            errTxt += "The command <twt-timeshift> must be set if output parameter <seismic-timeshift-segy> is given.\n";
        } else {
            modelSettings_->SetOutputTimeshiftSegy(value);
        }
    }
    if (ParseBool(root, "seismic-depth-segy", value, errTxt) == true) {
        modelSettings_->SetOutputDepthSegy(value);
    }

    if (ParseBool(root, "elastic-parameters-time-segy", value, errTxt) == true) {
        modelSettings_->SetOutputElasticParametersTimeSegy(value);
    }

    if (ParseBool(root, "elastic-parameters-depth-segy", value, errTxt) == true) {
        modelSettings_->SetOutputElasticParametersDepthSegy(value);
    }


    if (ParseBool(root, "extra-parameters-time-segy", value, errTxt) == true) {
        modelSettings_->SetOutputExtraParametersTimeSegy(value);
    }

    if (ParseBool(root, "extra-parameters-depth-segy", value, errTxt) == true) {
        modelSettings_->SetOutputExtraParametersDepthSegy(value);
    }

    ParseSeismicStack(root, errTxt);

    if (ParseBool(root, "seismic-prenmo-time-segy", value, errTxt) == true) {
      modelSettings_->SetOutputPrenmoTimeSegy(value);
    }

    CheckForJunk(root, errTxt, legalCommands);

    return true;
}

bool XmlModelFile::ParseSegyIndexes(TiXmlNode *node, std::string &errTxt, bool area_from_segy) {
    TiXmlNode *root = node->FirstChildElement("segy-indexes");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("inline-start");
    legalCommands.push_back("xline-start");
    legalCommands.push_back("inline-direction");
    legalCommands.push_back("inline-step");
    legalCommands.push_back("xline-step");

      int inline_start, xline_start;
      std::string inline_direction;
      if (ParseValue(root, "inline-start", inline_start, errTxt) && ParseValue(root, "xline-start", xline_start, errTxt) && ParseValue(root, "inline-direction", inline_direction, errTxt)) {
        if (area_from_segy == false) {
          modelSettings_->SetSegyInlineStart(inline_start);
          modelSettings_->SetSegyXlineStart(xline_start);
          if (NRLib::Uppercase(inline_direction) == "X" || NRLib::Uppercase(inline_direction) == "Y") {
              modelSettings_->SetSegyInlineDirection(inline_direction);
          } else
              errTxt += "\"" + inline_direction + "\" is not a legal operation under keyword <inline-direction>, under command <" + root->ValueStr() + ">  on line " + NRLib::ToString(root->Row()) + ", column "
                      + NRLib::ToString(root->Column()) + ".\n";
        }
      }
      else
        if (area_from_segy == false)
          errTxt += "One or more keyword under command <" + root->ValueStr() + ">  on line " + NRLib::ToString(root->Row()) + ", column "
                  + NRLib::ToString(root->Column()) + " is not legal or missing. Three keywords; <inline-start>, <xline-start> and <inline-direction> are required.\n";

      int inline_step, xline_step;
      if (ParseValue(root, "inline-step", inline_step, errTxt) && area_from_segy == false) {
          modelSettings_->SetSegyInlineStep(inline_step);
      }
      if (ParseValue(root, "xline-step", xline_step, errTxt) && area_from_segy == false) {
          modelSettings_->SetSegyXlineStep(xline_step);
      }

    if (area_from_segy)  {
      printf("WARNING:: Both <area-from-segy> and <segy-indexes> are given. Indexes are taken from segy file, and values given in <segy-indexes> are not used.\n\n");
    }

    CheckForJunk(root, errTxt, legalCommands);

    return true;
}

bool XmlModelFile::ParseSeismicStack(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("seismic-stack");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("time-storm");
    legalCommands.push_back("timeshift-storm");
    legalCommands.push_back("depth-storm");
    legalCommands.push_back("time-segy");
    legalCommands.push_back("timeshift-segy");
    legalCommands.push_back("depth-segy");

    bool val;
    if (ParseBool(root, "time-storm", val, errTxt) == true) {
        modelSettings_->SetOutputSeismicStackTimeStorm(val);
    }

    if (ParseBool(root, "timeshift-storm", val, errTxt) == true) if (modelSettings_->GetTwtFileName() == "") {
        errTxt += "The command <twt-timeshift> must be set if output parameter <seismic-stack><time-shift-storm> is given.\n";
    } else {
        modelSettings_->SetOutputSeismicStackTimeShiftStorm(val);
    }

    if (ParseBool(root, "depth-storm", val, errTxt) == true) {
        modelSettings_->SetOutputSeismicStackDepthStorm(val);
    }

    if (ParseBool(root, "time-segy", val, errTxt) == true) {
        modelSettings_->SetOutputSeismicStackTimeSegy(val);
    }

    if (ParseBool(root, "timeshift-segy", val, errTxt) == true) if (modelSettings_->GetTwtFileName() == "") {
        errTxt += "The command <twt-timeshift> must be set if output parameter <seismic-stack><time-shift-segy> is given.\n";
    } else {
        modelSettings_->SetOutputSeismicStackTimeShiftSegy(val);
    }

    if (ParseBool(root, "depth-segy", val, errTxt) == true) {
        modelSettings_->SetOutputSeismicStackDepthSegy(val);
    }

    CheckForJunk(root, errTxt, legalCommands);

    return true;
}

bool XmlModelFile::ParseWhiteNoise(TiXmlNode *node, std::string &errTxt) {
    TiXmlNode *root = node->FirstChildElement("white-noise");
    if (root == 0) {
        return (false);
    }

    std::vector<std::string> legalCommands;
    legalCommands.push_back("standard-deviation");
    legalCommands.push_back("seed");

    double std_dev;
    if (ParseValue(root, "standard-deviation", std_dev, errTxt)) {
        modelSettings_->SetStandardDeviation(std_dev);
    }

    double seed;
    if (ParseValue(root, "seed", seed, errTxt)) {
        modelSettings_->SetSeed(seed);
    }


    CheckForJunk(root, errTxt, legalCommands);

    return true;
}


bool XmlModelFile::ParseBool(TiXmlNode *node, const std::string &keyword, bool &value, std::string &errTxt, bool allowDuplicates) {
    std::string tmpVal;
    std::string tmpErr = "";
    if (ParseValue(node, keyword, tmpVal, tmpErr, allowDuplicates) == false) {
        return (false);
    }

    //Keyword is found.
    if (tmpErr == "") {
        if (tmpVal == "yes") {
            value = true;
        } else if (tmpVal == "no") {
            value = false;
        } else {
            tmpErr = "Found '" + tmpVal + "' under keyword <" + keyword + ">, expected 'yes' or 'no'. This happened in command <" +
                    node->ValueStr() + "> on line " + NRLib::ToString(node->Row()) + ", column " + NRLib::ToString(node->Column()) + ".\n";
        }
    }

    //No junk-clearing call, done in parseValue.
    errTxt += tmpErr;
    return (true);
}


void XmlModelFile::CheckForJunk(TiXmlNode *root, std::string &errTxt, const std::vector<std::string> &legalCommands, bool allowDuplicates) {
    TiXmlNode *child = root->FirstChild();
    unsigned int startLength = static_cast<unsigned int>(errTxt.size());
    while (child != NULL) {
        switch (child->Type()) {
            case TiXmlNode::COMMENT :
            case TiXmlNode::DECLARATION :
                break;
            case TiXmlNode::TEXT :
                errTxt = errTxt + "Unexpected value '" + child->Value() + "' is not part of command <" + root->Value() +
                        "> on line " + NRLib::ToString(child->Row()) + ", column " + NRLib::ToString(child->Column()) + ".\n";
                break;
            case TiXmlNode::ELEMENT :
                errTxt = errTxt + "Unexpected command <" + child->Value() + "> is not part of command <" + root->Value() +
                        "> on line " + NRLib::ToString(child->Row()) + ", column " + NRLib::ToString(child->Column()) + ".\n";
                break;
            default :
                errTxt = errTxt + "Unexpected text '" + child->Value() + "' is not part of command <" + root->Value() +
                        "> on line " + NRLib::ToString(child->Row()) + ", column " + NRLib::ToString(child->Column()) + ".\n";
                break;
        }
        root->RemoveChild(child);
        child = root->FirstChild();
    }

    if (startLength < errTxt.size()) {
        errTxt = errTxt + "Legal commands are:\n";
        for (unsigned int i = 0; i < legalCommands.size(); i++) {
            errTxt = errTxt + "  <" + legalCommands[i] + ">\n";
        }
    }

    TiXmlNode *parent = root->Parent();

    if (parent != NULL) {
        std::string cmd = root->ValueStr();
        parent->RemoveChild(root);

        if (allowDuplicates == false) {
            root = parent->FirstChildElement(cmd);
            if (root != NULL) {
                int n = 0;
                while (root != NULL) {
                    n++;
                    parent->RemoveChild(root);
                    root = parent->FirstChildElement(cmd);
                }
                errTxt += "Found " + NRLib::ToString(n) + " extra occurences of command <" + cmd + "> under command <" + parent->Value() +
                        "> on line " + NRLib::ToString(parent->Row()) + ", column " + NRLib::ToString(parent->Column()) + ".\n";
            }
        }
    }
}

