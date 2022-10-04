// $Id: xmlmodelfile.cpp 41 2014-03-28 09:42:21Z vigsnes $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification
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

#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/iotools/fileio.hpp"
#include "nrlib/iotools/logkit.hpp"

#include "nrlib/exception/exception.hpp"

#include "nrlib/segy/traceheader.hpp"

#include "modelsettings.hpp"
#include "xmlmodelfile.hpp"
#include "tasklist.hpp"

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>

XmlModelFile::XmlModelFile(const std::string  & fileName)
{
  modelSettings_ = new ModelSettings();
  failed_        = false;

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

  if (doc.Error()) {
    NRLib::LogKit::WriteHeader("Invalid XML file");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\n%s is not a valid XML file. %s In line %d, column %d.",
                                fileName.c_str(), doc.ErrorDesc(), doc.ErrorRow(), doc.ErrorCol());
    if (doc.ErrorId() == 9) {  // Not very robust check, but a start
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\nPossible cause: Mis-spelled or forgotten end tag.");
    }
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\nAborting\n");
    failed_ = true;
  }
  else {
    std::string errTxt = "";
    if (ParseSeismicForward(&doc, errTxt) == false) {
      errTxt = "'" + std::string(fileName) + "' is not a SeismicForward model file (lacks the <seismic-forward> keyword.)\n";
    }

    if (errTxt == "") {
      modelSettings_->CheckConsistency(errTxt);
      if (errTxt == "") {
        modelSettings_->SetDerivedVariables();
      }
      else {
        NRLib::LogKit::LogMessage(NRLib::LogKit::Error, errTxt);
        NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\nAborting\n");
        failed_ = true;
      }
    }
    else {
      NRLib::LogKit::WriteHeader("Invalid model file");
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\n%s is not a valid model file:\n", fileName.c_str());
      NRLib::LogKit::LogMessage(NRLib::LogKit::Error, errTxt);
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Error, "\nAborting\n");
      failed_ = true;
    }
  }
}

XmlModelFile::~XmlModelFile()
{
}

bool XmlModelFile::ParseSeismicForward(TiXmlNode *node, std::string &errTxt)
{
  TiXmlNode *root = node->FirstChildElement("seismic-forward");
  if (root == 0) {
    return (false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("angle");
  legalCommands.push_back("default-underburden");
  legalCommands.push_back("elastic-parameters");    // PAALADDED
  legalCommands.push_back("elastic-param");
  legalCommands.push_back("max-threads");           // Deprecated
  legalCommands.push_back("nmo-stretch");
  legalCommands.push_back("output-grid");
  legalCommands.push_back("output-parameters");
  legalCommands.push_back("project-settings");      // PAALADDED
  legalCommands.push_back("ps-seismic");
  legalCommands.push_back("timeshift-twt");
  legalCommands.push_back("traces-in-memory");      // Deprecated
  legalCommands.push_back("wavelet");
  legalCommands.push_back("white-noise");

  if (ParseNMOStretch(root, errTxt)){
    modelSettings_->SetNMOCorr(true);
  }
  ParseAngle(root, errTxt);
  ParseOutputGrid(root, errTxt);
  ParseElasticParam(root, errTxt);
  ParseWavelet(root, errTxt);
  ParseProjectSettings(root, errTxt);

  if (ParseWhiteNoise(root, errTxt)) {
    modelSettings_->SetAddWhiteNoise();
  }
  if (ParseReflCoefNoise(root, errTxt)) {
    modelSettings_->SetAddNoiseToReflCoef();
  }

  std::string file_name;
  if (ParseValue(root, "timeshift-twt", file_name, errTxt)) {
    modelSettings_->SetTwtFileName(file_name);
  }

  bool bolval;
  if (ParseBool(root, "ps-seismic", bolval, errTxt)) {
    modelSettings_->SetPSSeismic(bolval);
  }

  bool bolval2;
  if (ParseBool(root, "default-underburden", bolval2, errTxt)) {
    modelSettings_->SetDefaultUnderburden(bolval2);
  }

  ParseOutputParameters(root, errTxt);

  //  ------ START Moved to new section project setting ----------------

  double number;
  if (ParseValue(root, "traces-in-memory", number, errTxt)) {
    modelSettings_->SetTracesInMemory(static_cast<size_t>(number));
    TaskList::AddTask("Keyword <traces-in-memory> has been made a sub-element of section <project-settings>. Current\n    placement is deprecated.");
  }

  double n_threads;
  if (ParseValue(root, "max-threads", n_threads, errTxt)) {
    modelSettings_->SetMaxThreads(static_cast<size_t>(n_threads));
    TaskList::AddTask("Keyword <max-threads> has been made a sub-element of section <project-settings>. Current\n    placement is deprecated.");
  }

  //  ------ END Moved to new section project setting ----------------


  CheckForJunk(root, errTxt, legalCommands);
  return (true);
}

//------------------------------------------------------------
bool XmlModelFile::ParseProjectSettings(TiXmlNode   * node,
                                        std::string & errTxt)
//------------------------------------------------------------
{
  TiXmlNode *root = node->FirstChildElement("project-settings");
  if (root == 0) {
    return (false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("max-threads");
  legalCommands.push_back("traces-in-memory");


  double number;
  if (ParseValue(root, "traces-in-memory", number, errTxt)) {
    modelSettings_->SetTracesInMemory(static_cast<size_t>(number));
  }

  double n_threads;
  if (ParseValue(root, "max-threads", n_threads, errTxt)) {
    modelSettings_->SetMaxThreads(static_cast<size_t>(n_threads));
  }

  std::string level;
  if(ParseValue(root, "log-level", level, errTxt) == true) {
    int log_level = NRLib::LogKit::Error;
    if(level=="error")
      log_level = NRLib::LogKit::L_Error;
    else if(level=="warning")
      log_level = NRLib::LogKit::L_Warning;
    else if(level=="low")
      log_level = NRLib::LogKit::L_Low;
    else if(level=="medium")
      log_level = NRLib::LogKit::L_Medium;
    else if(level=="high")
      log_level = NRLib::LogKit::L_High;
    else if(level=="debuglow")
      log_level = NRLib::LogKit::L_DebugLow;
    else if(level=="debughigh")
      log_level = NRLib::LogKit::L_DebugHigh;
    else {
      errTxt += "Unknown log level " + level + " in command <log-level>. ";
      errTxt += "Choose from: error, warning, low, medium, and high\n";
      CheckForJunk(root, errTxt, legalCommands);
      return(false);
    }
    modelSettings_->SetLogLevel(log_level);
  }

  CheckForJunk(root, errTxt, legalCommands);
  return true;
}

//--------------------------------------------------------
bool XmlModelFile::ParseElasticParam(TiXmlNode   * node,
                                     std::string & errTxt)
//--------------------------------------------------------
{
  TiXmlNode *root  = node->FirstChildElement("elastic-param");
  TiXmlNode *root2 = node->FirstChildElement("elastic-paramameters");
  if (root == 0 || root2) {
    return (false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("cornerpt-interpolation-in-depth");
  legalCommands.push_back("cornerpt-interpolation-at-faults");
  legalCommands.push_back("vertical-interpolation-of-undefined-cells");
  legalCommands.push_back("default-values");
  legalCommands.push_back("eclipse-file");
  legalCommands.push_back("extra-parameters");
  legalCommands.push_back("parameter-names");
  legalCommands.push_back("remove-negative-delta-z");
  legalCommands.push_back("resampl-param-to-segy-with-interpol");
  legalCommands.push_back("zero-thickness-limit");

  std::string value;
  if (ParseValue(root, "eclipse-file", value, errTxt)) {
    modelSettings_->SetEclipseGrid(value);
  } else {
    errTxt += " An Eclipse file name must be given\n";
  }

  ParseDefaultValues(root, errTxt);
  ParseParameterNames(root, errTxt);

  double val;
  if (ParseValue(root, "zero-thickness-limit", val, errTxt)) {
    modelSettings_->SetZeroThicknessLimit(val);
  }

  bool bolval;
  if (ParseBool(root, "cornerpt-interpolation-in-depth", bolval, errTxt)) {
    modelSettings_->SetUseCornerpointInterpol(bolval);
  }

  if (ParseBool(root, "cornerpt-interpolation-at-faults", bolval, errTxt)) {
    modelSettings_->SetCornerpointInterpolationAtFaults(bolval);
  }

  if (ParseBool(root, "vertical-interpolation-of-undefined-cells", bolval, errTxt)) {
    modelSettings_->SetUseVerticalInterpolation(bolval);
  }

  if (ParseBool(root, "remove-negative-delta-z", bolval, errTxt)) {
    modelSettings_->SetRemoveNegativeDeltaZ(bolval);
  }

  if (ParseBool(root, "resampl-param-to-segy-with-interpol", bolval, errTxt)) {
    modelSettings_->SetResamplParamToSegyInterpol(bolval);
  }

  size_t i = 0;
  while (ParseExtraParameters(root, errTxt)) {
    ++i;
  }

  CheckForJunk(root, errTxt, legalCommands);
  return true;
}

bool XmlModelFile::ParseDefaultValues(TiXmlNode   * node,
                                      std::string & errTxt)
{
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
  if (ParseValue(root, "vp-top", value, errTxt)) {
    modelSettings_->SetVpTop(value);
  } else {
    errTxt += "Value for Vp above reservoir is not given.\n";
  }

  if (ParseValue(root, "vp-mid", value, errTxt)) {
    modelSettings_->SetVpMid(value);
  } else {
    errTxt += "Value for Vp in missing cells is not given.\n";
  }

  if (ParseValue(root, "vp-bot", value, errTxt)) {
    modelSettings_->SetVpBot(value);
  } else {
    errTxt += "Value for Vp below reservoir is not given\n";
  }

  if (ParseValue(root, "vs-top", value, errTxt)) {
        modelSettings_->SetVsTop(value);
  } else {
    errTxt += "Value for Vs above reservoir is not given.\n";
  }

  if (ParseValue(root, "vs-mid", value, errTxt)) {
    modelSettings_->SetVsMid(value);
  } else {
    errTxt += "Value for Vs in missing cells is not given.\n";
  }

  if (ParseValue(root, "vs-bot", value, errTxt)) {
    modelSettings_->SetVsBot(value);
  } else {
    errTxt += "Value for Vs below reservoir is not given\n";
  }

  if (ParseValue(root, "rho-top", value, errTxt)) {
    modelSettings_->SetRhoTop(value);
  } else {
    errTxt += "Value for Rho above reservoir is not given.\n";
  }

  if (ParseValue(root, "rho-mid", value, errTxt)) {
    modelSettings_->SetRhoMid(value);
  } else {
    errTxt += "Value for rho in missing cells is not given.\n";
  }

  if (ParseValue(root, "rho-bot", value, errTxt)) {
    modelSettings_->SetRhoBot(value);
  } else {
    errTxt += "Value for rho below reservoir is not given\n";
  }

  CheckForJunk(root, errTxt, legalCommands);
  return true;
}

bool XmlModelFile::ParseParameterNames(TiXmlNode   * node,
                                       std::string & errTxt)
{
  TiXmlNode *root = node->FirstChildElement("parameter-names");
  if (root == 0) {
    return (false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("vp");
  legalCommands.push_back("vs");
  legalCommands.push_back("rho");

  std::string value;
  if (ParseValue(root, "vp", value, errTxt)) {
    modelSettings_->SetVpName(value);
  } else {
    errTxt += "Name for Vp parameter in Eclipse grid is not given.\n";
  }

  if (ParseValue(root, "vs", value, errTxt)) {
    modelSettings_->SetVsName(value);
  } else {
    errTxt += "Name for Vs parameter in Eclipse grid is not given.\n";
  }
  if (ParseValue(root, "rho", value, errTxt)) {
    modelSettings_->SetRhoName(value);
  } else {
    errTxt += "Name for rho parameter in Eclipse grid is not given.\n";
  }

  CheckForJunk(root, errTxt, legalCommands);
  return true;

}

bool XmlModelFile::ParseExtraParameters(TiXmlNode   * node,
                                        std::string & errTxt)
{
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

//------------------------------------------------------
bool XmlModelFile::ParseNMOStretch(TiXmlNode   * node,
                                   std::string & errTxt)
//------------------------------------------------------
{
  TiXmlNode *root = node->FirstChildElement("nmo-stretch");
  if (root == 0) {
    return (false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("extrapol-constant");
  legalCommands.push_back("offset");
  legalCommands.push_back("offset-without-stretch");
  legalCommands.push_back("seafloor-depth");
  legalCommands.push_back("velocity-water");

  double value;
  if (ParseValue(root, "seafloor-depth", value, errTxt)) {
    modelSettings_->SetZw(value);
  } else {
    errTxt += "Value for seafloor depth is not given.\n";
  }
  if (ParseValue(root, "velocity-water", value, errTxt)) {
    modelSettings_->SetVw(value);
  } else {
    errTxt += "Value for velocity in water is not given.\n";
  }
  if (ParseValue(root, "extrapol-constant", value, errTxt)) {
    modelSettings_->SetZExtrapolFactor(value);
  }
  ParseOffset(root, errTxt);

  bool val;
  if (ParseBool(root, "offset-without-stretch", val, errTxt)) {
    modelSettings_->SetOffsetWithoutStretch(val);
  }

  CheckForJunk(root, errTxt, legalCommands, true);
  return true;
}

//--------------------------------------------------
bool XmlModelFile::ParseOffset(TiXmlNode   * node,
                               std::string & errTxt)
//--------------------------------------------------
{
  TiXmlNode *root = node->FirstChildElement("offset");
  if (root == 0) {
    return (false);
  }
  std::vector<std::string> legalCommands;
  legalCommands.push_back("offset-0");
  legalCommands.push_back("doffset");
  legalCommands.push_back("offset-max");

  double value;
  if (ParseValue(root, "offset-0", value, errTxt)) {
    modelSettings_->SetOffset0(value);
  } else {
    errTxt += "Value for minimum offset is not given.\n";
  }
  if (ParseValue(root, "doffset", value, errTxt)) {
    modelSettings_->SetDOffset(value);
  } else {
    errTxt += "Value for offset increment is not given.\n";
  }
  if (ParseValue(root, "offset-max", value, errTxt)) {
    modelSettings_->SetOffsetMax(value);
  } else {
    errTxt += "Value for maximum offset is not given.\n";
  }

  CheckForJunk(root, errTxt, legalCommands, true);
  return true;
}

//-------------------------------------------------
bool XmlModelFile::ParseAngle(TiXmlNode   * node,
                              std::string & errTxt)
//-------------------------------------------------
{
  TiXmlNode *root = node->FirstChildElement("angle");
  if (root == 0) {
    return (false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("dtheta");
  legalCommands.push_back("theta-0");
  legalCommands.push_back("theta-max");

  double value;
  if (ParseValue(root, "theta-0", value, errTxt)) {
    modelSettings_->SetTheta0(value);
  } else {
    errTxt += "Value for minimum angle is not given.\n";
  }

  if (ParseValue(root, "dtheta", value, errTxt)) {
    modelSettings_->SetDTheta(value);
  } else {
    errTxt += "Value for angle increment is not given.\n";
  }

  if (ParseValue(root, "theta-max", value, errTxt)) {
    modelSettings_->SetThetaMax(value);
  } else {
    errTxt += "Value for maximum angle is not given.\n";
  }

  CheckForJunk(root, errTxt, legalCommands);
  return true;
}

//---------------------------------------------------
bool XmlModelFile::ParseWavelet(TiXmlNode   * node,
                                std::string & errTxt)
//---------------------------------------------------
{
  TiXmlNode *root = node->FirstChildElement("wavelet");
  if (root == 0) {
    return (false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("from-file");
  legalCommands.push_back("ricker");
  legalCommands.push_back("scale");
  legalCommands.push_back("length");

  if (ParseRicker(root, errTxt)) {
    modelSettings_->SetRicker(true);
  } else {
    modelSettings_->SetRicker(false);
  }

  double value;
  if (ParseValue(root, "scale", value, errTxt)) {
    modelSettings_->SetWaveletScale(value);
  }

  if (ParseWaveletFromFile(root, errTxt) == false && modelSettings_->GetRicker() == false) {
    errTxt += "No wavelet is given. Should be given as either <ricker> or <from-file>.\n";
  }

  if (ParseValue(root, "length", value, errTxt)) {
    modelSettings_->SetWaveletLength(value);
  }

  CheckForJunk(root, errTxt, legalCommands);
  return true;
}

//--------------------------------------------------
bool XmlModelFile::ParseRicker(TiXmlNode   * node,
                               std::string & errTxt)
//--------------------------------------------------
{
  TiXmlNode *root = node->FirstChildElement("ricker");
  if (root == 0) {
    return (false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("peak-frequency");

  double value;
  if (ParseValue(root, "peak-frequency", value, errTxt)) {
    modelSettings_->SetPeakF(value);
  } else {
    errTxt += "Value for peak frequency in Ricker wavelet is not given.\n";
  }

  CheckForJunk(root, errTxt, legalCommands);
  return true;
}

//-----------------------------------------------------------
bool XmlModelFile::ParseWaveletFromFile(TiXmlNode   * node,
                                        std::string & errTxt)
//-----------------------------------------------------------
{
  TiXmlNode *root = node->FirstChildElement("from-file");
  if (root == 0) {
    return (false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("file-name");
  legalCommands.push_back("format");

  std::string format;
  std::string file_name;

  if (ParseValue(root, "format", format, errTxt))
    modelSettings_->SetWaveletFileFormat(format);
  else
    modelSettings_->SetWaveletFileFormat("LANDMARK");


  if (ParseValue(root, "file-name", file_name, errTxt))
    modelSettings_->SetWaveletFileName(file_name);
  else
    errTxt += "One or more keyword under command <" + root->ValueStr() + ">  on line " + NRLib::ToString(root->Row()) + ", column "
      + NRLib::ToString(root->Column()) + " is not legal or missing. <file-name> is required.\n";

  CheckForJunk(root, errTxt, legalCommands);
  return true;
}

//------------------------------------------------------
bool XmlModelFile::ParseOutputGrid(TiXmlNode   * node,
                                   std::string & errTxt)
//------------------------------------------------------
{
  TiXmlNode *root = node->FirstChildElement("output-grid");
  if (root == 0) {
    return (false);
  }
  std::vector<std::string> legalCommands;
  legalCommands.push_back("area");
  legalCommands.push_back("area-from-segy");
  legalCommands.push_back("area-from-surface");
  legalCommands.push_back("cell-size");
  legalCommands.push_back("depth");
  legalCommands.push_back("depth-window");
  legalCommands.push_back("segy-indexes");
  legalCommands.push_back("segy-file-format");
  legalCommands.push_back("time-window");
  legalCommands.push_back("top-time");
  legalCommands.push_back("utm-precision");

  if (ParseArea(root, errTxt)) {
    modelSettings_->SetAreaGiven(true);
  }
  std::string value;
  if (ParseValue(root, "area-from-surface", value, errTxt)) {
    modelSettings_->SetAreaFromSurface(value);
  }
  bool area_from_segy = false;
  if (ParseAreaFromSegy(root, errTxt)) {
    area_from_segy = true;
  }

  if (modelSettings_->GetAreaGiven() && modelSettings_->GetAreaFromSurface() != "" && modelSettings_->GetAreaFromSegy() != "") {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Warning, "WARNING: Area defined in three different ways. The area specified by the <area-from-segy> command is used.\n");
    TaskList::AddTask("Inconsistent XML model file specified. See beginning of log file.");
  }
  else if ((modelSettings_->GetAreaFromSurface() != "" && modelSettings_->GetAreaFromSegy() != "") ||
           (modelSettings_->GetAreaGiven()             && modelSettings_->GetAreaFromSegy() != "")) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Warning, "WARNING: Area defined in two different ways. The area specified by the <area-from-segy> command is used.\n");
    TaskList::AddTask("Inconsistent XML model file specified. See beginning of log file.");
  }
  else if (modelSettings_->GetAreaGiven() && modelSettings_->GetAreaFromSurface() != "") {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Warning, "WARNING: Area defined in two different ways. The area specified by the <area> command is used.\n");
    TaskList::AddTask("Inconsistent XML model file specified. See beginning of log file.");
  }

  if (ParseTopTime(root, errTxt, "top-time") == false) {
    if (ParseTopTime(root, errTxt, "depth")) {
      TaskList::AddTask("The command <depth> has changed its name to <top-time>. Your xml-file should be updated.");
    }
  }
  else {
    //check for "depth" anyway in order to not get it in the junk
    //The option with "depth" should be removed in a while (per 15.12.15)
    ParseDummyTopTime(root, errTxt);
  }

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
  std::string format;
  if (ParseValue(root, "segy-file-format", format, errTxt)) {
    if (format == "seisworks")
      modelSettings_->SetOutputSegyFileFormat(NRLib::TraceHeaderFormat::SEISWORKS);
    else if (format == "charisma")
      modelSettings_->SetOutputSegyFileFormat(NRLib::TraceHeaderFormat::CHARISMA);
    else if (format == "sip")
      modelSettings_->SetOutputSegyFileFormat(NRLib::TraceHeaderFormat::SIP);
    else
      errTxt += "Segy file format "+format+" is unknown. Please choose between 'seisworks', 'charisma' and 'sip'\n";
  }

  ParseTimeWindow(root, errTxt);
  ParseDepthWindow(root, errTxt);

  CheckForJunk(root, errTxt, legalCommands);
  return true;
}

bool XmlModelFile::ParseAreaFromSegy(TiXmlNode   * node,
                                     std::string & errTxt)
{
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
  legalCommands.push_back("il0-loc");
  legalCommands.push_back("xl0-loc");
  legalCommands.push_back("utmx-loc");
  legalCommands.push_back("utmy-loc");
  legalCommands.push_back("scalco-loc");
  legalCommands.push_back("start-time-loc");

  int value;

  // START - deprecated
  if (ParseValue(root, "il0", value, errTxt)) {
    modelSettings_->SetIL0Loc(value);
    TaskList::AddTask("Keyword <il0> has been made deprecated. Please use <il0-loc>");
  }
  if (ParseValue(root, "xl0", value, errTxt)) {
    modelSettings_->SetXL0Loc(value);
    TaskList::AddTask("Keyword <xl0> has been made deprecated. Please use <xl0-loc>");
  }
  if (ParseValue(root, "utmxLoc", value, errTxt)) {
    modelSettings_->SetUtmxLoc(value);
    TaskList::AddTask("Keyword <utmxLoc> has been made deprecated. Please use <utmx-loc>");
  }
  if (ParseValue(root, "utmyLoc", value, errTxt)) {
    modelSettings_->SetUtmyLoc(value);
    TaskList::AddTask("Keyword <utmyLoc> has been made deprecated. Please use <utmy-loc>");
  }
  // END - deprecated

  if (ParseValue(root, "il0-loc", value, errTxt)) {
    modelSettings_->SetIL0Loc(value);
  }
  if (ParseValue(root, "xl0-loc", value, errTxt)) {
    modelSettings_->SetXL0Loc(value);
  }
  if (ParseValue(root, "utmx-loc", value, errTxt)) {
    modelSettings_->SetUtmxLoc(value);
  }
  if (ParseValue(root, "utmy-loc", value, errTxt)) {
    modelSettings_->SetUtmyLoc(value);
  }
  if (ParseValue(root, "scalco-loc", value, errTxt)) {
    modelSettings_->SetScalcoLoc(value);
  }
  if (ParseValue(root, "start-time-loc", value, errTxt)) {
    modelSettings_->SetStartTimeLoc(value);
  }

  std::string filename;
  if (ParseValue(root, "filename", filename, errTxt)) {
    modelSettings_->SetAreaFromSegy(filename);
  }
  CheckForJunk(root, errTxt, legalCommands);
  return true;
}

bool XmlModelFile::ParseArea(TiXmlNode   * node,
                             std::string & errTxt)
{
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
  if (ParseValue(root, "x0", value, errTxt)) {
    modelSettings_->SetX0(value);
  } else {
    errTxt += "Value for x0 in area is not given.\n";
  }

  if (ParseValue(root, "y0", value, errTxt)) {
    modelSettings_->SetY0(value);
  } else {
    errTxt += "Value for y0 in area is not given.\n";
  }

  if (ParseValue(root, "lx", value, errTxt)) {
    modelSettings_->SetLx(value);
  } else {
    errTxt += "Value for lx in area is not given.\n";
  }

  if (ParseValue(root, "ly", value, errTxt)) {
    modelSettings_->SetLy(value);
  } else {
    errTxt += "Value for ly in area is not given.\n";
  }

  if (ParseValue(root, "angle", value, errTxt)) {
    modelSettings_->SetAreaAngle(value);
  } else {
    errTxt += "Value for angle in area is not given.\n";
  }

  CheckForJunk(root, errTxt, legalCommands);
  return true;
}

bool XmlModelFile::ParseTopTime(TiXmlNode   * node,
                                std::string & errTxt,
                                std::string   cname)
{
  TiXmlNode *root = node->FirstChildElement(cname);
  if (root == 0) {
    return (false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("top-time-surface");
  legalCommands.push_back("top-time-constant");

  std::string value;

  if (ParseValue(root, "top-time-surface", value, errTxt)) {
    modelSettings_->SetTopTimeSurface(value);
  }
  double val;
  if (ParseValue(root, "top-time-constant", val, errTxt)) {
    modelSettings_->SetTopTimeConstant(val);
  }

  CheckForJunk(root, errTxt, legalCommands);
  return true;
}

bool XmlModelFile::ParseDummyTopTime(TiXmlNode *node, std::string &errTxt)
{
  TiXmlNode *root = node->FirstChildElement("depth");
  if (root == 0) {
    return (false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("top-time-surface");
  legalCommands.push_back("top-time-constant");

  std::string value;
  ParseValue(root, "top-time-surface", value, errTxt);
  ParseValue(root, "top-time-constant", value, errTxt);
  CheckForJunk(root, errTxt, legalCommands);
  return true;
}

bool XmlModelFile::ParseCellSize(TiXmlNode   * node,
                                 std::string & errTxt,
                                 bool        & area_from_segy)
{
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
  if (ParseValue(root, "dx", value, errTxt)) {
    if (area_from_segy) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Warning, "WARNING: Both <dx> and <area-from-segy> is specified, and dx will be taken from segy file.\n");
      TaskList::AddTask("Inconsistent XML model file specified. See beginning of log file.");
    }
    else {
      modelSettings_->SetDx(value);
    }
  }
  // else
  //   errTxt += "Value for dx is not given.\n";

  if (ParseValue(root, "dy", value, errTxt)) {
    if (area_from_segy) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Warning, "WARNING: Both <dy> and <area-from-segy> is specified, and dy will be taken from segy file.\n");
      TaskList::AddTask("Inconsistent XML model file specified. See beginning of log file.");
   }
    else {
      modelSettings_->SetDy(value);
    }
  }
  // else
  //   errTxt += "Value for dy is not given.\n";

  if (ParseValue(root, "dz", value, errTxt)) {
    modelSettings_->SetDz(value);
  }
  // else
  //   errTxt += "Value for dz is not given.\n";

  if (ParseValue(root, "dt", value, errTxt)) {
    modelSettings_->SetDt(value);
  }
  // else
  //   errTxt += "Value for dt is not given.\n";


  CheckForJunk(root, errTxt, legalCommands);

  return true;
}

bool XmlModelFile::ParseTimeWindow(TiXmlNode   * node,
                                   std::string & errTxt)
{
  TiXmlNode *root = node->FirstChildElement("time-window");
  if (root == 0) {
    return (false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("top");
  legalCommands.push_back("bot");

  double value;
  bool top_value = false;
  if (ParseValue(root, "top", value, errTxt)) {
    modelSettings_->SetTopTimeWindow(value);
    modelSettings_->SetTimeWindowSpecified(true);
    top_value = true;
  }
  if (ParseValue(root, "bot", value, errTxt)) {
    modelSettings_->SetBotTimeWindow(value);
    modelSettings_->SetTimeWindowSpecified(true);
    if (top_value == false) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Warning, "WARNING: Value for <top> under <time-window> is not given. Is set to top reservoir.\n");
      TaskList::AddTask("<time-window> not properly specified. See beginning of log file");
    }
  }
  else if (top_value) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Warning, "WARNING: Value for <bot> under <time-window> is not given. Is set to bottom reservoir.\n");
    TaskList::AddTask("<time-window> not properly specified. See beginning of log file");
  }

  CheckForJunk(root, errTxt, legalCommands);
  return true;
}

bool XmlModelFile::ParseDepthWindow(TiXmlNode   * node,
                                    std::string & errTxt)
{
  TiXmlNode *root = node->FirstChildElement("depth-window");
  if (root == 0) {
    return (false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("top");
  legalCommands.push_back("bot");

  double value;
  bool top_value = false;
  if (ParseValue(root, "top", value, errTxt)) {
    modelSettings_->SetTopDepthWindow(value);
    modelSettings_->SetDepthWindowSpecified(true);
    top_value = true;
  }
  if (ParseValue(root, "bot", value, errTxt)) {
    modelSettings_->SetBotDepthWindow(value);
    modelSettings_->SetDepthWindowSpecified(true);
    if (top_value == false) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Warning, "WARNING: Value for <top> under <depth-window> is not given. Is set to top reservoir.\n");
      TaskList::AddTask("<depth-window> not properly specified. See beginning of log file");
    }
  }
  else if (top_value) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Warning, "WARNING: Value for <bot> under <depth-window> is not given. Is set to bottom reservoir.\n");
    TaskList::AddTask("<depth-window> not properly specified. See beginning of log file");
  }
  CheckForJunk(root, errTxt, legalCommands);
  return true;
}

bool XmlModelFile::ParseOutputParameters(TiXmlNode   * node,
                                         std::string & errTxt)
{
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
  legalCommands.push_back("wavelet");
  legalCommands.push_back("twt");
  legalCommands.push_back("vrms");
  legalCommands.push_back("twt-offset-segy");
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
  legalCommands.push_back("seismic-time-prenmo-segy");

  bool value;
  if (ParseBool(root, "elastic-parameters", value, errTxt)) {
    modelSettings_->SetOutputVp(value);
  }

  if (ParseBool(root, "zvalues", value, errTxt)) {
    modelSettings_->SetOutputZvalues(value);
  }

  if (ParseBool(root, "reflections", value, errTxt)) {
    modelSettings_->SetOutputReflections(value);
  }

  if (ParseBool(root, "seismic-time", value, errTxt)) {
        modelSettings_->SetOutputSeismicTime(value);
  }

  if (ParseBool(root, "seismic-timeshift", value, errTxt)) {
    if (value && modelSettings_->GetTwtFileName() == "") {
      errTxt += "The command <twt-timeshift> must be set if output parameter <seismic-timeshift> is given.\n";
    } else {
      modelSettings_->SetOutputSeismicTimeshift(value);
    }
  }

  if (ParseBool(root, "seismic-depth", value, errTxt)) {
    modelSettings_->SetOutputSeismicDepth(value);
  }

  if (ParseBool(root, "time-surfaces", value, errTxt)) {
    modelSettings_->SetOutputTimeSurfaces(value);
  }

  if (ParseBool(root, "depth-surfaces", value, errTxt)) {
    modelSettings_->SetOutputDepthSurfaces(value);
  }

  if (ParseBool(root, "wavelet", value, errTxt)) {
    modelSettings_->SetOutputWavelet(value);
  }

  if (ParseBool(root, "twt", value, errTxt)) {
    modelSettings_->SetOutputTwt(value);
  }

  if (ParseBool(root, "vrms", value, errTxt)) {
    modelSettings_->SetOutputVrms(value);
  }

  if (ParseBool(root, "twt-offset-segy", value, errTxt)) {
    modelSettings_->SetOutputTwtOffset(value);
  }

  std::string val;
  if (ParseValue(root, "prefix", val, errTxt)) {
    modelSettings_->SetPrefix(val);
  }

  if (ParseValue(root, "suffix", val, errTxt)) {
    modelSettings_->SetSuffix(val);
  }

  if (ParseBool(root, "seismic-time-segy", value, errTxt)) {
    modelSettings_->SetOutputTimeSegy(value);
  }

  if (ParseBool(root, "seismic-timeshift-segy", value, errTxt)) {
    if (value && modelSettings_->GetTwtFileName() == "") {
      errTxt += "The command <twt-timeshift> must be set if output parameter <seismic-timeshift-segy> is given.\n";
    } else {
      modelSettings_->SetOutputTimeshiftSegy(value);
    }
  }
  if (ParseBool(root, "seismic-depth-segy", value, errTxt)) {
    modelSettings_->SetOutputDepthSegy(value);
  }

  if (ParseBool(root, "elastic-parameters-time-segy", value, errTxt)) {
    modelSettings_->SetOutputElasticParametersTimeSegy(value);
  }

  if (ParseBool(root, "elastic-parameters-depth-segy", value, errTxt)) {
    modelSettings_->SetOutputElasticParametersDepthSegy(value);
  }


  if (ParseBool(root, "extra-parameters-time-segy", value, errTxt)) {
    modelSettings_->SetOutputExtraParametersTimeSegy(value);
  }

  if (ParseBool(root, "extra-parameters-depth-segy", value, errTxt)) {
    modelSettings_->SetOutputExtraParametersDepthSegy(value);
  }

  ParseSeismicStack(root, errTxt);

  if (ParseBool(root, "seismic-time-prenmo-segy", value, errTxt)) {
    modelSettings_->SetOutputPrenmoTimeSegy(value);
    if (modelSettings_->GetOffsetWithoutStretch() && value) {
      modelSettings_->SetOutputPrenmoTimeSegy(false);
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Warning, "WARNING: <seismic-time-prenmo-segy> under <output-parameters> is specified. For seismic specified with <nmo-stretch> and <offset-withouth-stretch>, pre-nmo seismic will not be given.\n");
    }
  }

  CheckForJunk(root, errTxt, legalCommands);

  return true;
}

bool XmlModelFile::ParseSegyIndexes(TiXmlNode   * node,
                                    std::string & errTxt,
                                    bool          area_from_segy)
{
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
  else {
    if (area_from_segy == false)
      errTxt += "One or more keyword under command <" + root->ValueStr() + ">  on line " + NRLib::ToString(root->Row()) + ", column "
                + NRLib::ToString(root->Column()) + " is not legal or missing. Three keywords; <inline-start>, <xline-start> and <inline-direction> are required.\n";
  }
  int inline_step, xline_step;
  if (ParseValue(root, "inline-step", inline_step, errTxt) && area_from_segy == false) {
    modelSettings_->SetSegyInlineStep(inline_step);
  }
  if (ParseValue(root, "xline-step", xline_step, errTxt) && area_from_segy == false) {
    modelSettings_->SetSegyXlineStep(xline_step);
  }

  if (area_from_segy)  {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Warning, "WARNING: Both <segy-indexes> <area-from-segy> and are given. Indexes are taken from segy file, and values given in <segy-indexes> are not used.\n");
    TaskList::AddTask("Inconsistent model file specified. See beginning of log file.");
  }

  CheckForJunk(root, errTxt, legalCommands);

  return true;
}

bool XmlModelFile::ParseSeismicStack(TiXmlNode   * node,
                                     std::string & errTxt)
{
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
  if (ParseBool(root, "time-storm", val, errTxt)) {
    modelSettings_->SetOutputSeismicStackTimeStorm(val);
  }

  if (ParseBool(root, "timeshift-storm", val, errTxt)) {
    if (modelSettings_->GetTwtFileName() == "" && val) {
      errTxt += "The command <twt-timeshift> must be set if output parameter <seismic-stack><timeshift-storm> is given.\n";
      }
    else {
      modelSettings_->SetOutputSeismicStackTimeShiftStorm(val);
    }
  }

  if (ParseBool(root, "depth-storm", val, errTxt)) {
    modelSettings_->SetOutputSeismicStackDepthStorm(val);
  }

  if (ParseBool(root, "time-segy", val, errTxt)) {
    modelSettings_->SetOutputSeismicStackTimeSegy(val);
  }

  if (ParseBool(root, "timeshift-segy", val, errTxt)) {
    if (modelSettings_->GetTwtFileName() == "" && val) {
      errTxt += "The command <twt-timeshift> must be set if output parameter <seismic-stack><timeshift-segy> is given.\n";
    }
    else {
      modelSettings_->SetOutputSeismicStackTimeShiftSegy(val);
    }
  }

  if (ParseBool(root, "depth-segy", val, errTxt)) {
    modelSettings_->SetOutputSeismicStackDepthSegy(val);
  }

  CheckForJunk(root, errTxt, legalCommands);

  return true;
}

//-------------------------------------------------------
bool XmlModelFile::ParseWhiteNoise(TiXmlNode   * node,
                                   std::string & errTxt)
//-------------------------------------------------------
{
  TiXmlNode *root = node->FirstChildElement("white-noise");
  if (root == 0) {
    return (false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("standard-deviation");
  legalCommands.push_back("seed");

  double std_dev;
  if (ParseValue(root, "standard-deviation", std_dev, errTxt)) {
    modelSettings_->SetStandardDeviation1(std_dev);
  }

  double seed;
  if (ParseValue(root, "seed", seed, errTxt)) {
    modelSettings_->SetSeed1(seed);
  }

  CheckForJunk(root, errTxt, legalCommands);

  return true;
}

//---------------------------------------------------------
bool XmlModelFile::ParseReflCoefNoise(TiXmlNode   * node,
                                      std::string & errTxt)
//---------------------------------------------------------
{
  TiXmlNode *root = node->FirstChildElement("add-noise-to-refl-coef");
  if (root == 0) {
    return (false);
  }

  std::vector<std::string> legalCommands;
  legalCommands.push_back("standard-deviation");
  legalCommands.push_back("seed");

  double std_dev;
  if (ParseValue(root, "standard-deviation", std_dev, errTxt)) {
    modelSettings_->SetStandardDeviation2(std_dev);
  }

  double seed;
  if (ParseValue(root, "seed", seed, errTxt)) {
    modelSettings_->SetSeed2(seed);
  }

  CheckForJunk(root, errTxt, legalCommands);

  return true;
}

//---------------------------------------------------------------
bool XmlModelFile::ParseBool(TiXmlNode         * node,
                             const std::string & keyword,
                             bool              & value,
                             std::string       & errTxt,
                             bool                allowDuplicates)
//---------------------------------------------------------------
{
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


void XmlModelFile::CheckForJunk(TiXmlNode                      * root,
                                std::string                    & errTxt,
                                const std::vector<std::string> & legalCommands,
                                bool                             allowDuplicates)
{
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
