// $Id: xmlmodelfile.hpp 39 2014-02-06 08:53:13Z vigsnes $

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

#ifndef XMLMODELFILE_HPP
#define XMLMODELFILE_HPP

#include <nrlib/iotools/stringtools.hpp>
#include <nrlib/tinyxml/tinyxml.h>

#include <stdio.h>

class ModelSettings;

class XmlModelFile
{
public:

  XmlModelFile(const std::string & fileName);

    ~XmlModelFile(void);

    ModelSettings * getModelSettings(void) const { return modelSettings_ ;}
    bool            getParsingFailed(void) const { return failed_        ;}

private:

  bool ParseSeismicForward(TiXmlNode   * node,
                           std::string & errTxt);

  bool ParseAngle(TiXmlNode   * node,
                  std::string & errTxt);

  bool ParseProjectSettings(TiXmlNode   * node,
                            std::string & errTxt);

  bool ParseNMOStretch(TiXmlNode *node,
                       std::string &errTxt);
  bool    ParseOffset(TiXmlNode   * node,
                      std::string & errTxt);

  bool ParseElasticParam(TiXmlNode   * node,
                         std::string & errTxt);
  bool    ParseDefaultValues(TiXmlNode   * node,
                             std::string & errTxt);

  bool ParseOutputGrid(TiXmlNode   * node,
                       std::string & errTxt);
  bool    ParseArea(TiXmlNode   * node,
                    std::string & errTxt);
  bool    ParseAreaFromSegy(TiXmlNode   * node,
                            std::string & errTxt);
  bool    ParseTopTime(TiXmlNode   * node,
                      std::string & errTxt,
                     std::string   cname);
  bool    ParseTimeWindow(TiXmlNode   * node,
                          std::string & errTxt);
  bool    ParseDepthWindow(TiXmlNode   * node,
                           std::string & errTxt);
  bool    ParseCellSize(TiXmlNode   * node,
                        std::string & errTxt,
                        bool        & area_from_segy);

  bool ParseWavelet(TiXmlNode   * node,
                    std::string & errTxt);
  bool    ParseRicker(TiXmlNode   * node,
                      std::string & errTxt);
  bool    ParseWaveletFromFile(TiXmlNode   * node,
                               std::string & errTxt);

  bool ParseOutputParameters(TiXmlNode   * node,
                             std::string & errTxt);

  bool ParseSegyIndexes(TiXmlNode   * node,
                        std::string & errTxt,
                        bool          area_from_segy);

  bool ParseWhiteNoise(TiXmlNode   * node,
                       std::string & errTxt);

  bool ParseZValueExtrapolation(TiXmlNode   * node,
                                std::string & errTxt);

  bool ParseSeismicStack(TiXmlNode   * node,
                         std::string & errTxt);

  bool ParseDummyTopTime(TiXmlNode   * node,
                         std::string & errTxt);

  template<typename T>
  bool ParseValue(TiXmlNode         * node,
                  const std::string & keyword,
                  T                 & value,
                  std::string       & errTxt,
                  bool                allowDuplicates = false);

  bool ParseBool(TiXmlNode         * node,
                 const std::string & keyword,
                 bool              & value,
                 std::string       & errTxt,
                 bool                allowDuplicates = false);

  bool ParseParameterNames(TiXmlNode   * node,
                           std::string & errTxt);

  bool ParseExtraParameters(TiXmlNode   * node,
                            std::string & errTxt);

  void CheckForJunk(TiXmlNode                      * root,
                    std::string                    & errTxt,
                    const std::vector<std::string> & legalCommands,
                    bool                             allowDuplicates = false);

  void SetMissing(int    & value)      { value = -99999    ;}
  void SetMissing(float  & value)      { value = -99999.0f ;}
  void SetMissing(double & value)      { value = -99999.0  ;}
  void SetMissing(std::string & value) { value = ""        ;}

  ModelSettings * modelSettings_;
  bool            failed_;
};

template<typename T>
bool XmlModelFile::ParseValue(TiXmlNode         * node,
                              const std::string & keyword,
                              T                 & value,
                              std::string       & errTxt,
                              bool                allowDuplicates)
{
  SetMissing(value);
  TiXmlNode *root = node->FirstChild(keyword);
  if (root == NULL) {
    return (false);
  }

  std::vector<std::string> legalCommands(1);

  TiXmlNode *valNode = root->FirstChild();
  while (valNode != NULL && valNode->Type() != TiXmlNode::TEXT) {
    valNode = valNode->NextSibling();
  }

  std::string tmpErr = "";
  if (valNode == NULL)
    tmpErr = "Error: No value found under keyword '" + keyword +
      " on line " + NRLib::ToString(root->Row()) + ", column " + NRLib::ToString(root->Column()) + ".\n";
  else {
    try {
      value = NRLib::ParseType<T>(valNode->ValueStr());
    } catch (NRLib::Exception &e) {
      SetMissing(value);
      tmpErr = "Error: " + std::string(e.what()) + " on line " + NRLib::ToString(valNode->Row()) + ", column " + NRLib::ToString(valNode->Column()) + ".\n";
    }
    root->RemoveChild(valNode);
  }

  CheckForJunk(root, tmpErr, legalCommands, allowDuplicates);

  errTxt += tmpErr;
  return (true);
}


#endif
