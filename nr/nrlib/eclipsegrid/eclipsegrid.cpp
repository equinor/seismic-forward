// $Id: eclipsegrid.cpp 1121 2012-11-21 08:27:28Z hgolsen $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// •  Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// •  Redistributions in binary form must reproduce the above copyright notice, this list of
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

#include "../geometry/bilinearsurface.hpp"
#include "../geometry/line.hpp"
#include "../geometry/point.hpp"
#include "../geometry/transformation.hpp"
#include "../iotools/fileio.hpp"
#include "../iotools/stringtools.hpp"

#include "eclipsecarfin.hpp"
#include "eclipsecarfingroup.hpp"
#include "eclipsegrid.hpp"
#include "eclipsegeometry.hpp"
#include "eclipsetools.hpp"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <istream>
#include <list>
#include <map>
#include <ostream>
#include <string>
#include <set>
#include <vector>

namespace NRLib {

EclipseGrid::EclipseGrid()
{}


EclipseGrid::EclipseGrid(size_t nx, size_t ny, size_t nz)
  : geometry_(nx, ny, nz)
{}


EclipseGrid::EclipseGrid(const std::string& file_name)
{
  ReadFromFile(file_name);
}


EclipseGrid::EclipseGrid(const std::vector<std::string>& file_names)
{
  ReadFromFiles(file_names);
}


EclipseGrid::~EclipseGrid()
{
  ClearRefinementGroups();

  for (size_t i = 0; i < grid_refinements_.size(); ++i) {
    delete grid_refinements_[i];
  }
  grid_refinements_.clear();
}


void EclipseGrid::Resize(size_t nx, size_t ny, size_t nz)
{
  geometry_.Resize(nx, ny, nz);
  ParameterIterator p_it = ParametersBegin();
  for ( ; p_it != ParametersEnd(); ++p_it) {
    (*p_it).second.Resize(nx, ny, nz);
  }
}

void EclipseGrid::SetGeometry(const EclipseGeometry& geometry)
{
  geometry_ = geometry;
}


void EclipseGrid::AddTransMult(const EclipseTransMult& trans_mult)
{
  trans_mult_.push_back(trans_mult);
}


void EclipseGrid::AddParameter(const std::string& parameter_name, double val)
{
  assert (!HasParameter(parameter_name));
  continuous_parameters_[parameter_name].Resize(GetNI(), GetNJ(), GetNK(), val);
}


void EclipseGrid::AddGridRefinement(const EclipseCarFin& refinement)
{
  grid_refinements_.push_back(new EclipseCarFin(refinement));
}


void EclipseGrid::AddRefinementGroup(const EclipseCarFinGroup& refinement_group)
{
  refinement_groups_.push_back(new EclipseCarFinGroup(refinement_group));
}


void EclipseGrid::ClearRefinementGroups()
{
  for (size_t i = 0; i < refinement_groups_.size(); ++i) {
    delete refinement_groups_[i];
  }
  refinement_groups_.clear();
}



void EclipseGrid::ReadContinuousParameter(ifstream& in_file, const std::string& map_key )
{
  std::string buffer = ReadParameterBuffer(in_file);
  continuous_parameters_[map_key].Resize(geometry_.GetNI(), geometry_.GetNJ(), geometry_.GetNK(), 0);
  ParseAsciiArrayFast( buffer, continuous_parameters_[map_key].begin(), geometry_.GetN());
}


void EclipseGrid::ReadFaults(std::ifstream& in_file)
{
  std::string in_line;
  std::vector<std::string> line_tokens;

  bool be_in_loop = true;
  while (be_in_loop == true){
    std::getline(in_file, in_line);
    line_tokens = GetTokens(in_line);
    if ( line_tokens.empty() ){
    }
    else if (line_tokens[0] == "--"){ // -- is a comment
    }
    else if (line_tokens[0] == "/"){
      be_in_loop = false;
    }
    else{
      if (line_tokens.size() == 9)
      {
        line_tokens[0].erase(line_tokens[0].size() - 1, 1);
        line_tokens[0].erase(0, 1);
        faults_[line_tokens[0]].SetName(line_tokens[0]);
        faults_[line_tokens[0]].ReadSegments(line_tokens);
      }
      else{
        throw Exception("Not a segment");
      }
    }
  }
}


void EclipseGrid::ReadMapUnits( std::ifstream& in_file, std::string& map_units )
{
  int line = 0; //line not used
  SkipComments(in_file);
  ReadNextQuoted(in_file, '\'', map_units, line);
  map_units = NRLib::Chomp(map_units);
  units_ = map_units;
  SkipComments(in_file);

  std::string token = ReadNext<std::string>(in_file, line);
  if (token != "/")
    throw Exception("Expected a /");

}


void EclipseGrid::ReadMapAxes( std::ifstream& in_file, double& map_ax_x1, double& map_ax_y1,
                                                       double& map_ax_x2, double& map_ax_y2,
                                                       double& map_ax_x3, double& map_ax_y3 )
{
  int line = 0; //line not used
  SkipComments(in_file);

  map_ax_x1 = ReadNext<double>( in_file, line );
  map_ax_y1 = ReadNext<double>( in_file, line );
  map_ax_x2 = ReadNext<double>( in_file, line );
  map_ax_y2 = ReadNext<double>( in_file, line );
  map_ax_x3 = ReadNext<double>( in_file, line );
  map_ax_y3 = ReadNext<double>( in_file, line );

  SkipComments(in_file);
  std::string token = ReadNext<std::string>(in_file, line);
  if (token != "/")
    throw Exception("Expected a /");
}


void EclipseGrid::ReadGridUnit( std::ifstream& in_file, std::string& grid_unit, bool& local_coord )
{
  int line = 0; //line not used
  SkipComments(in_file);
  ReadNextQuoted(in_file, '\'', grid_unit, line);
  grid_unit = NRLib::Chomp(grid_unit);
  units_ = grid_unit;

  std::string local_coordinates;
  ReadNextQuoted(in_file, '\'', local_coordinates, line);
  local_coordinates = NRLib::Chomp(local_coordinates);
  if ( local_coordinates == "" ){
    local_coord = true;
  }
  else if( local_coordinates == "MAP" ) {
    local_coord = false;
  }
  else{
    throw Exception("Wrong input");
  }

  SkipComments(in_file);

  std::string token = ReadNext<std::string>(in_file, line);
  if (token != "/")
    throw Exception("Expected a /");
}


void EclipseGrid::ReadFromFiles(const std::vector<std::string>& file_names)
{
  size_t i;
  for(i = 0; i < file_names.size(); i++){
    ReadFromFile( file_names[i] );
  }
}


void EclipseGrid::ReadFromFile(const std::string& file_name)
{
  ifstream in_file;
  OpenRead(in_file, file_name, std::ios::in | std::ios::binary);

  std::string token = "a";
  int line = 0;

  double map_ax_x1 = 0, map_ax_y1 = 1,
         map_ax_x2 = 0, map_ax_y2 = 0,
         map_ax_x3 = 1, map_ax_y3 = 0;

  bool local_coord = false;
  bool has_map_axes = false;
  std::string map_units;
  std::string grid_unit;

  while (ReadNextToken(in_file, token, line)) {
    if (token == "MAPAXES") {
      has_map_axes = true;
      ReadMapAxes(in_file, map_ax_x1, map_ax_y1, map_ax_x2, map_ax_y2, map_ax_x3, map_ax_y3);
    }
    else if (token == "MAPUNITS") {
      ReadMapUnits(in_file, map_units);
    }
    else if (token == "GRIDUNIT") {
      ReadGridUnit(in_file, grid_unit, local_coord);
    }
    else {
      ReadKeyword(in_file, token);
    }
  }

  if (local_coord == true) {
    if (has_map_axes == true) {
      NRLib::Transformation transformation;
      double vector_length = sqrt( pow((map_ax_x3 - map_ax_x2), 2) + pow( (map_ax_y3 - map_ax_y2), 2) );
      double theta = acos( (map_ax_x3 - map_ax_x2)/vector_length );

      transformation.RotateZ(theta);

      double cross_prod =(map_ax_y1 - map_ax_y2)*(map_ax_x3 - map_ax_x2) - (map_ax_x1 - map_ax_x2)*(map_ax_y3 - map_ax_y2);
      if (cross_prod <0){
        transformation.MirrorXZ();
      }
      transformation.TranslateX(map_ax_x2);
      transformation.TranslateY(map_ax_y2);

      size_t i, j;
      for (j = 0; j <= geometry_.GetNJ(); j++) {
        for (i = 0; i <= geometry_.GetNI(); i++) {
          Point p1 = geometry_.GetPillar(i,j).GetPt1();
          Point p2 = geometry_.GetPillar(i,j).GetPt2();

          transformation.Transform(p1);
          transformation.Transform(p2);
          Line line(p1, p2);
          geometry_.SetCoordLine(i, j, line);
        }
      }
    }
    else {
      throw Exception("Do not have map axes");
    }
  }
  geometry_.FindMinAndMaxZValueAndSetPolygons();
}


void EclipseGrid::ReadKeyword(std::ifstream& in_file, const std::string& token)
{
  // -- comment line
  if (token == "--" || token == "NOECHO" || token == "ECHO"){
    std::string s;
    std::getline(in_file, s);
  }
  else if (token == "ZCORN") {
    geometry_.ReadZCorn(in_file);
  }
  else if (token == "COORD") {
    geometry_.ReadCoord(in_file);
  }
  else if (token == "ACTNUM") {
    geometry_.ReadActNum(in_file);
  }
  else if (token == "SPECGRID") {
    geometry_.ReadSpecGrid(in_file);
  }
  else if (token == "FAULTS") {
    ReadFaults(in_file);
  }
  else if (token == "CARFIN") {
    EclipseCarFin* refinement = new EclipseCarFin(*this);
    bool more_refinements = refinement->ReadCarFin(in_file);
    grid_refinements_.push_back(refinement);
    while (more_refinements) {
      // Read more refinements
      refinement = new EclipseCarFin(*this);
      more_refinements = refinement->ReadCarFin(in_file);
      grid_refinements_.push_back(refinement);
    }
  }
  else if (token == "MULTIPLY") {
    cout << "Warning, reading MULTIPLY is not supported" << endl;
    string s;
    bool cont_read = true;
    while (cont_read){
      std::getline(in_file, s);
      std::vector<std::string> line_tokens;
      line_tokens = GetTokens(s);
      if ( line_tokens.empty() ){
      }
      else if ( line_tokens[0] == "/" ){
        cont_read = false;
      }
    }
  }
  else if (token == "GDORIENT") {
    cout << "Warning, reading GDORIENT is not supported" << endl;
    string s;
    bool cont_read = true;
    while (cont_read){
      std::getline(in_file, s);
      std::vector<std::string> line_tokens;
      line_tokens = GetTokens(s);
      if (!line_tokens.empty()) {
        size_t last = line_tokens.size() - 1;
        if ( line_tokens[last] == "/" )
          cont_read = false;
      }
    }
  }
  else if(token != "") {
    ReadContinuousParameter(in_file, token);
  }
}


void EclipseGrid::WriteToFile(const std::string& file_name,
                              bool write_param) const
{
  std::ofstream out_file;
  OpenWrite(out_file, file_name);

  geometry_.WriteGeometry(out_file);
  if(write_param){
    ParameterConstIterator p_it;
    for (p_it = ParametersBegin(); p_it != ParametersEnd(); ++p_it) {
      WriteContParameter(out_file, p_it->first, p_it->second);
    }
  }
  WriteFaults(out_file);

  RefinementConstIterator r_it;
  for (r_it = RefinementsBegin(); r_it != RefinementsEnd(); ++r_it) {
    (*r_it)->WriteCarFin(out_file);
  }

  WriteRefinementGroups(out_file);

  WriteEditNNC(out_file);
  WriteMultiply(out_file);
}


void EclipseGrid::WriteToFiles(const std::vector<std::string>             & file_names,
                               const std::vector<std::list<std::string> > & keywords ) const
{
  assert(file_names.size() == keywords.size());

  std::set<std::string> parameters_written;
  for (size_t j = 0; j < keywords.size(); j++) {
    if (keywords[j].size() == 0) {
      // Remainig keywords can only be written to last file.
      assert(j == keywords.size() - 1);
      WriteRemainingParameters(file_names[j], parameters_written);
    }
    else {
      std::ofstream out_file;
      OpenWrite(out_file, file_names[j]);
      std::list<std::string>::const_iterator it;
      for (it = keywords[j].begin(); it != keywords[j].end(); it++){
        WriteKeyword(out_file, *it);
        parameters_written.insert(*it);
      }
    }
  }
}


void EclipseGrid::WriteKeyword(std::ofstream& out_file, const std::string& keyword) const
{
  if (keyword == "SPECGRID" ){
    geometry_.WriteSpecGrid(out_file);
  }
  else if( keyword == "COORD" ){
    geometry_.WriteCoord(out_file);
  }
  else if( keyword == "ZCORN" ){
    geometry_.WriteZCorn(out_file);
  }
  else if( keyword == "ACTNUM" ){
    geometry_.WriteActNum(out_file);
  }
  else if (keyword == "FAULTS") {
    WriteFaults(out_file);
  }
  else if (keyword == "CARFIN") {
    RefinementConstIterator r_it;
    for (r_it = RefinementsBegin(); r_it != RefinementsEnd(); ++r_it) {
      (*r_it)->WriteCarFin(out_file);
    }
  }
  else if (keyword == "AMALGAM") {
    WriteRefinementGroups(out_file);
  }
  else {
    ParameterConstIterator param_it = continuous_parameters_.find(keyword);
    if (param_it != continuous_parameters_.end()) {
      WriteContParameter(out_file, param_it->first, param_it->second);
    }
    else{
      throw Exception("Keyword \'" + keyword + "\' does not exist in grid.");
    }
  }
}


void EclipseGrid::WriteContinousParameters(std::ofstream& out_file) const
{
  ParameterConstIterator param_it = continuous_parameters_.begin();
  for ( ; param_it != continuous_parameters_.end(); ++param_it) {
    WriteContParameter(out_file, param_it->first, param_it->second);
  }
}


void EclipseGrid::WriteRemainingParameters(const std::string           & file_name,
                                           const std::set<std::string> & parameters_written) const
{
  std::ofstream out_file;
  OpenWrite(out_file, file_name);

  if (parameters_written.find("SPECGRID") == parameters_written.end())
    geometry_.WriteSpecGrid(out_file);

  if (parameters_written.find("COORD") == parameters_written.end())
    geometry_.WriteCoord(out_file);

  if (parameters_written.find("ZCORN") == parameters_written.end())
    geometry_.WriteZCorn(out_file);

  if (parameters_written.find("ACTNUM") == parameters_written.end())
    geometry_.WriteActNum(out_file);

  if (parameters_written.find("FAULTS") == parameters_written.end())
    WriteFaults(out_file);

  ParameterConstIterator it;
  for (it = ParametersBegin() ; it != ParametersEnd(); ++it) {
    if (parameters_written.find(it->first) == parameters_written.end())
      WriteContParameter(out_file, it->first, it->second);
  }

  if (parameters_written.find("CARFIN") == parameters_written.end())
    WriteKeyword(out_file, "CARFIN");

  if (parameters_written.find("AMALGAM") == parameters_written.end())
    WriteKeyword(out_file, "AMALGAM");
}


void EclipseGrid::WriteContParameter(std::ofstream             & out_file,
                                     const std::string         & name,
                                     const NRLib::Grid<double> & parameter) const
{
  out_file.precision(5);      // this value might be changed
  out_file.setf( ios_base::fixed );

  out_file << name <<"\n";

  size_t i;
  for (i = 0; i < geometry_.GetN(); i++){
    if(i%6 == 0){
      out_file << "\n";
    }
    out_file << std::setw(15) << parameter(i);
  }
  out_file << "\n" << "  /\n\n";
}


void EclipseGrid::WriteFaults(std::ofstream& out_file) const
{
  if ( !faults_.empty() ){
    out_file << "FAULTS\n";
    out_file << "-- NAME" << std::setw(15) << "IX1" << std::setw(8) << "IX2"
                          << std::setw(10) << "IY1" << std::setw(8) << "IY2"
                          << std::setw(10) << "IZ1" << std::setw(8) << "IZ2"
                          << std::setw(10) << "FACE\n\n";

    std::map < std::string, EclipseFault >::const_iterator it;

    for (it = faults_.begin(); it != faults_.end(); it++){
      it->second.WriteFault(out_file, it->first);
    }
    out_file << "\n" << "  /\n\n";
  }
}

void EclipseGrid::WriteEditNNC(std::ofstream& out_file) const
{
  if (trans_mult_.empty())
    return;

  out_file << "EDITNNC\n";

  out_file << "--" << std::setw(4) << "IX" << std::setw(7) << "IY" << std::setw(7)
           << "IZ" << std::setw(7) << "JX" << std::setw(7) << "JY" << std::setw(7)
           << "JZ" << std::setw(10) << "TRANM" << "\n\n";

  //index +1
  std::list<EclipseTransMult>::const_iterator it;
  for (it = trans_mult_.begin(); it != trans_mult_.end(); it++){
    out_file << std::setw(5);
    out_file << it->GetCellFace().i_ + 1 << std::setw(7);
    out_file << it->GetCellFace().j_ + 1 << std::setw(7);
    out_file << it->GetCellFace().k_ + 1 << std::setw(7);

    if (it->GetCellFace().face_ == PosX){
      out_file << it->GetCellFace().i_ + 1 + 1 << std::setw(7);
      out_file << it->GetCellFace().j_ + 1 << std::setw(7);
    }
    else if (it->GetCellFace().face_ == NegX){
      out_file << it->GetCellFace().i_ - 1 + 1 << std::setw(7);
      out_file << it->GetCellFace().j_ + 1 << std::setw(7);
    }
    else if (it->GetCellFace().face_ == PosY){
      out_file << it->GetCellFace().i_ + 1 << std::setw(7);
      out_file << it->GetCellFace().j_ + 1 + 1 << std::setw(7);
    }
    else if (it->GetCellFace().face_ == NegY){
      out_file << it->GetCellFace().i_ + 1 << std::setw(7);
      out_file << it->GetCellFace().j_ - 1 + 1 << std::setw(7);
    }
    else if (it->GetCellFace().face_ == PosZ || it->GetCellFace().face_ == NegZ){
      out_file << it->GetCellFace().i_ + 1 << std::setw(7);
      out_file << it->GetCellFace().j_ + 1 << std::setw(7);
    }

    out_file << it->GetNeighbour() + 1 << "  " << it->GetMultiplier() << std::setw(5) << "/\n";
  }

  out_file << "/\n\n";

}


void EclipseGrid::WriteRefinementGroups(std::ofstream& out_file) const
{
  if (!refinement_groups_.empty()) {
    out_file << "AMALGAM\n";
    RefinementGroupConstIterator rg_it;
    for (rg_it = RefinementGroupsBegin(); rg_it != RefinementGroupsEnd(); ++rg_it) {
      (*rg_it)->WriteToFile(out_file);
    }
    out_file << "/\n\n";
  }
}


void EclipseGrid::WriteMultiply(std::ofstream& out_file) const
{
  if (!trans_mult_.empty()) {
    WriteTran(out_file, 0);
    WriteTran(out_file, 1);
    WriteTran(out_file, 2);
  }
}


void EclipseGrid::WriteTran(std::ofstream& out_file, int direction) const
{
  out_file << "MULTIPLY\n";

  std::list<EclipseTransMult>::const_iterator it;
  for (it = trans_mult_.begin(); it != trans_mult_.end(); it++){
    int i_out = it->GetCellFace().i_;
    int j_out = it->GetCellFace().j_;
    bool write = false;
    if (direction == 0){
      if (it->GetCellFace().face_ == PosX){
        out_file << "'TRANX'" << std::setw(10);
        write = true;
      }
      else if (it->GetCellFace().face_ == NegX){
        out_file << "'TRANX'" << std::setw(10);
        i_out--;
        write = true;
      }
    }
    else if (direction == 1){
      if (it->GetCellFace().face_ == PosY){
        out_file << "'TRANY'" << std::setw(10);
        write = true;
      }
      else if (it->GetCellFace().face_ == NegY){
        out_file << "'TRANY'" << std::setw(10);
        j_out--;
        write = true;
      }
    }
    else if (direction == 2){
      if (it->GetCellFace().face_ == PosZ || it->GetCellFace().face_ == NegZ){
        out_file << "'TRANZ'" << std::setw(10);
        write = true;
      }
    }
    if (write == true){
      //index +1
      out_file << it->GetMultiplier() << std::setw(5);
      out_file << i_out + 1 << std::setw(7) << i_out + 1 << std::setw(7)
               << j_out + 1 << std::setw(7) << j_out + 1 << std::setw(7)
               << it->GetCellFace().k_ + 1 << std::setw(7) << it->GetCellFace().k_ + 1
               << std::setw(5);

      out_file << "/\n";
    }

  }

  out_file << "/\n\n";
}


void EclipseGrid::GetNeighbourColumn(int i, int j, Face face, int& i_out, int& j_out) const
{
  if (face == PosX){
    i_out = i + 1;
    j_out = j;
  }
  else if (face == NegX){
    i_out = i - 1;
    j_out = j;
  }
  else if (face == PosY){
    i_out = i;
    j_out = j + 1;
  }
  else if (face == NegY){
    i_out = i;
    j_out = j - 1;
  }
  assert( i_out >= 0 && i_out < static_cast<int>(geometry_.GetNI()) );
  assert( j_out >= 0 && j_out < static_cast<int>(geometry_.GetNJ()) );
}


int EclipseGrid::GetNeighbourK(int i, int j, int k, Face face) const
{
  int k_out = k;
  bool search = true;
  if (face == PosZ){
    while(search == true){
      k_out++;
      if ( k_out >= static_cast<int>(geometry_.GetNK()) ){
        return k;
      }
      else{
        if (geometry_.IsActive(i, j, k_out) == true){
          bool same_corners;
          double err = pow(10.0, -4);
          if ( abs( geometry_.GetZCorner(i, j, k, 0, 0, 1) - geometry_.GetZCorner(i, j, k_out, 0, 0, 0) ) < err
            && abs( geometry_.GetZCorner(i, j, k, 1, 0, 1) - geometry_.GetZCorner(i, j, k_out, 1, 0, 0) ) < err
            && abs( geometry_.GetZCorner(i, j, k, 0, 1, 1) - geometry_.GetZCorner(i, j, k_out, 0, 1, 0) ) < err
            && abs( geometry_.GetZCorner(i, j, k, 1, 1, 1) - geometry_.GetZCorner(i, j, k_out, 1, 1, 0) ) < err )
          {
            same_corners = true;
          }
          else{
            same_corners = false;
          }
          if(same_corners == true){
            return k_out;
          }
          else{
            return k;
          }
        }


      }
    }// end while
  }
  else if (face == NegZ){
    while(search == true){
      k_out--;
      if ( k_out < 0 ){
        return k;
      }
      else{
        if (geometry_.IsActive(i, j, k_out) == true){
          bool same_corners;
          double err = pow(10.0, -4);
          if ( abs( geometry_.GetZCorner(i, j, k, 0, 0, 0) - geometry_.GetZCorner(i, j, k_out, 0, 0, 1) ) < err
            && abs( geometry_.GetZCorner(i, j, k, 1, 0, 0) - geometry_.GetZCorner(i, j, k_out, 1, 0, 1) ) < err
            && abs( geometry_.GetZCorner(i, j, k, 0, 1, 0) - geometry_.GetZCorner(i, j, k_out, 0, 1, 1) ) < err
            && abs( geometry_.GetZCorner(i, j, k, 1, 1, 0) - geometry_.GetZCorner(i, j, k_out, 1, 1, 1) ) < err )
          {
            same_corners = true;
          }
          else{
            same_corners = false;
          }
          if(same_corners == true){
            return k_out;
          }
          else{
            return k;
          }
        }


      }
    }// end while
  }
  else{
    throw Exception("Face must be PosZ or NegZ");
  }
  return(0);
}


std::vector<Point> EclipseGrid::GetFaceCornerPoints(size_t i_in, size_t j_in,
                                                    size_t k_in, Face face) const

{
  std::vector<Point> corners(4);

  switch(face){
  case PosX:
    corners[0] = geometry_.FindCornerPoint(i_in, j_in, k_in, 1, 0, 0);
    corners[1] = geometry_.FindCornerPoint(i_in, j_in, k_in, 1, 1, 0);
    corners[2] = geometry_.FindCornerPoint(i_in, j_in, k_in, 1, 1, 1);
    corners[3] = geometry_.FindCornerPoint(i_in, j_in, k_in, 1, 0, 1);
    break;

  case NegX:
    corners[0] = geometry_.FindCornerPoint(i_in, j_in, k_in, 0, 0, 0);
    corners[1] = geometry_.FindCornerPoint(i_in, j_in, k_in, 0, 1, 0);
    corners[2] = geometry_.FindCornerPoint(i_in, j_in, k_in, 0, 1, 1);
    corners[3] = geometry_.FindCornerPoint(i_in, j_in, k_in, 0, 0, 1);
    break;

  case PosY:
    corners[0] = geometry_.FindCornerPoint(i_in, j_in, k_in, 0, 1, 0);
    corners[1] = geometry_.FindCornerPoint(i_in, j_in, k_in, 1, 1, 0);
    corners[2] = geometry_.FindCornerPoint(i_in, j_in, k_in, 1, 1, 1);
    corners[3] = geometry_.FindCornerPoint(i_in, j_in, k_in, 0, 1, 1);
    break;

  case NegY:
    corners[0] = geometry_.FindCornerPoint(i_in, j_in, k_in, 0, 0, 0);
    corners[1] = geometry_.FindCornerPoint(i_in, j_in, k_in, 1, 0, 0);
    corners[2] = geometry_.FindCornerPoint(i_in, j_in, k_in, 1, 0, 1);
    corners[3] = geometry_.FindCornerPoint(i_in, j_in, k_in, 0, 0, 1);
    break;

  case PosZ:
    corners[0] = geometry_.FindCornerPoint(i_in, j_in, k_in, 0, 0, 1);
    corners[1] = geometry_.FindCornerPoint(i_in, j_in, k_in, 1, 0, 1);
    corners[2] = geometry_.FindCornerPoint(i_in, j_in, k_in, 1, 1, 1);
    corners[3] = geometry_.FindCornerPoint(i_in, j_in, k_in, 0, 1, 1);
    break;

  case NegZ:
    corners[0] = geometry_.FindCornerPoint(i_in, j_in, k_in, 0, 0, 0);
    corners[1] = geometry_.FindCornerPoint(i_in, j_in, k_in, 1, 0, 0);
    corners[2] = geometry_.FindCornerPoint(i_in, j_in, k_in, 1, 1, 0);
    corners[3] = geometry_.FindCornerPoint(i_in, j_in, k_in, 0, 1, 0);
  }
  return corners;
}


std::vector<Point> EclipseGrid::GetAllCornerPoints(size_t i_in,
                                                   size_t j_in,
                                                   size_t k_in) const
{
  std::vector<Point> corners(8);
  corners[0] = geometry_.FindCornerPoint(i_in, j_in, k_in, 0, 0, 0);
  corners[1] = geometry_.FindCornerPoint(i_in, j_in, k_in, 1, 0, 0);
  corners[2] = geometry_.FindCornerPoint(i_in, j_in, k_in, 1, 1, 0);
  corners[3] = geometry_.FindCornerPoint(i_in, j_in, k_in, 0, 1, 0);
  corners[4] = geometry_.FindCornerPoint(i_in, j_in, k_in, 0, 0, 1);
  corners[5] = geometry_.FindCornerPoint(i_in, j_in, k_in, 1, 0, 1);
  corners[6] = geometry_.FindCornerPoint(i_in, j_in, k_in, 1, 1, 1);
  corners[7] = geometry_.FindCornerPoint(i_in, j_in, k_in, 0, 1, 1);
  return corners;
}


void EclipseGrid::FindIntersectionAreaFaceXY(CellFace cf, int k_other_cell, double& a_x, double& a_y, double& a_z) const
{
  Point p11, p12, p13, p14,
        p21, p22, p23, p24,
        q11, q12, q13, q14,
        q21, q22, q23, q24;
  switch(cf.face_){
  case PosX:
    FindPointsPosX(p11, p12, p13, p14,
                   p21, p22, p23, p24,
                   q11, q12, q13, q14,
                   q21, q22, q23, q24, cf, k_other_cell);
    break;

  case NegX:
    FindPointsNegX(p11, p12, p13, p14,
                   p21, p22, p23, p24,
                   q11, q12, q13, q14,
                   q21, q22, q23, q24, cf, k_other_cell);
    break;

  case PosY:
    FindPointsPosY(p11, p12, p13, p14,
                   p21, p22, p23, p24,
                   q11, q12, q13, q14,
                   q21, q22, q23, q24, cf, k_other_cell);
    break;

  case NegY:
    FindPointsNegY(p11, p12, p13, p14,
                   p21, p22, p23, p24,
                   q11, q12, q13, q14,
                   q21, q22, q23, q24, cf, k_other_cell);
    break;
  case PosZ:
  case NegZ:
    assert(0); // Not implemented yet.
  }

  FindAreaZ(p11, p12, p13, p14,
            p21, p22, p23, p24,
            q11, q12, q13, q14,
            q21, q22, q23, q24, a_z);

  FindAreaX(p11, p12, p13, p14,
            p21, p22, p23, p24,
            q11, q12, q13, q14,
            q21, q22, q23, q24, a_x);

  FindAreaY(p11, p12, p13, p14,
            p21, p22, p23, p24,
            q11, q12, q13, q14,
            q21, q22, q23, q24, a_y);
}


void EclipseGrid::AddFault(const EclipseFault& fault_trace)
{
  faults_[fault_trace.GetName()] = fault_trace;
}


bool EclipseGrid::HasStairSteppedFaults() const
{
  std::map<std::string, EclipseFault>::const_iterator it = faults_.begin();
  for ( ; it != faults_.end(); ++it) {
    if (it->second.IsStairStepped())
      return true;
  }
  return false;
}


void EclipseGrid::FindStartAndEndForLines(double& kmin11, double& kmax11,
                                          double& kmin12, double& kmax12,
                                          double& kmin13, double& kmax13,
                                          double& kmin14, double& kmax14,
                                          double& kmin21, double& kmax21,
                                          double& kmin22, double& kmax22,
                                          double& kmin23, double& kmax23,
                                          double& kmin24, double& kmax24) const
{
  if (kmin11 > kmax11){
    double temp;
    temp = kmin11;
    kmin11 = kmax11;
    kmax11 = temp;
  }

  if (kmin12 > kmax12){
    double temp;
    temp = kmin12;
    kmin12 = kmax12;
    kmax12 = temp;
  }

  if (kmin13 > kmax13){
    double temp;
    temp = kmin13;
    kmin13 = kmax13;
    kmax13 = temp;
  }

  if (kmin14 > kmax14){
    double temp;
    temp = kmin14;
    kmin14 = kmax14;
    kmax14 = temp;
  }

  if (kmin21 > kmax21){
    double temp;
    temp = kmin21;
    kmin21 = kmax21;
    kmax21 = temp;
  }

  if (kmin22 > kmax22){
    double temp;
    temp = kmin22;
    kmin22 = kmax22;
    kmax22 = temp;
  }

  if (kmin23 > kmax23){
    double temp;
    temp = kmin23;
    kmin23 = kmax23;
    kmax23 = temp;
  }

  if (kmin24 > kmax24){
    double temp;
    temp = kmin24;
    kmin24 = kmax24;
    kmax24 = temp;
  }
}


void EclipseGrid::FindStartAndEndForAlgorithm(double kmin11, double kmin12,
                                              double kmin21, double kmin22,
                                              double kmax11, double kmax12,
                                              double kmax21, double kmax22,
                                              double& k_start, double& k_end) const
{
  double min_kmin1;
  if (kmin11 < kmin12){
    min_kmin1 = kmin11;
  }
  else{
    min_kmin1 = kmin12;
  }
  double min_kmin2;
  if (kmin21 < kmin22){
    min_kmin2 = kmin21;
  }
  else{
    min_kmin2 = kmin22;
  }

  if (min_kmin1 < min_kmin2){
    k_start = min_kmin2;
  }
  else{
    k_start = min_kmin1;
  }


  double max_kmax1;
  if (kmax11 > kmax12){
    max_kmax1 = kmax11;
  }
  else{
    max_kmax1 = kmax12;
  }

  double max_kmax2;
  if (kmax21 > kmax22){
    max_kmax2 = kmax21;
  }
  else{
    max_kmax2 = kmax22;
  }

  if (max_kmax1 < max_kmax2){
    k_end = max_kmax1;
  }
  else{
    k_end = max_kmax2;
  }
}


double EclipseGrid::FindIntersectYFromX(double x, const BilinearSurface& surface) const
{
  double a1 = surface.GetPoint(0,0).x + surface.GetPoint(1,1).x - surface.GetPoint(1,0).x - surface.GetPoint(0,1).x;
  double b1 = surface.GetPoint(1,0).x + surface.GetPoint(0,1).x - 2*surface.GetPoint(0,0).x;
  double c1 = surface.GetPoint(0,0).x - x;

  double a2 = surface.GetPoint(0,0).y + surface.GetPoint(1,1).y - surface.GetPoint(1,0).y - surface.GetPoint(0,1).y;
  double b2 = surface.GetPoint(1,0).y + surface.GetPoint(0,1).y - 2*surface.GetPoint(0,0).y;
  double c2 = surface.GetPoint(0,0).y;

  double u = FindUFromK(a1, b1, c1);
  double y = a2*u*u + b2*u + c2;
  return y;
}


double EclipseGrid::FindIntersectYFromZ(double z, const BilinearSurface& surface) const
{
  double a1 = surface.GetPoint(0,0).z + surface.GetPoint(1,1).z - surface.GetPoint(1,0).z - surface.GetPoint(0,1).z;
  double b1 = surface.GetPoint(1,0).z + surface.GetPoint(0,1).z - 2*surface.GetPoint(0,0).z;
  double c1 = surface.GetPoint(0,0).z - z;

  double a2 = surface.GetPoint(0,0).z + surface.GetPoint(1,1).z - surface.GetPoint(1,0).z - surface.GetPoint(0,1).z;
  double b2 = surface.GetPoint(1,0).z + surface.GetPoint(0,1).z - 2*surface.GetPoint(0,0).z;
  double c2 = surface.GetPoint(0,0).z;

  double u = FindUFromK(a1, b1, c1);
  double y = a2*u*u + b2*u + c2;
  return y;
}


double EclipseGrid::FindIntersectXFromZ(double z, const BilinearSurface& surface) const
{
  double a1 = surface.GetPoint(0,0).z + surface.GetPoint(1,1).z - surface.GetPoint(1,0).z - surface.GetPoint(0,1).z;
  double b1 = surface.GetPoint(1,0).z + surface.GetPoint(0,1).z - 2*surface.GetPoint(0,0).z;
  double c1 = surface.GetPoint(0,0).z - z;

  double a2 = surface.GetPoint(0,0).z + surface.GetPoint(1,1).z - surface.GetPoint(1,0).z - surface.GetPoint(0,1).z;
  double b2 = surface.GetPoint(1,0).z + surface.GetPoint(0,1).z - 2*surface.GetPoint(0,0).z;
  double c2 = surface.GetPoint(0,0).z;

  double u = FindUFromK(a1, b1, c1);
  double x = a2*u*u + b2*u + c2;
  return x;
}


double EclipseGrid::FindIntersectYFromXPillar(double x, const Line& line) const
{
  if ( abs( line.GetPt2().x - line.GetPt1().x ) > pow(10.0, -14) ){
    double u =( x - line.GetPt1().x ) / ( line.GetPt2().x - line.GetPt1().x );
    double y = ( line.GetPt2().y - line.GetPt1().y )*u + line.GetPt1().y;
    return y;
  }
  else{
    throw Exception("No intersection");
  }
}


double EclipseGrid::FindIntersectYFromZPillar(double z, const Line& line) const
{
  if ( abs( line.GetPt2().z - line.GetPt1().z ) > pow(10.0, -14) ){
    double u =( z - line.GetPt1().z ) / ( line.GetPt2().z - line.GetPt1().z );
    double y = ( line.GetPt2().y - line.GetPt1().y )*u + line.GetPt1().y;
    return y;
  }
  else{
    throw Exception("No intersection");
  }
}


double EclipseGrid::FindIntersectXFromZPillar(double z, const Line& line) const
{
  if ( abs( line.GetPt2().z - line.GetPt1().z ) > pow(10.0, -14) ){
    double u =( z - line.GetPt1().z ) / ( line.GetPt2().z - line.GetPt1().z );
    double x = ( line.GetPt2().x - line.GetPt1().x )*u + line.GetPt1().x;
    return x;
  }
  else{
    throw Exception("No intersection");
  }
}


double EclipseGrid::FindUFromK(double a, double b, double c) const
{
  if (a != 0) {
    double u1 = ( -b + sqrt(b*b - 4*a*c) )/(2*a);
    if (u1 < 0 || u1 > 1) {
      double u2 = ( -b - sqrt(b*b - 4*a*c) )/(2*a);
      assert(u2 >= 0 && u2 <= 1);
      return u2;
    }
    else {
      return u1;
    }
  }
  else {
    return (-c/b);
  }
}


bool EclipseGrid::FindHeight(double a_11, double a_12, double a_21, double a_22,
                             double& height) const
{
  double a1_big, a1_small;
  if (a_11 < a_12) {
    a1_big = a_12;
    a1_small = a_11;
  }
  else {
    a1_big = a_11;
    a1_small = a_12;
  }

  double a2_big, a2_small;
  if (a_21 < a_22) {
    a2_big = a_22;
    a2_small = a_21;
  }
  else {
    a2_big = a_21;
    a2_small = a_22;
  }

  if (a1_big < a2_big) {
    if (a1_big > a2_small) {
      if (a1_small > a2_small) {
        height = a1_big - a1_small;
        return true;
      }
      else {
        height = a1_big - a2_small;
        return true;
      }
    }
    else {
      return false;
    }
  }
  else {
    if (a2_big > a1_small) {
      if (a2_small > a1_small) {
        height = a2_big - a2_small;
        return true;
      }
      else {
        height = a2_big - a1_small;
        return true;
      }
    }
    else {
      return false;
    }
  }
}


void EclipseGrid::FindPointsPosX(Point& p11, Point& p12, Point& p13, Point& p14,
                                 Point& q11, Point& q12, Point& q13, Point& q14,
                                 Point& p21, Point& p22, Point& p23, Point& p24,
                                 Point& q21, Point& q22, Point& q23, Point& q24,
                                 CellFace cf, int k_other_cell) const
{
  //qab is an extra point with the same value as pab
  //in pab and qab, a is 1 or 2 for cell 1 (current cell) and 2
  //b is 1 - 4 for each of the four points for a face
  p11 = geometry_.FindCornerPoint(cf.i_, cf.j_, cf.k_, 1, 0, 0);
  p12 = geometry_.FindCornerPoint(cf.i_, cf.j_, cf.k_, 1, 1, 0);
  p13 = geometry_.FindCornerPoint(cf.i_, cf.j_, cf.k_, 1, 0, 1);
  p14 = geometry_.FindCornerPoint(cf.i_, cf.j_, cf.k_, 1, 1, 1);

  q11 = geometry_.FindPointAtPillar(cf.i_ + 1, cf.j_ + 1, p11.z);
  q12 = geometry_.FindPointAtPillar(cf.i_ + 1, cf.j_, p12.z);
  q13 = geometry_.FindPointAtPillar(cf.i_ + 1, cf.j_ + 1, p13.z);
  q14 = geometry_.FindPointAtPillar(cf.i_ + 1, cf.j_, p14.z);

  p21 = geometry_.FindCornerPoint(cf.i_ + 1, cf.j_, k_other_cell, 0, 0, 0);
  p22 = geometry_.FindCornerPoint(cf.i_ + 1, cf.j_, k_other_cell, 0, 1, 0);
  p23 = geometry_.FindCornerPoint(cf.i_ + 1, cf.j_, k_other_cell, 0, 0, 1);
  p24 = geometry_.FindCornerPoint(cf.i_ + 1, cf.j_, k_other_cell, 0, 1, 1);

  q21 = geometry_.FindPointAtPillar(cf.i_ + 1, cf.j_ + 1, p21.z);
  q22 = geometry_.FindPointAtPillar(cf.i_ + 1, cf.j_, p22.z);
  q23 = geometry_.FindPointAtPillar(cf.i_ + 1, cf.j_ + 1, p23.z);
  q24 = geometry_.FindPointAtPillar(cf.i_ + 1, cf.j_, p24.z);
}


void EclipseGrid::FindPointsNegX(Point& p11, Point& p12, Point& p13, Point& p14,
                                 Point& q11, Point& q12, Point& q13, Point& q14,
                                 Point& p21, Point& p22, Point& p23, Point& p24,
                                 Point& q21, Point& q22, Point& q23, Point& q24, CellFace cf, int k_other_cell) const
{
  //qab is an extra point with the same value as pab
  //in pab and qab, a is 1 or 2 for cell 1 (current cell) and 2
  //b is 1 - 4 for each of the four points for a face
  p11 = geometry_.FindCornerPoint(cf.i_, cf.j_, cf.k_, 0, 0, 0);
  p12 = geometry_.FindCornerPoint(cf.i_, cf.j_, cf.k_, 0, 1, 0);
  p13 = geometry_.FindCornerPoint(cf.i_, cf.j_, cf.k_, 0, 0, 1);
  p14 = geometry_.FindCornerPoint(cf.i_, cf.j_, cf.k_, 0, 1, 1);

  q11 = geometry_.FindPointAtPillar(cf.i_, cf.j_ + 1, p11.z);
  q12 = geometry_.FindPointAtPillar(cf.i_, cf.j_, p12.z);
  q13 = geometry_.FindPointAtPillar(cf.i_, cf.j_ + 1, p13.z);
  q14 = geometry_.FindPointAtPillar(cf.i_, cf.j_, p14.z);

  p21 = geometry_.FindCornerPoint(cf.i_ - 1, cf.j_, k_other_cell, 1, 0, 0);
  p22 = geometry_.FindCornerPoint(cf.i_ - 1, cf.j_, k_other_cell, 1, 1, 0);
  p23 = geometry_.FindCornerPoint(cf.i_ - 1, cf.j_, k_other_cell, 1, 0, 1);
  p24 = geometry_.FindCornerPoint(cf.i_ - 1, cf.j_, k_other_cell, 1, 1, 1);

  q21 = geometry_.FindPointAtPillar(cf.i_, cf.j_ + 1, p21.z);
  q22 = geometry_.FindPointAtPillar(cf.i_, cf.j_, p22.z);
  q23 = geometry_.FindPointAtPillar(cf.i_, cf.j_ + 1, p23.z);
  q24 = geometry_.FindPointAtPillar(cf.i_, cf.j_, p24.z);
}


void EclipseGrid::FindPointsPosY(Point& p11, Point& p12, Point& p13, Point& p14,
                                 Point& q11, Point& q12, Point& q13, Point& q14,
                                 Point& p21, Point& p22, Point& p23, Point& p24,
                                 Point& q21, Point& q22, Point& q23, Point& q24, CellFace cf, int k_other_cell) const
{
  //qab is an extra point with the same value as pab
  //in pab and qab, a is 1 or 2 for cell 1 (current cell) and 2
  //b is 1 - 4 for each of the four points for a face
  p11 = geometry_.FindCornerPoint(cf.i_, cf.j_, cf.k_, 0, 1, 0);
  p12 = geometry_.FindCornerPoint(cf.i_, cf.j_, cf.k_, 1, 1, 0);
  p13 = geometry_.FindCornerPoint(cf.i_, cf.j_, cf.k_, 0, 1, 1);
  p14 = geometry_.FindCornerPoint(cf.i_, cf.j_, cf.k_, 1, 1, 1);

  q11 = geometry_.FindPointAtPillar(cf.i_ + 1, cf.j_, p11.z);
  q12 = geometry_.FindPointAtPillar(cf.i_, cf.j_, p12.z);
  q13 = geometry_.FindPointAtPillar(cf.i_ + 1, cf.j_, p13.z);
  q14 = geometry_.FindPointAtPillar(cf.i_, cf.j_, p14.z);

  p21 = geometry_.FindCornerPoint(cf.i_, cf.j_ - 1, k_other_cell, 0, 0, 0);
  p22 = geometry_.FindCornerPoint(cf.i_, cf.j_ - 1, k_other_cell, 1, 0, 0);
  p23 = geometry_.FindCornerPoint(cf.i_, cf.j_ - 1, k_other_cell, 0, 0, 1);
  p24 = geometry_.FindCornerPoint(cf.i_, cf.j_ - 1, k_other_cell, 1, 0, 1);

  q21 = geometry_.FindPointAtPillar(cf.i_ + 1, cf.j_, p21.z);
  q22 = geometry_.FindPointAtPillar(cf.i_, cf.j_, p22.z);
  q23 = geometry_.FindPointAtPillar(cf.i_ + 1, cf.j_, p23.z);
  q24 = geometry_.FindPointAtPillar(cf.i_, cf.j_, p24.z);
}


void EclipseGrid::FindPointsNegY(Point& p11, Point& p12, Point& p13, Point& p14,
                                 Point& q11, Point& q12, Point& q13, Point& q14,
                                 Point& p21, Point& p22, Point& p23, Point& p24,
                                 Point& q21, Point& q22, Point& q23, Point& q24, CellFace cf, int k_other_cell) const
{
  //qab is an extra point with the same value as pab
  //in pab and qab, a is 1 or 2 for cell 1 (current cell) and 2
  //b is 1 - 4 for each of the four points for a face
  p11 = geometry_.FindCornerPoint(cf.i_, cf.j_, cf.k_, 0, 0, 0);
  p12 = geometry_.FindCornerPoint(cf.i_, cf.j_, cf.k_, 1, 0, 0);
  p13 = geometry_.FindCornerPoint(cf.i_, cf.j_, cf.k_, 0, 0, 1);
  p14 = geometry_.FindCornerPoint(cf.i_, cf.j_, cf.k_, 1, 0, 1);

  q11 = geometry_.FindPointAtPillar(cf.i_ + 1, cf.j_ + 1, p11.z);
  q12 = geometry_.FindPointAtPillar(cf.i_, cf.j_ + 1, p12.z);
  q13 = geometry_.FindPointAtPillar(cf.i_ + 1, cf.j_ + 1, p13.z);
  q14 = geometry_.FindPointAtPillar(cf.i_, cf.j_ + 1, p14.z);

  p21 = geometry_.FindCornerPoint(cf.i_, cf.j_ + 1, k_other_cell, 0, 1, 0);
  p22 = geometry_.FindCornerPoint(cf.i_, cf.j_ + 1, k_other_cell, 1, 1, 0);
  p23 = geometry_.FindCornerPoint(cf.i_, cf.j_ + 1, k_other_cell, 0, 1, 1);
  p24 = geometry_.FindCornerPoint(cf.i_, cf.j_ + 1, k_other_cell, 1, 1, 1);

  q21 = geometry_.FindPointAtPillar(cf.i_ + 1, cf.j_ + 1, p21.z);
  q22 = geometry_.FindPointAtPillar(cf.i_, cf.j_ + 1, p22.z);
  q23 = geometry_.FindPointAtPillar(cf.i_ + 1, cf.j_ + 1, p23.z);
  q24 = geometry_.FindPointAtPillar(cf.i_, cf.j_ + 1, p24.z);
}


void EclipseGrid::FindAreaZ(Point& p11, Point& p12, Point& p13, Point& p14,
                            Point& q11, Point& q12, Point& q13, Point& q14,
                            Point& p21, Point& p22, Point& p23, Point& p24,
                            Point& q21, Point& q22, Point& q23, Point& q24, double& a_z) const
{
  //Default values. Values may change in FindStartAndEndForLines
  double xmin11 = p11.x, xmax11 = p12.x;
  double xmin12 = p13.x, xmax12 = p14.x;
  double xmin13 = p11.x, xmax13 = p13.x;
  double xmin14 = p12.x, xmax14 = p14.x;

  double xmin21 = p21.x, xmax21 = p22.x;
  double xmin22 = p23.x, xmax22 = p24.x;
  double xmin23 = p21.x, xmax23 = p23.x;
  double xmin24 = p22.x, xmax24 = p24.x;

  FindStartAndEndForLines(xmin11, xmax11, xmin12, xmax12, xmin13, xmax13, xmin14, xmax14,
                          xmin21, xmax21, xmin22, xmax22, xmin23, xmax23, xmin24, xmax24);

  double x_start, x_end;
  FindStartAndEndForAlgorithm(xmin11, xmin12, xmin21, xmin22, xmax11, xmax12, xmax21, xmax22, x_start, x_end);
  double delta_x = (x_end - x_start) * 0.001;

  double x_cur = x_start;
  a_z = 0;

  BilinearSurface surf_11(p11, q11, q12, p12);
  BilinearSurface surf_12(p13, q13, q14, p14);
  BilinearSurface surf_21(p21, q21, q22, p22);
  BilinearSurface surf_22(p23, q23, q24, p24);
  Line line_13(p11, p13);
  Line line_14(p12, p14);
  Line line_23(p21, p23);
  Line line_24(p22, p24);

  while (x_cur < x_end){
    double y;
    std::vector<double> y_val;
    if (x_cur <= xmax11 && x_cur >= xmin11){
      y = FindIntersectYFromX(x_cur, surf_11);
      y_val.push_back(y);
    }
    if (x_cur < xmax12 && x_cur >= xmin12){
      y = FindIntersectYFromX(x_cur, surf_12);
      y_val.push_back(y);
    }
    if (x_cur < xmax13 && x_cur >= xmin13){
      y = FindIntersectYFromXPillar(x_cur, line_13);
      y_val.push_back(y);
    }
    if (x_cur < xmax14 && x_cur >= xmin14){
      y = FindIntersectYFromXPillar(x_cur, line_14);
      y_val.push_back(y);
    }


    if (x_cur < xmax21 && x_cur >= xmin21){
      y = FindIntersectYFromX(x_cur, surf_21);
      y_val.push_back(y);
    }
    if (x_cur < xmax22 && x_cur >= xmin22){
      y = FindIntersectYFromX(x_cur, surf_22);
      y_val.push_back(y);
    }
    if (x_cur < xmax23 && x_cur >= xmin23){
      y = FindIntersectYFromXPillar(x_cur, line_23);
      y_val.push_back(y);
    }
    if (x_cur < xmax24 && x_cur >= xmin24){
      y = FindIntersectYFromXPillar(x_cur, line_24);
      y_val.push_back(y);
    }


    double height;
    bool isHeight = FindHeight(y_val[0], y_val[1], y_val[2], y_val[3], height);
    if (isHeight == true){
      a_z += (height*delta_x);
    }

    x_cur+=delta_x;
  }
}


void EclipseGrid::FindAreaX(Point& p11, Point& p12, Point& p13, Point& p14,
                            Point& q11, Point& q12, Point& q13, Point& q14,
                            Point& p21, Point& p22, Point& p23, Point& p24,
                            Point& q21, Point& q22, Point& q23, Point& q24, double& a_x) const
{
  //Default values. Values may change in FindStartAndEndForLines
  double zmin11 = p11.z, zmax11 = p12.z;
  double zmin12 = p13.z, zmax12 = p14.z;
  double zmin13 = p11.z, zmax13 = p13.z;
  double zmin14 = p12.z, zmax14 = p14.z;

  double zmin21 = p21.z, zmax21 = p22.z;
  double zmin22 = p23.z, zmax22 = p24.z;
  double zmin23 = p21.z, zmax23 = p23.z;
  double zmin24 = p22.z, zmax24 = p24.z;

  FindStartAndEndForLines(zmin11, zmax11, zmin12, zmax12, zmin13, zmax13, zmin14, zmax14,
                          zmin21, zmax21, zmin22, zmax22, zmin23, zmax23, zmin24, zmax24);

  double z_start, z_end;
  FindStartAndEndForAlgorithm(zmin11, zmin12, zmin21, zmin22, zmax11, zmax12, zmax21, zmax22, z_start, z_end);
  double delta_z = (z_end - z_start) * 0.001;

  double z_cur = z_start;
  a_x = 0;

  BilinearSurface surf_11(p11, q11, q12, p12);
  BilinearSurface surf_12(p13, q13, q14, p14);
  BilinearSurface surf_21(p21, q21, q22, p22);
  BilinearSurface surf_22(p23, q23, q24, p24);
  Line line_13(p11, p13);
  Line line_14(p12, p14);
  Line line_23(p21, p23);
  Line line_24(p22, p24);

  while (z_cur < z_end){
    double y;
    std::vector<double> y_val;
    if (z_cur <= zmax11 && z_cur >= zmin11){
      y = FindIntersectYFromZ(z_cur, surf_11);
      y_val.push_back(y);
    }
    if (z_cur < zmax12 && z_cur >= zmin12){
      y = FindIntersectYFromZ(z_cur, surf_12);
      y_val.push_back(y);
    }
    if (z_cur < zmax13 && z_cur >= zmin13){
      y = FindIntersectYFromZPillar(z_cur, line_13);
      y_val.push_back(y);
    }
    if (z_cur < zmax14 && z_cur >= zmin14){
      y = FindIntersectYFromZPillar(z_cur, line_14);
      y_val.push_back(y);
    }


    if (z_cur < zmax21 && z_cur >= zmin21){
      y = FindIntersectYFromZ(z_cur, surf_21);
      y_val.push_back(y);
    }
    if (z_cur < zmax22 && z_cur >= zmin22){
      y = FindIntersectYFromZ(z_cur, surf_22);
      y_val.push_back(y);
    }
    if (z_cur < zmax23 && z_cur >= zmin23){
      y = FindIntersectYFromZPillar(z_cur, line_23);
      y_val.push_back(y);
    }
    if (z_cur < zmax24 && z_cur >= zmin24){
      y = FindIntersectYFromZPillar(z_cur, line_24);
      y_val.push_back(y);
    }


    double height;
    bool isHeight = FindHeight(y_val[0], y_val[1], y_val[2], y_val[3], height);
    if (isHeight == true){
      a_x += (height*delta_z);
    }

    z_cur+=delta_z;
  }

}


void EclipseGrid::FindAreaY(Point& p11, Point& p12, Point& p13, Point& p14,
                            Point& q11, Point& q12, Point& q13, Point& q14,
                            Point& p21, Point& p22, Point& p23, Point& p24,
                            Point& q21, Point& q22, Point& q23, Point& q24, double& a_y) const
{
  //Default values. Values may change in FindStartAndEndForLines
  double zmin11 = p11.z, zmax11 = p12.z;
  double zmin12 = p13.z, zmax12 = p14.z;
  double zmin13 = p11.z, zmax13 = p13.z;
  double zmin14 = p12.z, zmax14 = p14.z;

  double zmin21 = p21.z, zmax21 = p22.z;
  double zmin22 = p23.z, zmax22 = p24.z;
  double zmin23 = p21.z, zmax23 = p23.z;
  double zmin24 = p22.z, zmax24 = p24.z;

  FindStartAndEndForLines(zmin11, zmax11, zmin12, zmax12, zmin13, zmax13, zmin14, zmax14,
                          zmin21, zmax21, zmin22, zmax22, zmin23, zmax23, zmin24, zmax24);

  double z_start, z_end;
  FindStartAndEndForAlgorithm(zmin11, zmin12, zmin21, zmin22, zmax11, zmax12, zmax21, zmax22, z_start, z_end);
  double delta_z = (z_end - z_start) * 0.001;

  double z_cur = z_start;
  a_y = 0;


  BilinearSurface surf_11(p11, q11, q12, p12);
  BilinearSurface surf_12(p13, q13, q14, p14);
  BilinearSurface surf_21(p21, q21, q22, p22);
  BilinearSurface surf_22(p23, q23, q24, p24);
  Line line_13(p11, p13);
  Line line_14(p12, p14);
  Line line_23(p21, p23);
  Line line_24(p22, p24);

  while (z_cur < z_end){
    double x;
    std::vector<double> x_val;
    if (z_cur <= zmax11 && z_cur >= zmin11){
      x = FindIntersectXFromZ(z_cur, surf_11);
      x_val.push_back(x);
    }
    if (z_cur < zmax12 && z_cur >= zmin12){
      x = FindIntersectXFromZ(z_cur, surf_12);
      x_val.push_back(x);
    }
    if (z_cur < zmax13 && z_cur >= zmin13){
      x = FindIntersectXFromZPillar(z_cur, line_13);
      x_val.push_back(x);
    }
    if (z_cur < zmax14 && z_cur >= zmin14){
      x = FindIntersectXFromZPillar(z_cur, line_14);
      x_val.push_back(x);
    }


    if (z_cur < zmax21 && z_cur >= zmin21){
      x = FindIntersectXFromZ(z_cur, surf_21);
      x_val.push_back(x);
    }
    if (z_cur < zmax22 && z_cur >= zmin22){
      x = FindIntersectXFromZ(z_cur, surf_22);
      x_val.push_back(x);
    }
    if (z_cur < zmax23 && z_cur >= zmin23){
      x = FindIntersectXFromZPillar(z_cur, line_23);
      x_val.push_back(x);
    }
    if (z_cur < zmax24 && z_cur >= zmin24){
      x = FindIntersectXFromZPillar(z_cur, line_24);
      x_val.push_back(x);
    }


    double height;
    bool isHeight = FindHeight(x_val[0], x_val[1], x_val[2], x_val[3], height);
    if (isHeight == true){
      a_y += (height*delta_z);
    }

    z_cur+=delta_z;
  }
}


EclipseGrid EclipseGrid::ExtractSubGrid(size_t i_from, size_t i_to,
                                        size_t j_from, size_t j_to,
                                        size_t k_from, size_t k_to) const
{
  assert(i_from < GetNI() && i_to > i_from && i_to < GetNI() &&
         j_from < GetNJ() && j_to > j_from && j_to < GetNJ() &&
         k_from < GetNK() && k_to > k_from && k_to < GetNK());

  size_t ni = i_to - i_from + 1;
  size_t nj = j_to - j_from + 1;
  size_t nk = k_to - k_from + 1;

  EclipseGrid grid(ni, nj, nk);

  ParameterConstIterator it = ParametersBegin();
  for ( ; it != ParametersEnd(); ++it) {
    grid.AddParameter(it->first);
    NRLib::Grid<double>& grid_param = grid.GetParameter(it->first);
    for (size_t k = k_from; k <= k_to; ++k) {
      for (size_t j = j_from; j <= j_to; ++j) {
        for (size_t i = i_from; i <= i_to; ++i) {
          grid_param(i - i_from, j - j_from, k - k_from) = (it->second)(i, j, k);
        }
      }
    }
  }

  return grid;
}

} // namespace NRLib
