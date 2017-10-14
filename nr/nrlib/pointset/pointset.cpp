// $Id: pointset.cpp 1188 2013-06-26 07:33:27Z perroe $

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

#include <cassert>
#include <iostream>
#include <fstream>

#include "pointset.hpp"
#include "../exception/exception.hpp"
#include "../geometry/polygon.hpp"
#include "../iotools/fileio.hpp"

using namespace NRLib;

PointSet::PointSet()
{}


PointSet::PointSet(const std::vector<NRLib::Point> & points)
  : points_(points)
{}


PointSet::PointSet(size_t n, const Point & init_val)
  : points_(n, init_val)
{}


PointSet::PointSet(const std::string & filename,
                   bool                filter_inactive_horizon_points)
{
  PointSet::FileFormat format = FindFileType(filename);

  if (format == PointSet::RoxarText)
    ReadRoxarText(filename);
  else if (format == PointSet::RmsInternalPoints)
    ReadRMSInternalPoints(filename, filter_inactive_horizon_points);
  else if (format == PointSet::JasonXYZ)
    ReadJasonXYZ(filename);
  else if (format == PointSet::GeneralXYZ)
    ReadGeneralXYZ(filename);
  else
    throw FileFormatError(filename  + " is not a valid point-set file.");
}


PointSet::PointSet(const NRLib::RegularSurface<double> & storm_surface)
{
  for (size_t j = 0; j < storm_surface.GetNJ(); j++) {
    for (size_t i = 0; i < storm_surface.GetNI(); i++) {
      double x, y;
      storm_surface.GetXY(i, j, x, y);
      double z = storm_surface(i,j);
      if (!storm_surface.IsMissing(z))
        points_.push_back(Point(x,y,z));
    }
  }
}


PointSet::FileFormat
PointSet::FindFileType(const std::string & filename)
{
  if (CheckFileIsInternalPoints(filename))
    return RmsInternalPoints;
  else if (CheckFileIsRoxarText(filename))
    return RoxarText;
  else if (CheckFileIsJasonXYZ(filename))
    return JasonXYZ;
  else if (CheckFileIsGeneralXYZ(filename))
    return GeneralXYZ;
  else if (CheckFileIsOpenWorksXYZ(filename))
    return OpenWorksXYZ;
  else
    return PointSet::UnknownFormat;
}


void
PointSet::AddPoint(const NRLib::Point& point)
{
  points_.push_back(point);

  std::map<std::string, std::vector<int> >::iterator iter;
  for(iter = integer_tag_.begin(); iter != integer_tag_.end(); ++iter)
    iter->second.push_back(0);

  std::map<std::string, std::vector<double> >::iterator iter2;
  for(iter2 = float_tag_.begin(); iter2 != float_tag_.end(); ++iter2)
    iter2->second.push_back(0.0);

  std::map<std::string, std::vector<std::string> >::iterator iter3;
  for(iter3 = string_tag_.begin(); iter3 != string_tag_.end(); ++iter3)
    iter3->second.push_back("");
}


void
PointSet::AddDiscreteParameter(const std::string& name, std::vector<int> &values)
{
  assert(values.size() == GetSize());
  integer_tag_[name] = values;
}


void
PointSet::AddContParameter(const std::string& name, std::vector<double> &values)
{
  assert(values.size() == GetSize());
  float_tag_[name] = values;
}


void
PointSet::AddStringParameter(const std::string& name, std::vector<std::string> &values)
{
  assert(values.size() == GetSize());
  string_tag_[name] = values;
}


void
PointSet::WriteToFile(const std::string& filename, FileFormat format) const
{
  if (format == RmsInternalPoints)
    WriteInternalPoints(filename);
  else if (format == RoxarText)
    WriteRoxarText(filename);
  else if (format ==  GeneralXYZ)
    WriteGeneralXYZ(filename);
  else
    assert(0); // Unknown file format.
}



bool PointSet::CheckFileIsRoxarText(const std::string& filename)
{
  std::ifstream in_file;
  OpenRead(in_file, filename);

  std::string first_line;
  std::getline(in_file, first_line);

  std::vector<std::string> tokens = GetTokens(first_line);
  if (   tokens.size() == 3 && IsType<double>(tokens[0])
      && IsType<double>(tokens[1]) && IsType<double>(tokens[2])) {

    // Check that last line is -999 -999 -999
    std::string last_nonempty_line = FindLastNonEmptyLine(in_file, 100);

    tokens = GetTokens(last_nonempty_line);
    if (tokens.size() == 3) {
      std::vector<double> xyz(3);
      NRLib::ParseAsciiArrayFast(last_nonempty_line, xyz.begin(), 3);
      if (IsRoxarFileEnder(xyz))
        return true;
    }
  }

  return false;
}


bool PointSet::CheckFileIsGeneralXYZ(const std::string& filename)
{
  std::ifstream in_file;
  OpenRead(in_file, filename);

  std::string first_line;
  std::getline(in_file, first_line);

  std::vector<std::string> tokens = GetTokens(first_line);
  if (   tokens.size() == 3 && IsType<double>(tokens[0])
      && IsType<double>(tokens[1]) && IsType<double>(tokens[2])) {
    return true;
  }

  return false;
}


bool PointSet::CheckFileIsInternalPoints(const std::string& filename)
{
  std::ifstream in_file;
  OpenRead(in_file, filename);

  std::string first_token;
  if (!(in_file >> first_token)) {
    // Empty file.
    return false;
  }

  if (first_token == "Discrete" || first_token == "Float" || first_token == "String")
    return true;

  return false;
}


bool PointSet::CheckFileIsJasonXYZ(const std::string& filename)
{
  std::ifstream in_file;
  OpenRead(in_file, filename);

  int line_number = 0;
  SkipComments(in_file, '*', line_number);

  std::string token;
  in_file >> token;

  if (token == "JGWFILE")
    return true;

  return false;
}


bool PointSet::CheckFileIsOpenWorksXYZ(const std::string& filename)
{
  std::ifstream in_file;
  OpenRead(in_file, filename);

  int line_number = 0;
  SkipComments(in_file, '!', line_number);

  std::string token;
  in_file >> token;

  if (token == "New")
    return true;

  return false;
}


void
PointSet::ReadRoxarText(const std::string &filename)
{
  std::ifstream in_file;
  OpenRead(in_file, filename);
  std::vector<double> xyz(3);
  std::string temp;

  while (!NRLib::CheckEndOfFile(in_file)) {
    std::getline(in_file, temp);
    NRLib::ParseAsciiArrayFast(temp, xyz.begin(), 3);
    while (!IsRoxarFileEnder(xyz)) {
      NRLib::Point pt(xyz[0], xyz[1], xyz[2]);
      points_.push_back(pt);
      std::getline(in_file, temp);
      NRLib::ParseAsciiArrayFast(temp, xyz.begin(), 3);
    }

    if (in_file.eof() && points_.size() > 0) {
      throw NRLib::Exception("Format error reading Roxar Text file " + filename);
    }
  }
}


void
PointSet::ReadRMSInternalPoints(const std::string & filename,
                                bool                filter_inactive_horizon_points)
{
  bool filter_inactive_points = false;
  size_t active_index(999); //Initialized to dummy value to prevent compilation warning on Linux
  std::vector<std::string> type;
  std::vector<std::string> names;
  std::ifstream file;
  NRLib::OpenRead(file, filename);
  std::vector<std::string> help;
  std::string line;

  std::getline(file, line);
  size_t first_white_pos = line.find_first_of(NRLib::Whitespace());
  std::string type_test = line.substr(0, first_white_pos);

  while (!IsType<double>(type_test)){
    type.push_back(type_test);
    size_t name_pos = line.find_first_not_of(NRLib::Whitespace(), first_white_pos);
    size_t end_pos_string = line.find_last_not_of(NRLib::Whitespace());
    std::string name_test = line.substr(name_pos, (end_pos_string - name_pos + 1));
    names.push_back(name_test);
    if (name_test == "Active" && filter_inactive_horizon_points){
      filter_inactive_points = true;
      active_index = type.size();
    }
    if (NRLib::CheckEndOfFile(file))
      break;
    std::getline(file, line);
    first_white_pos = line.find_first_of(NRLib::Whitespace());
    type_test = line.substr(0, first_white_pos);
  }
  size_t n_tags = type.size();
  std::vector<std::vector<int> > discrete(n_tags);
  std::vector<std::vector<double> > floats(n_tags);
  std::vector<std::vector<std::string> > strings(n_tags);
  std::vector<double> xyz(3);
  double infloat;
  size_t n_token;

  bool first_point = true;
  bool active;
  while(!NRLib::CheckEndOfFile(file)) {
    active = false;
    if (!first_point)
      std::getline(file, line);
    first_point = false;
    help = NRLib::GetQuotedTokens(line);
    n_token = help.size();
    for(size_t i = n_token; i < 3 + n_tags; i++)
      help.push_back("UNDEF");

    if (filter_inactive_points){
      if (ParseType<int>(help[active_index - 1 + 3]) == 1)
        active = true;
    }

    if (!filter_inactive_points || (filter_inactive_points && active)){
      xyz[0] = ParseType<double>(help[0]);
      xyz[1] = ParseType<double>(help[1]);
      xyz[2] = ParseType<double>(help[2]);
      // NRLib::ReadAsciiArrayFast(file, xyz.begin(), 3);
      NRLib::Point pt(xyz[0], xyz[1], xyz[2]);
      points_.push_back(pt);
      for(size_t i = 0; i < n_tags; i++){
        if(type[i] == "Discrete"){
          if(IsType<int>(help[i + 3]))
            discrete[i].push_back(ParseType<int>(help[i + 3]));
          else if(IsType<std::string>(help[i + 3])){
            if(help[i + 3] == "UNDEF")
              discrete[i].push_back(-999);
            else
              throw NRLib::Exception("Error in Discrete tag in RMSInternalPoints file " + filename);
          }
        }
        else if(type[i] == "Float"){
          if(IsType<double>(help[i + 3]))
            infloat = ParseType<double>(help[i + 3]);
          else { //if(IsType<std::string>(help[i + 3]))
            if(help[i + 3] == "UNDEF")
              infloat = -999.0;
            else
              throw NRLib::Exception("Error in Float tag in RMSInternalPoints file " + filename);
          }
          floats[i].push_back(infloat);
        }
        else if(type[i] == "String"){
          strings[i].push_back(help[i + 3]);
        }
        else
          throw NRLib::Exception("Wrong tag specification in RMSInternalPoints file " + filename);
      }
    }
  }
  for(size_t i = 0; i< n_tags; i++){
    if(type[i] == "Discrete")
      AddDiscreteParameter(names[i], discrete[i]);
    else if(type[i] == "Float")
      AddContParameter(names[i], floats[i]);
    else if(type[i] == "String")
      AddStringParameter(names[i], strings[i]);
  }
}


void PointSet::ReadJasonXYZ(const std::string &filename)
{
  std::ifstream in_file;
  OpenRead(in_file, filename);

  int line_number = 0;
  SkipComments(in_file, '*', line_number);

  std::string line;
  line_number++;
  if (!std::getline(in_file, line)) {
    throw FileFormatError("End of file not expected.");
  }
  std::vector<std::string> tokens = GetTokens(line);

  if (tokens[0] != "JGWFILE") {
    throw FileFormatError("File is not a valid Jason horizon file.");
  }

  size_t n_traces = NRLib::ParseType<size_t>(tokens[tokens.size() - 1]);

  if (!std::getline(in_file, line)) {
    throw FileFormatError("End of file not expected.");
  }
  line_number++;

  Resize(n_traces);

  std::vector<int> il(n_traces);
  std::vector<int> xl(n_traces);

  // On each line: IL, XL, X, Y, Z and optionally survey name.
  for (size_t i = 0; i < n_traces; ++i) {
    std::getline(in_file, line);
    tokens = GetTokens(line);
    if (tokens.size() < 5 || tokens.size() > 6) {
      throw FileFormatError("Error in file, wrong number of tokens on line "
                            + NRLib::ToString(line_number));
    }
    il[i] = NRLibPrivate::UnsafeParser<int>::ParseType(tokens[0].c_str());
    xl[i] = NRLibPrivate::UnsafeParser<int>::ParseType(tokens[1].c_str());
    points_[i].x = NRLibPrivate::UnsafeParser<int>::ParseType(tokens[2].c_str());
    points_[i].y = NRLibPrivate::UnsafeParser<int>::ParseType(tokens[3].c_str());
    points_[i].z = NRLibPrivate::UnsafeParser<int>::ParseType(tokens[4].c_str());
  }

  AddDiscreteParameter("IL", il);
  AddDiscreteParameter("XL", xl);
}

void PointSet::ReadOpenWorksXYZ(const std::string& filename)
{
  std::ifstream in_file;
  OpenRead(in_file, filename);

  int line_number = 0;
  SkipComments(in_file, '!', line_number);

  std::string line;
  line_number++;
  for (size_t i=0; i<8; i++) {
    if (!std::getline(in_file, line)) {
      throw FileFormatError("End of file not expected.");
    }
    line_number++;
  }

  std::vector<double> xyz(3);
  while (!NRLib::CheckEndOfFile(in_file)) {
    std::getline(in_file, line);
    NRLib::ParseAsciiArrayFast(line, xyz.begin(), 3);
    NRLib::Point pt(xyz[0], xyz[1], xyz[2]);
    points_.push_back(pt);
    line_number++;
  }
}

void PointSet::ReadGeneralXYZ(const std::string& filename)
{
  std::ifstream in_file;
  OpenRead(in_file, filename);
  std::vector<double> xyz(3);
  std::string temp;

  while (!NRLib::CheckEndOfFile(in_file)) {
    std::getline(in_file, temp);
    NRLib::ParseAsciiArrayFast(temp, xyz.begin(), 3);
    NRLib::Point pt(xyz[0], xyz[1], xyz[2]);
    points_.push_back(pt);
  }
}


void
PointSet::WriteInternalPoints(const std::string& filename) const
{
  std::ofstream file;
  NRLib::OpenWrite(file, filename);

  file << std::fixed;
  std::map<std::string, std::vector<int> >::const_iterator iter;
  for(iter = integer_tag_.begin(); iter != integer_tag_.end(); ++iter)
    file << "Discrete   " + iter->first << "\n";

  std::map<std::string, std::vector<double> >::const_iterator iter2;
  for(iter2 = float_tag_.begin(); iter2 != float_tag_.end(); ++iter2)
    file << "Float   " + iter2->first << "\n";

  std::map<std::string, std::vector<std::string> >::const_iterator iter3;
  for(iter3 = string_tag_.begin(); iter3 != string_tag_.end(); ++iter3)
    file << "String   " + iter3->first << "\n";


  for(size_t i = 0; i < points_.size(); i++){
    NRLib::Point pt = points_[i];
    file << std::left << std::setprecision(3) << std::setw(15) << pt.x << " "
           << std::setw(15) << pt.y << " " << std::setw(15) << pt.z << " ";
      for(iter = integer_tag_.begin(); iter != integer_tag_.end(); ++iter)
        file << std::setw(15) << iter->second[i] << " ";
      for(iter2 = float_tag_.begin(); iter2 != float_tag_.end(); ++iter2) {
        if (iter2->second[i] == -999.0)
          file << std::setw(15) << "UNDEF";
        else
          file << std::setw(15) << iter2->second[i];
        file << " ";
      }
      for(iter3 = string_tag_.begin(); iter3 != string_tag_.end(); ++iter3){
        if (iter3->second[i].find(" ") != std::string::npos)
          file << "\"" << iter3->second[i] << "\"";
        else
          file << std::setw(15) << iter3->second[i];
        file << " ";
      }
    file << "\n";
  }

}


void
PointSet::WriteRoxarText(const std::string& filename) const
{
  std::ofstream file;
  NRLib::OpenWrite(file, filename);
  std::vector<double> xyz(3);
  file << std::fixed;
  for(size_t i = 0; i < points_.size(); i++){
    xyz[0] = points_[i].x;
    xyz[1] = points_[i].y;
    xyz[2] = points_[i].z;
    NRLib::WriteAsciiArray(file, xyz.begin(), xyz.end(), 3);
  }

  xyz[0] = 999.000000;
  xyz[1] = 999.000000;
  xyz[2] = 999.000000;
  NRLib::WriteAsciiArray(file, xyz.begin(), xyz.end(), 3);

}


void
PointSet::WriteGeneralXYZ(const std::string& filename) const
{
  std::ofstream file;
  NRLib::OpenWrite(file, filename);
  std::vector<double> xyz(3);
  file << std::fixed;
  for(size_t i = 0; i < points_.size(); i++){
    xyz[0] = points_[i].x;
    xyz[1] = points_[i].y;
    xyz[2] = points_[i].z;
    NRLib::WriteAsciiArray(file, xyz.begin(), xyz.end(), 3);
  }
}


bool PointSet::IsRoxarFileEnder(const std::vector<double> &xyz)
{
  return(xyz[0] == 999.000000 && xyz[1] == 999.000000 && xyz[2] == 999.000000);
}


const std::vector<int>&
PointSet::GetDiscreteParameter(const std::string& name) const
{
  std::map<std::string, std::vector<int> >::const_iterator iter = integer_tag_.find(name);
  if(iter != integer_tag_.end())
    return iter->second;
  else
    throw NRLib::Exception("No discrete tag with name " + name);
}

void
PointSet::SetDiscreteParameter(const std::string& name,
                               size_t             point_nr,
                               int                new_parameter)
{
  std::map<std::string, std::vector<int> >::const_iterator iter = integer_tag_.find(name);

  std::vector<int> parameters;
  if(iter != integer_tag_.end()){
    parameters = iter->second;
    parameters[point_nr] = new_parameter;
    integer_tag_[name] = parameters;
  }
  else
    throw NRLib::Exception("No discrete tag with name " + name);
}


bool
PointSet::HasDiscreteParameter(const std::string& name) const
{
  std::map<std::string, std::vector<int> >::const_iterator iter = integer_tag_.find(name);

  if(iter != integer_tag_.end())
    return true;
  else
    return false;
}


const std::vector<double>&
PointSet::GetContParameter(const std::string& name) const
{
  std::map<std::string, std::vector<double> >::const_iterator iter = float_tag_.find(name);

  if(iter != float_tag_.end())
    return iter->second;
  else
    throw NRLib::Exception("No float tag with name " + name);
}


bool
PointSet::HasContParameter(const std::string& name) const
{
  std::map<std::string, std::vector<double> >::const_iterator iter = float_tag_.find(name);

  if(iter != float_tag_.end())
    return true;
  else
    return false;
}


const std::vector<std::string>&
PointSet::GetStringParameter(const std::string& name) const
{
  std::map<std::string, std::vector<std::string> >::const_iterator iter = string_tag_.find(name);

  if(iter != string_tag_.end())
    return iter->second;
  else
    throw NRLib::Exception("No string tag with name " + name);
}


bool
PointSet::HasStringParameter(const std::string& name) const
{
  std::map<std::string, std::vector<std::string> >::const_iterator iter = string_tag_.find(name);

  if(iter != string_tag_.end())
    return true;
  else
    return false;
}


bool
PointSet::HasTags() const
{
  if(integer_tag_.size() > 0 || float_tag_.size() > 0 || string_tag_.size() > 0)
    return true;
  else
    return false;

}


void
PointSet::Resize(size_t new_size)
{
  points_.resize(new_size);

  std::map<std::string, std::vector<int> >::iterator iter;
  for(iter = integer_tag_.begin(); iter != integer_tag_.end(); ++iter)
    iter->second.resize(new_size);

  std::map<std::string, std::vector<double> >::iterator iter2;
  for(iter2 = float_tag_.begin(); iter2 != float_tag_.end(); ++iter2)
    iter2->second.resize(new_size);

  std::map<std::string, std::vector<std::string> >::iterator iter3;
  for(iter3 = string_tag_.begin(); iter3 != string_tag_.end(); ++iter3)
    iter3->second.resize(new_size);
}


void
PointSet::Erase(size_t index)
{
  assert(index < GetSize());

  points_.erase(points_.begin() + index);

  std::map<std::string, std::vector<int> >::iterator iter;
  for(iter = integer_tag_.begin(); iter != integer_tag_.end(); ++iter)
    iter->second.erase(iter->second.begin() + index);

  std::map<std::string, std::vector<double> >::iterator iter2;
  for(iter2 = float_tag_.begin(); iter2 != float_tag_.end(); ++iter2)
    iter2->second.erase(iter2->second.begin() + index);

  std::map<std::string, std::vector<std::string> >::iterator iter3;
  for(iter3 = string_tag_.begin(); iter3 != string_tag_.end(); ++iter3)
    iter3->second.erase(iter3->second.begin() + index);

}


Polygon PointSet::GetConvexHull()
{
  int minpointindex=0; //The index of the point with the lowest y-coordinate
  std::vector<int> indexes;
  indexes.reserve(points_.size());
  for (unsigned int i = 0; i < points_.size(); i++) {
    indexes.push_back(i);
    if (points_[i].y<points_[minpointindex].y) {
      minpointindex=i;
    }
  }
  indexes[minpointindex]=0;
  indexes[0]=minpointindex; //Switching, so that the element corresponding to lowest y-coordinate is first
  AngleSort(indexes, minpointindex, 1, static_cast<int>(indexes.size())-1); //Sorting from second element to end
  RemoveEqualAngles(indexes, minpointindex); //If two points has the same angle wrt p=points_[minpointindex], the one closest to p will be deleted
  //debug
 // std::vector<NRLib::Point> points;
//  for(unsigned int i = 0; i < indexes.size(); i++){
//    points.push_back(points_[indexes[i]]);
 // }
 // NRLib::PointSet polygon(points);
 // polygon.WriteToFile("reducedpoints.dat", NRLib::PointSet::RoxarText);

  std::vector<int> candidates; //Storing indexes of points in the convex hull for the first i (as in the variable i in the for-loop) indexes in the vector indexes
  candidates.push_back(indexes[0]);
  candidates.push_back(indexes[1]);
  candidates.push_back(indexes[2]); //The point with the lowest y-coordinate, the point with the smallest angle and the point with the greatest angle is always in the convex hull
  int top=2;
  int next_to_top=1;
  Point vec1, vec2;
  double cross;
  for (unsigned int i=3; i<indexes.size(); i++) {
    vec1=points_[candidates[top]]-points_[candidates[next_to_top]];
    vec2=points_[indexes[i]]-points_[candidates[top]];
    cross=vec1.x*vec2.y-vec1.y*vec2.x;
    while( cross<=0.0){ //Adding points_[indexes[i]] implies that points_[candidates[top]] is not in the convex hull
      if(candidates.size()==3){
        candidates.pop_back();
        cross = 1.0;
        top--;
        next_to_top--;
      }
      else
      {
        candidates.pop_back();
        top--;
        next_to_top--;
        vec1=points_[candidates[top]]-points_[candidates[next_to_top]];
        vec2=points_[indexes[i]]-points_[candidates[top]];
        cross=vec1.x*vec2.y-vec1.y*vec2.x;
      }
    }
    candidates.push_back(indexes[i]);
    top++;
    next_to_top++;
  }
  Polygon convex_hull;
  for (unsigned int i=0; i<candidates.size(); i++) {
    convex_hull.AddPoint(points_[candidates[i]]);
  }
  return convex_hull;
}



void PointSet::AngleSort(std::vector<int> &indexes, int p_index, int start, int end)
{
  if (start<end) {
    std::vector<int> sorting;
    sorting.reserve(end-start+1);
    int midpoint=int((start+end)/2)+1;
    AngleSort(indexes, p_index, start, midpoint-1); //Sorting left half
    AngleSort(indexes, p_index, midpoint, end); //Sorting right half
    int iter1=start;
    int iter2=midpoint;
    Point vec1, vec2;
    Point vec0;
    vec0.x = 1.0;
    vec0.y = 0.0;
    vec0.z = 0.0;
   // Point helppt;
    //helppt.x = points_[p_index].x + 1000;
    //helppt.y = points_[p_index].y;
    //helppt.z = 0.0;
    //vec0 = helppt-points_[p_index];
    //vec0 = vec0/sqrt(vec0.x*vec0.x+vec0.y*vec0.y);
    //double cross;
    double dot1, dot2;
    while (iter1<midpoint && iter2<=end) { //Merging
      vec1=points_[indexes[iter1]]-points_[p_index];
      if(vec1.x != 0.0 || vec1.y !=0.0)
        vec1 = vec1/sqrt(vec1.x*vec1.x+vec1.y*vec1.y);
      vec2=points_[indexes[iter2]]-points_[p_index];
      if(vec2.x != 0.0 || vec2.y !=0.0)
        vec2 = vec2/sqrt(vec2.x*vec2.x+vec2.y*vec2.y);
      //cross=vec1.x*vec2.y-vec1.y*vec2.x;
      dot1 = vec1.x*vec0.x + vec1.y*vec0.y;
      dot2 = vec0.x*vec2.x + vec0.y*vec2.y;
      //if (cross<=0.0) { //the point at indexes[iter2] has less than or equal angle
      if (dot2 >= dot1) { //the point at indexes[iter2] has less than or equal angle
        sorting.push_back(indexes[iter2]);
        iter2++;
      }
      //if (cross>=0.0) { //the point at indexes[iter2] has less than or equal angle
      if (dot1 >= dot2) { //the point at indexes[iter2] has less than or equal angle
        sorting.push_back(indexes[iter1]);
        iter1++;
      }
    }
    while (iter1<midpoint) { //Filling in the rest of the lsh.
    //[If iter2<=end, we used all the elements in the lhs first, which means that the rest of the elements on the rhs are in the correct cell]
      sorting.push_back(indexes[iter1]);
      iter1++;
    }
    for (unsigned int i=0; i<sorting.size(); i++) {
      indexes[start+i]=sorting[i];
    }
  }
}


void PointSet::RemoveEqualAngles(std::vector<int> &indexes, int p_index)
{
  unsigned int iter1=1;
  unsigned int iter2=2;
  //double cross;
  Point vec1, vec2, vec0;
  vec0.x = 1.0;
  vec0.y = 0.0;
  vec0.z = 0.0;
  double dot1, dot2;
  double norm1, norm2;
  while (iter2<indexes.size()) {
    vec1=points_[indexes[iter1]]-points_[p_index];
    norm1=sqrt(vec1.x*vec1.x+vec1.y*vec1.y);
    vec1.x = vec1.x/norm1;
    vec1.y = vec1.y/norm1;
    vec2=points_[indexes[iter2]]-points_[p_index];
    norm2=sqrt(vec2.x*vec2.x+vec2.y*vec2.y);
    vec2.x = vec2.x/norm2;
    vec2.y = vec2.y/norm2;
    dot1 = vec1.x*vec0.x + vec1.y*vec0.y;
    dot2 = vec0.x*vec2.x + vec0.y*vec2.y;
    //cross=vec1.x*vec2.y-vec1.y*vec2.x;
   // if (fabs(cross) < 0.001) { //the points has the same angle wrt point p
    if(fabs(dot1-dot2) < 1e-10){
      if (norm2>norm1) {
        indexes[iter1]=indexes[iter2];
      }
      iter2++; //The index corresponding to the greatest distance is now stored at indexes[iter1], so the one at iter2 can be deleted later
    }
    else {
      iter1++;
      indexes[iter1]=indexes[iter2];
      iter2++;
    }
  }
  indexes.resize(iter1+1);
}
