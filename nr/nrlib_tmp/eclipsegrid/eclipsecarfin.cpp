// $Id: eclipsecarfin.cpp 882 2011-09-23 13:10:16Z perroe $

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

#include "eclipsecarfin.hpp"

#include "eclipsegrid.hpp"
#include "../iotools/fileio.hpp"

#include <fstream>

namespace NRLib {

EclipseCarFin::EclipseCarFin(const EclipseGrid & mother_grid)
  : mother_grid_(mother_grid),
    nwmax_(-1)
{}


EclipseCarFin::EclipseCarFin(const std::string & name,
                             const EclipseGrid & mother_grid,
                             size_t              i0,
                             size_t              i1,
                             size_t              j0,
                             size_t              j1,
                             size_t              k0,
                             size_t              k1,
                             size_t              nx,
                             size_t              ny,
                             size_t              nz,
                             int                 nwmax)
  : eclipse_grid_(nx, ny, nz),
    mother_grid_(mother_grid),
    name_(name),
    i0_(i0), i1_(i1),
    j0_(j0), j1_(j1),
    k0_(k0), k1_(k1),
    nx_(nx), ny_(ny),
    nz_(nz), nwmax_(nwmax)
{
  // There must be a match between global and local cells:
  assert (nx_%(i1_-i0_+1) == 0 && ny_%(j1_-j0_+1) == 0 && nz_%(k1_-k0_+1) == 0);

  for (size_t index = 0; index < eclipse_grid_.GetN(); ++index) {
    eclipse_grid_.SetActive(index, true);
  }
}


EclipseCarFin::EclipseCarFin(const std::string & name,
                             const EclipseGrid & grid,
                             const EclipseGrid & mother_grid,
                             size_t              i0,
                             size_t              i1,
                             size_t              j0,
                             size_t              j1,
                             size_t              k0,
                             size_t              k1,
                             int nwmax)
  : eclipse_grid_(grid),
    mother_grid_(mother_grid),
    name_(name),
    i0_(i0),
    i1_(i1),
    j0_(j0),
    j1_(j1),
    k0_(k0),
    k1_(k1),
    nwmax_(nwmax)
{
  nx_ = eclipse_grid_.GetNI();
  ny_ = eclipse_grid_.GetNJ();
  nz_ = eclipse_grid_.GetNK();
}


void EclipseCarFin::GetLocation(size_t& i0, size_t& i1,
                                size_t& j0, size_t& j1,
                                size_t& k0, size_t& k1) const
{
  i0 = i0_;
  i1 = i1_;
  j0 = j0_;
  j1 = j1_;
  k0 = k0_;
  k1 = k1_;
}


void EclipseCarFin::GetRefinements(size_t& ref_x,
                                   size_t& ref_y,
                                   size_t& ref_z) const
{
  assert (nx_%(i1_-i0_+1) == 0 && ny_%(j1_-j0_+1) == 0 && nz_%(k1_-k0_+1) == 0);

  ref_x = nx_ / (i1_ - i0_ + 1);
  ref_y = ny_ / (j1_ - j0_ + 1);
  ref_z = nz_ / (k1_ - k0_ + 1);
}


//const ECLIPSE& EclipseCarFin::getGrid() const
//{
//  return *eclipseGrid_;
//}


//void EclipseCarFin::readProperties(const ECLIPSE& grid) {
//  out.message() << "Reading properties for " << name_ << std::endl;
//  EclipseKeywords key;
//  unsigned int nPix = nx_ * ny_ * nz_;
//  for (int c = 0; c < EclipseKeywords::n_comm; ++c) {
//    if ( keywordAllowed(c) && grid.hasKeyword(c) ) {
//      switch (c)
//  {
//  case EclipseKeywords::actnum: {
//    eclipseGrid_->setCommand(EclipseKeywords::actnum, true);
//    for (unsigned int p = 0; p < nPix; ++p) {
//      bool status(grid.isActive(motherCell(p, grid)));
//      eclipseGrid_->setActive(p, status);
//    }
//
//    break;
//  }
//  case EclipseKeywords::permx: {
//    float* array = eclipseGrid_->getFloatArray("PERMX");
//    for (unsigned int p = 0; p < nPix; ++p) {
//      array[p] = grid.getPERMX(motherCell(p, grid));
//    }
//
//    break;
//  }
//  case EclipseKeywords::permy: {
//    float* array = eclipseGrid_->getFloatArray("PERMY");
//    for (unsigned int p = 0; p < nPix; ++p) {
//      array[p] = grid.getPERMY(motherCell(p, grid));
//    }
//    break;
//  }
//  case EclipseKeywords::permz: {
//    float* array = eclipseGrid_->getFloatArray("PERMZ");
//    for (unsigned int p = 0; p < nPix; ++p) {
//      array[p] = grid.getPERMZ(motherCell(p, grid));
//    }
//    break;
//  }
//  case EclipseKeywords::poro: {
//    float* array = eclipseGrid_->getFloatArray("PORO");
//    for (unsigned int p = 0; p < nPix; ++p) {
//      array[p] = grid.getPORO(motherCell(p, grid));
//    }
//    break;
//  }
//  case EclipseKeywords::ntg: {
//    float* array = eclipseGrid_->getFloatArray("NTG");
//    for (unsigned int p = 0; p < nPix; ++p) {
//      array[p] = grid.getNTG(motherCell(p, grid));
//    }
//    break;
//  }
//       case EclipseKeywords::facies: {
//    short int* array = eclipseGrid_->getShortIntArray("FACIES");
//    for (unsigned int p = 0; p < nPix; ++p) {
//      array[p] = grid.getFACIES(motherCell(p, grid));
//    }
//    break;
//  }
//       case EclipseKeywords::strain: {
//    float* array = eclipseGrid_->getFloatArray("STRAIN");
//          for (unsigned int p = 0; p < nPix; ++p) {
//      array[p] = grid.getSTRAIN(motherCell(p, grid));
//    }
//    break;
//  }
//       case EclipseKeywords::lower: {
//    float* array = eclipseGrid_->getFloatArray("LOWER");
//          for (unsigned int p = 0; p < nPix; ++p) {
//      array[p] = grid.getLOWER(motherCell(p, grid));
//    }
//    break;
//  }
//       case EclipseKeywords::upper: {
//    float* array = eclipseGrid_->getFloatArray("UPPER");
//          for (unsigned int p = 0; p < nPix; ++p) {
//      array[p] = grid.getUPPER(motherCell(p, grid));
//    }
//    break;
//  }
//      }
//      out.message() << "Reading " << key.getKeyword(c) << std::endl;
//    }
//  }
//}
//


// TODO: We only support CARFIN boxes that end with ENDFIN, not
// CARFIN boxes that are ended by a superseding CARFIN box.

bool EclipseCarFin::ReadCarFin(std::ifstream& file)
{
  std::string token;
  int dummy_line = 0; // not used.
  try {
    name_ = NRLib::ReadNext<std::string>(file, dummy_line);
    // The indices in the file format are 1 based.
    i0_ = NRLib::ReadNext<size_t>(file, dummy_line) - 1;
    i1_ = NRLib::ReadNext<size_t>(file, dummy_line) - 1;
    j0_ = NRLib::ReadNext<size_t>(file, dummy_line) - 1;
    j1_ = NRLib::ReadNext<size_t>(file, dummy_line) - 1;
    k0_ = NRLib::ReadNext<size_t>(file, dummy_line) - 1;
    k1_ = NRLib::ReadNext<size_t>(file, dummy_line) - 1;
    nx_ = NRLib::ReadNext<size_t>(file, dummy_line);
    ny_ = NRLib::ReadNext<size_t>(file, dummy_line);
    nz_ = NRLib::ReadNext<size_t>(file, dummy_line);

    token = NRLib::ReadNext<std::string>(file, dummy_line);
    if (NRLib::IsType<int>(token)) {
      nwmax_ = NRLib::ParseType<int>(token);
      token = NRLib::ReadNext<std::string>(file, dummy_line);
    }
    else {
      nwmax_ = -1;
    }
  }
  catch (NRLib::EndOfFile& ) {
    throw NRLib::FileFormatError("Unexpected end of file while reading CARFIN named \'"
                                 + name_ + "\'");
  }
  catch (NRLib::Exception& e) {
    throw NRLib::FileFormatError("Error reading CARFIN named \'" + name_ + "\' :"
                                 + e.what());
  }

  if (token != "/") {
    throw NRLib::FileFormatError("Error reading CARFIN named \'" + name_ + "\' :"
                                 + " missing final /.");
  }

  try {
    token = NRLib::ReadNext<std::string>(file, dummy_line);

    while (token != "ENDFIN" && token != "CARFIN") {
      if (!IsKeywordAllowed(token)) {
        throw NRLib::FileFormatError( "Not supported " + token + " command "
            + " appears inside a CARFIN box named " + name_ + ".");
      } else {
        eclipse_grid_.ReadKeyword(file, token);
      }
      token = NRLib::ReadNext<std::string>(file, dummy_line);
    }
  }
  catch (NRLib::EndOfFile& ) {
    throw NRLib::FileFormatError("Unexpected end of file while reading parameters for "
                                 "CARFIN named \'" + name_ + "\'");
  }
  catch (NRLib::Exception& e) {
    throw NRLib::FileFormatError("Error reading parameters for CARFIN named \'"
                                 + name_ + "\' :" + e.what());
  }

  return (token == "CARFIN") ? true : false;
}


void EclipseCarFin::WriteCarFin(std::ofstream& of) const
{
  of << "CARFIN\n";
  // Adding 1 to get base 1.
  //  of << "-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ  NWMAX\n";
  of << "\'" << name_ << "\'  "
     << i0_ + 1 << "  " << i1_ + 1 << "  "
     << j0_ + 1 << "  " << j1_ + 1 << "  "
     << k0_ + 1 << "  " << k1_ + 1 << "  "
     << nx_ << "  " << ny_ << "  " << nz_;
  if (nwmax_ >= 0)
    of << "  " << nwmax_;
  of <<  "    /\n";

  // eclipse_grid_.WriteKeyword(of, "ACTNUM");
  // eclipse_grid_.WriteKeyword(of, "COORD");
  // eclipse_grid_.WriteKeyword(of, "ZCORN");
  eclipse_grid_.WriteContinousParameters(of);

  of << "ENDFIN\n\n";
}


bool EclipseCarFin::IsConnected(const EclipseCarFin& car_fin) const
{
  assert(&mother_grid_ == &(car_fin.mother_grid_));

  size_t o_i0, o_i1, o_j0, o_j1, o_k0, o_k1;
  car_fin.GetLocation(o_i0, o_i1, o_j0, o_j1, o_k0, o_k1);

  if (o_i0 > i1_ + 1 || o_i1 < (i0_ > 0 ? i0_ - 1 : i0_) ||
      o_j0 > j1_ + 1 || o_j1 < (j0_ > 0 ? j0_ - 1 : j0_) ||
      o_k0 > k1_ + 1 || o_k1 < (k0_ > 0 ? k0_ - 1 : k0_)    ) {
    return false;
  }

  return true;
}


//void EclipseCarFin::ConvertToGrid(EclipseGrid& grid_out) const
//{
//  grid_out =
//
//  if (!grid->hasKeyword(EclipseKeywords::coord)) {
//    grid->refineCornerpointGrid(motherGrid, i0_, i1_, j0_, j1_, k0_, k1_);
//    grid->setCommand(EclipseKeywords::coord, true);
//    grid->setCommand(EclipseKeywords::zcorn, true);
//  }
//
//  return grid;
//}


bool EclipseCarFin::IsKeywordAllowed(const std::string& keyword) const
{
  // Since we support general parameters, not supported by ECLIPSE, we test on
  // disallowed keywords.
  if (keyword == "AMALGAM"  ||
      keyword == "CARFIN"   ||
      keyword == "ECHO"     ||
      keyword == "FAULTS"   ||
      keyword == "GRIDUNIT" ||
      keyword == "MAPAXES"  ||
      keyword == "MAPUNITS" ||
      keyword == "NOECHO"   ||
      keyword == "SPECGRID"
      )
  {
    return false;
  }

  return true;
}


//size_t EclipseCarFin::FindMotherCell(size_t             index,
//                                     const EclipseGrid& mother_grid) const
//{
//  // LGR indices
//  unsigned int i, j, k;
//  eclipseGrid_->pixel2ijk(index, i, j, k);
//  // Mother indices
//  unsigned int I, J, K;
//  unsigned int rx = nx_/(i1_ - i0_ + 1);
//  unsigned int ry = ny_/(j1_ - j0_ + 1);
//  unsigned int rz = nz_/(k1_ - k0_ + 1);
//  I = i0_ + i/rx;
//  J = j0_ + j/ry;
//  K = k0_ + k/rz;
//  return grid.ijk2pixel(I, J, K);
//}

}
