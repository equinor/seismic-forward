// $Id: eclipsecarfin.hpp 882 2011-09-23 13:10:16Z perroe $

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

#ifndef NRLIB_ECLIPSEGRID_ECLIPSECARFIN_HPP
#define NRLIB_ECLIPSEGRID_ECLIPSECARFIN_HPP

#include <iosfwd>
#include <string>

#include "eclipsegrid.hpp"


namespace NRLib {

/// A local cartesian grid refinement in an Eclipse Grid
class EclipseCarFin {

public:
  /// Contructs a non-initialised EclipseCarFin object.
  EclipseCarFin(const EclipseGrid& mother_grid);

  /// Contructs a EclipseCarFin object initialised with the given input
  /// values. The indices start at 0. (0,0,0) is the first cell in the grid.
  /// \param name name of local grid refinement section
  /// \param grid parent grid.
  /// \param i0 index from in x direction (base 0)
  /// \param i1 index to in x direction (base 0)
  /// \param j0 index from in y direction (base 0)
  /// \param j1 index to in y direction (base 0)
  /// \param k0 index from in z direction (base 0)
  /// \param k1 index to in z direction (base 0)
  /// \param nx number of cells in the x direction
  /// \param ny number of cells in the y direction
  /// \param nz number of cells in the z direction
  /// \param nz number of cells in the z direction
  /// \param nwmax max number of wells inside refinement. Defaults to -1,
  ///              if the parameter is not given.
  EclipseCarFin(const std::string & name,
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
                int                 nwmax = -1);

  /// Constructs a EclipseCarFin-object based on a given grid.
  /// \param name name of local grid refinement section.
  /// \param grid grid with properties for the local grid. (Not copied.
  ///            Will be deleted by the destructor.)
  /// \param parent_grid parent grid.
  /// \param i0 index from in x direction (base 0)
  /// \param i1 index to in x direction (base 0)
  /// \param j0 index from in y direction (base 0)
  /// \param j1 index to in y direction (base 0)
  /// \param k0 index from in z direction (base 0)
  /// \param k1 index to in z direction (base 0)
  /// \param nwmax max number of wells inside refinement. Defaults to -1
  ///             (unspesified), if the parameter is not given.
  EclipseCarFin(const std::string & name,
                const EclipseGrid & grid,
                const EclipseGrid & parent_grid,
                size_t              i0,
                size_t              i1,
                size_t              j0,
                size_t              j1,
                size_t              k0,
                size_t              k1,
                int nwmax = -1);


  /// Returns the name of the local grid.
  /// \return name of grid
  const std::string& GetName() const      { return name_; }


  /// Sets the name of the local grid.
  /// \param name name of grid
  void SetName(const std::string& name)  { name_ = name; }


  /// Returns the location within the mother grid of the refinement.
  /// The returned indices start at 0 for the first cell.
  /// \param i0 from in i direction.
  /// \param i1 to in i direction.
  /// \param j0 from in j direction.
  /// \param j1 to in j direction.
  /// \param k0 from in k direction.
  /// \param k1 to in k direction.
  void GetLocation(size_t& i0, size_t& i1,
                   size_t& j0, size_t& j1,
                   size_t& k0, size_t& k1) const;


  /// Returns the grid refinement pr cell, provided that it is a whole number.
  /// If not the closest whole number is returned, and a warning is given.
  /// \param refX refinement in x direction.
  /// \param refY refinement in y direction.
  /// \param refZ refinement in z direction.
  void GetRefinements(size_t & ref_x,
                      size_t & ref_y,
                      size_t & ref_z) const;


  /// Returns a reference to the ECLIPSE grid with the local grid properties.
  /// \return reference to ECLIPSE grid.
  const EclipseGrid& GetGrid() const { return eclipse_grid_; };


  /// Read a CARFIN section from an ECLIPSE file.
  /// \param file Input file.
  /// \returns true if there are more refinements (CARFIN is encountered),
  ///          or false if there are no more refinements yet (ENDFIN is encountered).
  bool ReadCarFin(std::ifstream& file);


  /// Outputs the object to an ECLIPSE-file.
  /// \param file filestream to write to.
  void WriteCarFin(std::ofstream& file) const;


  /// Reads properties from an ECLIPSE mother grid.
  /// \param grid mother grid.
  // void readProperties(const ECLIPSE& grid);


  /// Checks if the given cartesian refinement is connected to this
  /// cartesian refinement.
  /// \param carFin  cartesian refinement.
  /// \return true if connected, else false.
  bool IsConnected(const EclipseCarFin& refinment) const;


  /// Exports the local grid as an independent ECLIPSE grid. The geometry of
  /// the mother grid is used for the grid if no local geometry is given.
  /// \param motherGrid mother grid.
  /// \return pointer to generated ECLIPSE grid. Should be released by caller.
  //void ConvertToGrid(EclipseGrid& grid_out) const;

private:
  /// Mother cell.
  /// \param index lgr cell index.
  /// \param grid mother grid.
  /// \return mother grid index.
  //size_t FindMotherCell(size_t              index,
  //                      const EclipseGrid & mother_grid) const;

  /// Checks if a keyword is allowed within a local grid
  /// \param i index of keyword as defined in EclipseKeywords.
  /// \return true if keyword is allowed, else false.
  bool IsKeywordAllowed(const std::string& keyword) const;

  /// No assignment operator possible due to constant grid reference.
  // In the future we
  // might check that the grid references matches, or make the reference non-constant.
  EclipseCarFin& operator=(const EclipseCarFin& rhs);

  /// Eclipse object holding local grid properties.
  EclipseGrid eclipse_grid_;

  const EclipseGrid& mother_grid_;

  /// Name of refinement
  std::string name_;

  /// Index from in x direction
  size_t      i0_;

  /// Index to in x direction
  size_t      i1_;

  /// Index from in y direction
  size_t      j0_;

  /// Index to in y direction
  size_t      j1_;

  /// Index from in z direction
  size_t      k0_;

  /// Index to in z direction
  size_t      k1_;

  /// Number of blocks inside local refinement in x direction
  size_t      nx_;

  /// Number of blocks inside local refinement in y direction
  size_t      ny_;

  /// Number of blocks inside local refinement in z direction
  size_t      nz_;

  /// Max number of wells. -1 if unspesified.
  int         nwmax_;
};

} // namespace NRLib

#endif // NRLIB_ECLIPSEGRID_ECLIPSECARFIN_HPP
