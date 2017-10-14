// $Id: eclipsecarfingroup.hpp 882 2011-09-23 13:10:16Z perroe $

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

#ifndef NRLIB_ECLIPSEGRID_ECLIPSECARFINGROUP_HPP
#define NRLIB_ECLIPSEGRID_ECLIPSECARFINGROUP_HPP

#include <string>
#include <vector>

namespace NRLib {
  class EclipseCarFin;
  class EclipseGrid;
}


namespace NRLib {

/// Class handling a group of cartesian local grids.
class EclipseCarFinGroup {
public:
  typedef std::vector<const EclipseCarFin*>::iterator        RefinementIterator;
  typedef std::vector<const EclipseCarFin*>::const_iterator  RefinementConstIterator;

  /// Constructs a CarFinGroup with a given name.
  /// \param name Name of the group of refinements.
  EclipseCarFinGroup(const std::string& name);

  /// Returns the name of the local grid.
  /// \return name of grid
  const std::string& GetName() const { return name_; }

  RefinementIterator      RefinementsBegin()       { return cartesian_refinements_.begin(); }
  RefinementConstIterator RefinementsBegin() const { return cartesian_refinements_.begin(); }
  RefinementIterator      RefinementsEnd()         { return cartesian_refinements_.end(); }
  RefinementConstIterator RefinementsEnd()   const { return cartesian_refinements_.end(); }

  /// Returns the location within the mother grid of the smallest box that can
  /// contaub the whole local grid. The parameters are used for returnvalues.
  /// The first cell in the grid has indices (0,0,0).
  /// \param i_from from in i direction.
  /// \param i_to to in i direction.
  /// \param j_from from in j direction.
  /// \param j_to to in j direction.
  /// \param k_from from in k direction.
  /// \param k_to to in k direction.
  void GetLocation(size_t &i_from, size_t &i_to,
                   size_t &j_from, size_t &j_to,
                   size_t &k_from, size_t &k_to) const;

  /// Outputs the object to an ECLIPSE-file using the ECLIPSE AMALGAM keyword.
  /// \param file filestream to write to.
  void WriteToFile(std::ofstream& file) const;

  /// Adds a local grid to the group.
  /// \param cartesian_refinement grid to add.
  void AddCarFin(const EclipseCarFin* cartesian_refinement);

  /// Exports the local grid as an independent ECLIPSE grid. The geometry of
  /// the mother grid is used for the grid if no local geometry is given.
  /// \param[in]  mother_grid mother grid.
  /// \param[out] grid_out    grid consisting of all the local grid.
  void ConvertToEclipseGrid(const EclipseGrid& mother_grid, EclipseGrid& grid_out) const;

  /// Checks if the given cartesian refinement is connected to any
  /// of the refinements already in the group.
  /// \param cartesian_refinement  cartesian refinement.
  /// \return true if connected, else false.
  bool IsConnected(const EclipseCarFin& cartesian_refinement) const;

private:
  /// Vector of local grids
  std::vector<const EclipseCarFin*> cartesian_refinements_;

  /// Name of local grid
  std::string name_;

  /// Refinement of each cell in x direction
  size_t refinement_x_;

  /// Refinement of each cell in x direction
  size_t refinement_y_;

  /// Refinement of each cell in x direction
  size_t refinement_z_;

  static void RefineGrid(const EclipseGrid & mother_grid,
                         size_t              i_from,
                         size_t              i_to,
                         size_t              j_from,
                         size_t              j_to,
                         size_t              k_from,
                         size_t              k_to,
                         EclipseGrid       & grid_out);
};

} // namespace NRLib

#endif // NRLIB_ECLIPSEGRID_ECLIPSECARFINGROUP_HPP
