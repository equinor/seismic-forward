// $Id: eclipsegrid.hpp 1187 2013-06-14 13:34:45Z perroe $

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

#ifndef NRLIB_ECLIPSEGRID_ECLIPSEGRID_HPP
#define NRLIB_ECLIPSEGRID_ECLIPSEGRID_HPP

#include "eclipsegeometry.hpp"
#include "eclipsefault.hpp"
#include "eclipsetransmult.hpp"

#include "../grid/grid.hpp"

#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>


namespace NRLib {

class BilinearSurface;
class EclipseCarFin;
class EclipseCarFinGroup;

class EclipseGrid{

public:
  typedef std::map<std::string, EclipseFault>::iterator               FaultIterator;
  typedef std::map<std::string, EclipseFault>::const_iterator         FaultConstIterator;

  typedef std::map<std::string, NRLib::Grid<double> >::iterator       ParameterIterator;
  typedef std::map<std::string, NRLib::Grid<double> >::const_iterator ParameterConstIterator;

  typedef std::vector<EclipseCarFin*>::iterator                       RefinementIterator;
  typedef std::vector<EclipseCarFin*>::const_iterator                 RefinementConstIterator;

  typedef std::vector<EclipseCarFinGroup*>::iterator                  RefinementGroupIterator;
  typedef std::vector<EclipseCarFinGroup*>::const_iterator            RefinementGroupConstIterator;

  //constructors
  EclipseGrid();
  EclipseGrid(size_t nx, size_t ny, size_t nz);
  EclipseGrid(const std::string& file_name);
  EclipseGrid(const std::vector<std::string>& file_names);

  ~EclipseGrid();

  // ====================== Grid geometry methods ==============================

  size_t GetNI() const    { return geometry_.GetNI(); }
  size_t GetNJ() const    { return geometry_.GetNJ(); }
  size_t GetNK() const    { return geometry_.GetNK(); }
  size_t GetN()  const    { return geometry_.GetN();  }

  inline void   GetIJK(size_t index, size_t& i, size_t& j, size_t& k) const   { geometry_.GetIJK(index, i ,j, k); }
  inline size_t GetIndex(size_t i, size_t j, size_t k) const                  { return geometry_.GetIndex(i, j, k); }

  const Line& GetPillar(size_t i, size_t j) const {return geometry_.GetPillar(i, j); }

  void SetPillar(size_t i, size_t j, const Line & pillar)
  { geometry_.SetCoordLine(i, j, pillar); }

  double GetZCorner(size_t i, size_t j, size_t k, size_t a, size_t b, size_t c) const
  { return geometry_.GetZCorner(i, j, k, a, b, c); }

  void SetZCorner(size_t i, size_t j, size_t k, size_t a, size_t b, size_t c, double z)
  { geometry_.SetZCorner(i, j, k, a, b, c, z); }

  /// Resizes the grid. Everything, parameters, z-corners, etc., is nulled out.
  void Resize(size_t nx, size_t ny, size_t nz);

  /// Check if pillar is active, i.e. if there is at least one active cell connected
  /// to this pillar.
  bool IsPillarActive(size_t i, size_t j) const
  { return geometry_.IsPillarActive(i, j); }

  /// Is cell active
  bool IsActive(size_t i, size_t j, size_t k) const
  { return geometry_.IsActive(i, j, k); }

  /// Set cell to be active.
  void SetActive(size_t i, size_t j, size_t k, bool is_active = true)
  { geometry_.SetActive(i, j, k, is_active); }

  /// Set cell to be active.
  void SetActive(size_t index, bool is_active = true)
  { geometry_.SetActive(index, is_active); }

  const EclipseGeometry& GetGeometry() const    { return geometry_;}

  void SetGeometry(const EclipseGeometry& geometry);

  // Gives k-value of top-most active cell in column returns nk if all cells are inactive
  size_t FindTopCell(size_t i, size_t j) const { return geometry_.FindTopCell(i, j); }

  // Gives k-value of lowest active cell in column.
  size_t FindBottomCell(size_t i, size_t j) const { return geometry_.FindBottomCell(i, j); }

  NRLib::Point FindPointAtPillar(size_t i, size_t j, double z) const
  { return geometry_.FindPointAtPillar(i, j, z); }

  NRLib::Point FindPointAtPillarInsideGrid(size_t i, size_t j, double z, bool & found) const
  { return geometry_.FindPointAtPillarInsideGrid(i, j, z, found); }

  //double  FindZTopAtPillar(size_t i, size_t j) const
  //{ return geometry_.FindZTopAtPillar(i,j); }

  //double  FindZBotAtPillar(size_t i, size_t j) const
  //{ return geometry_.FindZBotAtPillar(i,j); }

  //Gives the neighbour column in the face direction
  void GetNeighbourColumn(int i, int j, Face face, int& i_out, int& j_out) const;

  //Gives the k value for the neighbour if active, otherwise it returns the input k
  int GetNeighbourK(int i, int j, int k, Face face) const;

  std::vector<Point> GetFaceCornerPoints(size_t i_in, size_t j_in, size_t k_in, Face face) const;
  std::vector<Point> GetAllCornerPoints(size_t i_in, size_t j_in, size_t k_in) const;

  /// \return false if the given point is not within the grid.
  bool FindIndex(double x_in, double y_in, double z_in,
                 size_t& i_out, size_t& j_out, size_t& k_out ) const
  { return geometry_.FindIndex(x_in, y_in, z_in, i_out, j_out, k_out); }

  // Find the intersection area for the current cell and the cell in the face direction
  // (not PosZ or NegZ)
  void FindIntersectionAreaFaceXY(CellFace   cf,
                                  int        k_other_cell,
                                  double   & a_x,
                                  double   & a_y,
                                  double   & a_z) const;

  EclipseGrid ExtractSubGrid(size_t i_from, size_t i_to,
                             size_t j_from, size_t j_to,
                             size_t k_from, size_t k_to) const;

  // ====================== Grid parameter methods ==============================

  /// Constant access operator - general parameter.
  /// Fails with assert if paramter does not exist - check with HasParameter first!
  inline const Grid<double>& operator()(const std::string& parameter_name) const;

  inline bool HasParameter(const std::string& parameter_name) const;

  /// Get general parameter - same as constant variant of access operator.
  /// Fails with assert if paramter does not exist - check with HasParameter first!
  const Grid<double>& GetParameter(const std::string& parameter_name) const
  { return (*this)(parameter_name); }

  /// Get general parameter
  /// Fails with assert if paramter does not exist - check with HasParameter first!
  Grid<double>& GetParameter(const std::string& parameter_name)
  { return continuous_parameters_[parameter_name]; }

  ParameterIterator       ParametersBegin()       { return continuous_parameters_.begin(); }
  ParameterConstIterator  ParametersBegin() const { return continuous_parameters_.begin(); }
  ParameterIterator       ParametersEnd()         { return continuous_parameters_.end(); }
  ParameterConstIterator  ParametersEnd()   const { return continuous_parameters_.end(); }

  void AddParameter(const std::string& parameter_name, double val = 0.0);


  // ====================== Fault accessors ==============================

  const std::map<std::string, EclipseFault>& GetFaults() const { return faults_; }

  FaultIterator      FaultsBegin()       { return faults_.begin(); }
  FaultConstIterator FaultsBegin() const { return faults_.begin(); }
  FaultIterator      FaultsEnd()         { return faults_.end(); }
  FaultConstIterator FaultsEnd()   const { return faults_.end(); }

  /// \return true if at least one of the fault traces in the grid is stair-stepped.
  bool HasStairSteppedFaults() const;

  void AddFault(const EclipseFault& fault_trace);

  // ====================== Local grid refinements ========================

  RefinementIterator      RefinementsBegin()       { return grid_refinements_.begin(); }
  RefinementConstIterator RefinementsBegin() const { return grid_refinements_.begin(); }
  RefinementIterator      RefinementsEnd()         { return grid_refinements_.end(); }
  RefinementConstIterator RefinementsEnd()   const { return grid_refinements_.end(); }

  void AddGridRefinement(const EclipseCarFin& refinement);


  RefinementGroupIterator      RefinementGroupsBegin()       { return refinement_groups_.begin(); }
  RefinementGroupConstIterator RefinementGroupsBegin() const { return refinement_groups_.begin(); }
  RefinementGroupIterator      RefinementGroupsEnd()         { return refinement_groups_.end(); }
  RefinementGroupConstIterator RefinementGroupsEnd()   const { return refinement_groups_.end(); }

  void AddRefinementGroup(const EclipseCarFinGroup& refinement_group);

  void ClearRefinementGroups();

  // ======================= Transmissibilites ==========================

  /// \note This code is still not complete, at least not completely tested!!
  void AddTransMult(const EclipseTransMult& trans_mult);

  // ===========================  File IO  ===============================

  //Reads multiple files
  void ReadFromFiles(const std::vector<std::string>& file_names);
  //Reads one file
  void ReadFromFile(const std::string& file_name);

  void ReadKeyword(std::ifstream& file, const std::string& keyword);

  /// Writes to one file
  void WriteToFile(const std::string& file_name,
                   bool write_param = true) const;

  /// Writes to several files (number of files is the size of the vectors)
  /// \param file_names Names of output eclipse files.
  /// \param keywords   Lists of keywords, one for each file. If the list of keywords
  ///                   is empty for the last file, all remaining parameters are written to
  ///                   this file.
  /// \throws Exception
  void WriteToFiles(const std::vector<std::string>             & file_names,
                    const std::vector<std::list<std::string> > & keywords ) const;

  //Returns the 4 corner points for in the face direction for cell i_in, j_in, k_in

  void WriteKeyword(std::ofstream& file, const std::string& keyword) const;

  void WriteContinousParameters(std::ofstream& file) const;

private:

  // ================== PRIVATE MEMBER FUNCTIONS ========================

  //Reads a continuous parameter. One double for each cell.
  void ReadContinuousParameter(std::ifstream& in_file, const std::string& map_key);

  void ReadFaults(std::ifstream& in_file);
  void ReadMapUnits(std::ifstream& in_file, std::string& map_units);
  void ReadMapAxes(std::ifstream& in_file, double& map_ax_x1, double& map_ax_y1,
                                           double& map_ax_x2, double& map_ax_y2,
                                           double& map_ax_x3, double& map_ax_y3 );
  void ReadGridUnit(std::ifstream& in_file, std::string& grid_unit, bool& local_coordinates);
  void ReadEditNNC(std::ifstream& in_file);

  void WriteRemainingParameters(const std::string           & file_name,
                                const std::set<std::string> & parameters_written) const;
  void WriteContParameter(std::ofstream& out_file,
                          const std::string         & name,
                          const NRLib::Grid<double> & parameter) const;
  void WriteFaults(std::ofstream& out_file) const;
  void WriteRefinementGroups(std::ofstream& out_file) const;
  void WriteEditNNC(std::ofstream& out_file) const;
  void WriteMultiply(std::ofstream& out_file) const;

  //0 is x, 1 is y and 2 is z (for the direction)
  void WriteTran(std::ofstream& out_file, int direction) const;

  // "k" is used to indicate that it can be used for x, y and z
  // Finds the start and end point for each of the 8 lines
  void FindStartAndEndForLines(double& kmin11, double& kmax11, double& kmin12, double& kmax12,
                               double& kmin13, double& kmax13, double& kmin14, double& kmax14,
                               double& kmin21, double& kmax21, double& kmin22, double& kmax22,
                               double& kmin23, double& kmax23, double& kmin24, double& kmax24) const;

  // Finds the start and end point for the while loop
  void FindStartAndEndForAlgorithm(double kmin11, double kmin12, double kmin21,
                                   double kmin22, double kmax11, double kmax12, double kmax21,
                                   double kmax22, double& k_start, double& k_end) const;

  // FindIntersectMFromN finds the M value from the intersection between the surface diagonal
  // and a constant N value
  double FindIntersectYFromX( double x, const BilinearSurface& surface) const;
  double FindIntersectYFromZ( double z, const BilinearSurface& surface) const;
  double FindIntersectXFromZ( double z, const BilinearSurface& surface) const;

  // FindIntersectMFromNPillar finds the M value from the intersection between the line and
  // a constant N value
  double FindIntersectYFromXPillar( double x, const Line& line) const;
  double FindIntersectYFromZPillar( double x, const Line& line) const;
  double FindIntersectXFromZPillar( double x, const Line& line) const;

  // Finds the u value used in the FindIntersection functions
  double FindUFromK(double a, double b, double c) const;

  // Finds the height in the intersecting area
  bool FindHeight(double a_11, double a_12, double a_21, double a_22, double& height) const;

  //FindPointsN gives the 16 values involved in the area calculations in the N direction
  void FindPointsPosX(Point& p11, Point& p12, Point& p13, Point& p14,
                      Point& q11, Point& q12, Point& q13, Point& q14,
                      Point& p21, Point& p22, Point& p23, Point& p24,
                      Point& q21, Point& q22, Point& q23, Point& q24,
                      CellFace cf, int k_other_cell) const;

  void FindPointsNegX(Point& p11, Point& p12, Point& p13, Point& p14,
                      Point& q11, Point& q12, Point& q13, Point& q14,
                      Point& p21, Point& p22, Point& p23, Point& p24,
                      Point& q21, Point& q22, Point& q23, Point& q24,
                      CellFace cf, int k_other_cell) const;

  void FindPointsPosY(Point& p11, Point& p12, Point& p13, Point& p14,
                      Point& q11, Point& q12, Point& q13, Point& q14,
                      Point& p21, Point& p22, Point& p23, Point& p24,
                      Point& q21, Point& q22, Point& q23, Point& q24,
                      CellFace cf, int k_other_cell) const;

  void FindPointsNegY(Point& p11, Point& p12, Point& p13, Point& p14,
                      Point& q11, Point& q12, Point& q13, Point& q14,
                      Point& p21, Point& p22, Point& p23, Point& p24,
                      Point& q21, Point& q22, Point& q23, Point& q24,
                      CellFace cf, int k_other_cell) const;

  // FindAreaZ finds the area of the projection in the xy plane
  void FindAreaZ(Point& p11, Point& p12, Point& p13, Point& p14,
                 Point& q11, Point& q12, Point& q13, Point& q14,
                 Point& p21, Point& p22, Point& p23, Point& p24,
                 Point& q21, Point& q22, Point& q23, Point& q24, double& a_z) const;

  void FindAreaX(Point& p11, Point& p12, Point& p13, Point& p14,
                 Point& q11, Point& q12, Point& q13, Point& q14,
                 Point& p21, Point& p22, Point& p23, Point& p24,
                 Point& q21, Point& q22, Point& q23, Point& q24, double& a_x) const;

  void FindAreaY(Point& p11, Point& p12, Point& p13, Point& p14,
                 Point& q11, Point& q12, Point& q13, Point& q14,
                 Point& p21, Point& p22, Point& p23, Point& p24,
                 Point& q21, Point& q22, Point& q23, Point& q24, double& a_y) const;

  // ======================= MEMBER VARIABLES =========================

  /// The grid geometry.
  EclipseGeometry                      geometry_;

  /// Grid parameters. One double value per cell.
  std::map<std::string, Grid<double> > continuous_parameters_;

  /// Faults.
  std::map<std::string, EclipseFault>  faults_;

  /// Cartesian local grid refinements.
  std::vector<EclipseCarFin*>          grid_refinements_;

  /// Groups of local grid refinements.
  std::vector<EclipseCarFinGroup*>     refinement_groups_;

  /// Grid units description. \todo Move to geometry?
  std::string                          units_;

  /// Transmissibility multipliers. Currently not completely supported.
  std::list<EclipseTransMult>          trans_mult_;
};


// ================ IMPLEMENTATION OF INLINE FUNCTIONS =====================

const Grid<double>& EclipseGrid::operator()(const std::string& parameter_name) const
{
  std::map<std::string, Grid<double> >::const_iterator
    it = continuous_parameters_.find(parameter_name);
  assert(it != continuous_parameters_.end());

  return it->second;
}


bool EclipseGrid::HasParameter(const std::string& parameter_name) const
{
  std::map<std::string, Grid<double> >::const_iterator
    it = continuous_parameters_.find(parameter_name);

  if (it == continuous_parameters_.end())
    return false;
  return true;
}

} // namespace NRLib

#endif
