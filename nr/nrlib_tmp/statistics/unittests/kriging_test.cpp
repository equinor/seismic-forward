// $Id: kriging_test.cpp 1306 2015-07-08 09:55:41Z veralh $

/// \file Unit tests for the FileIO functions in the NRLib IOTools library.

#include <nrlib/statistics/kriging.hpp>
#include <nrlib/variogram/variogram.hpp>
#include <nrlib/variogram/variogramtypes.hpp>
#include <nrlib/variogram/covgrid2d.hpp>
#include <nrlib/statistics/nodalkrigingdata2d.hpp>
#include <nrlib/grid/grid2d.hpp>
#include <nrlib/surface/regularsurface.hpp>
#include <nrlib/surface/surfaceio.hpp>

#include <nrlib/random/random.hpp>
#include <nrlib/random/randomgenerator.hpp>
#include <boost/test/unit_test.hpp>

using namespace NRLib;

BOOST_AUTO_TEST_CASE( Kriging1DOneDataTest )
{

  std::vector<double> field(3, -1.0);
  std::vector<bool>   is_known(3, false);
  is_known[1] = true;
  std::vector<double> obs(1, 1.0);
  double dx = 1.0;
  ExpVario vario(50);

  Kriging::Krig1D(field, is_known, obs, dx, vario);

  //double answer [] = { 0.88, 1.0, 0.88};

  BOOST_CHECK_CLOSE(field[0], 1.0, 14);
  BOOST_CHECK_EQUAL(field[1], 1.0);
  BOOST_CHECK_CLOSE(field[2], 1.0, 14);

}

BOOST_AUTO_TEST_CASE( Kriging1DTwoDataTest )
{

  std::vector<double> field(3, -1.0);
  std::vector<bool>   is_known(3, false);
  is_known[0] = true;
  is_known[2] = true;
  std::vector<double> obs(2, 1.0);
  double dx = 1.0;
  ExpVario vario(50);

  Kriging::Krig1D(field, is_known, obs, dx, vario);

  //double answer [] = { 0.88, 1.0, 0.88};

  BOOST_CHECK_CLOSE(field[0], 1.0, 14);
  BOOST_CHECK_CLOSE(field[1], 1.0, 14);
  BOOST_CHECK_CLOSE(field[2], 1.0, 14);

}

BOOST_AUTO_TEST_CASE( Kriging1DMoreDataTest )
{

  std::vector<double> field(10, -1.0);
  std::vector<bool>   is_known(10, false);
  is_known[0] = true;
  is_known[2] = true;
  is_known[7] = true;
  is_known[8] = true;
  std::vector<double> obs(4, 1.0);

  double dx = 1.0;
  ExpVario vario(50);

  Kriging::Krig1D(field, is_known, obs, dx, vario);
  //double answer [] = { 0.88, 1.0, 0.88, 0.88, 1.0, 0.88, 0.88, 1.0, 0.88, 1.0};

  std::vector<double> answer(10, 1.0);

  //CHECK_CLOSE_COLLECTION(field, answer, 15);

  for (int i = 0; i < 10; ++i){
    BOOST_CHECK_CLOSE(field[i], answer[i], 15);
  }

}

BOOST_AUTO_TEST_CASE( Kriging2DFirstTest )
{
  int nx = 100;
  int ny = nx;
  double range = 10.0;
  double dx = 1.0;
  double dy = dx;
  ExpVario vario(range, range);
  CovGrid2D cov(vario, nx, ny, dx, dy);
  Grid2D<double> trend(nx, ny, 20);
  //Grid2D<double> prediction_error(nx, ny, 0);

  int ndata = 1000;
  NodalKrigingData2D krigingdata(ndata);
  NRLib::Random::Initialize(635435);
  for (int i=0; i < ndata; ++i){
    int indexi = static_cast<int>(NRLib::Random::Unif01()*nx);
    int indexj = static_cast<int>(NRLib::Random::Unif01()*nx);
    if (indexi == nx) indexi = nx-1;
    if (indexj == nx) indexj = nx-1;
    krigingdata.AddData(indexi,indexj, static_cast<float>(i*0.0005+20));
  }
  krigingdata.FindMeanValues();

  Kriging::Krig2D(trend, krigingdata, cov);



  //WRITE DATA/SURFACES TO FILE
  //Grid2D<double> data(nx, ny, 0);
  //for (int i=0; i < krigingdata.GetNumberOfData(); ++i){
  //  size_t indexi = krigingdata.GetIndexI(i);
  //  size_t indexj = krigingdata.GetIndexJ(i);
  //  data(indexi, indexj) = krigingdata.GetData(i);
  //}
  //RegularSurface<double> surface_data(0, 0, nx, ny, data);
  //std::string filename_data = "C:/nrlib/data.irap";
  //WriteIrapClassicAsciiSurf(surface_data, 0.0, filename_data);

  RegularSurface<double> surface(0, 0, nx-1, ny-1, trend);
  std::string filename = "trend.irap";
  WriteIrapClassicAsciiSurf(surface, 0.0, filename);

  //RegularSurface<double> surface_pred(0, 0, nx, ny, prediction_error);
  //std::string filename0 = "C:/nrlib/prediction_error.irap";
  //WriteIrapClassicAsciiSurf(surface_pred, 0.0, filename0);

  BOOST_CHECK_CLOSE(trend(0,0), 20.090436, 1);
  BOOST_CHECK_CLOSE(trend(1,0), 20.159995, 1);
  BOOST_CHECK_CLOSE(trend(2,0), 20.226837, 1);
}


BOOST_AUTO_TEST_CASE( Kriging2DFirstTestContinuousInputIndexBasedKriging )
{
  int nx = 100;
  int ny = nx;
  double range = 10;
  ExpVario vario(range, range);
  RegularSurface<double> trend(0, 0, nx-1, ny-1, nx, ny, 20.0);

  // First make index-based data set, then convert it to a vector of points.
  int ndata = 1000;
  NodalKrigingData2D krigingdata(ndata);
  NRLib::Random::Initialize(635435);
  for (int i=0; i < ndata; ++i){
    int indexi = static_cast<int>(NRLib::Random::Unif01()*nx);
    int indexj = static_cast<int>(NRLib::Random::Unif01()*nx);
    if (indexi == nx) indexi = nx-1;
    if (indexj == nx) indexj = nx-1;
    krigingdata.AddData(indexi,indexj, static_cast<float>(i*0.0005+20));
  }
  krigingdata.FindMeanValues();

  std::vector<Point* > krigingdata_cont(krigingdata.GetNumberOfData());
  for(int i = 0; i < krigingdata.GetNumberOfData(); i++) {
    double x = static_cast<double>(krigingdata.GetIndexI(i));
    if (x > nx-1 ) x = nx-1;
    double y = static_cast<double>(krigingdata.GetIndexJ(i));
    if (y > ny-1 ) y = ny-1;
    krigingdata_cont[i] = new Point(x, y, static_cast<double>(krigingdata.GetData(i)));
  }

  bool use_index_grid = true;
  Kriging::Krig2D(trend, krigingdata_cont, vario, false, use_index_grid);

  std::string filename = "trend_continuousInput_indexbasedKriging.irap";
  WriteIrapClassicAsciiSurf(trend, 0.0, filename);

  // Release memory
  for(int i = 0; i < static_cast<int>(krigingdata_cont.size()); i++){
    delete krigingdata_cont[i];
  }

  BOOST_CHECK_CLOSE(trend(0,0), 20.090436, 0.0001);
  BOOST_CHECK_CLOSE(trend(1,0), 20.159995, 0.0001);
  BOOST_CHECK_CLOSE(trend(2,0), 20.226837, 0.0001);
}

BOOST_AUTO_TEST_CASE( Kriging2DFirstTestContinuousInputContinuousKriging )
{
  int nx = 100;
  int ny = nx;
  double range = 10.0;
  ExpVario vario(range, range);
  RegularSurface<double> trend(0, 0, nx-1, ny-1, nx, ny, 20.0);

  int ndata = 1000;
  NodalKrigingData2D krigingdata(ndata);
  NRLib::Random::Initialize(635435);
  for (int i=0; i < ndata; ++i){
    int indexi = static_cast<int>(NRLib::Random::Unif01()*nx);
    int indexj = static_cast<int>(NRLib::Random::Unif01()*nx);
    if (indexi == nx) indexi = nx-1;
    if (indexj == nx) indexj = nx-1;
    krigingdata.AddData(indexi,indexj, static_cast<float>(i*0.0005+20));
  }
  krigingdata.FindMeanValues();

  std::vector<Point* > krigingdata_cont(krigingdata.GetNumberOfData());
  for(int i = 0; i < krigingdata.GetNumberOfData(); i++) {
    double x = static_cast<double>(krigingdata.GetIndexI(i));
    if (x > static_cast<double>(nx-1) ) x = static_cast<double>(nx-1);
    double y = static_cast<double>(krigingdata.GetIndexJ(i));
    if (y > static_cast<double>(ny-1)) y = static_cast<double>(ny-1);
    krigingdata_cont[i] = new Point(x, y, static_cast<double>(krigingdata.GetData(i)));
  }

  bool use_index_grid = false;
  Kriging::Krig2D(trend, krigingdata_cont, vario, false, use_index_grid);

  std::string filename = "trend_continuousInput_continuousKriging.irap";
  WriteIrapClassicAsciiSurf(trend, 0.0, filename);

  // Release memory
  for(int i = 0; i < static_cast<int>(krigingdata_cont.size()); i++){
    delete krigingdata_cont[i];
  }


  BOOST_CHECK_CLOSE(trend(0,0), 20.090436, 0.0001);
  BOOST_CHECK_CLOSE(trend(1,0), 20.159995, 0.0001);
  BOOST_CHECK_CLOSE(trend(2,0), 20.226837, 0.0001);
}
