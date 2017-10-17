// $Id: segyio_test.cpp 1688 2017-09-13 08:48:53Z perroe $

/// \file Unit tests testing reading and writing of SEGY files.

#include "nrlib/segy/segy.hpp"
#include "nrlib/segy/segygeometry.hpp"
#include "unittests/util.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>


using namespace boost::filesystem;
using namespace NRLib;


BOOST_AUTO_TEST_CASE(SegyReadGeometryTest)
{
  path input_filename = GetTestDir() / "lnea32_mr_jura_IL_3461_4410_XL_1617_2378.segy";
  SegY in_file(input_filename.string(), 0);

  BOOST_CHECK_EQUAL(in_file.GetNz(), 542U);
  BOOST_CHECK_EQUAL(in_file.GetDz(), 4.0F);
  BOOST_CHECK_EQUAL(in_file.GetTop(), 0.0F);

  in_file.FindAndSetGridGeometry(false, true);

  BOOST_CHECK_EQUAL(in_file.GetNTraces(), 180975U);

  TraceHeaderFormat format = in_file.GetTraceHeaderFormat();
  BOOST_CHECK(format.GetStandardType());
  BOOST_CHECK_EQUAL(format.GetFormatName(), "SeisWorks");
  BOOST_CHECK_EQUAL(format.GetUtmxLoc(), 73);
  BOOST_CHECK_EQUAL(format.GetUtmyLoc(), 77);
  BOOST_CHECK_EQUAL(format.GetInlineLoc(), 9);
  BOOST_CHECK_EQUAL(format.GetCrosslineLoc(), 21);
  BOOST_CHECK_EQUAL(format.GetStartTimeLoc(), 109);
  BOOST_CHECK_EQUAL(format.GetScalCoLoc(), 71);
  BOOST_CHECK_EQUAL(format.GetBypassCoordScaling(), false);
  BOOST_CHECK_LT(format.GetOffsetLoc(), 0);
  BOOST_CHECK_EQUAL(format.GetCoordSys(), TraceHeaderFormat::UTM);
  BOOST_CHECK_NO_THROW(format.CheckFormat());

  const SegyGeometry* geometry = in_file.GetGeometry();
  BOOST_CHECK_EQUAL(geometry->GetNx(), 381U);
  BOOST_CHECK_EQUAL(geometry->GetNy(), 475U);
  BOOST_CHECK_CLOSE(geometry->GetDx(), 25.0, 0.01);
  BOOST_CHECK_CLOSE(geometry->GetDy(), 25.0, 0.01);
  BOOST_CHECK_CLOSE(geometry->GetX0(), 397823.17, 0.00001);
  BOOST_CHECK_CLOSE(geometry->GetY0(), 7230211.16, 0.00001);
  BOOST_CHECK_CLOSE(geometry->Getlx(), 9525.0, 0.001);
  BOOST_CHECK_CLOSE(geometry->Getly(), 11875.0, 0.001);
  BOOST_CHECK_CLOSE(geometry->GetAngle(), -0.41243, 0.001);
  BOOST_CHECK_CLOSE(geometry->GetCosRot(), 0.916148, 0.001);
  BOOST_CHECK_CLOSE(geometry->GetSinRot(), -0.400840, 0.001);
  BOOST_CHECK_CLOSE(geometry->GetInLine0(), 3461.0, 0.00001);
  BOOST_CHECK_CLOSE(geometry->GetCrossLine0(), 1617.0, 0.00001);
  BOOST_CHECK_CLOSE(geometry->GetILStepX(), 0.0320675, 0.001);
  BOOST_CHECK_CLOSE(geometry->GetILStepY(), 0.0732919, 0.001);
  BOOST_CHECK_CLOSE(geometry->GetXLStepX(), 0.0732919, 0.001);
  BOOST_CHECK_CLOSE(geometry->GetXLStepY(), -0.0320663, 0.001);
  BOOST_CHECK_EQUAL(geometry->GetMinIL(), 3462);
  BOOST_CHECK_EQUAL(geometry->GetMaxIL(), 4410);
  BOOST_CHECK_EQUAL(geometry->GetILStep(), 2);
  BOOST_CHECK_EQUAL(geometry->GetMinXL(), 1618);
  BOOST_CHECK_EQUAL(geometry->GetMaxXL(), 2378);
  BOOST_CHECK_EQUAL(geometry->GetXLStep(), 2);
  BOOST_CHECK_EQUAL(geometry->GetIL0(), 3462);
  BOOST_CHECK_EQUAL(geometry->GetXL0(), 1618);
  BOOST_CHECK(!geometry->GetFirstAxisIL());

  std::streampos pos1 = in_file.GetFilePos(3626, 1962);
  std::streampos pos2 = in_file.GetFilePos(402605.80, 7230376.84);

  BOOST_CHECK_EQUAL(pos1, pos2);

  std::vector<float> trace1, trace2, trace3;

  in_file.GetTraceData(pos1, trace1);
  in_file.GetTraceData(3626, 1962, trace2);
  in_file.GetTraceData(402605.80, 7230376.84, trace3);

  BOOST_CHECK_CLOSE(trace1[40], -2336.68, 0.01);
  BOOST_CHECK_EQUAL(trace1.size(), trace2.size());
  BOOST_CHECK_EQUAL(trace1.size(), trace3.size());

  for (size_t i = 0; i < trace1.size(); ++i) {
    BOOST_CHECK_EQUAL(trace1[i], trace2[i]);
    BOOST_CHECK_EQUAL(trace1[i], trace3[i]);
  }

  TraceHeader trace_header = in_file.GetTraceHeader(3626, 1962);
  BOOST_CHECK_EQUAL(trace_header.GetUtmx(), 402600.80);
  BOOST_CHECK_EQUAL(trace_header.GetUtmy(), 7230371.84);
  BOOST_CHECK_EQUAL(trace_header.GetInline(), 3626);
  BOOST_CHECK_EQUAL(trace_header.GetCrossline(), 1962);
  BOOST_CHECK_EQUAL(trace_header.GetStartTime(), 0.0f);
  BOOST_CHECK_EQUAL(trace_header.GetCoord1(), 402600.80);
  BOOST_CHECK_EQUAL(trace_header.GetCoord2(), 7230371.84);
  BOOST_CHECK_EQUAL(trace_header.GetDt(), 4000);
  BOOST_CHECK_EQUAL(trace_header.GetStatus(), 0);

}


BOOST_AUTO_TEST_CASE(SegyReadGeometryMinimalTraceheaderTest)
{
  path input_filename = GetTestDir() / "lnea32_mr_jura_IL_3461_4410_XL_1617_2378.segy";

  TraceHeaderFormat thf(-1, -1, -1, -1, -1, -1, TraceHeaderFormat::ILXL);

  SegY in_file(input_filename.string(), 0, thf);

  SegYTrace* trace = in_file.GetNextTrace();

  BOOST_CHECK_CLOSE(trace->GetTrace()[15], 6764.14, 0.01);

  // For some reason the constructor above constructs a standard header type...
  // BOOST_CHECK(!thf.GetStandardType());
}


BOOST_AUTO_TEST_CASE(SegyReadAllTracesTest)
{
  path input_filename = GetTestDir() / "lnea32_mr_jura_IL_3461_4410_XL_1617_2378.segy";

  SegY in_file(input_filename.string(), 0);

  in_file.ReadAllTraces(NULL, 0);

  in_file.CreateRegularGrid(true);

  std::vector<float> trace;
  bool missing = true;
  double z0_data = 0;
  in_file.GetNearestTrace(trace, missing, z0_data, 402600.80, 7230371.84);
  BOOST_CHECK_CLOSE(trace[40], -2336.68, 0.01);

  double x, y;
  in_file.GetGeometry()->FindXYFromILXL(3626, 1962, x, y);
  float val = in_file.GetValue(x, y, 162.0);

  BOOST_CHECK_CLOSE(val, -2336.68, 0.03);
}


BOOST_AUTO_TEST_CASE(SegyReadWriteTest)
{
  path input_filename = GetTestDir() / "lnea32_mr_jura_IL_3461_4410_XL_1617_2378.segy";
  path out_path = temp_directory_path() / "test_write_trace_by_trace.sgy";

  SegY in_file_orig(input_filename.string(), 0);

  TextualHeader test_header;

  {
    // Reduce scope, to ensure that out_file1 is closed before we read from it.
    SegY out_file1(out_path.string(), 0.0F, in_file_orig.GetNz(), in_file_orig.GetDz(),
                   test_header, in_file_orig.GetTraceHeaderFormat());

    SegYTrace* trace = in_file_orig.GetNextTrace();
    while (trace != NULL) {
      out_file1.WriteTrace(trace->GetTraceHeader(), trace->GetTrace(), NULL);
      delete trace;
      trace = in_file_orig.GetNextTrace();
    }
  }

  {
    SegY in_file(out_path.string(), 0);

    BOOST_CHECK_EQUAL(in_file.GetNz(), 542U);
    BOOST_CHECK_EQUAL(in_file.GetDz(), 4.0F);
    BOOST_CHECK_EQUAL(in_file.GetTop(), 0.0F);

    in_file.FindAndSetGridGeometry(false, true);

    BOOST_CHECK_EQUAL(in_file.GetNTraces(), 180975U);

    TraceHeaderFormat format = in_file.GetTraceHeaderFormat();
    BOOST_CHECK(format.GetStandardType());
    BOOST_CHECK_EQUAL(format.GetFormatName(), "SeisWorks");
    BOOST_CHECK_EQUAL(format.GetUtmxLoc(), 73);
    BOOST_CHECK_EQUAL(format.GetUtmyLoc(), 77);
    BOOST_CHECK_EQUAL(format.GetInlineLoc(), 9);
    BOOST_CHECK_EQUAL(format.GetCrosslineLoc(), 21);
    BOOST_CHECK_EQUAL(format.GetStartTimeLoc(), 109);
    BOOST_CHECK_EQUAL(format.GetScalCoLoc(), 71);
    BOOST_CHECK_EQUAL(format.GetBypassCoordScaling(), false);
    BOOST_CHECK_LT(format.GetOffsetLoc(), 0);
    BOOST_CHECK_EQUAL(format.GetCoordSys(), TraceHeaderFormat::UTM);
    BOOST_CHECK_NO_THROW(format.CheckFormat());


    // NB: We loose precision after writing and reading again.
    // See commented-out tests below.
    // This should be examined!

    const SegyGeometry* geometry = in_file.GetGeometry();
    BOOST_CHECK_EQUAL(geometry->GetNx(), 381U);
    BOOST_CHECK_EQUAL(geometry->GetNy(), 475U);
    BOOST_CHECK_CLOSE(geometry->GetDx(), 25.0, 0.01);
    BOOST_CHECK_CLOSE(geometry->GetDy(), 25.0, 0.01);
    BOOST_CHECK_CLOSE(geometry->GetX0(), 397823.17, 0.00001);
    BOOST_CHECK_CLOSE(geometry->GetY0(), 7230211.16, 0.00001);
    BOOST_CHECK_CLOSE(geometry->Getlx(), 9525.0, 0.001);
    BOOST_CHECK_CLOSE(geometry->Getly(), 11875.0, 0.001);
    BOOST_CHECK_CLOSE(geometry->GetAngle(), -0.41243, 0.001);
    BOOST_CHECK_CLOSE(geometry->GetCosRot(), 0.916148, 0.001);
    BOOST_CHECK_CLOSE(geometry->GetSinRot(), -0.400840, 0.001);
    BOOST_CHECK_CLOSE(geometry->GetInLine0(), 3461.0, 0.00001);
    BOOST_CHECK_CLOSE(geometry->GetCrossLine0(), 1617.0, 0.00001);
    BOOST_CHECK_CLOSE(geometry->GetILStepX(), 0.0320675, 0.001);
    BOOST_CHECK_CLOSE(geometry->GetILStepY(), 0.0732919, 0.001);
    BOOST_CHECK_CLOSE(geometry->GetXLStepX(), 0.0732919, 0.001);
    BOOST_CHECK_CLOSE(geometry->GetXLStepY(), -0.0320663, 0.001);
    BOOST_CHECK_EQUAL(geometry->GetMinIL(), 3462);
    BOOST_CHECK_EQUAL(geometry->GetMaxIL(), 4410);
    BOOST_CHECK_EQUAL(geometry->GetILStep(), 2);
    BOOST_CHECK_EQUAL(geometry->GetMinXL(), 1618);
    BOOST_CHECK_EQUAL(geometry->GetMaxXL(), 2378);
    BOOST_CHECK_EQUAL(geometry->GetXLStep(), 2);
    BOOST_CHECK_EQUAL(geometry->GetIL0(), 3462);
    BOOST_CHECK_EQUAL(geometry->GetXL0(), 1618);
    BOOST_CHECK(!geometry->GetFirstAxisIL());

    std::streampos pos1 = in_file.GetFilePos(3626, 1962);
    std::streampos pos2 = in_file.GetFilePos(402605.80, 7230376.84);

    BOOST_CHECK_EQUAL(pos1, pos2);

    std::vector<float> trace1, trace2, trace3;

    in_file.GetTraceData(pos1, trace1);
    in_file.GetTraceData(3626, 1962, trace2);
    in_file.GetTraceData(402605.80, 7230376.84, trace3);

    BOOST_CHECK_EQUAL(trace1.size(), 542U);
    BOOST_CHECK_CLOSE(trace1[0], -1929.11, 0.01);
    BOOST_CHECK_CLOSE(trace1[40], -2336.68, 0.01);
    BOOST_CHECK_CLOSE(trace1[540], 751.809, 0.01);
    BOOST_CHECK_CLOSE(trace1[541], 0.0, 0.01);

    BOOST_CHECK_EQUAL(trace1.size(), trace2.size());
    BOOST_CHECK_EQUAL(trace1.size(), trace3.size());

    for (size_t i = 0; i < trace1.size(); ++i) {
      BOOST_CHECK_EQUAL(trace1[i], trace2[i]);
      BOOST_CHECK_EQUAL(trace1[i], trace3[i]);
    }

    TraceHeader trace_header = in_file.GetTraceHeader(3626, 1962);
    BOOST_CHECK_EQUAL(trace_header.GetUtmx(), 402600.80);
    BOOST_CHECK_EQUAL(trace_header.GetUtmy(), 7230371.84);
    BOOST_CHECK_EQUAL(trace_header.GetInline(), 3626);
    BOOST_CHECK_EQUAL(trace_header.GetCrossline(), 1962);
    BOOST_CHECK_EQUAL(trace_header.GetStartTime(), 0.0f);
    BOOST_CHECK_EQUAL(trace_header.GetCoord1(), 402600.80);
    BOOST_CHECK_EQUAL(trace_header.GetCoord2(), 7230371.84);
    BOOST_CHECK_EQUAL(trace_header.GetDt(), 4000);
    BOOST_CHECK_EQUAL(trace_header.GetStatus(), 0);
  }

  remove(out_path);
}