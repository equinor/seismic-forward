// $Id: stringtools_test.cpp 1496 2017-05-23 12:02:22Z perroe $

/// \file Unit tests for the Stringtools functions in the NRLib IOTools library.

#include <nrlib/iotools/stringtools.hpp>
#include <boost/test/unit_test.hpp>

using namespace NRLib;

BOOST_AUTO_TEST_CASE(ExtensionsTest)
{

    std::string file = "f";
    file = AddExtension(file, "ex1");
    BOOST_CHECK_EQUAL(file, "f.ex1");

    file = ReplaceExtension(file, "ex2");
    BOOST_CHECK_EQUAL(file, "f.ex2");

}
