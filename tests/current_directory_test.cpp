#include <cutest.h>
#include "iostream"
#include <sys/param.h>
#include <sstream>

#define BOOST_FILESYSTEM_VERSION 2
#include <boost/filesystem.hpp>

std::string getWorkingPath() {
    char temp[PATH_MAX];

    if (getcwd(temp, PATH_MAX) != 0) {
        return std::string(temp);
    }

    int error = errno;

    switch (error) {
        // EINVAL can't happen - size argument > 0

        // PATH_MAX includes the terminating nul,
        // so ERANGE should not be returned

        case EACCES:
            throw std::runtime_error("Access denied");

        case ENOMEM:
            // I'm not sure whether this can happen or not
            throw std::runtime_error("Insufficient storage");

        default: {
            std::ostringstream str;
            str << "Unrecognised error" << error;
            throw std::runtime_error(str.str());
        }
    }
}


void test_boost_current_path() {
    std::string current_path = getWorkingPath();
    boost::filesystem::path my_current_path(current_path);

    TEST_CHECK(boost::filesystem::is_directory(my_current_path));

    boost::filesystem::path test_file("nr.xml");

    TEST_CHECK(boost::filesystem::is_regular_file(test_file));

    TEST_CHECK(boost::filesystem::exists(test_file));

    TEST_CHECK(test_file.extension() == ".xml");

    boost::filesystem::path same_test_file = my_current_path / "nr.xml";

    TEST_CHECK(boost::filesystem::equivalent(test_file, same_test_file));

}


TEST_LIST = {
        {"test_boost_current_path", test_boost_current_path},
        {0}
};

