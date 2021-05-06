cmake_minimum_required(VERSION 2.6)

set(C_FILES
    boost_1_65_1/libs/filesystem/src/operations.cpp
    boost_1_65_1/libs/filesystem/src/path.cpp
    boost_1_65_1/libs/filesystem/src/portability.cpp
    boost_1_65_1/libs/system/src/error_code.cpp
)

set(H_FILES
    boost_1_65_1/boost/filesystem/config.hpp
    boost_1_65_1/boost/filesystem/convenience.hpp
    boost_1_65_1/boost/filesystem/fstream.hpp
    boost_1_65_1/boost/filesystem/operations.hpp
    boost_1_65_1/boost/filesystem/path.hpp
    boost_1_65_1/boost/system/config.hpp
    boost_1_65_1/boost/system/cygwin_error.hpp
    boost_1_65_1/boost/system/error_code.hpp
    boost_1_65_1/boost/system/linux_error.hpp
    boost_1_65_1/boost/system/system_error.hpp
    boost_1_65_1/boost/system/windows_error.hpp
    boost_1_65_1/boost/cerrno.hpp
    boost_1_65_1/boost/filesystem.hpp
)


append_source_files_to_property(${CMAKE_CURRENT_SOURCE_DIR} "${C_FILES}")
append_source_files_to_property(${CMAKE_CURRENT_SOURCE_DIR} "${H_FILES}")
add_include_root(${CMAKE_CURRENT_SOURCE_DIR})
