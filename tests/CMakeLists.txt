cmake_minimum_required(VERSION 2.6)

project(Testing)

set(SOURCE_FILES
    main.cpp
	compare_traces.cpp
    test_test.cpp)

enable_testing()

# Turn of Boost auto-linking, due to link trouble.
set (Boost_USE_STATIC_LIBS TRUE)
set (Boost_NO_BOOST_CMAKE TRUE) # BOOST_CMAKE gives trouble on Redhat6.
# add_definitions( -DBOOST_ALL_NO_LIB )

# Custom Boost Library dir on Windows
if(WIN32)
    if(NOT BOOST_LIBRARYDIR)
        set(BOOST_LIBRARYDIR $ENV{BOOST_ROOT}/lib/x64)
    endif(NOT BOOST_LIBRARYDIR)
endif(WIN32)

find_package(Boost 1.35.0 COMPONENTS unit_test_framework REQUIRED)
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIRS})
# link_directories(${Boost_LIBRARY_DIRS})

add_executable(test_geo2seis ${SOURCE_FILES})
target_link_libraries(test_geo2seis geo2seis_lib ${MATH_LIBRARIES} ${Boost_LIBRARIES})

set (DATA_DIR ${CMAKE_SOURCE_DIR}/tests/data)
add_test(NAME test_geo2seis WORKING_DIRECTORY ${DATA_DIR} COMMAND test_geo2seis)
