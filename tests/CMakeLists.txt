cmake_minimum_required(VERSION 3.26)
project(Testing LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(SOURCE_FILES
    main.cpp
    compare_traces.cpp
    test_test.cpp)

enable_testing()

# Turn of Boost auto-linking, due to link trouble.
set (Boost_USE_STATIC_LIBS TRUE)

# Custom Boost Library dir on Windows
if(WIN32)
    if(NOT BOOST_LIBRARYDIR)
        set(BOOST_LIBRARYDIR $ENV{BOOST_ROOT}/lib/x64)
    endif(NOT BOOST_LIBRARYDIR)
endif(WIN32)

find_package(Boost COMPONENTS unit_test_framework REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIRS})

add_executable(test_seismic_forward ${SOURCE_FILES})
target_link_libraries(test_seismic_forward seismic-forward-lib ${MATH_LIBRARIES} ${Boost_LIBRARIES})

set (DATA_DIR ${CMAKE_SOURCE_DIR}/tests/data)
add_test(NAME test_seismic_forward WORKING_DIRECTORY ${DATA_DIR} COMMAND test_seismic_forward)
