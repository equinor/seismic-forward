cmake_minimum_required(VERSION 2.8)

project(Testing)

enable_testing()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

function(build_and_test TEST_NAME SOURCE_FILES INCLUDES LINK_LIBRARIES)
    add_executable(${TEST_NAME}_exe ${SOURCE_FILES} ${INCLUDES})
    target_link_libraries(${TEST_NAME}_exe ${LINK_LIBRARIES})

    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME}_exe WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../testdata)

endfunction()

set(LINK_LIBRARIES ${MATH_LIBRARIES})

build_and_test(cutest test_cutest.cpp "" "${LINK_LIBRARIES}")
build_and_test(cwd current_directory_test.cpp "$<TARGET_OBJECTS:geo2seis_object>" "${LINK_LIBRARIES}")

