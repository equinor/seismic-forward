cmake_minimum_required(VERSION 2.6)

append_source_files_to_property(${CMAKE_CURRENT_SOURCE_DIR} "${C_FILES}")
append_source_files_to_property(${CMAKE_CURRENT_SOURCE_DIR} "${H_FILES}")
add_include_root(${CMAKE_CURRENT_SOURCE_DIR})
