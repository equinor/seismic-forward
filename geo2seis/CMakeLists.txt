cmake_minimum_required(VERSION 2.8)

set(CPP_FILES
    modelsettings.cpp
    seismic_forward.cpp
    seismic_geometry.cpp
    seismic_output.cpp
    seismic_parameters.cpp
    seismic_regridding.cpp
    xmlmodelfile.cpp

    physics/wavelet.cpp
    physics/zoeppritz.cpp
    physics/zoeppritz_pp.cpp
    physics/zoeppritz_ps.cpp

    utils/segy_writer.cpp
    utils/storm_writer.cpp
)

append_source_files_to_property(${CMAKE_CURRENT_SOURCE_DIR} "${CPP_FILES}")
add_include_root(${CMAKE_CURRENT_SOURCE_DIR})