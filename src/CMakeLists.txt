cmake_minimum_required(VERSION 3.5)

set(CPP_FILES
    modelsettings.cpp
    seismic_forward.cpp
    seismic_geometry.cpp
    seismic_output.cpp
    seismic_parameters.cpp
    seismic_regridding.cpp
    tasklist.cpp
    wavelet.cpp
    xmlmodelfile.cpp

    utils/storm_writer.cpp
    utils/output.cpp
    utils/result_trace.cpp
    utils/trace.cpp
    utils/gen_seis_trace_params.cpp
    utils/gen_resampl_param.cpp
    utils/resampl_trace.cpp
    utils/resampl_output.cpp
)

set(H_FILES
    modelsettings.hpp
    seismic_forward.hpp
    seismic_geometry.hpp
    seismic_output.hpp
    seismic_parameters.hpp
    seismic_regridding.hpp
    tasklist.hpp
    wavelet.hpp
    xmlmodelfile.hpp

    utils/storm_writer.hpp
    utils/output.hpp
    utils/result_trace.hpp
    utils/trace.hpp
    utils/gen_seis_trace_params.hpp
    utils/gen_resampl_param.hpp
    utils/resampl_trace.hpp
    utils/resampl_output.hpp
)


append_source_files_to_property(${CMAKE_CURRENT_SOURCE_DIR} "${CPP_FILES}")
append_source_files_to_property(${CMAKE_CURRENT_SOURCE_DIR} "${H_FILES}")
add_include_root(${CMAKE_CURRENT_SOURCE_DIR})
