cmake_minimum_required(VERSION 2.6)

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
	utils/output.cpp
	utils/result_trace.cpp
	utils/trace.cpp
	utils/gen_seis_trace_params.cpp
)

set(H_FILES
    modelsettings.hpp
    seismic_forward.hpp
    seismic_geometry.hpp
    seismic_output.hpp
    seismic_parameters.hpp
    seismic_regridding.hpp
    xmlmodelfile.hpp

    physics/wavelet.hpp
    physics/zoeppritz.hpp
    physics/zoeppritz_pp.hpp
    physics/zoeppritz_ps.hpp

    utils/segy_writer.hpp
    utils/storm_writer.hpp
	utils/output.hpp
	utils/result_trace.hpp
	utils/trace.hpp
	utils/gen_seis_trace_params.hpp
)


append_source_files_to_property(${CMAKE_CURRENT_SOURCE_DIR} "${CPP_FILES}")
append_source_files_to_property(${CMAKE_CURRENT_SOURCE_DIR} "${H_FILES}")
add_include_root(${CMAKE_CURRENT_SOURCE_DIR})