cmake_minimum_required(VERSION 2.8.10) # Due to $<CONFIG:>

project(mex)

find_package(MatlabMex REQUIRED)
execute_process(COMMAND ${MATLAB_MEX}ext OUTPUT_VARIABLE MEX_EXTENSION)

# -----------------------------------------------------------------
# CLION IDE AUTO-COMPLETE SUPPORT
# -----------------------------------------------------------------
set(SOURCE_FILES
    class_handle.hpp
    g2s_model.cpp
    surface_container.cpp
)

include_directories(${MATLAB_ROOT}/extern/include)
add_library(ide_mex_library ${SOURCE_FILES})
# -----------------------------------------------------------------

function(mex_compile MEX_SOURCE_FILE LINK_LIBRARIES DEPENDENCY)
    if(NOT WIN32)
        set(MEX_CFLAGS -fPIC)
        set(MEX_CXXFLAGS -fPIC)
    endif(NOT WIN32)
    set(MEX_LDFLAGS)
    set(MEX_INCLUDES)

    get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
    foreach(dir ${dirs})
      set(MEX_INCLUDES ${MEX_INCLUDES} -I${dir})
    endforeach()

    if(APPLE)
        set(MACOSX_DEPLOYMENT_TARGET 10.5)
        set(ARCHS x86_64)
        set(MEX_CFLAGS ${MEX_CFLAGS} -fno-common -fexceptions -arch ${ARCHS} -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET})
        set(MEX_CXXFLAGS ${MEX_CXXFLAGS} -fno-common -fexceptions -arch ${ARCHS} -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET})
        set(MEX_LDFLAGS -Wl,-twolevel_namespace -undefined error -arch ${ARCHS} -bundle -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET} -lstdc++ "${LINK_LIBRARIES}")

        add_custom_target(mex_compile_${MEX_SOURCE_FILE} ALL
                          COMMAND
                                  ${MATLAB_MEX}
                                  CC="${CMAKE_C_COMPILER}"
                                  CXX="${CMAKE_CXX_COMPILER}"
                                  LD="${CMAKE_CXX_COMPILER}"
                                  CFLAGS="${MEX_CFLAGS}"
                                  CXXFLAGS="${MEX_CXXFLAGS}"
                                  LDFLAGS="${MEX_LDFLAGS}"
                                  ${MEX_INCLUDES}
                                  -outdir .
                                  ${CMAKE_CURRENT_LIST_DIR}/${MEX_SOURCE_FILE}.cpp
                          DEPENDS
                                  ${CMAKE_CURRENT_LIST_DIR}/${MEX_SOURCE_FILE}.cpp)
    elseif(WIN32)
        # Does probably not work with older versions of Matlab.
        set(MEX_OPT_PATH ${MATLAB_ROOT}/bin/win64/mexopts)
        if(MSVC90)
            set(MEX_OPT_FILE ${MEX_OPT_PATH}/msvcpp2008.xml)
        elseif(MSVC10)
            set(MEX_OPT_FILE ${MEX_OPT_PATH}/msvcpp2010.xml)
        elseif(MSVC11)
            set(MEX_OPT_FILE ${MEX_OPT_PATH}/msvcpp2012.xml)
        elseif(MSVC12)
            set(MEX_OPT_FILE ${MEX_OPT_PATH}/msvcpp2013.xml)
        endif(MSVC90)

        set(MEX_OPT_FLAG)
        if()
            set(MEX_OPT_FLAG "-g")
        endif()
        
        add_custom_target(mex_compile_${MEX_SOURCE_FILE} ALL
                          COMMAND
                                  ${MATLAB_MEX}
                                  -f ${MEX_OPT_FILE}
                                  $<$<CONFIG:Debug>:-g>
                                  $<$<CONFIG:Debug>:-D_DEBUG>
                                  ${MEX_INCLUDES}
                                  -outdir .
                                  ${CMAKE_CURRENT_LIST_DIR}/${MEX_SOURCE_FILE}.cpp
                                  ${LINK_LIBRARIES}
                          DEPENDS
                                  ${CMAKE_CURRENT_LIST_DIR}/${MEX_SOURCE_FILE}.cpp)
    else() #Linux
        set(MEX_LDFLAGS -shared)
        add_custom_target(mex_compile_${MEX_SOURCE_FILE} ALL
                          COMMAND
                                  ${MATLAB_MEX}
                                  CC="${CMAKE_C_COMPILER}"
                                  CXX="${CMAKE_CXX_COMPILER}"
                                  LD="${CMAKE_CXX_COMPILER}"
                                  CFLAGS="${MEX_CFLAGS}"
                                  CXXFLAGS="${MEX_CXXFLAGS}"
                                  LDFLAGS="${MEX_LDFLAGS}"
                                  ${MEX_INCLUDES}
                                  -outdir .
                                  ${CMAKE_CURRENT_LIST_DIR}/${MEX_SOURCE_FILE}.cpp
                                  ${LINK_LIBRARIES}
                          DEPENDS
                                  ${CMAKE_CURRENT_LIST_DIR}/${MEX_SOURCE_FILE}.cpp)
    endif(APPLE)

    add_dependencies(mex_compile_${MEX_SOURCE_FILE} ${DEPENDENCY})
endfunction()

set(MEX_LINK_LIBRARIES
    $<TARGET_FILE:geo2seis_lib>
    ${FFT_MATLAB_LIBRARIES}
    $<$<CONFIG:Debug>:${TBB_LIBRARIES_DEBUG}>
    $<$<NOT:$<CONFIG:Debug>>:${TBB_LIBRARIES_RELEASE}>
)

mex_compile(g2s_model "${MEX_LINK_LIBRARIES}" geo2seis_lib)
mex_compile(surface_container "${MEX_LINK_LIBRARIES}" geo2seis_lib)

set(MATLAB_SOURCE_FILES G2SModel.m ElasticParameters.m Angle.m Wavelet.m WhiteNoise.m OutputGrid.m SurfaceContainer.m)

add_custom_target(copy_m_files ALL)
foreach(MATLAB_SOURCE_FILE ${MATLAB_SOURCE_FILES})
    add_custom_command(TARGET copy_m_files PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_LIST_DIR}/${MATLAB_SOURCE_FILE} .)
endforeach()


