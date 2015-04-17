function(append_source_files_to_property CURRENT_DIRECTORY FILE_LIST)
    set(FILES)
    foreach (filename ${FILE_LIST})
        list(APPEND FILES "${CURRENT_DIRECTORY}/${filename}")
    endforeach()
    set_property(GLOBAL APPEND PROPERTY SOURCE_FILES ${FILES})
endfunction()


function(add_include_root INCLUDE_ROOT)
    set_property(GLOBAL APPEND PROPERTY INCLUDE_ROOTS ${INCLUDE_ROOT})
endfunction()
