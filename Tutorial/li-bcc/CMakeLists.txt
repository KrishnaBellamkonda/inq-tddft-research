add_executable(li_bcc li_bcc.cpp)
target_link_libraries(li_bcc PRIVATE inq)
set_target_properties(li_bcc PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

if(ENABLE_GPU)
    set_source_files_properties(li_bcc.cpp PROPERTIES LANGUAGE ${GPU_LANGUAGE})
endif()
