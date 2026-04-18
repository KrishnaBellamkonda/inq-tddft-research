add_executable(n2 n2.cpp)
target_link_libraries(n2 PRIVATE inq)
set_target_properties(n2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

if(ENABLE_GPU)
    set_source_files_properties(n2.cpp PROPERTIES LANGUAGE ${GPU_LANGUAGE})
endif()
