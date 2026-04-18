cmake_minimum_required(VERSION 3.21)

option(ENABLE_CUDA "Enable CUDA GPU support" OFF)

if(ENABLE_CUDA)
    project(inq_tutorials CXX Fortran CUDA)
    set(ENABLE_GPU ON)
    set(GPU_LANGUAGE CUDA)
else()
    project(inq_tutorials CXX Fortran)
endif()

# Import INQ once — all tutorial subdirectories share this target.
add_subdirectory(../inq inq EXCLUDE_FROM_ALL)

add_subdirectory(hello-world)
add_subdirectory(n2)
add_subdirectory(li-bcc)
