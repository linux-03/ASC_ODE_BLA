cmake_minimum_required(VERSION 3.15...3.27)
project(ASC_ODE_BLA)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

find_package(pybind11)

include_directories(bla/src)
include_directories(bla/HPC/src)
include_directories(bla/HPC/concurrentqueue)

#include_directories(bla/src)
#include_directories(bla/HPC/src)
#include_directories(bla/HPC/concurrentqueue)

# Include subdirectories
add_subdirectory(helper)
add_subdirectory(Elements)
add_subdirectory(nonlinfunc)
add_subdirectory(src)
include_directories(helper)

# Add the main executable
add_executable(main main.cc)

# Link all necessary libraries to the main executable
target_link_libraries(main PRIVATE helper Elements nonlinfunc src)

pybind11_add_module(my_module bind.cc)

target_link_libraries(my_module PUBLIC Elements nonlinfunc src)


# if invoked using scikit-build-core
#if(DEFINED SKBUILD)
#    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
#    find_package(pybind11 CONFIG REQUIRED)

    #rigid_body_fem
#    python_add_library(rigid_body_FEM MODULE rigid_body_FEM/bind_rigid_body_FEM.cc WITH_SOABI)
#    target_link_libraries(rigid_body_FEM PRIVATE pybind11::headers)

#    target_compile_definitions(rigid_body_FEM PRIVATE VERSION_INFO=${PROJECT_VERSION})


#    install(TARGETS rigid_body_FEM DESTINATION lib_rigid_body)

#else()
set (CXXFLAGS "-g -O0")

#add_subdirectory (src)

#add_subdirectory (rigid_body_FEM)
#endif()