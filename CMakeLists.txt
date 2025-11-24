cmake_minimum_required(VERSION 3.18)
project(ASC_ODE_BLA LANGUAGES CXX)

# --- Skip host-Python probing when cross-compiling to Emscripten ---
if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten" OR CMAKE_CROSSCOMPILING)
  message(STATUS "Cross-compiling detected: disabling PYBIND11_FINDPYTHON to avoid host-Python vs cross-compiler checks.")
  set(PYBIND11_FINDPYTHON OFF CACHE BOOL "Do not probe host python when cross-compiling")
endif()
# --- end guard ---

include(FetchContent)

# Fetch pybind11 for cross-builds or fallback
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
  message(STATUS "pybind11 not found via find_package — fetching with FetchContent")
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.10.4
  )
  FetchContent_MakeAvailable(pybind11)
endif()


# Allow caller to provide Python_EXECUTABLE on the cmake command line
#if(NOT DEFINED Python_EXECUTABLE)
#  find_program(Python_EXECUTABLE NAMES python3 python)
#endif()

# Detect cross-compilation for Emscripten
if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
  set(BUILDING_FOR_EMSCRIPTEN TRUE)
else()
  set(BUILDING_FOR_EMSCRIPTEN FALSE)
endif()

if(NOT BUILDING_FOR_EMSCRIPTEN)
  message(STATUS "Native build: locating Python interpreter + development files")
  # Native build: require interpreter + development (headers / embedding)
  find_package(Python COMPONENTS Interpreter Development REQUIRED)
else()
  message(STATUS "Cross build (Emscripten): using Python executable: ${Python_EXECUTABLE}")
  if(NOT Python_EXECUTABLE)
    message(FATAL_ERROR "Please provide -DPython_EXECUTABLE=/path/to/python when cross-compiling")
  endif()
  # Avoid attempting to link host Python libraries into wasm target:
  set(Python_EXECUTABLE ${Python_EXECUTABLE} CACHE FILEPATH "Host python executable")
  set(Python_INCLUDE_DIRS "" CACHE STRING "Empty when targeting Emscripten")
  set(Python_LIBRARIES "" CACHE STRING "Empty when targeting Emscripten")
endif()

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# Prefer vendored pybind11 if available; otherwise try find_package.
if(EXISTS "${CMAKE_SOURCE_DIR}/third_party/pybind11/CMakeLists.txt")
  add_subdirectory(${CMAKE_SOURCE_DIR}/third_party/pybind11)
else()
  find_package(pybind11 CONFIG QUIET)
  if(NOT pybind11_FOUND)
    message(STATUS "pybind11 not found with find_package; consider adding pybind11 as a submodule (third_party/pybind11) or installing the pybind11 python package in the build venv.")
  endif()
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions(-DDEBUG_MODE)
endif()

include_directories(bla/src)
include_directories(bla/HPC/src)
include_directories(bla/HPC/concurrentqueue)

# Include subdirectories
add_subdirectory(helper)
add_subdirectory(Elements)
add_subdirectory(nonlinfunc)
add_subdirectory(src)
include_directories(helper)

# Add the main executable
add_executable(main main.cc)

# Link libraries to the main executable
target_link_libraries(main PRIVATE helper Elements nonlinfunc src)

# pybind11 module
pybind11_add_module(rigid_body_FEM bind.cc)
target_link_libraries(rigid_body_FEM PUBLIC Elements nonlinfunc src)

# (Optional) any extra packaging / install logic follows...
