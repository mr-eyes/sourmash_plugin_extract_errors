project(extract_errors)
cmake_minimum_required(VERSION 3.15...3.27)
find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)

# Warn if the user invokes CMake directly
if (NOT SKBUILD)
message(WARNING "\
This CMake file is meant to be executed using 'scikit-build-core'.
Running it directly will almost certainly not produce the desired
result. If you are a user trying to install this package, use the
command below, which will install all necessary build dependencies,
compile the package in an isolated environment, and then install it.
=====================================================================
$ pip install .
=====================================================================
If you are a software developer, and this is your own package, then
it is usually much more efficient to install the build dependencies
in your environment once and use the following command that avoids
a costly creation of a new virtual environment at every compilation:
=====================================================================
$ pip install nanobind scikit-build-core[pyproject]
$ pip install --no-build-isolation -ve .
=====================================================================
You may optionally add -Ceditable.rebuild=true to auto-rebuild when
the package is imported. Otherwise, you need to rerun the above
after editing C++ files.")
endif()
cmake_minimum_required(VERSION 3.15...3.26)

project(nanobind_example LANGUAGES CXX)

if (NOT SKBUILD)
  message(WARNING "\
  This CMake file is meant to be executed using 'scikit-build'. Running
  it directly will almost certainly not produce the desired result. If
  you are a user trying to install this package, please use the command
  below, which will install all necessary build dependencies, compile
  the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ pip install --no-build-isolation -ve .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to re-run the above
  after editing C++ files.")
endif()



find_package(nanobind CONFIG REQUIRED)


set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -lstdc++fs -fPIC -lgomp -lrt -fopenmp -Ofast")

add_subdirectory(lib/json)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/simdjson)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/parallel-hashmap)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/json/include)


nanobind_add_module(
  _extract_errors_impl

  STABLE_ABI

  NB_STATIC

  src/extract_errors.cpp 
  lib/simdjson/simdjson.cpp
)

target_link_libraries(_extract_errors_impl PRIVATE nlohmann_json::nlohmann_json)

install(TARGETS _extract_errors_impl LIBRARY DESTINATION extract_errors)


