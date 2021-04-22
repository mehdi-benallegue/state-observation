# Copyright 2013, CNRS
#
# Author: Mehdi Benallegue
#
# This file is part of state-observation.
# state-observation is free software: you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# state-observation is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Lesser Public License for more details.  You should have
# received a copy of the GNU Lesser General Public License along with
# state-observation. If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.1.4)
if(POLICY CMP0063)
  cmake_policy(SET CMP0063 NEW)
endif()

set(PROJECT_NAME state-observation)
set(PROJECT_DESCRIPTION "General implementation of observers.")
set(PROJECT_URL "https://github.com/jrl-umi3218/state-observation")
set(PROJECT_VERSION 1.4.1)
set(CMAKE_CXX_STANDARD 11)
set(CXX_DISABLE_WERROR TRUE)
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set_property(GLOBAL PROPERTY USE_FOLDERS OFF)

include(cmake/base.cmake)
#INCLUDE(cmake/eigen.cmake)
include(cmake/msvc-specific.cmake)
list(APPEND CMAKE_MODULE_PATH "C:/devel/install/share/eigen3/cmake")
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

project(
  ${PROJECT_NAME}
  LANGUAGES CXX
  VERSION ${PROJECT_VERSION})

set(DOXYGEN_USE_MATHJAX "YES")
set(DOXYGEN_EXAMPLE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/doc")

option(BUILD_STATE_OBSERVATION_TOOLS "Build state-observation tools" ON)

# For MSVC, set local environment variable to enable finding the built dll
# of the main library when launching ctest with RUN_TESTS
if(MSVC)
  set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
endif(MSVC)

add_project_dependency(Boost REQUIRED COMPONENTS timer)
add_project_dependency(Eigen3 QUIET NO_CMAKE_PACKAGE_REGISTRY)
if(Eigen3_FOUND)
  add_project_dependency(Eigen3 REQUIRED NO_CMAKE_PACKAGE_REGISTRY)
else()
  add_project_dependency(Eigen3 MODULE REQUIRED)
endif()

file(
  GENERATE
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/include/state-observation/config.h"
  INPUT "${CMAKE_CURRENT_SOURCE_DIR}/include/state-observation/config.h.in")
install(FILES "${CMAKE_BINARY_DIR}/$<CONFIG>/include/state-observation/config.h"
        DESTINATION ${INCLUDE_INSTALL_DESTINATION})

add_subdirectory(src)

set(TEST_COMPUTATION_TIME false)
if(BUILD_TESTING)
  add_subdirectory(unit-testings)
endif()
if(BUILD_STATE_OBSERVATION_TOOLS)
  add_subdirectory(misc-tools)
endif()
