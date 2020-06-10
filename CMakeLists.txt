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

CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)

SET(PROJECT_NAME state-observation)
SET(PROJECT_DESCRIPTION "General implementation of observers.")
SET(PROJECT_URL "")
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(PROJECT_VERSION 1.2.1)
SET(CXX_DISABLE_WERROR TRUE)
set(PROJECT_USE_CMAKE_EXPORT TRUE)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/eigen.cmake)

PROJECT(${PROJECT_NAME} LANGUAGES CXX)

add_project_dependency(Boost COMPONENTS filesystem system thread program_options unit_test_framework timer)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
add_project_dependency(Eigen3 QUIET NO_CMAKE_PACKAGE_REGISTRY)
if(Eigen3_FOUND)
  add_project_dependency(Eigen3 REQUIRED NO_CMAKE_PACKAGE_REGISTRY)
else()
  add_project_dependency(Eigen3 MODULE REQUIRED)
endif()

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(unit-testings)
ADD_SUBDIRECTORY(misc-tools)
