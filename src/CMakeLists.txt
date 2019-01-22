# Copyright 2008, 2009, 2010, 2013 Olivier Stasse, Mehdi Benallegue JRL, CNRS/AIST
#
# This file is part of state-observation.
# state-observation is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# state-observation is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Lesser Public License for more details.
# You should have received a copy of the GNU Lesser General Public License
# along with state-observation.  If not, see <http://www.gnu.org/licenses/>.

# prefix and suffix each element of list by ${prefix}elemnt${suffix}
macro(ADDPREFIX newlist prefix list_name)
    # create empty list - necessary?
    SET(${newlist})

     # prefix and suffix elements
    foreach(l ${${list_name}})
      list(APPEND ${newlist} ${prefix}${l} )
    endforeach()

endmacro(ADDPREFIX)

if(${_Eigen_VERSION} VERSION_LESS "3.2.0")
    add_definitions(-DEIGEN_VERSION_LESS_THAN_3_2)
endif()

ADDPREFIX(${PROJECT_NAME}_ABSOLUTE_HEADERS  "${CMAKE_SOURCE_DIR}/"  ${PROJECT_NAME}_HEADERS)

SET(LIBRARY_NAME ${PROJECT_NAME})

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

ADD_LIBRARY(${LIBRARY_NAME}
  SHARED
  tilt-estimator.cpp
  kalman-filter-base.cpp
  observer-base.cpp
  extended-kalman-filter.cpp
  linear-kalman-filter.cpp
  zero-delay-observer.cpp
  sensor-base.cpp
  algebraic-sensor.cpp
  rotation-velocity.cpp
  linear-acceleration.cpp
  magnetic-field.cpp
  definitions.cpp
  logger.cpp
  accelerometer-gyrometer.cpp
  accelerometer-gyrometer-magnetometer.cpp
  probability-law-simulation.cpp
  gaussian-white-noise.cpp
  dynamical-system-functor-base.cpp
  dynamical-system-simulator.cpp
  imu-dynamical-system.cpp
  imu-mltpctive-dynamical-system.cpp
  imu-magnetometer-dynamical-system.cpp
  flexibility-estimator-base.cpp
  ekf-flexibility-estimator-base.cpp
  fixed-contact-ekf-flex-estimator-imu.cpp
  model-base-ekf-flex-estimator-imu.cpp
  imu-fixed-contact-dynamical-system.cpp
  stable-imu-fixed-contact-dynamical-system.cpp
  imu-elastic-local-frame-dynamical-system.cpp
  bidim-elastic-inv-pendulum-dyn-sys.cpp
  kinetics-observer.cpp
  kinetics-dynamical-system.cpp
  ${${PROJECT_NAME}_ABSOLUTE_HEADERS}
  )

SET_TARGET_PROPERTIES(${LIBRARY_NAME}
  PROPERTIES
  SOVERSION ${PROJECT_VERSION}
  INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)

INSTALL(TARGETS ${LIBRARY_NAME}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

ADD_REQUIRED_DEPENDENCY("eigen3 >= 3.0.0")
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIRS})
SET(LINEAR_ALGEBRA_LIB_FOUND 1)
SET(_EIGEN_MATRIX_FOUND_ 1)


