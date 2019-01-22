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

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/eigen.cmake)

SET(PROJECT_NAME state-observation)
SET(PROJECT_DESCRIPTION
  "General implementation of observers.")
SET(PROJECT_URL "")

SET(CUSTOM_HEADER_DIR "${PROJECT_NAME}")

SET(${PROJECT_NAME}_HEADERS
  include/state-observation/observer/tilt-estimator.hpp
  include/state-observation/observer/extended-kalman-filter.hpp
  include/state-observation/observer/kalman-filter-base.hpp
  include/state-observation/observer/linear-kalman-filter.hpp
  include/state-observation/observer/observer-base.hpp
  include/state-observation/observer/zero-delay-observer.hpp
  include/state-observation/observer/compile-time/compile-time-extended-kalman-filter.hpp
  include/state-observation/observer/compile-time/compile-time-extended-kalman-filter.hxx
  include/state-observation/observer/compile-time/compile-time-kalman-filter-base.hpp
  include/state-observation/observer/compile-time/compile-time-kalman-filter-base.hxx
  include/state-observation/observer/compile-time/compile-time-kalman-filter.hpp
  include/state-observation/observer/compile-time/compile-time-kalman-filter.hxx
  include/state-observation/observer/compile-time/compile-time-observer-base.hpp
  include/state-observation/observer/compile-time/compile-time-observer-base.hxx
  include/state-observation/observer/compile-time/compile-time-zero-delay-observer.hpp
  include/state-observation/observer/compile-time/compile-time-zero-delay-observer.hxx
  include/state-observation/sensors-simulation/sensor-base.hpp
  include/state-observation/sensors-simulation/algebraic-sensor.hpp
  include/state-observation/sensors-simulation/accelerometer-gyrometer.hpp
  include/state-observation/sensors-simulation/accelerometer-gyrometer-magnetometer.hpp
  include/state-observation/noise/noise-base.hpp
  include/state-observation/noise/gaussian-white-noise.hpp
  include/state-observation/sensors-simulation/algorithm/rotation-velocity.hpp
  include/state-observation/sensors-simulation/algorithm/linear-acceleration.hpp
  include/state-observation/sensors-simulation/algorithm/magnetic-field.hpp
  include/state-observation/tools/definitions.hpp
  include/state-observation/tools/definitions.hxx
  include/state-observation/tools/hrp2.hpp
  include/state-observation/tools/probability-law-simulation.hpp
  include/state-observation/tools/miscellaneous-algorithms.hpp
  include/state-observation/tools/rigid-body-kinematics.hpp
  include/state-observation/tools/rigid-body-kinematics.hxx
  include/state-observation/tools/logger.hpp
  include/state-observation/tools/logger.hxx
  include/state-observation/dynamical-system/dynamical-system-functor-base.hpp
  include/state-observation/dynamical-system/dynamical-system-simulator.hpp
  include/state-observation/dynamical-system/imu-dynamical-system.hpp
  include/state-observation/dynamical-system/imu-mltpctive-dynamical-system.hpp
  include/state-observation/dynamical-system/imu-magnetometer-dynamical-system.hpp
  include/state-observation/dynamical-system/bidim-elastic-inv-pendulum-dyn-sys.hpp
  include/state-observation/flexibility-estimation/flexibility-estimator-base.hpp
  include/state-observation/flexibility-estimation/ekf-flexibility-estimator-base.hpp
  include/state-observation/flexibility-estimation/fixed-contact-ekf-flex-estimator-imu.hpp
  include/state-observation/flexibility-estimation/model-base-ekf-flex-estimator-imu.hpp
  include/state-observation/flexibility-estimation/imu-fixed-contact-dynamical-system.hpp
  include/state-observation/flexibility-estimation/stable-imu-fixed-contact-dynamical-system.hpp
  include/state-observation/flexibility-estimation/imu-elastic-local-frame-dynamical-system.hpp
  include/state-observation/dynamic-estimators/kinetics-observer.hpp
  include/state-observation/dynamic-estimators/kinetics-dynamical-system.hpp
  include/state-observation/examples/imu-attitude-trajectory-reconstruction.hpp
  include/state-observation/examples/imu-attitude-trajectory-reconstruction.hxx
  include/state-observation/examples/imu-multiplicative-attitude-reconstruction.hpp
  include/state-observation/examples/imu-multiplicative-attitude-reconstruction.hxx
  include/state-observation/examples/offline-ekf-flexibility-estimation.hpp
  include/state-observation/examples/offline-ekf-flexibility-estimation.hxx
  include/state-observation/examples/offline-model-base-flex-estimation.hpp
  include/state-observation/examples/offline-model-base-flex-estimation.hxx
)
SET(CXX_DISABLE_WERROR TRUE)

SETUP_PROJECT()

SET(BOOST_COMPONENTS filesystem system thread program_options unit_test_framework timer)

SEARCH_FOR_BOOST()
SEARCH_FOR_EIGEN()

PKG_CONFIG_APPEND_LIBS(${PROJECT_NAME})

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(unit-testings)
ADD_SUBDIRECTORY(misc-tools)
SETUP_PROJECT_FINALIZE()
