# Copyright 2008, 2009, 2010, 2013 Olivier Stasse, Mehdi Benallegue JRL, CNRS/AIST
#
# This file is part of state-observation.
# state-observation is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# state-observation is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Lesser Public License for more details.
# You should have received a copy of the GNU Lesser General Public License
# along with state-observation.  If not, see <http://www.gnu.org/licenses/>.

set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME})

set(HEADERS
    ${PROJECT_INCLUDE_DIR}/observer/tilt-estimator.hpp
    ${PROJECT_INCLUDE_DIR}/observer/extended-kalman-filter.hpp
    ${PROJECT_INCLUDE_DIR}/observer/kalman-filter-base.hpp
    ${PROJECT_INCLUDE_DIR}/observer/linear-kalman-filter.hpp
    ${PROJECT_INCLUDE_DIR}/observer/observer-base.hpp
    ${PROJECT_INCLUDE_DIR}/observer/zero-delay-observer.hpp
    ${PROJECT_INCLUDE_DIR}/sensors-simulation/sensor-base.hpp
    ${PROJECT_INCLUDE_DIR}/sensors-simulation/algebraic-sensor.hpp
    ${PROJECT_INCLUDE_DIR}/sensors-simulation/accelerometer-gyrometer.hpp
    ${PROJECT_INCLUDE_DIR}/sensors-simulation/accelerometer-gyrometer-magnetometer.hpp
    ${PROJECT_INCLUDE_DIR}/noise/noise-base.hpp
    ${PROJECT_INCLUDE_DIR}/noise/gaussian-white-noise.hpp
    ${PROJECT_INCLUDE_DIR}/sensors-simulation/algorithm/rotation-velocity.hpp
    ${PROJECT_INCLUDE_DIR}/sensors-simulation/algorithm/linear-acceleration.hpp
    ${PROJECT_INCLUDE_DIR}/sensors-simulation/algorithm/magnetic-field.hpp
    ${PROJECT_INCLUDE_DIR}/tools/definitions.hpp
    ${PROJECT_INCLUDE_DIR}/tools/definitions.hxx
    ${PROJECT_INCLUDE_DIR}/tools/hrp2.hpp
    ${PROJECT_INCLUDE_DIR}/tools/probability-law-simulation.hpp
    ${PROJECT_INCLUDE_DIR}/tools/miscellaneous-algorithms.hpp
    ${PROJECT_INCLUDE_DIR}/tools/rigid-body-kinematics.hpp
    ${PROJECT_INCLUDE_DIR}/tools/rigid-body-kinematics.hxx
    ${PROJECT_INCLUDE_DIR}/tools/logger.hpp
    ${PROJECT_INCLUDE_DIR}/tools/logger.hxx
    ${PROJECT_INCLUDE_DIR}/tools/state-vector-arithmetics.hpp
    ${PROJECT_INCLUDE_DIR}/dynamical-system/dynamical-system-functor-base.hpp
    ${PROJECT_INCLUDE_DIR}/dynamical-system/dynamical-system-simulator.hpp
    ${PROJECT_INCLUDE_DIR}/dynamical-system/imu-dynamical-system.hpp
    ${PROJECT_INCLUDE_DIR}/dynamical-system/imu-mltpctive-dynamical-system.hpp
    ${PROJECT_INCLUDE_DIR}/dynamical-system/imu-magnetometer-dynamical-system.hpp
    ${PROJECT_INCLUDE_DIR}/dynamical-system/bidim-elastic-inv-pendulum-dyn-sys.hpp
    ${PROJECT_INCLUDE_DIR}/flexibility-estimation/flexibility-estimator-base.hpp
    ${PROJECT_INCLUDE_DIR}/flexibility-estimation/ekf-flexibility-estimator-base.hpp
    ${PROJECT_INCLUDE_DIR}/flexibility-estimation/fixed-contact-ekf-flex-estimator-imu.hpp
    ${PROJECT_INCLUDE_DIR}/flexibility-estimation/model-base-ekf-flex-estimator-imu.hpp
    ${PROJECT_INCLUDE_DIR}/flexibility-estimation/imu-fixed-contact-dynamical-system.hpp
    ${PROJECT_INCLUDE_DIR}/flexibility-estimation/stable-imu-fixed-contact-dynamical-system.hpp
    ${PROJECT_INCLUDE_DIR}/flexibility-estimation/imu-elastic-local-frame-dynamical-system.hpp
    ${PROJECT_INCLUDE_DIR}/dynamics-estimators/kinetics-observer.hpp
    ${PROJECT_INCLUDE_DIR}/dynamics-estimators/kinetics-observer.hxx
    ${PROJECT_INCLUDE_DIR}/dynamics-estimators/lipm-dcm-estimator.hpp
    ${PROJECT_INCLUDE_DIR}/dynamics-estimators/unidim-lipm-dcm-estimator.hpp
    ${PROJECT_INCLUDE_DIR}/examples/imu-attitude-trajectory-reconstruction.hpp
    ${PROJECT_INCLUDE_DIR}/examples/imu-attitude-trajectory-reconstruction.hxx
    ${PROJECT_INCLUDE_DIR}/examples/imu-multiplicative-attitude-reconstruction.hpp
    ${PROJECT_INCLUDE_DIR}/examples/imu-multiplicative-attitude-reconstruction.hxx
    ${PROJECT_INCLUDE_DIR}/examples/offline-ekf-flexibility-estimation.hpp
    ${PROJECT_INCLUDE_DIR}/examples/offline-ekf-flexibility-estimation.hxx
    ${PROJECT_INCLUDE_DIR}/examples/offline-model-base-flex-estimation.hpp
    ${PROJECT_INCLUDE_DIR}/examples/offline-model-base-flex-estimation.hxx)

set(SRC
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
    state-vector-arithmetics.cpp
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
    lipm-dcm-estimator.cpp
    unidim-lipm-dcm-estimator.cpp)

set(LIBRARY_NAME ${PROJECT_NAME})
add_library(${LIBRARY_NAME} SHARED ${SRC} ${HEADERS})
set_target_properties(${LIBRARY_NAME} PROPERTIES COMPILE_FLAGS "-DSTATE_OBSERVATION_EXPORTS")
set_target_properties(${LIBRARY_NAME} PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR} VERSION ${PROJECT_VERSION})
target_include_directories(
  ${LIBRARY_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/include
                         $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>/$<CONFIG>/include $<INSTALL_INTERFACE:include>)

target_link_libraries(${LIBRARY_NAME} PUBLIC Boost::timer Boost::disable_autolinking)
if(TARGET Eigen3::Eigen)
  target_link_libraries(${LIBRARY_NAME} PUBLIC Eigen3::Eigen)
else()
  target_include_directories(${LIBRARY_NAME} PUBLIC "${EIGEN3_INCLUDE_DIR}")
endif()
set_target_properties(${LIBRARY_NAME} PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR} VERSION ${PROJECT_VERSION})

install(
  TARGETS ${LIBRARY_NAME}
  EXPORT "${TARGETS_EXPORT_NAME}"
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin)

install(DIRECTORY ${PROJECT_INCLUDE_DIR}/ DESTINATION ${INCLUDE_INSTALL_DESTINATION})

pkg_config_append_libs(${PROJECT_NAME})
