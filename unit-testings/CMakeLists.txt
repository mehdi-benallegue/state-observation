
# provide path to dependency libraries


MACRO(ADD_STATE_ESTIMATOR_TEST NAME)
  add_executable(${NAME} ${NAME}.cpp)
  target_link_libraries(${NAME} PUBLIC state-observation Boost::unit_test_framework Eigen3::Eigen) 
  target_compile_definitions(${NAME} PRIVATE -DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)
  add_test(${NAME} ${NAME})
ENDMACRO()

ADD_STATE_ESTIMATOR_TEST(test-kalman-filter)
ADD_STATE_ESTIMATOR_TEST(test-kalman-filter-template)
ADD_STATE_ESTIMATOR_TEST(imu-test)
ADD_STATE_ESTIMATOR_TEST(imu-multiplicative-test)
ADD_STATE_ESTIMATOR_TEST(HRP2LAAS-IMU-Observation)
ADD_STATE_ESTIMATOR_TEST(flex-estimation-test)
ADD_STATE_ESTIMATOR_TEST(simple-flex-estimation-test)
ADD_STATE_ESTIMATOR_TEST(real-offline-flexibility-estimator)
ADD_STATE_ESTIMATOR_TEST(other-tests)
ADD_STATE_ESTIMATOR_TEST(test_acceleration_stabilization)
ADD_STATE_ESTIMATOR_TEST(test_imu-elastic-local-frame-dynamical-system)
ADD_STATE_ESTIMATOR_TEST(test_model-base-ekf-flex-estimator-imu)
ADD_STATE_ESTIMATOR_TEST(test_model-base-ekf-flex-estimator-imu_withComBias_withForcesSensors)
ADD_STATE_ESTIMATOR_TEST(zmpEstimation)
ADD_STATE_ESTIMATOR_TEST(test-kinetics-observer)
