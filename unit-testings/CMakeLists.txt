MACRO(ADD_STATE_ESTIMATOR_TEST NAME)
  add_executable(${NAME} ${NAME}.cpp)
  target_link_libraries(${NAME} PUBLIC state-observation Boost::unit_test_framework) 
  target_compile_definitions(${NAME} PRIVATE -DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)
  add_test(${NAME} ${NAME})
  # Adding a project configuration file (for MSVC only)
  generate_msvc_dot_user_file(${NAME})
ENDMACRO()

ADD_STATE_ESTIMATOR_TEST(test-kalman-filter)
ADD_STATE_ESTIMATOR_TEST(imu-test)
ADD_STATE_ESTIMATOR_TEST(imu-multiplicative-test)
ADD_STATE_ESTIMATOR_TEST(flex-estimation-test)
ADD_STATE_ESTIMATOR_TEST(simple-flex-estimation-test)
ADD_STATE_ESTIMATOR_TEST(other-tests)
ADD_STATE_ESTIMATOR_TEST(test_acceleration_stabilization)
ADD_STATE_ESTIMATOR_TEST(test_model-base-ekf-flex-estimator-imu)
ADD_STATE_ESTIMATOR_TEST(test-kinetics-observer)

add_custom_command(TARGET test_model-base-ekf-flex-estimator-imu POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                       ${CMAKE_SOURCE_DIR}/unit-testings/inputFiles/ $<TARGET_FILE_DIR:test_model-base-ekf-flex-estimator-imu>/inputFiles/)
