MACRO(ADD_SO_TEST NAME)
  add_executable(${NAME} ${NAME}.cpp)
  target_link_libraries(${NAME} PUBLIC state-observation ${ARGN}) 

  add_test(NAME ${NAME} 
	  COMMAND $<TARGET_FILE:${NAME}>
	  #CONFIGURATIONS $<CONFIGURATION>
	  WORKING_DIRECTORY $<TARGET_FILE_DIR:${NAME}>)
  # Adding a project configuration file (for MSVC only)
  generate_msvc_dot_user_file(${NAME})
ENDMACRO()

ADD_SO_TEST(test-kalman-filter Boost::boost)
ADD_SO_TEST(imu-test)
ADD_SO_TEST(imu-multiplicative-test)
ADD_SO_TEST(flex-estimation-test)
ADD_SO_TEST(simple-flex-estimation-test)
ADD_SO_TEST(other-tests)
ADD_SO_TEST(test_acceleration_stabilization)
ADD_SO_TEST(test_model-base-ekf-flex-estimator-imu)
ADD_SO_TEST(test-kinetics-observer)

add_custom_command(TARGET test_model-base-ekf-flex-estimator-imu POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                       ${CMAKE_SOURCE_DIR}/unit-testings/inputFiles/ $<TARGET_FILE_DIR:test_model-base-ekf-flex-estimator-imu>/inputFiles/)
