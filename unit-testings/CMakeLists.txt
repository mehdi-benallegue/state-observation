macro(ADD_SO_TEST name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PUBLIC state-observation ${ARGN})

  add_test(
    NAME ${name}
    COMMAND $<TARGET_FILE:${name}>
    #CONFIGURATIONS $<CONFIGURATION>
    WORKING_DIRECTORY $<TARGET_FILE_DIR:${name}>)
  # Adding a project configuration file (for MSVC only)
  generate_msvc_dot_user_file(${name})
endmacro()

add_so_test(test-kalman-filter Boost::boost)
add_so_test(imu-test)
add_so_test(imu-multiplicative-test)
add_so_test(flex-estimation-test)
add_so_test(simple-flex-estimation-test)
add_so_test(other-tests)
add_so_test(test_acceleration_stabilization)
add_so_test(test_model-base-ekf-flex-estimator-imu)
add_so_test(test-kinetics-observer)

add_custom_command(
  TARGET test_model-base-ekf-flex-estimator-imu
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/unit-testings/inputFiles/
          $<TARGET_FILE_DIR:test_model-base-ekf-flex-estimator-imu>/inputFiles/)
