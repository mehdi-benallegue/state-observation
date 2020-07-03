
MACRO(ADD_STATE_ESTIMATOR_TOOL NAME)
add_executable(${NAME} ${NAME}.cpp)
target_link_libraries(${NAME} PUBLIC state-observation Boost::unit_test_framework) 
ENDMACRO()


ADD_STATE_ESTIMATOR_TOOL(model-base-flex-estimator-file-tester)
ADD_STATE_ESTIMATOR_TOOL(file-splitter)
ADD_STATE_ESTIMATOR_TOOL(real-offline-flexibility-estimator)
ADD_STATE_ESTIMATOR_TOOL(test_model-base-ekf-flex-estimator-imu_withComBias_withForcesSensors)
ADD_STATE_ESTIMATOR_TOOL(zmpEstimation)
ADD_STATE_ESTIMATOR_TOOL(IMU-Observation-tester)

