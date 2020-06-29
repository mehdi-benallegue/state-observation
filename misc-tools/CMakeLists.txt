ADD_EXECUTABLE(flexEstimatorTester model-base-flex-estimator-file-tester.cpp)
TARGET_LINK_LIBRARIES(flexEstimatorTester PUBLIC ${PROJECT_NAME})

ADD_EXECUTABLE(fileSplitter file-splitter.cpp)
TARGET_LINK_LIBRARIES(fileSplitter PUBLIC ${PROJECT_NAME})

ADD_EXECUTABLE(real-offline-flexibility-estimator real-offline-flexibility-estimator.cpp)
TARGET_LINK_LIBRARIES(real-offline-flexibility-estimator PUBLIC ${PROJECT_NAME})


