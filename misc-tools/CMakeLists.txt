ADD_EXECUTABLE(flexEstimatorTester model-base-flex-estimator-file-tester.cpp)
TARGET_LINK_LIBRARIES(flexEstimatorTester PUBLIC ${PROJECT_NAME})

ADD_EXECUTABLE(fileSplitter file-splitter.cpp)
TARGET_LINK_LIBRARIES(fileSplitter PUBLIC ${PROJECT_NAME})
