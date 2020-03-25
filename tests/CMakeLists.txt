cmake_minimum_required(VERSION 3.12)
project(testEvergreen LANGUAGES C CXX)

add_custom_target(test_all ${CMAKE_CTEST_COMMAND} -V)

# Find the dependency package
find_package(AnEnIO REQUIRED)

# Add tests
add_library("testEvergreen" "${CMAKE_CURRENT_SOURCE_DIR}/testEvergreen.cpp" ${EXTRA_CPP})
target_link_libraries("testEvergreen" PUBLIC ssc soltrack AnEnIO::AnEnIO ${CPPUNIT_LIBRARIES})
target_include_directories("testEvergreen" PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../evergreen>)

add_executable("runEvergreen" "${CMAKE_CURRENT_SOURCE_DIR}/runnerEvergreen.cpp")
target_link_libraries("runEvergreen" "testEvergreen")

add_test(NAME TestEvergreen COMMAND "runEvergreen")
install(TARGETS "runEvergreen" RUNTIME DESTINATION test)

