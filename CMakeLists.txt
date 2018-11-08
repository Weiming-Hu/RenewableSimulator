cmake_minimum_required(VERSION 3.0)
project(RenewableSimulator)

enable_language(CXX)

set(COMMON_OUTPUT_DIR "${CMAKE_SOURCE_DIR}/output")
if (EXISTS ${COMMON_OUTPUT_DIR})
    file(REMOVE_RECURSE ${COMMON_OUTPUT_DIR})
endif(EXISTS ${COMMON_OUTPUT_DIR})

file(MAKE_DIRECTORY ${COMMON_OUTPUT_DIR})
file(MAKE_DIRECTORY "${COMMON_OUTPUT_DIR}/bin")
file(MAKE_DIRECTORY "${COMMON_OUTPUT_DIR}/lib")

# Adding requirement for C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define path to SAM package
if (DEFINED SAM_PATH)
    message(STATUS "SAM_PATH is set to ${SAM_PATH}")
else (DEFINED SAM_PATH)
    message(ERROR "Please define SAM_PATH to the directory that includes sscapi.h.")
endif (DEFINED SAM_PATH)

# Import SAM library
include_directories(${SAM_PATH})
add_library(ssc SHARED IMPORTED)
set_target_properties(ssc PROPERTIES IMPORTED_LOCATION "${SAM_PATH}/linux64/ssc.so")

# Build the simulator
message(STATUS "build application renewableSimulator")
add_executable(${PROJECT_NAME} "renewableSimulator.cpp")

# Add dependency to ssc
add_dependencies(${PROJECT_NAME} ssc)
target_link_libraries(${PROJECT_NAME} ssc)

# Additional linking packages for the SSC
target_link_libraries(${PROJECT_NAME} dl)

# Set the output directory
set_target_properties(${PROJECT_NAME}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${COMMON_OUTPUT_DIR}/bin")
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
