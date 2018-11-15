cmake_minimum_required(VERSION 3.0)
project(RenewableSimulator)

enable_language(CXX)

# The NEW behavior for this policy is to link libraries by full
# path even if they are in implicit link directories.
cmake_policy(SET CMP0060 NEW)

# The paths listed there are taken as "implicit" search
# paths and full library paths are replaced accordingly
message("CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES: ${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")

set(COMMON_OUTPUT_DIR "${CMAKE_SOURCE_DIR}/output")

if (EXISTS ${COMMON_OUTPUT_DIR})
    file(REMOVE_RECURSE ${COMMON_OUTPUT_DIR})
endif(EXISTS ${COMMON_OUTPUT_DIR})

file(MAKE_DIRECTORY ${COMMON_OUTPUT_DIR})
file(MAKE_DIRECTORY "${COMMON_OUTPUT_DIR}/bin")

# Define RPATH for MacOS
option(CMAKE_MACOSX_RPATH "Whether RPATH should be used on MacOS"  ON)
option(CMAKE_INSTALL_RPATH_USE_LINK_PATH "Whether to use link path during installation" ON)

SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# Adding requirement for C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define path to SAM package
if (DEFINED SAM_PATH)
    message(STATUS "SAM_PATH is set to ${SAM_PATH}")
else (DEFINED SAM_PATH)
    message(ERROR "Please define SAM_PATH to the directory that includes sscapi.h.")
endif (DEFINED SAM_PATH)

# Build the simulator and import SAM library
message(STATUS "build application renewableSimulator")
include_directories("${SAM_PATH}")
add_executable(${PROJECT_NAME} "renewableSimulator.cpp")
if (APPLE)
    message(STATUS "Link ssc in oxs64")
    target_link_libraries(${PROJECT_NAME} "${SAM_PATH}/osx64/ssc.dylib")
elseif (UNIX)
    message(STATUS "Link ssc in linux64")
    target_link_libraries(${PROJECT_NAME} "${SAM_PATH}/linux64/ssc.so")
else (APPLE)
endif (APPLE)

# Additional linking packages for the SSC
target_link_libraries(${PROJECT_NAME} dl)

# Set the output directory
set_target_properties(${PROJECT_NAME}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${COMMON_OUTPUT_DIR}/bin")
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
