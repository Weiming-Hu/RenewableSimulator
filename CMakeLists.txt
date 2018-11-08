cmake_minimum_required(VERSION 3.0)
project(RenewableSimulator)

enable_language(CXX)

set(COMMON_OUTPUT_DIR "${CMAKE_SOURCE_DIR}/output")
if (EXISTS ${COMMON_OUTPUT_DIR})
    file(REMOVE_RECURSE ${COMMON_OUTPUT_DIR})
endif(EXISTS ${COMMON_OUTPUT_DIR})

file(MAKE_DIRECTORY ${COMMON_OUTPUT_DIR})
file(MAKE_DIRECTORY "${COMMON_OUTPUT_DIR}/bin")

# Define RPATH for MacOS
option(CMAKE_MACOSX_RPATH "Whether RPATH should be used on MacOS" OFF)
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

# Import SAM library
add_library(ssc SHARED IMPORTED GLOBAL)
set_target_properties(ssc PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${SAM_PATH}")
if (APPLE)
    message(STATUS "System: APPLE")
    set_target_properties(ssc PROPERTIES
        IMPORTED_LOCATION "${SAM_PATH}/osx64/ssc.dylib")
    
    if (CMAKE_INSTALL_RPATH)
        set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${SAM_PATH}/osx64/")
    else (CMAKE_INSTALL_RPATH)
        set(CMAKE_INSTALL_RPATH "${SAM_PATH}/osx64/")
    endif (CMAKE_INSTALL_RPATH)

elseif (UNIX)
    message(STATUS "System: UNIX")
    set_target_properties(ssc PROPERTIES IMPORTED_LOCATION "${SAM_PATH}/linux64/ssc.so")
else (APPLE)
    message(ERROR "I don't know your system.")
endif (APPLE)

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
