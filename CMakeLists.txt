###################################################################################
# Author: Weiming Hu <weiming@psu.edu>                                            #
#         Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu) #
#         Department of Geography                                                 #
#         Institute for Computational and Data Science                            #
#         The Pennsylvania State University                                       #
###################################################################################
#
# This project builds the renewable simulator -- evergreen
#

cmake_minimum_required(VERSION 3.12)
project(evergreen LANGUAGES C CXX)
message(STATUS "
***********************************************************************
*                                                                     *
*      evergreen -- our pursuit for a more sustainable future         *
*                                                                     *
***********************************************************************
")

# I'm skipping the following components when building SAM SSC to remove dependencies
option(skip_tools "Skip building SDKtools and TCSconsole" ON)
option(skip_tests "Skip building SAM tests" ON)

option(BUILD_SSC "Build SSC from source. This might take a while because SSC is large." OFF)

# Only supporting Linux
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    # Linus is supported
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    # Mac OS X is supported
else(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    # Other OS detected but currently they are not supported
    message(FATAL_ERROR "Evergreen currently only support OS X and Linux")
endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")

# This is where our extra cmake scripts are
set(EXTRA_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH ${EXTRA_MODULE_PATH})

# Include our macros
include("${EXTRA_MODULE_PATH}/macros.cmake")

# Add library soltrack
add_soltrack()

# Add library SSC
add_ssc()

# Add a folder to the search path for possible AnEn and AnEnIO locations
# This folder is added simply because of MY PERSONAL HABIT.
# I usually install my packages to this user space rather than the system location
#
set(CMAKE_PREFIX_PATH "$ENV{HOME}/packages/release/lib" CACHE STRING "locations to search for libraries")
message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")

# Add library AnEnIO
find_package(OpenMP)
find_package(Boost 1.58.0 REQUIRED COMPONENTS date_time serialization regex date_time program_options filesystem)
find_package(AnEn REQUIRED)
find_package(AnEnIO REQUIRED)

# Create target
add_executable(${PROJECT_NAME}
    ${CMAKE_CURRENT_SOURCE_DIR}/evergreen/evergreen.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/evergreen/AnEnContainer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/evergreen/Scenarios.cpp)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

target_link_libraries(${PROJECT_NAME} PRIVATE ssc)
target_link_libraries(${PROJECT_NAME} PRIVATE soltrack)
target_link_libraries(${PROJECT_NAME} PRIVATE AnEnIO::AnEnIO Boost::program_options)

# Export the executable
install(TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}Targets
    RUNTIME DESTINATION bin)

