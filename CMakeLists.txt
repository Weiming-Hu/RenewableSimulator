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
if(UNIX AND NOT APPLE)
    message(STATUS "Linux detected")
else(UNIX AND NOT APPLE)
    message(FATAL_ERROR "Only supporting Linux")
endif(UNIX AND NOT APPLE)

# This is where our extra cmake scripts are
set(EXTRA_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH ${EXTRA_MODULE_PATH})

# Include our macros
include("${EXTRA_MODULE_PATH}/macros.cmake")

# Add library SSC
add_ssc()

# Add library AnEnIO
find_package(OpenMP)
find_package(Boost 1.58.0 REQUIRED COMPONENTS date_time serialization regex date_time filesystem)
find_package(AnEn REQUIRED)
find_package(AnEnIO REQUIRED)

# Create target
add_executable(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/evergreen.cpp)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

target_link_libraries(${PROJECT_NAME} PUBLIC ssc)
target_link_libraries(${PROJECT_NAME} PUBLIC AnEnIO::AnEnIO)

# Export the executable
install(TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}Targets
    RUNTIME DESTINATION bin)

