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
*      evergreen -- our pursuit for a more sustainable society        *
*                                                                     *
***********************************************************************
")

# I'm skipping the following components when building SAM SSC to remove dependencies
option(skip_tools "Skip building SDKtools and TCSconsole" ON)
option(skip_tests "Skip building SAM tests" ON)

# Check whether SSC has been downloaded already
check_ssc_download()

# Add SSC to our current project
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/ssc)

# Create target
add_executable(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/evergreen.cpp)

# Configure the properties of this target
target_link_libraries(${PROJECT_NAME} PUBLIC sam_simulation_core)

# Export the executable
install(TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}Targets
    RUNTIME DESTINATION bin)

