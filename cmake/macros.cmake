###################################################################################
# Author: Weiming Hu <weiming@psu.edu>                                            #
#         Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu) #
#         Department of Geography                                                 #
#         Institute for Computational and Data Science                            #
#         The Pennsylvania State University                                       #
###################################################################################
#
# This file includes some macros for the project
#

# Check whether SSC source files have been downloaded
macro(check_ssc_download)

    if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/ssc/CMakeLists.txt)
        find_package(Git 1.9 REQUIRED QUIET)
        message(STATUS "Downloading source files for SSC ...")
        execute_process(COMMAND ${GIT_EXECUTABLE} "clone"
            "https://github.com/NREL/ssc.git"
            "${CMAKE_CURRENT_BINARY_DIR}/ssc")
        message(STATUS "SSC source files downloaded at ${CMAKE_CURRENT_BINARY_DIR}/ssc")
    endif(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/ssc/CMakeLists.txt)

endmacro(check_ssc_download)

