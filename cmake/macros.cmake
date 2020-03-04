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

# This macro takes care of creating a target of ssc
macro(add_ssc)

    # Check whether SSC has been downloaded already
    if(BUILD_SSC)
        check_ssc_download()

        # Add SSC to our current project
        add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/ssc)

        # This property was not set in the sub project. Take care of setting
        # this here to save troubles elsewhere.
        #
        set_target_properties(ssc PROPERTIES 
            INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_BINARY_DIR}/ssc/ssc")

    else(BUILD_SSC)
        add_ssc_library()
    endif(BUILD_SSC)

endmacro(add_ssc)


# Add a prebuilt version of SSC. This only supports Linux platforms.
#
# You don't need to call this function explicitly. Use add_ssc.
#
macro(add_ssc_library)

    # These are required files
    set(SSC_INCLUDE_FILE "${CMAKE_SOURCE_DIR}/ssc/sscapi.h")
    set(SSC_LIB_FILE "${CMAKE_SOURCE_DIR}/ssc/ssc.so")

    message(STATUS "Use pre-built SSC library")

    # Add target
    add_library(ssc SHARED IMPORTED)

    set_target_properties(ssc PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/ssc"
        IMPORTED_LOCATION ${SSC_LIB_FILE}
        PUBLIC_HEADER ${SSC_INCLUDE_FILE})

    target_link_libraries(ssc INTERFACE dl)

endmacro(add_ssc_library)


# Check whether SSC source files have been downloaded. If not, download the source
# files to the current binary source directory
#
# You don't need to call this function explicitly. Use add_ssc.
#
macro(check_ssc_download)

    if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/ssc/CMakeLists.txt)
        find_package(Git 1.9 REQUIRED QUIET)

        message(STATUS "Downloading source files for SSC")
        message(STATUS "The SSC github repository is roughly 1.5 GB")
        message(STATUS "This might take a while depending on your internet connection")

        execute_process(COMMAND ${GIT_EXECUTABLE} "clone"
            "https://github.com/NREL/ssc.git"
            "${CMAKE_CURRENT_BINARY_DIR}/ssc")

        message(STATUS "SSC source files downloaded at ${CMAKE_CURRENT_BINARY_DIR}/ssc")

    else(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/ssc/CMakeLists.txt)
        message(STATUS "SSC source files have been downloaded")
    endif(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/ssc/CMakeLists.txt)

endmacro(check_ssc_download)

