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

# This macro imports a library soltrack
macro(add_soltrack)

    # These are required files
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(SOLTRACK_LIB_FILENAME "libSolTrack.so")
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        set(SOLTRACK_LIB_FILENAME "libSolTrack.dylib")
    else(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        message(FATAL_ERROR "Currently only supports Linux and OSX")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")

    set(SOLTRACK_ROOT "$ENV{HOME}/packages/release/" CACHE STRING "root directory for soltrack")

    set(SOLTRACK_LIB_DIR "${SOLTRACK_ROOT}/lib")
    set(SOLTRACK_INCLUDE_DIR "${SOLTRACK_ROOT}/include")
    set(SOLTRACK_LIB_FILE "${SOLTRACK_LIB_DIR}/${SOLTRACK_LIB_FILENAME}")

    if(NOT EXISTS ${SOLTRACK_LIB_FILE})
        message(FATAL_ERROR "${SOLTRACK_LIB_FILE} is not found. \
        Download and install SolTrack from http://soltrack.sourceforge.net/ and \
        specify the installation root directory in SOLTRACK_ROOT")
    endif(NOT EXISTS ${SOLTRACK_LIB_FILE})

    if(NOT EXISTS ${SOLTRACK_INCLUDE_DIR})
        message(FATAL_ERROR "${SOLTRACK_INCLUDE_DIR} is not found")
    endif(NOT EXISTS ${SOLTRACK_INCLUDE_DIR})

    # Add target
    add_library(soltrack INTERFACE IMPORTED)

    set_target_properties(soltrack PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES ${SOLTRACK_INCLUDE_DIR})

    target_link_libraries(soltrack INTERFACE ${SOLTRACK_LIB_FILE})

    # SolTrack is found from the CMAKE_PREFIX_PATH, so I need to add this directory to RPATH
    list(APPEND CMAKE_INSTALL_RPATH ${CMAKE_PREFIX_PATH})

endmacro(add_soltrack)


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
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(SSC_LIB_FILENAME "libssc.so")
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        set(SSC_LIB_FILENAME "libssc.dylib")
    else(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        message(FATAL_ERROR "Currently only compiled librares for Linux and OSX have been uploaded")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")

    set(LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ssc")
    set(SSC_INCLUDE_FILE "${LIB_DIR}/sscapi.h")
    set(SSC_LIB_FILE_BUILD "${LIB_DIR}/${SSC_LIB_FILENAME}")
    set(SSC_LIB_FILE_INSTALL "lib/${SSC_LIB_FILENAME}")

    # SSC is found from the current github repo and will be installed to the installation path,
    # so I need to add this directory to RPATH
    #
    list(APPEND CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)

    message(STATUS "Pre-built SSC library: ${SSC_LIB_FILE_BUILD}")
    message(STATUS "SSC inlude header file: ${SSC_INCLUDE_FILE}")

    # Add target
    add_library(ssc INTERFACE IMPORTED)

    set_target_properties(ssc PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES ${LIB_DIR}
        PUBLIC_HEADER ${SSC_INCLUDE_FILE})

    target_link_libraries(ssc INTERFACE
        $<BUILD_INTERFACE:${SSC_LIB_FILE_BUILD}>
        $<INSTALL_INTERFACE:${SSC_LIB_FILE_INSTALL}>)

    target_link_libraries(ssc INTERFACE dl)

    install(FILES ${SSC_LIB_FILE_BUILD} DESTINATION lib)

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

