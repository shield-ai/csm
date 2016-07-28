cmake_minimum_required(VERSION 2.6)

PROJECT(csm_eigen)

find_package(catkin REQUIRED)
find_package(cmake_modules REQUIRED)
find_package(Eigen REQUIRED)

catkin_package(
    INCLUDE_DIRS include ${PROJECT_SOURCE_DIR}/src
    LIBRARIES csm_eigen
    DEPENDS Eigen
    )

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_RELEASE "-O3")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

if(WIN32)
    add_definitions(-DWINDOWS)
    if(MSVC)
        # SSE2 optimizations
        ADD_DEFINITIONS("/arch:SSE2")
    endif()
endif()

#SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_BUILD_TYPE RelWithDebInfo)
#SET(CMAKE_BUILD_TYPE Release)

#SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../bin )
#SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../lib )
#SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../lib )

OPTION(ENABLE_C11 "Flag to enable C++11" ON)
OPTION(ENABLE_OPTIMIZATION "Flag to enable optimizations" ON)

IF (ENABLE_C11)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ENDIF()

CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/src/csm/flags.h.in"
    "${PROJECT_SOURCE_DIR}/src/csm/flags.h")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

INCLUDE_DIRECTORIES(
    include
    ${EIGEN_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/src
    )

IF (ENABLE_OPTIMIZATION)
    SET(CSM_SRC
        src/csm/laser_data.cpp
        src/csm/math_utils.cpp
        src/csm/math_utils_gsl.cpp
        src/csm/logging.cpp
        src/csm/orientation.cpp
        src/csm/clustering.cpp)
ELSE()
    SET(CSM_SRC
        src/csm/laser_data.cpp
        src/csm/math_utils.cpp
        src/csm/math_utils_gsl.cpp
        src/csm/utils.cpp
        src/csm/logging.cpp
        src/csm/sm_options.cpp
        src/csm/orientation.cpp
        src/csm/clustering.cpp)
ENDIF()

SET(EGSL_SRC
    src/egsl/egsl.cpp
    src/egsl/egsl_conversions.cpp
    src/egsl/egsl_misc.cpp
    src/egsl/egsl_ops.cpp)

SET(GPC_SRC
    src/gpc/gpc.cpp
    src/gpc/gpc_utils.cpp)

SET(ICP_SRC
    src/icp/icp.cpp
    src/icp/icp_corr_dumb.cpp
    src/icp/icp_corr_tricks.cpp
    src/icp/icp_covariance.cpp
    src/icp/icp_debug.cpp
    src/icp/icp_loop.cpp
    src/icp/icp_outliers.cpp)

IF (NOT ENABLE_OPTIMIZATION)
    SET(OPTIONS_SRC
        src/options/options.cpp
        src/options/options_interface.cpp)
ENDIF()

IF (ENABLE_OPTIMIZATION)
    ADD_LIBRARY(csm_eigen
        ${CSM_SRC}
        ${EGSL_SRC}
        ${GPC_SRC}
        ${ICP_SRC})
ELSE()
    ADD_LIBRARY(csm_eigen
        ${CSM_SRC}
        ${EGSL_SRC}
        ${GPC_SRC}
        ${OPTIONS_SRC}
        ${ICP_SRC})
ENDIF()
#install(TARGETS csm_eigen
#  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
#)
