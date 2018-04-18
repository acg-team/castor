cmake_minimum_required(VERSION 3.5)

set(MINIJATI_SOFTWARENAME "miniJATI")
set(MINIJATI_DESCRIPTION "Joint Alignment Tree Inference")

set(MINIJATI_MAJOR_VERSION 0)
set(MINIJATI_MINOR_VERSION 1)
set(MINIJATI_PATCH_VERSION 1)
set(MINIJATI_VERSION ${MINIJATI_MAJOR_VERSION}.${MINIJATI_MINOR_VERSION}.${MINIJATI_PATCH_VERSION})

project(${MINIJATI_SOFTWARENAME}
        VERSION ${MINIJATI_VERSION}
        DESCRIPTION ${MINIJATI_DESCRIPTION}
        LANGUAGES "CXX")

## Store the git hash of the current head
if(EXISTS "${PROJECT_SOURCE_DIR}/.git/HEAD")
    file(READ "${PROJECT_SOURCE_DIR}/.git/HEAD"
            PROJECT_SOURCE_VERSION)
    if("${PROJECT_SOURCE_VERSION}" MATCHES "^ref:")
        string(REGEX REPLACE "^ref: *([^ \n\r]*).*" "\\1"
                PROJECT_GIT_REF "${PROJECT_SOURCE_VERSION}")
        file(READ "${PROJECT_SOURCE_DIR}/.git/${PROJECT_GIT_REF}"
                PROJECT_SOURCE_VERSION)
    endif()
    string(STRIP "${PROJECT_SOURCE_VERSION}"
            PROJECT_SOURCE_VERSION)

endif()

# Store the build date
if(WIN32)
    execute_process(COMMAND "cmd" " /c date /t"
            OUTPUT_VARIABLE DATE)
    string(REGEX REPLACE "[^0-9]*(..).*" "\\1" MONTH "${DATE}")
    set(MONTHS ""
            "Jan" "Feb" "Mar" "Apr" "May" "Jun"
            "Jul" "Aug" "Sep" "Oct" "Nov" "Dec")
    list(GET MONTHS "${MONTH}" MONTH)
    string(REGEX REPLACE "[^/]*/(..)/(....).*" "\\1 ${MONTH} \\2"
            PROJECT_BUILD_DATE "${DATE}")
    execute_process(COMMAND "cmd" " /c echo %TIME%"
            OUTPUT_VARIABLE TIME)
    string(REGEX REPLACE "[^0-9]*(..:..:..).*" "\\1"
            PROJECT_BUILD_TIME "${TIME}")
else()
    execute_process(COMMAND "date" "+%d %b %Y/%H:%M:%S"
            OUTPUT_VARIABLE DATE_TIME)
    string(REGEX REPLACE "([^/]*)/.*" "\\1"
            PROJECT_BUILD_DATE "${DATE_TIME}")
    string(REGEX REPLACE "[^/]*/([0-9:]*).*" "\\1"
            PROJECT_BUILD_TIME "${DATE_TIME}")
endif()

set(MINIJATI_RELDATE ${PROJECT_BUILD_DATE})
set(MINIJATI_RELTIME ${PROJECT_BUILD_TIME})


add_definitions(-DPRJ_GITBRANCH=\"${PROJECT_GIT_REF}\"
        -DPRJ_VERSION=\"${MINIJATI_VERSION}\"
        -DPRJ_GITREF=\"${PROJECT_SOURCE_VERSION}\"
        -DPRJ_DESC=\"${MINIJATI_DESCRIPTION}\"
        -DPRJ_NAME=\"${MINIJATI_SOFTWARENAME}\"
        -DPRJ_DATE=\"${MINIJATI_RELDATE}\"
        -DPRJ_TIME=\"${MINIJATI_RELTIME}\" )

set(CMAKE_CXX_STANDARD 14)

SET(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wshadow -Weffc++ -O0 -g -DNDEBUG=1")
SET(CMAKE_C_FLAGS_DEBUG "-Wall -Wshadow -O0 -g -DNDEBUG=1")

SET(CMAKE_CXX_FLAGS_RELEASE "-O3  -g0 -DNDEBUG=0")
SET(CMAKE_C_FLAGS_RELEASE "-O3  -g0 -DNDEBUG=0")

SET(CMAKE_CXX_FLAGS_INTEL "-O3  -g0 -DNDEBUG=0")
SET(CMAKE_C_FLAGS_INTEL "-O3  -g0 -DNDEBUG=0")

string(FIND "${CMAKE_BUILD_TYPE}" "static" STATIC)

message("-- Compilation will be performed with the following release of the software:
\tbranch   ${PROJECT_GIT_REF}
\tref      ${PROJECT_SOURCE_VERSION}
\ttime     ${MINIJATI_RELDATE} ${MINIJATI_RELTIME}")

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(WINDOWS TRUE)
    message(STATUS "Compilation will be performed under Windows")

elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(LINUX TRUE)
    message(STATUS "Compilation will be performed under Linux")
    find_package(glog 0.3.5 REQUIRED)
    
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    if (STATIC GREATER 0)
        message(FATAL_ERROR "You can not build this project statically on Mac OS. Ask Apple why! CMake will exit.")
    endif ()
    message(STATUS "Compilation will be performed under Apple MacOS")
    find_package(glog 0.3.5 REQUIRED)

endif ()

# Dependencies not covered by find package should be found in the following directories
if (${CMAKE_PREFIX_PATH})
    include_directories("${CMAKE_PREFIX_PATH}/include")
    LINK_DIRECTORIES("${CMAKE_PREFIX_PATH}}/lib")
    LINK_DIRECTORIES("${CMAKE_PREFIX_PATH}/lib64")
    message(STATUS "Looking for libraries in the following directory: ${CMAKE_PREFIX_PATH}/lib")
    message(STATUS "Looking for libraries in the following directory: ${CMAKE_PREFIX_PATH}/lib64")
    message(STATUS "Looking for headers in the following directory: ${CMAKE_PREFIX_PATH}/include")
endif ()

find_package(Boost REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package(Eigen3 REQUIRED)
if (EIGEN3_FOUND)
    message("-- Found Eigen3 in: ${EIGEN3_INCLUDE_DIR}")
    include_directories(${EIGEN3_INCLUDE_DIR})
endif ()

find_package(TshLib REQUIRED)
if (TSHLIB_FOUND)
    message("-- Found TSHLIB in ${TSHLIB_INCLUDE_DIRS}")
    include_directories(${TSHLIB_INCLUDE_DIRS})
endif ()

find_package(bpp-core REQUIRED)
find_package(bpp-seq REQUIRED)
find_package(bpp-phyl REQUIRED)


set(SOURCES
        main.cpp
        Utilities.cpp
        Utilities.hpp
        PIP.cpp
        PIP.hpp
        Version.hpp
        RHomogeneousTreeLikelihood_PIP.cpp
        RHomogeneousTreeLikelihood_PIP.hpp
        ExtendedAlphabet.cpp
        ExtendedAlphabet.hpp
        TSHHomogeneousTreeLikelihood.cpp
        TSHHomogeneousTreeLikelihood.hpp
        TSHSearchable.hpp
        TSHTopologySearch.cpp
        TSHTopologySearch.hpp
        JATIApplication.cpp
        JATIApplication.hpp
        pPIP.cpp
        pPIP.hpp
        Optimizators.cpp
        Optimizators.hpp
        UnifiedTSHomogeneousTreeLikelihood.cpp
        UnifiedTSHomogeneousTreeLikelihood.hpp
        UnifiedTSHomogeneousTreeLikelihood_PIP.cpp
        UnifiedTSHomogeneousTreeLikelihood_PIP.hpp
        UnifiedTSHSearchable.cpp UnifiedTSHSearchable.hpp)


if (STATIC GREATER 0)
    message(STATUS "Compilation will produce a static executable")
    add_executable(miniJATI ${SOURCES})
    set_target_properties(miniJATI PROPERTIES LINK_FLAGS "-static")
    target_link_libraries(miniJATI tshlib-static)
    target_link_libraries(miniJATI glog)
    target_link_libraries(miniJATI ${BPP_LIBS_STATIC})
    target_link_libraries(miniJATI ${Boost_LIBRARIES})
else ()
    add_executable(miniJATI ${SOURCES})
    target_link_libraries(miniJATI tshlib)
    target_link_libraries(miniJATI glog)
    target_link_libraries(miniJATI ${BPP_LIBS_SHARED})
    target_link_libraries(miniJATI ${Boost_LIBRARIES})
endif ()