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


message("-- Compilation will be performed with the following release of the software:
branch   ${PROJECT_GIT_REF}
ref      ${PROJECT_SOURCE_VERSION}
time     ${MINIJATI_RELDATE} ${MINIJATI_RELTIME}")

set(CMAKE_CXX_STANDARD 14)

SET(CMAKE_CXX_FLAGS_DEBUG  "-Wall -Wshadow -Weffc++ -O0 -g")
SET(CMAKE_C_FLAGS_DEBUG "-Wall -Wshadow -O0 -g")

SET(CMAKE_CXX_FLAGS_RELEASE  "-O3")
SET(CMAKE_C_FLAGS_RELEASE "-O3")

SET(CMAKE_CXX_FLAGS_INTEL  "-O3")
SET(CMAKE_C_FLAGS_INTEL "-O3")


if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(WINDOWS TRUE)
    message(STATUS "Compilation will be performed under Windows")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(LINUX TRUE)
    message(STATUS "Compilation will be performed under Linux")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lstdc++ -ldl -pthread ")
    if (GOOGLE_INCLUDE_PATH)
        message(STATUS "Google Path is ${GOOGLE_INCLUDE_PATH}")
        include_directories("${GOOGLE_INCLUDE_PATH}/glog")
        include_directories("${GOOGLE_INCLUDE_PATH}/gflags")
    else ()
        include_directories("/usr/local/include/glog")
        include_directories("/usr/local/include/gflags")
    endif ()
    #find_library(gflags gflags)
    #message(STATUS "Found Gflags at ${gflags}")
    #find_package(glog REQUIRED)
    #find_package(gflags 2.2 REQUIRED)

elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    message(STATUS "Compilation will be performed under Apple MacOS")
    find_package(glog 0.3.5 REQUIRED)
    find_package(gflags REQUIRED)
endif ()

find_package(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()



find_package("Eigen3")
if (EIGEN3_FOUND)
    message("Found Eigen3 in: ${EIGEN3_INCLUDE_DIR}")
endif (EIGEN3_FOUND)
include_directories(${EIGEN3_INCLUDE_DIR})


find_package("TshLib")
if (TSHLIB_FOUND)
    message("Found TSHLIB in ${TSHLIB_INCLUDE_DIRS}")
endif(TSHLIB_FOUND)
include_directories(${TSHLIB_INCLUDE_DIRS})


find_package(bpp-core)
find_package(bpp-seq)
find_package(bpp-phyl)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
#find_package(OpenMP)
#if (OPENMP_FOUND)
#    message("OpenMP found")
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

add_definitions(-DOMPENABLED=1)

set(SOURCES
        main.cpp
        Utilities.cpp
        Utilities.hpp
        PIP.cpp
        PIP.hpp
        CommandLineFlags.hpp
        Version.hpp
        progressivePIP.hpp
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
        JATIApplication.hpp)

if (CMAKE_BUILD_TYPE EQUAL "Release")
    set_target_properties(miniJATI PROPERTIES LINK_FLAGS "-static")
endif ()

add_executable(miniJATI ${SOURCES})

target_link_libraries(miniJATI tshlib)
target_link_libraries(miniJATI glog)
target_link_libraries(miniJATI gflags)
target_link_libraries(miniJATI ${BPP_LIBS_SHARED})
target_link_libraries(miniJATI ${Boost_LIBRARIES})


SET(SOURCES_EXTBPP  Jati.cpp PIP.cpp PIP.hpp Version.hpp Jati.cpp)

add_executable(miniJATI_onBPP ${SOURCES_EXTBPP})

target_link_libraries(miniJATI_onBPP tshlib)
target_link_libraries(miniJATI_onBPP glog)
target_link_libraries(miniJATI_onBPP gflags)
target_link_libraries(miniJATI_onBPP ${BPP_LIBS_SHARED})

