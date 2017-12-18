cmake_minimum_required(VERSION 3.9)
project(miniJATI)

set(CMAKE_CXX_STANDARD 11)

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
    include_directories("/usr/local/include/glog")

elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    message(STATUS "Compilation will be performed under Apple MacOS")
    find_package(glog 0.3.5 REQUIRED)
    find_package(gflags REQUIRED)
endif ()

if (UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lstdc++ -ldl -pthread ")
endif (UNIX)


find_package("Eigen3")
if (EIGEN3_FOUND)
    message("Found Eigen3 in: ${EIGEN3_INCLUDE_DIR}")
endif (EIGEN3_FOUND)
include_directories(${EIGEN3_INCLUDE_DIR})

add_subdirectory("../tshlib/src" ${CMAKE_CACHEFILE_DIR}/lib/tshlib)
include_directories("../tshlib/src")

set(SOURCES
        main.cpp)


add_executable(miniJATI ${SOURCES})

target_link_libraries(miniJATI tshlib)
target_link_libraries(miniJATI glog)
target_link_libraries(miniJATI gflags)