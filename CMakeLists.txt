cmake_minimum_required(VERSION 3.5)

set(THSLIB_MAJOR_VERSION 1)
set(THSLIB_MINOR_VERSION 0)
set(THSLIB_PATCH_VERSION 0)
set(THSLIB_VERSION ${THSLIB_MAJOR_VERSION}.${THSLIB_MINOR_VERSION}.${THSLIB_PATCH_VERSION})

project(miniJATI VERSION ${THSLIB_VERSION}
        DESCRIPTION "Joint Alignment Tree Inferecne -- Minimal Version"
        LANGUAGES "CXX")


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
    include_directories("/usr/local/include/glog")

elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    message(STATUS "Compilation will be performed under Apple MacOS")
    find_package(glog 0.3.5 REQUIRED)
    find_package(gflags REQUIRED)
endif ()




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



set(SOURCES
        main.cpp utils.cpp utils.hpp PIP.cpp PIP.hpp cli_parser.hpp Version.hpp)


add_executable(miniJATI ${SOURCES})

target_link_libraries(miniJATI tshlib)
target_link_libraries(miniJATI glog)
target_link_libraries(miniJATI gflags)
target_link_libraries(miniJATI ${BPP_LIBS_SHARED})


SET(SOURCES_EXTBPP  Jati.cpp PIP.cpp PIP.hpp Version.hpp Jati.cpp)

add_executable(miniJATI_onBPP ${SOURCES_EXTBPP})

target_link_libraries(miniJATI_onBPP tshlib)
target_link_libraries(miniJATI_onBPP glog)
target_link_libraries(miniJATI_onBPP gflags)
target_link_libraries(miniJATI_onBPP ${BPP_LIBS_SHARED})
