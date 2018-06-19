cmake_minimum_required(VERSION 3.5)

set(MINIJATI_SOFTWARENAME "miniJATI")
set(MINIJATI_DESCRIPTION "Joint Alignment Tree Inference")

project(${MINIJATI_SOFTWARENAME}
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
    execute_process(COMMAND git describe --abbrev=0 --tags
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_TAG_VERSION)
    string(REGEX REPLACE "\n" "" GIT_TAG_VERSION_STRIPPED "${GIT_TAG_VERSION}")
else()
    execute_process(COMMAND "date" "+%d %b %Y/%H:%M:%S"
            OUTPUT_VARIABLE DATE_TIME)
    string(REGEX REPLACE "([^/]*)/.*" "\\1"
            PROJECT_BUILD_DATE "${DATE_TIME}")
    string(REGEX REPLACE "[^/]*/([0-9:]*).*" "\\1"
            PROJECT_BUILD_TIME "${DATE_TIME}")
    execute_process(COMMAND git describe --abbrev=0 --tags
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_TAG_VERSION)
    string(REGEX REPLACE "\n" "" GIT_TAG_VERSION_STRIPPED "${GIT_TAG_VERSION}")
endif()

set(MINIJATI_RELTIME ${PROJECT_BUILD_TIME})
set(MINIJATI_RELDATE ${PROJECT_BUILD_DATE})
set(MINIJATI_VERSION ${GIT_TAG_VERSION_STRIPPED})

project(${MINIJATI_SOFTWARENAME} VERSION ${GIT_TAG_VERSION_STRIPPED})


add_definitions(-DPRJ_GITBRANCH=\"${PROJECT_GIT_REF}\"
        -DPRJ_VERSION=\"${MINIJATI_VERSION}\"
        -DPRJ_GITREF=\"${PROJECT_SOURCE_VERSION}\"
        -DPRJ_DESC=\"${MINIJATI_DESCRIPTION}\"
        -DPRJ_NAME=\"${MINIJATI_SOFTWARENAME}\"
        -DPRJ_DATE=\"${MINIJATI_RELDATE}\"
        -DPRJ_TIME=\"${MINIJATI_RELTIME}\" )

set(CMAKE_CXX_STANDARD 14)

SET(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wshadow -Weffc++ -O0 -g")
SET(CMAKE_C_FLAGS_DEBUG "-Wall -Wshadow -O0 -g")

SET(CMAKE_CXX_FLAGS_DEBUGINTEL "-O0 -g -xHost")
SET(CMAKE_C_FLAGS_DEBUGINTEL "-Wall -Wshadow -O0 -g -xHost")

SET(CMAKE_CXX_FLAGS_RELWITHINFOINTEL "-O3 -g -xHost -DNDEBUG")
SET(CMAKE_C_FLAGS_RELWITHINFOINTEL "-O3 -g -xHost -DNDEBUG")

SET(CMAKE_CXX_FLAGS_RELEASE "-O3  -g0 -DNDEBUG")
SET(CMAKE_C_FLAGS_RELEASE "-O3  -g0 -DNDEBUG")

SET(CMAKE_CXX_FLAGS_INTEL "-O3  -g0 -xHost -DNDEBUG")
SET(CMAKE_C_FLAGS_INTEL "-O3  -g0 -xHost -DNDEBUG")

string(FIND "${CMAKE_BUILD_TYPE}" "static" STATIC)

message("-- Compilation will be performed with the following release of the software:
\tbranch   ${PROJECT_GIT_REF}
\tref      ${PROJECT_SOURCE_VERSION}
\ttime     ${MINIJATI_RELDATE} ${MINIJATI_RELTIME}
\tcurrent  ${MINIJATI_VERSION} (latest version)")

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

if (STATIC GREATER 0)
set(Boost_USE_STATIC_LIBS        ON) # only find static libs
else ()
set(Boost_USE_STATIC_LIBS        OFF) # only find dynamic libs
endif()
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost REQUIRED COMPONENTS system filesystem)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    message(STATUS "  static lib:  ${Boost_LIBRARY_DIRS}" )
    message(STATUS "  include dir: ${Boost_INCLUDE_DIRS}" )
else ()

if (STATIC GREATER 0)
   message(STATUS "You have requested to link the static version of the Boost libraries.")
   message(STATUS "However, no static version of the Boost libraries has been found. " )
   message(STATUS "If you downloaded Boost source files, run the following command: " )
   message(STATUS "./b2 release --link=static install " )
endif()

endif()

find_package(TshLib 2.0.0 REQUIRED)
if (TSHLIB_FOUND)
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
        #TSHHomogeneousTreeLikelihood.cpp
        #TSHHomogeneousTreeLikelihood.hpp
        #TSHSearchable.hpp
        UnifiedTSHTopologySearch.cpp
        UnifiedTSHTopologySearch.hpp
        JATIApplication.cpp
        JATIApplication.hpp
        pPIP.cpp
        pPIP.hpp
        Optimizators.cpp
        Optimizators.hpp
        UnifiedTSHomogeneousTreeLikelihood_Generic.cpp
        UnifiedTSHomogeneousTreeLikelihood_Generic.hpp
        UnifiedTSHomogeneousTreeLikelihood_PIP.cpp
        UnifiedTSHomogeneousTreeLikelihood_PIP.hpp
        UnifiedTSHSearchable.cpp
        UnifiedTSHSearchable.hpp
        RHomogeneousTreeLikelihood_Generic.cpp
        RHomogeneousTreeLikelihood_Generic.hpp
        SupportMeasures.cpp
        SupportMeasures.hpp
        progressivePIP.hpp
        progressivePIP.cpp)


if (STATIC GREATER 0)
    message(STATUS "Compilation will produce a static executable")
    add_executable(miniJATI ${SOURCES})
    set_target_properties(miniJATI PROPERTIES LINK_FLAGS "-static")
    target_link_libraries(miniJATI tshlib-static)
    target_link_libraries(miniJATI glog::glog)
    target_link_libraries(miniJATI ${BPP_LIBS_STATIC})
    target_link_libraries(miniJATI ${Boost_LIBRARIES})
else ()
    add_executable(miniJATI ${SOURCES})
    target_link_libraries(miniJATI tshlib)
    target_link_libraries(miniJATI glog::glog)
    target_link_libraries(miniJATI ${BPP_LIBS_SHARED})
    target_link_libraries(miniJATI ${Boost_LIBRARIES})
endif ()

# Test target
#enable_testing()
#add_subdirectory (tests)
