cmake_minimum_required (VERSION 2.8.5 FATAL_ERROR)
project (JETJet)
set(JET_VERSION_MAJOR 0)
set(JET_VERSION_MINOR 2)
set(JET_VERSION "${JET_VERSION_MAJOR}.${JET_VERSION_MINOR}" )


# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

option(JET_ENABLE_FASTJET "Build FastJet Plugin" ON)
option(JET_ENABLE_PYTHON "Build Python Module" ON)

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)

CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")

list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" _isSystemPlatformLibDir)
list(FIND CMAKE_C_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" _isSystemCLibDir)
list(FIND CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" _isSystemCxxLibDir)
if("${_isSystemPlatformLibDir}" STREQUAL "-1" AND "${_isSystemCLibDir}" STREQUAL "-1" AND "${_isSystemCxxLibDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
endif("${_isSystemPlatformLibDir}" STREQUAL "-1" AND "${_isSystemCLibDir}" STREQUAL "-1" AND "${_isSystemCxxLibDir}" STREQUAL "-1")

set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(JET_ENABLE_FASTJET)
    add_subdirectory(fastjet)
endif(JET_ENABLE_FASTJET)
if(JET_ENABLE_PYTHON)
    add_subdirectory(python)
endif(JET_ENABLE_PYTHON)
#add_definitions(-std=c++11 -DDEBUG -g -ggdb)
#add_definitions(-std=c++11 -O2 -g -ggdb -DDEBUG -D_USE_MATH_DEFINES -march=native)
add_definitions(-D_USE_MATH_DEFINES)

add_library(JETJet SHARED
    src/JetDefinition.cpp
    src/JetFinder.cpp
    src/Vector.cpp)

add_executable(jetjet_example
    src/example.cpp)

target_link_libraries(jetjet_example JETJet)

install(TARGETS JETJet DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(FILES src/Vector.h src/JetDefinition.h src/JetFinder.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/jetjet)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
add_custom_target(uninstall
    "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")
