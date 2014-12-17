cmake_minimum_required (VERSION 2.6)
project (JETJet)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
#add_definitions(-std=c++11 -DDEBUG -g -ggdb)
#add_definitions(-std=c++11 -O2 -g -ggdb -DDEBUG -D_USE_MATH_DEFINES -march=native)
add_definitions(-D_USE_MATH_DEFINES)

add_library(jetjet SHARED
    src/JetDefinition.cpp
    src/JetFinder.cpp
    src/Vector.cpp)

add_executable(jetjet_example
    src/example.cpp)

target_link_libraries(jetjet_example jetjet)
