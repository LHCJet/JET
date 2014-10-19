cmake_minimum_required (VERSION 2.6)
project (SlowJet)
#add_definitions(-std=c++11 -DDEBUG -g -ggdb)
add_definitions(-std=c++11 -O2 -D_USE_MATH_DEFINES -march=native)
add_executable(slowjet
    src/example.cpp
    src/JetDefinition.cpp
    src/JetFinder.cpp
    src/Vector.cpp)
