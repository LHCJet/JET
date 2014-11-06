cmake_minimum_required (VERSION 2.6)
project (JETJet)
#add_definitions(-std=c++11 -DDEBUG -g -ggdb)
#add_definitions(-std=c++11 -O2 -g -ggdb -DDEBUG -D_USE_MATH_DEFINES -march=native)
add_definitions(-std=c++11 -O2 -D_USE_MATH_DEFINES -march=native)
add_executable(jetjet
    src/example.cpp
    src/JetDefinition.cpp
    src/JetFinder.cpp
    src/Vector.cpp)
