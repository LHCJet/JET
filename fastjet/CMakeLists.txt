set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")
include_directories("${CMAKE_SOURCE_DIR}/src")

find_package(FastJet)

if(FASTJET_FOUND)
    include_directories(FASTJET_INCLUDE_DIRS)
    add_library(JETJetPlugin SHARED
        JETJetPlugin.cpp)
    add_executable(fastjet_plugin_example
        example.cpp)
    target_link_libraries(fastjet_plugin_example JETJetPlugin JETJet ${FASTJET_LIBRARIES})
else(FASTJET_FOUND)
    message(STATUS "FastJet not found, plugin will not be built!")
endif(FASTJET_FOUND)

