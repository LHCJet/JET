include_directories("${CMAKE_SOURCE_DIR}/src")
include_directories("${CMAKE_SOURCE_DIR}/fastjet")

add_executable(jetjet_example
    jetjet_example.cpp)

target_link_libraries(jetjet_example JETJet)

if(FASTJET_FOUND)
    include_directories(FASTJET_INCLUDE_DIRS)
    add_executable(fastjet_plugin_example
        fastjet_plugin_example.cpp)
    target_link_libraries(fastjet_plugin_example JETJetPlugin JETJet ${FASTJET_LIBRARIES})
endif()
