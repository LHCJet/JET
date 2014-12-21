include_directories("${CMAKE_SOURCE_DIR}/src")

include_directories(FASTJET_INCLUDE_DIRS)
add_library(JETJetPlugin SHARED
    JETJetPlugin.cpp)
add_executable(fastjet_plugin_example
    example.cpp)
target_link_libraries(fastjet_plugin_example JETJetPlugin JETJet ${FASTJET_LIBRARIES})
install(TARGETS JETJetPlugin DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES JETJetPlugin.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/jetjet)
