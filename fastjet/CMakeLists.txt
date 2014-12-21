include_directories("${CMAKE_SOURCE_DIR}/src")

include_directories(FASTJET_INCLUDE_DIRS)
add_library(JETJetPlugin SHARED
    JETJetPlugin.cpp)
install(TARGETS JETJetPlugin DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES JETJetPlugin.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/jetjet)
