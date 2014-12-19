set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")
include_directories("${CMAKE_SOURCE_DIR}/src")

find_package(Cython)

if(CYTHON_FOUND)
    include(UseCython)
    execute_process ( COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "Python module will be installed to:
    ${PYTHON_SITE_DIR}")
    set_source_files_properties(jetjet.pyx PROPERTIES CYTHON_IS_CXX TRUE)
    cython_add_module(jetjet jetjet.pyx)
    target_link_libraries(jetjet JETJet)
    install(TARGETS jetjet DESTINATION ${PYTHON_SITE_DIR})
endif(CYTHON_FOUND)

