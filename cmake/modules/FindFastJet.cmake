find_path(FASTJET_INCLUDE_DIR NAMES fastjet/JetDefinition.hh)
mark_as_advanced(FASTJET_INCLUDE_DIR)

# Look for the library (sorted from most current/relevant entry to least).
find_library(FASTJET_LIBRARY NAMES fastjet)

mark_as_advanced(FASTJET_LIBRARY)

# handle the QUIETLY and REQUIRED arguments and set FASTJET_FOUND to TRUE if
# all listed variables are TRUE

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FASTJET REQUIRED_VARS FASTJET_LIBRARY FASTJET_INCLUDE_DIR)

if(FASTJET_FOUND)
  set(FASTJET_LIBRARIES ${FASTJET_LIBRARY})
  set(FASTJET_INCLUDE_DIRS ${FASTJET_INCLUDE_DIR})
endif()
