project(bayestracking)
cmake_minimum_required(VERSION 2.6)
set(CPACK_PACKAGE_NAME "bayestracking")
set(VERSION "1.0.0")

include_directories($ENV{HOME}/local/include)
link_directories($ENV{HOME}/local/lib /usr/local/lib)

add_subdirectory(src)

OPTION (BAYESTRACKING_BUILD_EXAMPLES "Enables compilation of the examples" OFF)
if(BAYESTRACKING_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif(BAYESTRACKING_BUILD_EXAMPLES)

OPTION (BAYESTRACKING_BUILD_DOC "Generates API documentation" OFF)
if(BAYESTRACKING_BUILD_DOC)
  # check if doxygen is even installed
  find_package(Doxygen)
  if (DOXYGEN_FOUND STREQUAL "NO")
      message(FATAL_ERROR "Doxygen not found. Please get a copy http://www.doxygen.org")
  endif (DOXYGEN_FOUND STREQUAL "NO")

  # prepare doxygen configuration file
  configure_file(${CMAKE_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

  # add doxygen as target
  add_custom_target(doc ALL ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  COMMENT "Generating API documentation with Doxygen" VERBATIM)

  # cleanup $build/doc on "make clean"
  set_property(DIRECTORY APPEND PROPERTY
        ADDITIONAL_MAKE_CLEAN_FILES doc)

  # install HTML API documentation and manual pages
  set(DOC_PATH "share/doc/${CPACK_PACKAGE_NAME}-${VERSION}")

  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html
    DESTINATION ${DOC_PATH}
  )

endif(BAYESTRACKING_BUILD_DOC)
