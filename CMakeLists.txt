project(bayestracking)
cmake_minimum_required(VERSION 2.6)
set(CPACK_PACKAGE_NAME "bayestracking")
set(VERSION "1.0.0")

include_directories(${CMAKE_SOURCE_DIR}/src)

add_subdirectory(src)

OPTION (BAYESTRACKING_CREATE_PKGCONFIG "Create pkg-config file" ON)
if(BAYESTRACKING_CREATE_PKGCONFIG)
  message(STATUS "Created pkg-config file: ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig/bayestracking.pc")
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/config/bayestracking.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/config/bayestracking.pc
    @ONLY
  )
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/config/bayestracking.pc
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig/
  )
endif(BAYESTRACKING_CREATE_PKGCONFIG)

OPTION (BAYESTRACKING_CREATE_CMAKE_CONFIG "Create cmake config files" ON)
if(BAYESTRACKING_CREATE_CMAKE_CONFIG)
  message(STATUS "Created cmake config files: ${CMAKE_INSTALL_PREFIX}/lib/cmake/bayestracking")
  export(TARGETS bayestracking
    FILE "${CMAKE_CURRENT_BINARY_DIR}/config/BayestrackingTargets.cmake")
  export(PACKAGE Bayestracking)
  set(CONF_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/include/bayestracking")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/BayestrackingConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/config/BayestrackingConfig.cmake" @ONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/BayestrackingConfigVersion.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/config/BayestrackingConfigVersion.cmake" @ONLY)
  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/config/BayestrackingConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/config/BayestrackingConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/bayestracking" COMPONENT dev)
  install(EXPORT BayestrackingTargets DESTINATION
    "${CMAKE_INSTALL_PREFIX}/lib/cmake/bayestracking" COMPONENT dev)
endif(BAYESTRACKING_CREATE_CMAKE_CONFIG)

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
