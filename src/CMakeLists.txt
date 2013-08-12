add_library(bayestracking STATIC associationmatrix.cpp ekfilter.cpp ukfilter.cpp pfilter.cpp trackwin.cpp models.cpp)

install(TARGETS bayestracking
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${files} DESTINATION include/bayestracking)
