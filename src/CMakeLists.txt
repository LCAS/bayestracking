add_library(bayestracking STATIC associationmatrix.cpp ekfilter.cpp ukfilter.cpp pfilter.cpp trackwin.cpp models.cpp
BayesFilter/bayesFltAlg.cpp
BayesFilter/bayesFlt.cpp
BayesFilter/CIFlt.cpp
BayesFilter/covFlt.cpp
BayesFilter/infFlt.cpp
BayesFilter/infRtFlt.cpp
BayesFilter/itrFlt.cpp
BayesFilter/matSup.cpp
BayesFilter/SIRFlt.cpp
BayesFilter/UDFlt.cpp
BayesFilter/UdU.cpp
BayesFilter/unsFlt.cpp

)

install(TARGETS bayestracking
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${files} DESTINATION include/bayestracking)

file(GLOB bffiles "${CMAKE_CURRENT_SOURCE_DIR}/BayesFilter/*.hpp")
install(FILES ${bffiles} DESTINATION include/bayestracking/BayesFilter)
