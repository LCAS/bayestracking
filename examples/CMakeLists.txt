include_directories(../src)

add_executable(simple_tracking simple_tracking.cpp)


find_package( OpenCV REQUIRED )

target_link_libraries(simple_tracking bayestracking ${OpenCV_LIBS})
