add_library(Elements STATIC Beam.cc Connector.cc Spring.cc rigid_body.cc)
target_include_directories(Elements PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Elements PRIVATE helper)