# CMakeLists.txt for helper

# Create a library for helper (actual code)
add_library(helper STATIC helper.cc)

# Specify include directories for helper
target_include_directories(helper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}) # helper's own headers
    #${CMAKE_SOURCE_DIR}/bla      # bla's headers


# Link the bla library
#target_link_libraries(helper PUBLIC bla)

