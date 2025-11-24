# CMakeLists.txt for nonlinfunc

# Create a library for nonlinfunc (actual code)
add_library(nonlinfunc nonlinfunc.cc Newton.cc)

# Specify include directories for nonlinfunc
target_include_directories(nonlinfunc PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}  # nonlinfunc's own headers
    ${CMAKE_SOURCE_DIR}/bla      # bla's headers
    ${CMAKE_SOURCE_DIR}/helper   # helper's headers
)

# Link dependencies
target_link_libraries(nonlinfunc PUBLIC
    #bla
    helper
)
