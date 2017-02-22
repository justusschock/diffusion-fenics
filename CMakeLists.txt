
cmake_minimum_required(VERSION 3.1)

set(build_dir "${CMAKE_CURRENT_LIST_DIR}/build/")

SET(_flags_all -Wall -Wextra -pedantic -std=c++14)
add_subdirectory(src)

file(COPY test_mesh DESTINATION ${build_dir})

