cmake_minimum_required(VERSION 3.10)

# set the project name
project(tft_test VERSION 0.1)


# add the executable
add_executable(tft_test main.cpp)
target_include_directories(tft_test PUBLIC "../include/tft")
target_link_libraries (tft_test tft)

