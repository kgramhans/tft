cmake_minimum_required(VERSION 3.14)
project(tft_gtest)

# GoogleTest requires at least C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  googletest
  DOWNLOAD_EXTRACT_TIMESTAMP
  URL https://github.com/google/googletest/archive/refs/tags/v1.17.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

#FetchContent_Declare(
#  tft
#  GIT_REPOSITORY https://github.com/kgramhans/tft.git
#  GIT_TAG main
#)
#FetchContent_MakeAvailable(tft)

enable_testing()

add_executable(
  tft_gtest
  test_polygonregion.cpp
  test_dyadicfilter.cpp
  test_waveletvoice_unbuffered.cpp
  test_waveletvoice_buffered.cpp
  test_wavelettransformer.cpp
)
target_include_directories(tft_gtest PUBLIC "../include")
target_link_libraries(
  tft_gtest
  tft
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(tft_gtest)
