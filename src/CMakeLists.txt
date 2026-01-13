cmake_minimum_required(VERSION 3.10)

# set the project name
project(tft VERSION 0.4)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(version.h.in version.h)
set(SOURCE_FILES
WaveletCalculator.cpp 
StftCalculator.cpp 
DyadicFilter.cpp 
WaveletVoice.cpp
WaveletVoiceUnbuffered.cpp
ConfinedGaussianWaveletVoice.cpp
RealFFTf.cpp

)
add_library(tft ${SOURCE_FILES})
target_include_directories(tft
PUBLIC
  ${CMAKE_CURRENT_LIST_DIR}/include
PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/include/tft
  ${CMAKE_CURRENT_BINARY_DIR}/src
)
