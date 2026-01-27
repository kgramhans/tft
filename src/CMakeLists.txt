cmake_minimum_required(VERSION 3.10)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


set(SOURCE_FILES
WaveletCalculator.cpp
WaveletTransformer.cpp
WaveletContainer.cpp
StftCalculator.cpp
DyadicFilter.cpp
WaveletVoice.cpp
ConfinedGaussianWaveletVoice.cpp
RealFFTf.cpp
)
add_library(tft ${SOURCE_FILES})
target_include_directories(tft PUBLIC "../include")

