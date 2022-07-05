cmake_minimum_required(VERSION 3.14)

include(FetchContent)

## Project-wide setup
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)

FetchContent_Declare(googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG main)
        
FetchContent_MakeAvailable(
        googlebenchmark)

# Add CUDA benchmarks
add_subdirectory(cuda)
