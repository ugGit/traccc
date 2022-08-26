cmake_minimum_required(VERSION 3.14)

include(FetchContent)

FetchContent_Declare(googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG main)
        
FetchContent_MakeAvailable(
        googlebenchmark)
