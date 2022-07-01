cmake_minimum_required(VERSION 3.14)

enable_testing()

find_package(CUDA REQUIRED)
enable_language(CUDA)
        
add_executable(benchmark_cuda bench.cu)
target_link_libraries(benchmark_cuda benchmark::benchmark)
