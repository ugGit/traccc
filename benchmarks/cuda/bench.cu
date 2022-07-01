#include <benchmark/benchmark.h>
#include "matrix_sums.cu"

const size_t MATRIX_SIZE_MIN = 2<<4; // min matrix side dimension
const size_t MATRIX_SIZE_MAX = 2<<13;// max matrix side dimension
const int BLOCK_SIZE_MIN = 1;
const int BLOCK_SIZE_MAX = 2 << 9;  // CUDA maximum is 1024 (i.e. 2<<9)

static void BM_CudaMatrixSum(benchmark::State& state){
  for (auto _ : state){
    cuda::matrix_sums(state.range(0), state.range(1));
  }
}

// Register the function as a benchmark
BENCHMARK(BM_CudaMatrixSum)
  ->ArgsProduct({
    // range for matrix arguments
    {MATRIX_SIZE_MIN, MATRIX_SIZE_MAX},
    // range for thread block size
    {benchmark::CreateRange(BLOCK_SIZE_MIN, BLOCK_SIZE_MAX, /*multi=*/2)}
  });

BENCHMARK_MAIN();
