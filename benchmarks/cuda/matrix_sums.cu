#include <stdio.h>
#include <iostream>
#include <chrono>

// error checking macro
#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)

namespace cuda{

// matrix row-sum kernel
__global__ void row_sums(const float *A, float *sums, size_t ds){

  int idx = threadIdx.x + blockIdx.x * blockDim.x; // create typical 1D thread index from built-in variables
  if (idx < ds){
    float sum = 0.0f;
    for (size_t i = 0; i < ds; i++)
      sum += A[idx*ds + i];         // write a for loop that will cause the thread to iterate across a row, keeeping a running sum, and write the result to sums
    sums[idx] = sum;
}}

// matrix column-sum kernel
__global__ void column_sums(const float *A, float *sums, size_t ds){

  int idx = threadIdx.x + blockIdx.x * blockDim.x; // create typical 1D thread index from built-in variables
  if (idx < ds){
    float sum = 0.0f;
    for (size_t i = 0; i < ds; i++)
      sum += A[idx + i*ds];         // write a for loop that will cause the thread to iterate down a column, keeeping a running sum, and write the result to sums
    sums[idx] = sum;
}}

bool validate(float *data, size_t sz){
  for (size_t i = 0; i < sz; i++)
    if (data[i] != (float)sz) {printf("results mismatch at %lu, was: %f, should be: %f\n", i, data[i], (float)sz); return false;}
    return true;
}

using std::cout;

int matrix_sums(int matrix_size, int block_size){
  // start crono
  // const auto t1 = std::chrono::high_resolution_clock::now();

  float *h_A, *h_sums, *d_A, *d_sums;
  h_A = new float[matrix_size*matrix_size];  // allocate space for data in host memory
  h_sums = new float[matrix_size]();
    
  for (int i = 0; i < matrix_size*matrix_size; i++)  // initialize matrix in host memory
    h_A[i] = 1.0f;
    
  cudaMalloc(&d_A, matrix_size*matrix_size*sizeof(float));  // allocate device space for A
  cudaMalloc(&d_sums, matrix_size*sizeof(float)); // allocate device space for vector d_sums
  cudaCheckErrors("cudaMalloc failure"); // error checking
    
  // copy matrix A to device:
  cudaMemcpy(d_A, h_A, matrix_size*matrix_size*sizeof(float), cudaMemcpyHostToDevice);
  cudaCheckErrors("cudaMemcpy H2D failure");
  
  cudaMemset(d_sums, 0, matrix_size*sizeof(float));
    
  column_sums<<<(matrix_size+block_size-1)/block_size, block_size>>>(d_A, d_sums, matrix_size);
  cudaCheckErrors("kernel launch failure");
  //cuda processing sequence step 2 is complete
    
  // copy vector sums from device to host:
  cudaMemcpy(h_sums, d_sums, matrix_size*sizeof(float), cudaMemcpyDeviceToHost);
  //cuda processing sequence step 3 is complete
  cudaCheckErrors("kernel execution failure or cudaMemcpy H2D failure");

  // stop crono
  // const auto t2 = std::chrono::high_resolution_clock::now();
  
  return 0;
}
  
}
