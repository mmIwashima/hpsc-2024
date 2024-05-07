#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void sort(int *key, int *bucket, int range, int n){
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  for (int i=0; i<range; i++) {
    bucket[i] = 0;
    __syncthreads();
  }
  for (int i=0; i<n; i++) {
    bucket[key[i]]++;
    __syncthreads();
  }
  for (int i=0, j=0; i<range; i++) {
    for (; bucket[i]>0; bucket[i]--) {
      key[j++] = i;
      __syncthreads();
    }
  }
}

int main() {
  int n = 50;
  int range = 5;
//  std::vector<int> key(n);
  int *key;
  cudaMallocManaged(&key, n*sizeof(int));

  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

//  std::vector<int> bucket(range);
  int *bucket; 
  cudaMallocManaged(&bucket, range*sizeof(int));
  
//    for (int i=0; i<range; i++) {
//    bucket[i] = 0;
//  }
//  for (int i=0; i<n; i++) {
//    bucket[key[i]]++;
//  }
//  for (int i=0, j=0; i<range; i++) {
//    for (; bucket[i]>0; bucket[i]--) {
//      key[j++] = i;
//    }
//  }

  sort<<<1,n>>>(key, bucket, range, n);
  cudaDeviceSynchronize();

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
  cudaFree(key);
  cudaFree(bucket);
}
