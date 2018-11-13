#include <stdlib.h>
#include <stdio.h>

#include "occa.hpp"

int main(int argc, char **argv){

  occa::device device;

  char setupString[BUFSIZ];
  sprintf(setupString, "mode: 'CUDA', device_id: 0");
  device.setup(setupString);
  
  occa::kernel atomicAdd = device.buildKernel("src/occaAtomicAdd.okl",
					      "occaAtomicAdd");

  occa::kernel atomicReduction = device.buildKernel("src/occaReduction.okl",
						    "occaReduction");


  int N = atoi(argv[1]);

  int *h_a = (int*) calloc(N, sizeof(int));
  for(int n=0;n<N;++n)
    h_a[n] = 1;

  int B = (N+255)/256;
  int *output = (int*) calloc(B, sizeof(int));
  occa::memory o_output = device.malloc(B*sizeof(int), output);

  occa::memory o_a = device.malloc(N*sizeof(int), h_a);

  //  atomicAdd(o_output);
  atomicReduction(N, o_a, o_output);

  o_output.copyTo(output);
  printf("output[0] = %d\n", output[0]);
  
  return 0;
}
