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

  int *output = (int*) calloc(1, sizeof(int));
  occa::memory o_output = device.malloc(1*sizeof(int), output);

  atomicAdd(o_output);

  o_output.copyTo(output);
  printf("output[0] = %d\n", output[0]);
  
  return 0;
}
