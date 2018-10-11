
#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"
#include "occa.hpp"

int main(int argc, char **argv){

  MPI_Init(&argc, &argv);

  // initialize DEVICE
  occa::device device("mode: 'CUDA', device_id: 0");

  for(int Nmil=4;Nmil<640;Nmil+=4){
    
    //  measure on DEVICE memcpy bandwidth
    size_t Nbytes = Nmil*1000000;
    
    occa::memory o_A = device.malloc(Nbytes);
    occa::memory o_B = device.malloc(Nbytes);
    
    occa::streamTag startCopy = device.tagStream();
    
    // read from B and write to A
    o_A.copyFrom(o_B);
    
    occa::streamTag endCopy = device.tagStream();

    device.finish(); 
    
    double elapsed = device.timeBetween(startCopy, endCopy);
    double bandwidth = 2*Nbytes/( elapsed*1.e9);
      
    printf("measured bwidth for %lu bytes is %lg GB/s in elapsed time %lg\n",
	   Nbytes, bandwidth, elapsed); // read + write 
    
    o_A.free();
    o_B.free();
  }
  
  MPI_Finalize();

  return 0;
}
