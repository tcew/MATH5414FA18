#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char **argv){
  
  // initialize MPI: always do this "early"
  MPI_Init(&argc, &argv);

  int rank;
  int size;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  printf(" This is process %d of %d total processes\n",
	 rank, size);

  // collective barrier
  MPI_Barrier(MPI_COMM_WORLD);

  if(rank==0) printf("hello from bob\n");

  // send a "message" from rank 0 to rank 1
  if(rank==0){
    int messageN = 10;
    int *messageOut = (int*) 
      calloc(messageN, sizeof(int));
    for(int i=0;i<messageN;++i){
      messageOut[i] = i;
    }
    int dest = 1;
    int tag = 666;
    MPI_Send(messageOut, 
	     messageN,
	     MPI_INT,
	     dest,
	     tag,
	     MPI_COMM_WORLD);

    free(messageOut);
  }

  // rank 1 receives a "message" from rank 0
  if(rank==1){
    int tag = 666;
    int messageN = 10;
    int *messageIn = (int*) 
      calloc(messageN, sizeof(int));

    int source = 0;
    MPI_Status status;

    MPI_Recv(messageIn, 
	     messageN,
	     MPI_INT,
	     source,
	     tag,
	     MPI_COMM_WORLD,
	     &status);

    for(int i=0;i<messageN;++i){
      printf("messageIn[%d] = %d\n",
	     i, messageIn[i]);
    }

    free(messageIn);

  }

  // tear down MPI (potentially block for all proceses)
  MPI_Finalize();

  return 0;
}
