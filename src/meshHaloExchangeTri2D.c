
#include "mesh.h"

void meshHaloExchangeTri2D(mesh_t *mesh,
			   unsigned char *q,
			   int bytesPerElement){

  int rank, size;
  
  int NhaloElements = mesh->NhaloElements;
  
  unsigned char *qout = (unsigned char)
    malloc(NhaloElements*bytesPerElement);

  for(int h=0;h<NhaloElements;++h){
    int e = mesh->haloElements[h];
    memcpy(qout+h*bytesPerElement,
	   q+e*bytesPerElement,
	   bytesPerElement);
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Request *sendRequests = (MPI_Request*)
    calloc(size, sizeof(MPI_Request));

  MPI_Request *recvRequests = (MPI_Request*)
    calloc(size, sizeof(MPI_Request));
  
  unsigned char *qin = q + bytesPerElement*mesh->Nelements;
  int cnt = 0;
  for(int r=0;r<size;++r){
    if(rank!=r){
      int Nexchange =
	bytesPerElement*mesh->NhaloExchangeElements[r];
      if(Nexchange>0){
	MPI_Isend(qout+cnt, Nexchange, MPI_CHAR,
		  r, tag, MPI_COMM_WORLD, sendRequests+r);
	MPI_Irecv(qin+cnt, Nexchange, MPI_CHAR,
		  r, tag, MPI_COMM_WORLD, recvRequests+r);
	cnt += Nexchange;
      }
    }
  }

  for(int r=0;r<size;++r){
    if(rank!=r){
      MPI_Status status;
      MPI_Wait(sendRequests+r, &status);
      MPI_Wait(recvRequests+r, &status);
    }
  }

  free(sendRequests);
  free(recvRequests);
  free(qout);
}
