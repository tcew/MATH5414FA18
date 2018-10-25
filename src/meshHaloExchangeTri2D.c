
#include "mesh.h"

void meshHaloExchangeTri2D(mesh_t *mesh,
			   void  *q,
			   int bytesPerElement){

  int rank, size;

  /* 
     NhaloElements
     haloElementIndices
     NhaloExchangeElements (number to send per rank)
  */
  
  int NhaloElements = mesh->NhaloElements;
  
  unsigned char *qout = (unsigned char*)
    malloc(NhaloElements*bytesPerElement);

  for(int h=0;h<NhaloElements;++h){
    int e = mesh->haloElementIndices[h];
    memcpy(qout+h*bytesPerElement,
	   (unsigned char*) q+e*bytesPerElement,
	   bytesPerElement);
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Request *sendRequests = (MPI_Request*)
    calloc(size, sizeof(MPI_Request));

  MPI_Request *recvRequests = (MPI_Request*)
    calloc(size, sizeof(MPI_Request));
  
  unsigned char *qin = (unsigned char*) q + bytesPerElement*mesh->Nelements;
  int cnt = 0;
  for(int r=0;r<size;++r){
    if(rank!=r){
      int Nexchange =
	bytesPerElement*mesh->NhaloExchangeElements[r];
      if(Nexchange>0){
	int tag = 999;
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
      int Nexchange =
	bytesPerElement*mesh->NhaloExchangeElements[r];
      if(Nexchange>0){
	MPI_Status status;
	MPI_Wait(sendRequests+r, &status);
	MPI_Wait(recvRequests+r, &status);
      }
    }
  }
  
  free(sendRequests);
  free(recvRequests);
  free(qout);
}


void meshHaloExtractTri2D(mesh_t *mesh,
			  int    NbytesPerElement,
			  occa::memory &o_q,
			  occa::memory &o_haloq,
			  dfloat *qout){

  int quarterNbytesPerElement = NbytesPerElement/4;

  printf("qNPE = %d, NHE=%d\n", quarterNbytesPerElement, mesh->NhaloElements);
  
  mesh->haloExtractKernel(mesh->NhaloElements,
			  quarterNbytesPerElement,
			  mesh->o_haloElementIndices,
			  o_q,
			  o_haloq);
  
  o_haloq.copyTo(qout);
  
}

void meshHaloInjectTri2D(mesh_t *mesh,
			 int    NbytesPerElement,
			 dfloat *qin,
			 occa::memory &o_q){
  
  int NhaloBytes = NbytesPerElement*mesh->NhaloElements;
  int offset = mesh->Nelements*NbytesPerElement;

  o_q.copyFrom(qin, NhaloBytes, offset);
  
}

void meshHaloExchangeStartTri2D(mesh_t *mesh,
				void   *haloq,
				void   *qin,
				int bytesPerElement,
				MPI_Request *sendRequests,
				MPI_Request *recvRequests){

  int rank, size;
  
  int NhaloElements = mesh->NhaloElements;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int cnt = 0;
  for(int r=0;r<size;++r){
    if(rank!=r){
      int Nexchange =
	bytesPerElement*mesh->NhaloExchangeElements[r];

      if(Nexchange>0){
	int tag = 999;
	MPI_Isend((unsigned char*)haloq+cnt, Nexchange, MPI_CHAR,
		  r, tag, MPI_COMM_WORLD, sendRequests+r);
	MPI_Irecv((unsigned char*)qin+cnt, Nexchange, MPI_CHAR,
		  r, tag, MPI_COMM_WORLD, recvRequests+r);
	cnt += Nexchange;
      }
    }
  }
}

void meshHaloExchangeFinishTri2D(mesh_t *mesh,
				 int bytesPerElement,
				 MPI_Request *sendRequests,
				 MPI_Request *recvRequests){

  int rank,size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  for(int r=0;r<size;++r){
    if(rank!=r){
      int Nexchange =
	bytesPerElement*mesh->NhaloExchangeElements[r];
      if(Nexchange>0){
	MPI_Status status;
	MPI_Wait(sendRequests+r, &status);
	MPI_Wait(recvRequests+r, &status);
      }
    }
  }
}

void meshHybridHaloExchangeTri2D(mesh_t *mesh,
				 int NbytesPerElement,
				 occa::memory &o_q,
				 occa::memory &o_haloq,
				 dfloat *haloqout,
				 dfloat *haloqin,
				 MPI_Request *sendRequests,
				 MPI_Request *recvRequests){



  if(mesh->NhaloElements){

    // extract halo data from DEVICE o_q to DEVICE o_haloq and copy to HOST haloq 
    meshHaloExtractTri2D(mesh,NbytesPerElement,o_q,o_haloq,haloqout);

    // start up requests to send data HOST<>HOST with MPI
    meshHaloExchangeStartTri2D(mesh, haloqout, haloqin, NbytesPerElement,sendRequests, recvRequests);
    
    // do something here
    
    // finalize MPI exchanges
    meshHaloExchangeFinishTri2D(mesh, NbytesPerElement, sendRequests,  recvRequests);
    
    // copy data from HOST qin to DEVICE o_q
    meshHaloInjectTri2D(mesh, NbytesPerElement, haloqin, o_q);

  }
}
