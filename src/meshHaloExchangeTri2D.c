
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
  
  mesh->haloExtractKernel(mesh->NhaloElements,
			  quarterNbytesPerElement,
			  o_q,
			  o_haloq);

  o_haloq.copyTo(qout);
  
}

void meshHaloInjectTri2D(mesh_t *mesh,
			 int    NbytesPerElement,
			 dfloat *qin,
			 occa::memory &o_q){
  
  int NhaloBytes = NbytesPerElement*mesh->NhaloElements;
  
  o_q.copyFrom(qin, NhaloBytes, offset);
  
}

void meshHaloExchangeStartTri2D(mesh_t *mesh,
				void   *q,
				void   *haloq,
				int bytesPerElement,
				MPI_Requests *sendRequests,
				MPI_Requests *recvRequests){

  int rank, size;
  
  int NhaloElements = mesh->NhaloElements;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  unsigned char *qin = haloq;
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
}

void meshHaloExchangeFinishTri2D(mesh_t *mesh,
				 int bytesPerElement,
				 MPI_Requests *sendRequests,
				 MPI_Requests *recvRequests){

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
				 dfloat *haloq,
				 dfloat *qin,
				 MPI_Requests *sendRequests,
				 MPI_Requests *recvRequests){

  meshHaloExtractTri2D(mesh,NbytesPerElement,o_q,o_haloq,haloq);

  meshHaloExchangeStartTri2D(mesh, haloq, NbytesPerElement,sendRequests,
			     recvRequests);

  // do something here
  
  meshHaloExchangeFinishTri2D(mesh,
			      NbytesPerElement,
			      sendRequests,
			      recvRequests);

  meshHaloInjectTri2D(mesh,
		      NbytesPerElement,
		      qin,
		      o_q);
  
}
