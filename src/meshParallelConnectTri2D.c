
#include "mesh.h"

typedef struct {

  int element;
  int face;
  int v1;
  int v2;
  int rank;

  int elementN; // neighbor element (in other rank local indices)
  int faceN;    // neighbor face
  int rankN;    // neighbor rank

} lonelyFace_t;

void meshParallelConnectTri2D(mesh_t *mesh){

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // serial connectivity (all local connections)
  meshConnectTri2D(mesh);

  // count the number of lonely faces being sent to each rank
  int *sendCounts = (int*) calloc(size, sizeof(int));

  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      if(mesh->EToE[e*mesh->Nfaces + f]==-1){
	int v1 = mesh->EToV[e*mesh->Nverts + f];
	int v2 = mesh->EToV[e*mesh->Nverts + (f+1)%mesh->Nverts];
	int targetRank = min(v1,v2)%size;
	sendCounts[targetRank] += sizeof(lonelyFace_t);
      }
    }
  }

  int *sendDispls = (int*) calloc(size+1, sizeof(int));
  sendDispls[0] = 0;
  for(int r=0;r<size;++r){
    sendDispls[r+1] = sendDispls[r]+sendCounts[r];
    sendCounts[r] = 0;
  }

  lonelyFace_t *sendBuf = (lonelyFace_t*) malloc(sendDispls[size]);

  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      if(mesh->EToE[e*mesh->Nfaces + f]==-1){
	int v1 = mesh->EToV[e*mesh->Nverts + f];
	int v2 = mesh->EToV[e*mesh->Nverts + (f+1)%mesh->Nverts];
	int targetRank = min(v1,v2)%size;

	lonelyFace_t lonely;

	lonely.element = e;
	lonely.face = f;
	lonely.v1 = v1;
	lonely.v2 = v2;
	lonely.rank =  rank;

	int id = sendDispls[targetRank] + sendCounts[targetRank];
	id /= sizeof(lonelyFace_t);
	lonelyFaces[id] = lonely;
		    
	sendCounts[targetRank] += sizeof(lonelyFace_t);
      }
    }
  }

  // find out how many things to receive
  int* recvCounts = (int*) calloc(size, sizeof(int));
  MPI_Alltoall(sendCounts, 1, MPI_INT,
	       recvCounts, 1, MPI_INT,
	       MPI_COMM_WORLD);
	       
  int *recvDispls = (int*) calloc(size+1, sizeof(int));
  recvDispls[0] = 0;
  for(int r=0;r<size;++r)
    recvDispls[r+1] = recvDispls[r]+recvCounts[r];

  int NlonelyFaces = recvDispls[size]/sizeof(lonelyFace_t);
  
  lonelyFace_t *recvBuf =
    (lonelyFace_t*) malloc(NlonelyFaces*sizeof(lonelyFace_t));

  MPI_Alltoallv(sendBuf, sendCounts, sendDispls, MPI_CHAR,
		recvBuf, recvCounts, recvDispls, MPI_CHAR,
		MPI_COMM_WORLD);
  
  qsort(recvBuf, NlonelyFaces,
	sizeof(lonelyFace_t),
	compareLonelyFaces);

  for(int n=0;n<NlonelyFaces-1;++n){
    if(compareLonelyFaces(recvBuf+n, recvBuf+n+1)==0){
      recvBuf[n].elementN = recvBuf[n+1].element;
      recvBuf[n].faceN = recvBuf[n+1].face;
      recvBuf[n].rankN = recvBuf[n+1].rank;
    }
  }
  
}
