
#include "mesh.h"

// lonely face: meta data associated with a disconnected face
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

int compareLonelyFaces(const void *a, const void *b){

  const lonelyFace_t *faceA = (lonelyFace_t*) a;
  const lonelyFace_t *faceB = (lonelyFace_t*) b;

    // assume that v1 is the 2nd digit
  if(faceA->v1 < faceB->v1)
    return -1;
  if(faceA->v1 > faceB->v1)
    return +1;

  if(faceA->v2 < faceB->v2)
    return -1;
  if(faceA->v2 > faceB->v2)
    return +1;

  return 0;
}


int compareRank(const void *a, const void *b){

  const lonelyFace_t *faceA = (lonelyFace_t*) a;
  const lonelyFace_t *faceB = (lonelyFace_t*) b;

    // assume that v1 is the 2nd digit
  if(faceA->rank < faceB->rank)
    return -1;
  if(faceA->rank > faceB->rank)
    return +1;

  return 0;
}
  

// find all element-to-element connections across all processes
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
	int targetRank = mymin(v1,v2)%size;
	// count in bytes
	sendCounts[targetRank] += sizeof(lonelyFace_t);
      }
    }
  }

  // construct the cumulative sum of sendCounts and reset sendCounts
  int *sendDispls = (int*) calloc(size+1, sizeof(int));
  sendDispls[0] = 0;
  for(int r=0;r<size;++r){
    sendDispls[r+1] = sendDispls[r]+sendCounts[r];
    sendCounts[r] = 0;
  }

  // collect the outgoing lonelyFaces 
  lonelyFace_t *sendBuf = (lonelyFace_t*) malloc(sendDispls[size]);

  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      if(mesh->EToE[e*mesh->Nfaces + f]==-1){
	int v1 = mesh->EToV[e*mesh->Nverts + f];
	int v2 = mesh->EToV[e*mesh->Nverts + (f+1)%mesh->Nverts];
	int targetRank = mymin(v1,v2)%size;

	// meta data for this lonelyFace
	lonelyFace_t lonely;
	lonely.element = e;
	lonely.face = f;
	lonely.v1 = mymin(v1,v2);
	lonely.v2 = mymax(v1,v2);
	lonely.rank =  rank;

	lonely.elementN = -111;
	lonely.faceN = -111;
	lonely.rankN = -111;
	
	// identify where this lonelyFace goes in the outgoing buffer
	int id = sendDispls[targetRank] + sendCounts[targetRank];
	id /= sizeof(lonelyFace_t); // adjust for bytes
	sendBuf[id] = lonely;

	// increment sendCounts again
	sendCounts[targetRank] += sizeof(lonelyFace_t);
      }
    }
  }

  // find out how many lonelyFaces to receive from each process
  int* recvCounts = (int*) calloc(size, sizeof(int));
  MPI_Alltoall(sendCounts, 1, MPI_INT,
	       recvCounts, 1, MPI_INT,
	       MPI_COMM_WORLD);

  // form cumulative sum of recvCounts
  int *recvDispls = (int*) calloc(size+1, sizeof(int));
  recvDispls[0] = 0;
  for(int r=0;r<size;++r)
    recvDispls[r+1] = recvDispls[r]+recvCounts[r];

  // count how many lonelyFaces you are going to receive
  int NlonelyFaces = recvDispls[size]/sizeof(lonelyFace_t);
  
  lonelyFace_t *recvBuf =
    (lonelyFace_t*) malloc(NlonelyFaces*sizeof(lonelyFace_t));

  // exchange lonelyFaces with all other ranks
  MPI_Alltoallv(sendBuf, sendCounts, sendDispls, MPI_CHAR,
		recvBuf, recvCounts, recvDispls, MPI_CHAR,
		MPI_COMM_WORLD);

  // sort based on [v1,v2]
  qsort(recvBuf, NlonelyFaces,
	sizeof(lonelyFace_t),
	compareLonelyFaces);

  // fill in meta data for any found match
  for(int n=0;n<NlonelyFaces-1;++n){
    if(compareLonelyFaces(recvBuf+n, recvBuf+n+1)==0){
      recvBuf[n].elementN = recvBuf[n+1].element;
      recvBuf[n].faceN = recvBuf[n+1].face;
      recvBuf[n].rankN = recvBuf[n+1].rank;
      recvBuf[n+1].elementN = recvBuf[n].element;
      recvBuf[n+1].faceN = recvBuf[n].face;
      recvBuf[n+1].rankN = recvBuf[n].rank;
    }
  }

  // reorder to return
  qsort(recvBuf, NlonelyFaces,
	sizeof(lonelyFace_t),
	compareRank);

  // exchange lonelyFaces with all other ranks
  MPI_Alltoallv(recvBuf, recvCounts, recvDispls, MPI_CHAR,
		sendBuf, sendCounts, sendDispls, MPI_CHAR,
		MPI_COMM_WORLD);

  mesh->EToP = (int*) calloc(mesh->Nelements*mesh->Nfaces, sizeof(int));
  for(int n=0;n<mesh->Nelements*mesh->Nfaces;++n){
    mesh->EToP[n] = -1;
  }

  for(int n=0;n<sendDispls[size]/sizeof(lonelyFace_t);++n){
    int e = sendBuf[n].element;
    int f = sendBuf[n].face;
    int eN = sendBuf[n].elementN;
    int fN = sendBuf[n].faceN;
    int rN = sendBuf[n].rankN;

    if(eN>=0 && fN>=0){
      mesh->EToE[e*mesh->Nfaces+f] = eN;
      mesh->EToF[e*mesh->Nfaces+f] = fN;
      mesh->EToP[e*mesh->Nfaces+f] = rN;
    }

    if(sendBuf[n].rank!=rank){
      printf("WRONG RANK!!!\n");
      printf("sendBuf[%d].rank =%d\n", n, sendBuf[n].rank);
    }
  }

#if 0
  // for debugging
  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      if(mesh->EToP[e*mesh->Nfaces+f]>=0){
      printf("conn: %d,%d,%d => %d,%d,%d\n",
	     e,f, rank,
	     mesh->EToE[e*mesh->Nfaces+f],
	     mesh->EToF[e*mesh->Nfaces+f],
	     mesh->EToP[e*mesh->Nfaces+f]);
      }
    }
  }
#endif
}
