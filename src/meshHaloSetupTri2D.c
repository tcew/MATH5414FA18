
#include "mesh.h"

typedef struct {

  int element;

  int elementN; // neighbor element (in other rank local indices)
  int faceN;    // neighbor face
  int rankN;    // neighbor rank

} haloFace_t;

int compareHaloFaces(const void *a, const void *b){

  const haloFace_t *faceA = (haloFace_t*) a;
  const haloFace_t *faceB = (haloFace_t*) b;

  if(faceA->rankN < faceB->rankN)
    return -1;
  if(faceA->rankN > faceB->rankN)
    return +1;

  if(faceA->elementN < faceB->elementN)
    return -1;
  if(faceA->elementN > faceB->elementN)
    return +1;

  if(faceA->faceN < faceB->faceN)
    return -1;
  if(faceA->faceN > faceB->faceN)
    return +1;

  return 0;
}

void meshHaloSetupTri2D(mesh_t *mesh){

  int rank, size;

  /* 
     NhaloElements
     haloElementIndices
     NhaloExchangeElements (number to send per rank)
  */

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int NhaloElements = 0;
  int *NhaloExchangeElements =
    (int*) calloc(size, sizeof(int));

  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      int r = mesh->EToP[e*mesh->Nfaces+f];
      if(r!=-1){
	++NhaloElements;
	++(NhaloExchangeElements[r]);
      }
    }
  }

  haloFace_t *haloFaces =
    (haloFace_t*) calloc(NhaloElements, sizeof(haloFace_t));

  int cnt = 0;

  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      int r = mesh->EToP[e*mesh->Nfaces+f];
      if(r!=-1){
	haloFaces[cnt].element = e;
	haloFaces[cnt].elementN = mesh->EToE[e*mesh->Nfaces+f];
	haloFaces[cnt].faceN =  mesh->EToF[e*mesh->Nfaces+f];
	haloFaces[cnt].rankN =  mesh->EToP[e*mesh->Nfaces+f];
	++cnt;
      }
    }
  }

  qsort(haloFaces, NhaloElements, sizeof(haloFace_t),
	compareHaloFaces);

  int *haloElementIndices
    = (int*) calloc(NhaloElements, sizeof(int));
  
  for(int n=0;n<NhaloElements;++n){
    
    fprintf(stdout, "haloFaces[%d].rankN/elementN/faceN = %d,%d,%d\n",
	    n, haloFaces[n].rankN,haloFaces[n].elementN, haloFaces[n].faceN);
    
    haloElementIndices[n] = haloFaces[n].element;
  }

  // shadow copy onto DEVICE
  mesh->o_haloElementIndices =
    mesh->device.malloc(NhaloElements*sizeof(int), haloElementIndices);
  
  mesh->NhaloElements = NhaloElements;
  mesh->NhaloExchangeElements = NhaloExchangeElements;
  mesh->haloElementIndices = haloElementIndices;
  
  /* hack EToE */
  cnt = mesh->Nelements;

  for(int r=0;r<size;++r){
    if(r!=rank){
      for(int e=0;e<mesh->Nelements;++e){
	for(int f=0;f<mesh->Nfaces;++f){
	  if(mesh->EToP[e*mesh->Nfaces+f]==r){
	    mesh->EToE[e*mesh->Nfaces+f] = cnt;
	    ++cnt;
	  }
	}
      }
    }
  }
}
