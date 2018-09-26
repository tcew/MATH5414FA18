
#include "mesh.h"

int main(int argc, char **argv){

  MPI_Init(&argc, &argv);
  
  mesh_t *mesh = meshParallelReaderTri2D(argv[1]);

  meshParallelConnectTri2D(mesh);

  meshHaloSetupTri2D(mesh);

  int *q =
    (int*) calloc(mesh->Nelements + mesh->NhaloElements,
		  sizeof(int));

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  for(int e=0;e<mesh->Nelements;++e){
    q[e] = rank;
  }

  meshHaloExchangeTri2D(mesh, q, sizeof(int));  
  
  mesh->EX = (dfloat*) realloc(mesh->EX,
			       (mesh->Nelements+mesh->NhaloElements)*
			       mesh->Nverts*sizeof(dfloat));

  mesh->EY = (dfloat*) realloc(mesh->EY,
			       (mesh->Nelements+mesh->NhaloElements)*
			       mesh->Nverts*sizeof(dfloat));

  meshHaloExchangeTri2D(mesh, mesh->EX, sizeof(dfloat)*mesh->Nverts);
  meshHaloExchangeTri2D(mesh, mesh->EY, sizeof(dfloat)*mesh->Nverts);
  
  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      int r = mesh->EToP[e*mesh->Nfaces+f];
      if(r!=-1){
	if(q[mesh->EToE[e*mesh->Nfaces+f]]!=r){
	  printf("Oh dear - got rank %d was expecting %d\n",
		 q[mesh->EToE[e*mesh->Nfaces+f]], r);
	}
      }
    }
  }
  
  meshVTUTri2D(mesh, "foo.vtu");
  
  MPI_Finalize();

  return 0;
}
