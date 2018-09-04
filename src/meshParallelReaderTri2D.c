#include "mesh.h"

mesh_t *meshParallelReaderTri2D(const char * fileName){

  int rank, size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh_t *localMesh = (mesh_t*) calloc(1, sizeof(mesh_t));
  localMesh->Nfaces = 3;
  localMesh->Nverts = 3;

  if(rank==0){
    mesh_t *globalMesh = meshReaderTri2D(fileName);

    // need to send Nelements, EX, EY, EToV, EToB for each rank
    int *rNelements = (int*) calloc(size, sizeof(int));
    for(int e=0;e<globalMesh->Nelements;++e){
      ++rNelements[e%size];
    }

    for(int r=0;r<size;++r){
      int Nel = rNelements[r];
      dfloat *EX = (dfloat*) calloc(Nel*globalMesh->Nverts, sizeof(dfloat));
      dfloat *EY = (dfloat*) calloc(Nel*globalMesh->Nverts, sizeof(dfloat));
      int *EToV = (int*) calloc(Nel*globalMesh->Nverts, sizeof(int));
      int *EToB = (int*) calloc(Nel*globalMesh->Nfaces, sizeof(int));

      int cnt = 0;
      for(int e=r;e<globalMesh->Nelements;e+=size){
	for(int v=0;v<globalMesh->Nverts;++v){
	  int idOut = cnt*globalMesh->Nverts+v;
	  EX[idOut] = globalMesh->VX[globalMesh->EToV[e*globalMesh->Nverts+v]];
	  EY[idOut] = globalMesh->VY[globalMesh->EToV[e*globalMesh->Nverts+v]];
	  EToV[idOut] = globalMesh->EToV[e*globalMesh->Nverts+v];
	}
	for(int f=0;f<globalMesh->Nfaces;++f){
	  int idOut = cnt*globalMesh->Nfaces+f;
	  EToB[idOut] = globalMesh->EToB[e*globalMesh->Nfaces+f];
	}
	++cnt;
      }

      if(r>0){
	MPI_Send(&Nel, 1, MPI_INT,
		 r, 999, MPI_COMM_WORLD);
	
	// send out element vertex coords
	MPI_Send(EX, Nel*globalMesh->Nverts, MPI_DFLOAT,
		 r, 999, MPI_COMM_WORLD);
	MPI_Send(EY, Nel*globalMesh->Nverts, MPI_DFLOAT,
		 r, 999, MPI_COMM_WORLD);
	MPI_Send(EToV, Nel*globalMesh->Nverts, MPI_INT,
		 r, 999, MPI_COMM_WORLD);
	MPI_Send(EToB, Nel*globalMesh->Nfaces, MPI_INT,
		 r, 999, MPI_COMM_WORLD);
	
	free(EX);
	free(EY);
	free(EToV);
	free(EToB);
      }
      else{
	localMesh->Nelements = Nel;
	localMesh->EX = EX;
	localMesh->EY = EY;
	localMesh->EToV = EToV;
	localMesh->EToB = EToB;
      }
    }
  }
  else{
    int Nel;
    int source = 0;
    MPI_Status status;

    MPI_Recv(&Nel, 1, MPI_INT,
	     source, 999, MPI_COMM_WORLD, &status);

    dfloat *EX = (dfloat*) calloc(Nel*localMesh->Nverts, sizeof(dfloat));
    dfloat *EY = (dfloat*) calloc(Nel*localMesh->Nverts, sizeof(dfloat));
    int *EToV = (int*) calloc(Nel*localMesh->Nverts, sizeof(int));
    int *EToB = (int*) calloc(Nel*localMesh->Nfaces, sizeof(int));
    
    // send out element vertex coords
    MPI_Recv(EX, Nel*localMesh->Nverts, MPI_DFLOAT,
	     source, 999, MPI_COMM_WORLD, &status);
    MPI_Recv(EY, Nel*localMesh->Nverts, MPI_DFLOAT,
	     source, 999, MPI_COMM_WORLD, &status);
    MPI_Recv(EToV, Nel*localMesh->Nverts, MPI_INT,
	     source, 999, MPI_COMM_WORLD, &status);
    MPI_Recv(EToB, Nel*localMesh->Nfaces, MPI_INT,
	     source, 999, MPI_COMM_WORLD, &status);
    
    localMesh->Nelements = Nel;
    localMesh->EX = EX;
    localMesh->EY = EY;
    localMesh->EToV = EToV;
    localMesh->EToB = EToB;
  }

  printf("rank %d got %d elements \n", rank, localMesh->Nelements);

  return localMesh;
}
