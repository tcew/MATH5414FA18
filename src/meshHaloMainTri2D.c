
#include "mesh.h"

int main(int argc, char **argv){

  MPI_Init(&argc, &argv);
  
  mesh_t *mesh = meshParallelReaderTri2D(argv[1]);

  // initialize DEVICE
  mesh->device.setup("mode: 'CUDA', device_id: 0");

  // load polynomial data for reference triangle
  int N = atoi(argv[2]); // use this for polynomial array ( use 2nd argument to define degree )
  meshLoadReferenceNodesTri2D(mesh, N);
  
  meshParallelConnectTri2D(mesh);

  meshHaloSetupTri2D(mesh);

  int *qrank = (int*) calloc(mesh->Nelements + mesh->NhaloElements, sizeof(int));

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  for(int e=0;e<mesh->Nelements;++e){
    qrank[e] = rank;
  }

  meshHaloExchangeTri2D(mesh, qrank, sizeof(int));  
  
  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      int r = mesh->EToP[e*mesh->Nfaces+f];
      if(r!=-1){
	if(qrank[mesh->EToE[e*mesh->Nfaces+f]]!=r){
	  printf("Oh dear - got rank %d was expecting %d\n",
		 qrank[mesh->EToE[e*mesh->Nfaces+f]], r);
	}
      }
    }
  }
  
  mesh->EX = (dfloat*) realloc(mesh->EX,
			       (mesh->Nelements+mesh->NhaloElements)*
			       mesh->Nverts*sizeof(dfloat));

  mesh->EY = (dfloat*) realloc(mesh->EY,
			       (mesh->Nelements+mesh->NhaloElements)*
			       mesh->Nverts*sizeof(dfloat));

  meshHaloExchangeTri2D(mesh, mesh->EX, sizeof(dfloat)*mesh->Nverts);
  meshHaloExchangeTri2D(mesh, mesh->EY, sizeof(dfloat)*mesh->Nverts);

  // compute geometric factors
  meshVolumeGeometricFactorsTri2D(mesh);

  // initialize HOST arrays
  dfloat *q     = (dfloat*) calloc(  mesh->Np*(mesh->Nelements+mesh->NhaloElements), sizeof(dfloat));

  for(int e=0;e<mesh->Nelements;++e){
    dfloat x0 = mesh->EX[e*mesh->Nverts+0];
    dfloat x1 = mesh->EX[e*mesh->Nverts+1];
    dfloat x2 = mesh->EX[e*mesh->Nverts+2];
    dfloat y0 = mesh->EY[e*mesh->Nverts+0];
    dfloat y1 = mesh->EY[e*mesh->Nverts+1];
    dfloat y2 = mesh->EY[e*mesh->Nverts+2];
    for(int n=0;n<mesh->Np;++n){
      dfloat rn =  mesh->r[n];
      dfloat sn =  mesh->s[n];
      q[e*mesh->Np+n] = -0.5*(rn+sn)*x0 + 0.5*(1+rn)*x1 + 0.5*(1+sn)*x2;
    }
  }

  // allocate array space on DEVICE
  occa::memory o_q = mesh->device.malloc(mesh->Np*(mesh->Nelements+mesh->NhaloElements)*sizeof(dfloat), q);

  occa::memory o_haloq = mesh->device.malloc(mesh->NhaloElements*mesh->Np*sizeof(dfloat));
  
  dfloat *haloqout = (dfloat*) calloc(mesh->NhaloElements*mesh->Np, sizeof(dfloat));
  dfloat *haloqin  = (dfloat*) calloc(mesh->NhaloElements*mesh->Np, sizeof(dfloat));

  MPI_Request *sendRequests = (MPI_Request*) calloc(size, sizeof(MPI_Request));
  MPI_Request *recvRequests = (MPI_Request*) calloc(size, sizeof(MPI_Request));
  
  meshHybridHaloExchangeTri2D(mesh, mesh->Np*sizeof(dfloat), o_q, o_haloq, haloqout, haloqin, sendRequests, recvRequests);
  
  meshVTUTri2D(mesh, "foo.vtu");
  
  MPI_Finalize();

  return 0;
}
