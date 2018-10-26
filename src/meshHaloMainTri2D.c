
#include "mesh.h"

int main(int argc, char **argv){

  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  mesh_t *mesh = meshParallelReaderTri2D(argv[1]);

  // initialize DEVICE
  char setupString[BUFSIZ];
  sprintf(setupString, "mode: 'CUDA', device_id: %d", rank%2);
  mesh->device.setup(setupString);
  
  //  mesh->device.setup("mode: 'Serial'");

  // load polynomial data for reference triangle
  int N = atoi(argv[2]); // use this for polynomial array ( use 2nd argument to define degree )
  meshLoadReferenceNodesTri2D(mesh, N);
  
  meshParallelConnectTri2D(mesh);

  // build kernel on DEVICE
  for(int r=0;r<size;++r){
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==r){
      mesh->haloExtractKernel =
	mesh->device.buildKernel("src/meshHaloExtract.okl",
				 "meshHaloExtract");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  meshHaloSetupTri2D(mesh);

  int *qrank = (int*) calloc(mesh->Nelements + mesh->NhaloElements, sizeof(int));

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

  dfloat *gradq = (dfloat*) calloc(2*mesh->Np*(mesh->Nelements+mesh->NhaloElements), sizeof(dfloat));

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
  printf("mesh->NhaloElements=%d\n", mesh->NhaloElements);
  occa::memory o_q = mesh->device.malloc(mesh->Np*(mesh->Nelements+mesh->NhaloElements)*sizeof(dfloat), q);

  occa::memory o_gradq =
    mesh->device.malloc(2*mesh->Nelements*mesh->Np*sizeof(dfloat), gradq);

  occa::memory o_vgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nvgeo*sizeof(dfloat), mesh->vgeo);

  occa::memory o_Dr =
    mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), mesh->Dr);

  occa::memory o_Ds =
    mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), mesh->Ds);

  // properties to be used for defining compiler variables
  occa::properties props;
  props["defines/p_RXID"]  = p_RXID;
  props["defines/p_RYID"]  = p_RYID;
  props["defines/p_SXID"]  = p_SXID;
  props["defines/p_SYID"]  = p_SYID;
  props["defines/p_JID"]   = p_JID;
  props["defines/p_Nvgeo"] = mesh->Nvgeo;
  props["defines/dfloat"]  = dfloatString;
  props["defines/p_Np"]    = mesh->Np;



  occa::memory o_haloq = mesh->device.malloc(mesh->NhaloElements*mesh->Np*sizeof(dfloat));
  
  dfloat *haloqout = (dfloat*) calloc(mesh->NhaloElements*mesh->Np, sizeof(dfloat));
  dfloat *haloqin  = (dfloat*) calloc(mesh->NhaloElements*mesh->Np, sizeof(dfloat));

  MPI_Request *sendRequests = (MPI_Request*) calloc(size, sizeof(MPI_Request));
  MPI_Request *recvRequests = (MPI_Request*) calloc(size, sizeof(MPI_Request));

  occa::kernel gradientKernel =
    mesh->device.buildKernel("src/meshGradientTri2D.okl",
			     "meshGradientTri2D_K08", props);
  
  occa::stream defaultStream = mesh->device.getStream();
  occa::stream dataStream    = mesh->device.createStream();
  occa::stream computeStream = mesh->device.createStream();

  
  mesh->device.setStream(dataStream);
    
  meshHybridHaloExchangeStartTri2D(mesh, mesh->Np*sizeof(dfloat), o_q, o_haloq, haloqout, haloqin, sendRequests, recvRequests);
  
  mesh->device.setStream(computeStream);

  gradientKernel(mesh->Nelements, mesh->Np, o_Dr, o_Ds, o_vgeo, o_q, o_gradq);

  mesh->device.setStream(dataStream);

  meshHybridHaloExchangeEndTri2D(mesh, mesh->Np*sizeof(dfloat), o_q, o_haloq, haloqout, haloqin, sendRequests, recvRequests);
  
  
  mesh->device.finish();
  //  meshVTUTri2D(mesh, "foo.vtu");

  o_gradq.copyTo(gradq);
  
  for(int n=0;n<1000;++n){
    printf("gradq[%d] = %g\n", n, gradq[n]);
  }
  

  
  MPI_Finalize();

  return 0;
}
