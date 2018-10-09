
#include "mesh.h"

int main(int argc, char **argv){

  MPI_Init(&argc, &argv);
  
  mesh_t *mesh = meshParallelReaderTri2D(argv[1]);

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

  // load polynomial data for reference triangle
  int N = atoi(argv[2]); // use this for polynomial array ( use 2nd argument to define degree )
  meshLoadReferenceNodesTri2D(mesh, N);

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
  // compute on HOST
  meshGradientTri2D(mesh, q, gradq);

#if 0
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      printf("gradq = %lg\n", gradq[e*2*mesh->Np+n]);
    }
  }
#endif
  
  // initialize DEVICE
  occa::device device("mode: 'CUDA', device_id: 0");

  // allocate array space on DEVICE
  occa::memory o_q = device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), q);
  occa::memory o_gradq = device.malloc(2*mesh->Nelements*mesh->Np*sizeof(dfloat), gradq);

  occa::memory o_vgeo = device.malloc(mesh->Nelements*mesh->Nvgeo*sizeof(dfloat), mesh->vgeo);

  occa::memory o_Dr = device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), mesh->Dr);
  occa::memory o_Ds = device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), mesh->Ds);

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

  int Nkernels = 7;
  char kernelName[BUFSIZ];
    
  for(int k=0;k<Nkernels;++k){

    sprintf(kernelName, "meshGradientTri2D_K%02d", k);
    
    // build kernel on DEVICE
    occa::kernel gradientKernel = device.buildKernel("src/meshGradientTri2D.okl",
						     kernelName, props);
    
    // compute on DEVICE
    gradientKernel(mesh->Nelements, mesh->Np, o_Dr, o_Ds, o_vgeo, o_q, o_gradq);

    // flush all DEVICE kernels
    device.finish(); 

    double tic = MPI_Wtime();
    
    int Ntests = 10;
    for(int test=0;test<Ntests;++test){
      gradientKernel(mesh->Nelements, mesh->Np, o_Dr, o_Ds, o_vgeo, o_q, o_gradq);
    }

    // flush all DEVICE kernels
    device.finish(); 
    
    double toc = MPI_Wtime();

    printf("average elapsed time for kernel %02d is %lg\n",
	   k, (toc-tic)/Ntests);
  }
  
  // copy data back to HOST
  o_q.copyTo(q);
  o_gradq.copyTo(gradq);

  printf("gradq[43] = %lg\n", gradq[43]);
  
  meshVTUTri2D(mesh, "foo.vtu");
  
  MPI_Finalize();

  return 0;
}
