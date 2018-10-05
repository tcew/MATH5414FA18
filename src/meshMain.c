
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
  int N = 5; // use this for polynomial array
  meshLoadReferenceNodesTri2D(mesh, N);

  // initialize HOST arrays
  dfloat *q     = (dfloat*) calloc(  mesh->Np*(mesh->Nelements+mesh->NhaloElements), sizeof(dfloat));
  dfloat *gradq = (dfloat*) calloc(2*mesh->Np*(mesh->Nelements+mesh->NhaloElements), sizeof(dfloat));

  // compute on HOST
  meshGradientTri2D(mesh, q, gradq);

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

  // build kernel on DEVICE
  occa::kernel gradientKernel = device.buildKernel("src/meshGradientTri2D.okl", "meshGradientTri2D", props);

  // compute on DEVICE
  gradientKernel(mesh->Nelements, mesh->Np, o_Dr, o_Ds, o_vgeo, o_q, o_gradq);

  // copy data back to HOST
  o_q.copyTo(q);
  o_gradq.copyTo(gradq);
  
  meshVTUTri2D(mesh, "foo.vtu");
  
  MPI_Finalize();

  return 0;
}
