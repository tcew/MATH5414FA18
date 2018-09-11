
#include "mesh.h"

int main(int argc, char **argv){

  MPI_Init(&argc, &argv);
  
  mesh_t *mesh = meshParallelReaderTri2D(argv[1]);

  meshParallelConnectTri2D(mesh);

  meshVTUTri2D(mesh, "foo.vtu");
  
  MPI_Finalize();

  return 0;
}
