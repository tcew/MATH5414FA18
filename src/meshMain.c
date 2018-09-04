
#include "mesh.h"

int main(int argc, char **argv){

  MPI_Init(&argc, &argv);
  
  mesh_t *mesh = meshParallelReaderTri2D(argv[1]);

  meshConnectTri2D(mesh);
  
  MPI_Finalize();

  return 0;
}
