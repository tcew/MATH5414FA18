
#include "mesh.h"

int main(int argc, char **argv){

  mesh_t *mesh = meshReaderTri2D(argv[1]);

  meshConnectTri2D(mesh);
  
  return 0;
}
