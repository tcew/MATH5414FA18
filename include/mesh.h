
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define dfloat double

typedef struct {
  
  int Nverts;    /* number of vertices per element */
  int Nfaces;    /* number of faces per element */
  
  int Nvertices; /* number of vertices in mesh */
  
  dfloat *VX; /* array of x-coordinates of vertices */
  dfloat *VY; /* array of y-coordinates of vertices */

  int Nelements; /* number of elements in mesh */
  
  int    *EToV; /* array of element-to-vertex conn 
		 -- Nelements x Nverts */

  int    *EToB; /* element to boundary face conn 
		 -- Nelements x Nfaces */
  
}mesh_t;


mesh_t *meshReaderTri2D(const char * fileName );
void meshConnectTri2D(mesh_t *mesh);
