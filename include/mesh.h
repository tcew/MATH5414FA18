
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#define dfloat double
#define MPI_DFLOAT MPI_DOUBLE

typedef struct {
  
  int Nverts;    /* number of vertices per element */
  int Nfaces;    /* number of faces per element */
  
  int Nvertices; /* number of vertices in mesh */
  
  dfloat *VX; /* array of x-coordinates of vertices */
  dfloat *VY; /* array of y-coordinates of vertices */

  dfloat *EX; /* array of x-coords of element vertices */
  dfloat *EY; /* array of y-coords of element vertices */

  int Nelements; /* number of elements in mesh */
  
  int    *EToV; /* array of element-to-vertex conn 
		 -- Nelements x Nverts */

  int    *EToB; /* element to boundary face conn 
		 -- Nelements x Nfaces */

  int    *EToE;
  int    *EToF;
  int    *EToP; // element to rank

  // halo stuff
  int NhaloElements; // total number of elements to send (also to recv)
  int *haloElementIndices; // sorted list of elements that need to be sent to other ranks
  int *NhaloExchangeElements; // number of elements to exchang with each other rank
}mesh_t;


mesh_t *meshReaderTri2D(const char * fileName );
void meshConnectTri2D(mesh_t *mesh);

mesh_t *meshParallelReaderTri2D(const char * fileName);

#define mymax(a,b) ( ((a)>(b)) ? (a):(b) )
#define mymin(a,b) ( ((a)<(b)) ? (a):(b) )

void meshVTUTri2D(mesh_t *msh, char *fileName);

void meshMortonOrderingTri2D(mesh_t *mesh);

void meshParallelConnectTri2D(mesh_t *mesh);

void meshHaloSetupTri2D(mesh_t *mesh);

void meshHaloExchangeTri2D(mesh_t *mesh,
			   unsigned char *q,
			   int bytesPerElement);
