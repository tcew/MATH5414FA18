
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"
#include "occa.hpp"

#define dfloat double
#define dfloatString "double"
#define dfloatFormat "%lf"
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

  // geometric factors
  int Nvgeo;
  dfloat *vgeo;

  int    N;  // polynomial degree of elements
  int    Np; // number of nodes per element

  // r,s coordinates of interpolation nodes
  dfloat *r, *s;
  
  // differentiation matrices
  dfloat *Dr, *Ds;

  // mass matrix
  dfloat *MM;

  // list of nodes on each face
  int *faceNodes;

  // LIFT matrix
  dfloat *LIFT;

  // OCCA INFO
  occa::device device;
  occa::kernel haloExtractKernel;
  occa::memory o_haloElementIndices; // sorted list of elements that need to be sent to other ranks
  
  
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
			   void *q,
			   int bytesPerElement);

void meshVolumeGeometricFactorsTri2D(mesh_t *mesh);

void meshGradientTri2D(mesh_t *mesh, dfloat *q, dfloat *gradq);

void meshLoadReferenceNodesTri2D(mesh_t *mesh, int N);

void readDfloatArray(FILE *fp, const char *label, dfloat **A, int *Nrows, int* Ncols);

void readIntArray(FILE *fp, const char *label, int **A, int *Nrows, int* Ncols);

void meshHybridHaloExchangeTri2D(mesh_t *mesh,
				 int NbytesPerElement,
				 occa::memory &o_q,
				 occa::memory &o_haloq,
				 dfloat *haloq,
				 dfloat *qin,
				 MPI_Request *sendRequests,
				 MPI_Request *recvRequests);


void meshHybridHaloExchangeStartTri2D(mesh_t *mesh,
				      int NbytesPerElement,
				      occa::memory &o_q,
				      occa::memory &o_haloq,
				      dfloat *haloqout,
				      dfloat *haloqin,
				      MPI_Request *sendRequests,
				      MPI_Request *recvRequests);

void meshHybridHaloExchangeEndTri2D(mesh_t *mesh,
				    int NbytesPerElement,
				    occa::memory &o_q,
				    occa::memory &o_haloq,
				    dfloat *haloqout,
				    dfloat *haloqin,
				    MPI_Request *sendRequests,
				    MPI_Request *recvRequests);


#define p_RXID 0
#define p_RYID 1
#define p_SXID 2
#define p_SYID 3
#define p_JID 4
