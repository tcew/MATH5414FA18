
#include "mesh.h"

void meshGradientTri2D(mesh_t *mesh, dfloat *q, dfloat *gradq){

  // loop over elements
  for(int e=0;e<mesh->Nelements;++e){

    // loop over nodes in element e
    for(int n=0;n<mesh->Np;++n){
      dfloat dqdr = 0;
      dfloat dqds = 0;

      // compute derivative at  node n of element e
      for(int m=0;m<mesh->Np;++m){
	dqdr += mesh->Dr[n*mesh->Np+m]*q[m+e*mesh->Np];
	dqds += mesh->Ds[n*mesh->Np+m]*q[m+e*mesh->Np];
      }

      dfloat drdx = mesh->vgeo[e*mesh->Nvgeo + p_RXID];
      dfloat dsdx = mesh->vgeo[e*mesh->Nvgeo + p_SXID];
      dfloat drdy = mesh->vgeo[e*mesh->Nvgeo + p_RYID];
      dfloat dsdy = mesh->vgeo[e*mesh->Nvgeo + p_SYID];

      dfloat dqdx = drdx*dqdr + dsdx*dqds;
      dfloat dqdy = drdy*dqdr + dsdy*dqds;

      // TW: explain this more later
      gradq[e*mesh->Np*2 + 0*mesh->Np + n] = dqdx;
      gradq[e*mesh->Np*2 + 1*mesh->Np + n] = dqdy;
      
    }
  }

}
