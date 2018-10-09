
#include "mesh.h"

void meshVolumeGeometricFactorsTri2D(mesh_t *mesh){

  //
  //
  //  x = -0.5*(r+s)*EX(0,e) + 0.5*(1+r)*EX(1,e) + 0.5*(1+s)*EX(2,e)
  //  y = -0.5*(r+s)*EY(0,e) + 0.5*(1+r)*EY(1,e) + 0.5*(1+s)*EY(2,e)

  // | drdx dsdx | . | dxdr dydr |  = | 1  0 | 
  // | drdy dsdy |   | dxds dyds |    | 0  1 |

  // | drdx dsdx | = 1 |  dyds -dydr |   where J = dxdr*dyds - dxds*dydr
  // | drdy dsdy |   J | -dxds  dxdr | 

  mesh->Nvgeo = 5;
  mesh->vgeo = (dfloat*) calloc(mesh->Nvgeo*mesh->Nelements,
				sizeof(dfloat));

  for(int e=0;e<mesh->Nelements;++e){
    dfloat x0 = mesh->EX[mesh->Nverts*e+0];
    dfloat x1 = mesh->EX[mesh->Nverts*e+1];
    dfloat x2 = mesh->EX[mesh->Nverts*e+2];
    dfloat y0 = mesh->EY[mesh->Nverts*e+0];
    dfloat y1 = mesh->EY[mesh->Nverts*e+1];
    dfloat y2 = mesh->EY[mesh->Nverts*e+2];

    dfloat dxdr = 0.5*(x1-x0);
    dfloat dxds = 0.5*(x2-x0);
    dfloat dydr = 0.5*(y1-y0);
    dfloat dyds = 0.5*(y2-y0);

    dfloat J = dxdr*dyds - dxds*dydr;

    dfloat drdx =  dyds/J;
    dfloat dsdx = -dydr/J;
    dfloat drdy = -dxds/J;
    dfloat dsdy =  dxdr/J;

    //    printf("vgeo = %g,%g,%g,%g,%g \n", drdx, dsdx, drdy, dsdy, J);
    
    mesh->vgeo[e*mesh->Nvgeo + p_RXID] = drdx;
    mesh->vgeo[e*mesh->Nvgeo + p_RYID] = drdy;
    mesh->vgeo[e*mesh->Nvgeo + p_SXID] = dsdx;
    mesh->vgeo[e*mesh->Nvgeo + p_SYID] = dsdy;
    mesh->vgeo[e*mesh->Nvgeo + p_JID]  = J;
  }
}
