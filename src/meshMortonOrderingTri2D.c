
#include "mesh.h"

#define ulli unsigned long long int

ulli mortonIndex(unsigned int x, unsigned int y){

  unsigned int xymask = 1;
  ulli mask = 1;

  ulli index = 0;

  for(int b=0;b<32;++b){

    index += mask * ( (x & xymask)!=0) ;
    mask <<= 1;

    //    printf("x: b =%d, index = %llu, mask = %llu\n",
    //	   b, index, mask);
    
    index += mask * ( (y & xymask)!=0) ;
    mask <<= 1;

    //    printf("y: b =%d, index = %llu, mask = %llu\n",
    //	   b, index, mask);
    
    xymask <<= 1;
  }

  return index;
  
}

typedef struct {
  int element;
  ulli mIndex;
  int v1, v2, v3;
}mortonNode_t;

int compareMortonNodes(const void *n1,
		       const void *n2){

  mortonNode_t *mn1 = (mortonNode_t*) n1;
  mortonNode_t *mn2 = (mortonNode_t*) n2;

  if(mn1->mIndex > mn2->mIndex) return 1;
  else if(mn1->mIndex < mn2->mIndex) return -1;

  return 0;

}

void meshMortonOrderingTri2D(mesh_t *mesh){

  int B = 1<<30;

  dfloat xmin = 1e9, xmax = -1e9;
  dfloat ymin = 1e9, ymax = -1e9;

  mortonNode_t *mortonNodes =
    (mortonNode_t*) calloc(mesh->Nelements, sizeof(mortonNode_t));
  
  for(int e=0;e<mesh->Nelements;++e){
    for(int v=0;v<mesh->Nverts;++v){
      xmin = mymin(xmin, mesh->VX[mesh->EToV[e*mesh->Nverts+v]]);
      ymin = mymin(ymin, mesh->VY[mesh->EToV[e*mesh->Nverts+v]]);
      xmax = mymax(xmax, mesh->VX[mesh->EToV[e*mesh->Nverts+v]]);
      ymax = mymax(ymax, mesh->VY[mesh->EToV[e*mesh->Nverts+v]]);

    }
  }
  
  for(int e=0;e<mesh->Nelements;++e){

    unsigned int ix, iy;

    dfloat cx = 0;
    dfloat cy = 0;

    for(int v=0;v<mesh->Nverts;++v){
      cx += mesh->VX[mesh->EToV[e*mesh->Nverts+v]];
      cy += mesh->VY[mesh->EToV[e*mesh->Nverts+v]];
    }
    
    cx /= mesh->Nverts;
    cy /= mesh->Nverts;

    ix = B*(cx-xmin)/(xmax-xmin);
    iy = B*(cy-ymin)/(ymax-ymin);
    
    mortonNodes[e].element = e;
    mortonNodes[e].mIndex  = mortonIndex(ix, iy);
    mortonNodes[e].v1 = mesh->EToV[e*mesh->Nverts+0];
    mortonNodes[e].v2 = mesh->EToV[e*mesh->Nverts+1];
    mortonNodes[e].v3 = mesh->EToV[e*mesh->Nverts+2];
  }

  qsort(mortonNodes, mesh->Nelements,
	sizeof(mortonNode_t), compareMortonNodes);

  for(int e=0;e<mesh->Nelements;++e){
    mesh->EToV[e*mesh->Nverts+0] =
      mortonNodes[e].v1;
    mesh->EToV[e*mesh->Nverts+1] =
      mortonNodes[e].v2;
    mesh->EToV[e*mesh->Nverts+2] =
      mortonNodes[e].v3;
  }

  free(mortonNodes);
}


#if 0
int main(int argc, char **argv){

  unsigned int x = 7, y = 3;

  ulli index = mortonIndex(x, y);
  printf("x=%d, y=%d => index=%llu\n", x, y, index);

  return 0;

}
#endif
