#include "mesh.h"

typedef struct {
  int element;
  int face;
  int v1;
  int v2;
}face_t;

int compareFaces(const void *a, const void*b){
  
  face_t *faceA = (face_t*) a;
  face_t *faceB = (face_t*) b;

  // assume that v1 is the 2nd digit
  if(faceA->v1 < faceB->v1)
    return -1;
  if(faceA->v1 > faceB->v1)
    return +1;

  if(faceA->v2 < faceB->v2)
    return -1;
  if(faceA->v2 > faceB->v2)
    return +1;
  
  return 0;

}

void meshConnectTri2D(mesh_t *mesh){

  int sortNfaces = mesh->Nelements*mesh->Nfaces;
  
  face_t *faces = (face_t*) 
    calloc(sortNfaces, sizeof(face_t));

  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      int id = e*mesh->Nfaces + f;
      faces[id].element = e;
      faces[id].face = f;
      faces[id].v1 = 
	mesh->EToV[e*mesh->Nverts + f];
      faces[id].v2 = 
	mesh->EToV[e*mesh->Nverts + (f+1)%mesh->Nverts];
      if(faces[id].v1 < faces[id].v2){
	int tmp = faces[id].v1;
	faces[id].v1 = faces[id].v2;
	faces[id].v2 = tmp;
      }
    }
  }

  // sort by double face vertices 
  qsort(faces,
	sortNfaces,
	sizeof(face_t),
	compareFaces);
  
  for(int f=0;f<mesh->Nfaces*mesh->Nelements;++f){
    printf("f: %d => %d,%d (%d,%d)\n", 
	   f, faces[f].v1, faces[f].v2,
	   faces[f].element, faces[f].face);
  }


}
