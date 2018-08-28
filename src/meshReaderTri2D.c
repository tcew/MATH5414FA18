#include "mesh.h"

mesh_t *meshReaderTri2D(const char * fileName ){
  
  mesh_t *mesh = (mesh_t*) calloc(1, sizeof(mesh_t));

  mesh->Nverts = 3;
  mesh->Nfaces = 3;
  
  char buf[BUFSIZ];
  
  FILE *fp = fopen(fileName, "r");
  if(!fp){
    printf("Mesh file: %s\n", fileName);
    exit(-1);
  }

  do{
    fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Nodes"));
  
  fscanf(fp, "%d", &(mesh->Nvertices));

  mesh->VX = (dfloat*) calloc(mesh->Nvertices, 
			      sizeof(dfloat));
  mesh->VY = (dfloat*) calloc(mesh->Nvertices, 
			      sizeof(dfloat));
  for(int v=0;v<mesh->Nvertices;++v){
    fscanf(fp, "%*d %lf %lf %*lf", 
	   mesh->VX+v,
	   mesh->VY+v);
    printf("%d %g %g \n", v, mesh->VX[v], mesh->VY[v]);
  }

  do{
    fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Elements"));

  fscanf(fp, "%d", &(mesh->Nelements));

  mesh->EToV = 
    (int*) calloc(mesh->Nelements*mesh->Nverts, 
		     sizeof(int));
  
  for(int e=0;e<mesh->Nelements;++e){
    int etype, Ntags;
    fscanf(fp, "%*d %d %d", 
	   &etype,
	   &Ntags);

    for(int t=0;t<Ntags;++t)
      fscanf(fp, "%*d");

    printf("EToV[e,:]= ");
    for(int v=0;v<mesh->Nverts;++v){
      fscanf(fp, "%d", mesh->EToV + e*mesh->Nverts+v);

      // change from 1-index to 0-index
      --(mesh->EToV[e*mesh->Nverts+v]);
      printf("%d ", mesh->EToV[e*mesh->Nverts+v]);
    }
    printf("\n");
  }

  fclose(fp);
  return mesh;
}
