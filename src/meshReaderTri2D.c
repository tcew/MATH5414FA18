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
  }

  do{
    fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Elements"));

  int Nelements;
  fscanf(fp, "%d", &Nelements);

  // set marker in file
  fpos_t pos;
  fgetpos(fp, &pos);

  // count number of triangles
  int Ntriangles = 0;
  for(int e=0;e<Nelements;++e){
    int etype;
    fscanf(fp, "%*d %d", 
	   &etype);
    if(etype==2) ++Ntriangles;
    fgets(buf, BUFSIZ, fp);
  }

  // rewind to element section
  fsetpos(fp, &pos);
  
  mesh->Nelements = Ntriangles;
  
  mesh->EToV = 
    (int*) calloc(mesh->Nelements*mesh->Nverts, sizeof(int));

  int triangle = 0;
  for(int e=0;e<Nelements;++e){
    int etype, Ntags;
    fscanf(fp, "%*d %d %d", &etype, &Ntags);
    
    if(etype==2){
      for(int t=0;t<Ntags;++t)
	fscanf(fp, "%*d");
      
      for(int v=0;v<mesh->Nverts;++v){
	fscanf(fp, "%d", mesh->EToV + triangle*mesh->Nverts+v);
	
	// change from 1-index to 0-index
	--(mesh->EToV[triangle*mesh->Nverts+v]);

      }

      ++triangle;
    }
    else{
      fgets(buf, BUFSIZ, fp);      
    }
  }

  fclose(fp);
  return mesh;
}
