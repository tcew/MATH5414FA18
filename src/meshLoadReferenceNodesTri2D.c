#include "mesh.h"

void meshLoadReferenceNodesTri2D(mesh_t *mesh, int N){

  char fname[BUFSIZ];
  sprintf(fname, "nodes/triangleN%02d.dat", N);

  FILE *fp = fopen(fname, "r");

  if (!fp) {
    printf("ERROR: Cannot open file: '%s'\n", fname);
    exit(-1);
  }

  mesh->N = N;
  mesh->Np = (N+1)*(N+2)/2;

  int Nrows, Ncols;

  /* Nodal Data */
  readDfloatArray(fp, "Nodal r-coordinates", &(mesh->r),&Nrows,&Ncols);
  readDfloatArray(fp, "Nodal s-coordinates", &(mesh->s),&Nrows,&Ncols);
  readDfloatArray(fp, "Nodal Dr differentiation matrix", &(mesh->Dr), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal Ds differentiation matrix", &(mesh->Ds), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal Mass Matrix", &(mesh->MM), &Nrows, &Ncols);
  readIntArray   (fp, "Nodal Face nodes", &(mesh->faceNodes), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal Lift Matrix", &(mesh->LIFT), &Nrows, &Ncols);

  fclose(fp);
}
