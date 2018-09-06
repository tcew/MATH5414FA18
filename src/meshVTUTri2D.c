#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "mesh.h"

void meshVTUTri2D(mesh_t *msh, char *fileName){
  
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  FILE *fp;
  int *allNelements = (int*) calloc(size, sizeof(int));
  int totalNelements = 0, maxNelements = 0;
  dfloat *tmpEX, *tmpEY;

  if(rank==0){
    fp = fopen(fileName, "w");
    printf("fp=%p\n for fileName=%s\n",fp, fileName);
  }

  // gather element counts to root
  MPI_Allgather(&(msh->Nelements), 1, MPI_INT, 
		allNelements, 1, MPI_INT, 
		MPI_COMM_WORLD);
  
  if(rank==0){
    
    for(int r=0;r<size;++r){
      totalNelements+=allNelements[r];
      maxNelements = mymax(maxNelements, allNelements[r]);
    }
    
    tmpEX = (dfloat*) calloc(maxNelements*msh->Nverts, sizeof(dfloat));
    tmpEY = (dfloat*) calloc(maxNelements*msh->Nverts, sizeof(dfloat));
    
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
    fprintf(fp, "  <UnstructuredGrid>\n");
    fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", totalNelements*msh->Nverts, totalNelements);

    // write out nodes
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");

    if(rank==0){  // root writes out its coordinates
      printf("printing local verts \n");
      for(int e=0;e<msh->Nelements;++e){
	fprintf(fp, "        ");
	for(int n=0;n<msh->Nverts;++n)
	  fprintf(fp, "%g %g 0.\n", msh->EX[e*msh->Nverts+n], msh->EY[e*msh->Nverts+n]);
      }
    }
  }
  
  for(int r=1;r<size;++r){
    
    if(rank==r){
      MPI_Send(msh->EX, msh->Nelements*msh->Nverts, 
	       MPI_DFLOAT, 0, 666, MPI_COMM_WORLD);
      MPI_Send(msh->EY, msh->Nelements*msh->Nverts, 
	       MPI_DFLOAT, 0, 666, MPI_COMM_WORLD);
    }
    
    if(rank==0){
      MPI_Status status;
      MPI_Recv(tmpEX, allNelements[r]*msh->Nverts, 
	       MPI_DFLOAT, r, 666, MPI_COMM_WORLD, &status);
      MPI_Recv(tmpEY, allNelements[r]*msh->Nverts, 
	       MPI_DFLOAT, r, 666, MPI_COMM_WORLD, &status);
      
      for(int e=0;e<allNelements[r];++e){
	fprintf(fp, "        ");
	for(int n=0;n<msh->Nverts;++n)
	  fprintf(fp, "%g %g 0\n", tmpEX[e*msh->Nverts+n], tmpEY[e*msh->Nverts+n]);
      }
    }
  }    
  if(rank==0){
    free(tmpEX);
    free(tmpEY);
    
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");
    
    // write out rank
    fprintf(fp, "      <CellData Scalars=\"scalars\">\n");
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"element_rank\" Format=\"ascii\">\n");
    
    for(int r=0;r<size;++r){
      for(int e=0;e<allNelements[r];++e)
	fprintf(fp, "         %d\n", r);
    }
    
    fprintf(fp, "       </DataArray>\n");
    fprintf(fp, "     </CellData>\n");
    
    fprintf(fp, "    <Cells>\n");
    fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");

    int cnt=0;
    for(int r=0;r<size;++r){
      for(int e=0;e<allNelements[r];++e){
	for(int n=0;n<msh->Nverts;++n){
	  fprintf(fp, "%d ", cnt);
	  ++cnt;
	}
	fprintf(fp, "\n");
      }
    }

    fprintf(fp, "        </DataArray>\n");
    
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
    for(int e=0;e<totalNelements;++e){
      if(e%10==0) fprintf(fp, "        ");
      fprintf(fp, "%d ", (e+1)*msh->Nverts);
      if(((e+1)%10==0) || (e==totalNelements-1))
	fprintf(fp, "\n");
    }
    fprintf(fp, "       </DataArray>\n");
    
    fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
    for(int e=0;e<totalNelements;++e){
      if(e%10==0) fprintf(fp, "        ");
      fprintf(fp, "5 ");
      if(((e+1)%10==0) || e==(totalNelements-1))
	fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Cells>\n");
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
  }

  MPI_Barrier(MPI_COMM_WORLD);
}
