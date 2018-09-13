# define variables
HDRDIR  = ./include

# set options for this machine
# specify which compilers to use for c and linking
CC	= mpicc
LD	= mpicc


# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I$(HDRDIR)  -g 


# link flags to be used 
LDFLAGS	=  -O3 

# libraries to be linked in
LIBS	=  -lm 

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c

# list of objects to be compiled
SOBJS    = \
src/meshParallelConnectTri2D.o \
src/meshConnectTri2D.o \
src/meshMain.o \
src/meshParallelReaderTri2D.o \
src/meshReaderTri2D.o \
src/meshVTUTri2D.o \
src/meshMortonOrderingTri2D.o 

all: meshMain

meshMain:$(SOBJS) 
	$(LD)  $(LDFLAGS) -o meshMain $(SOBJS) $(LIBS)

# what to do if user types "make clean"
clean :
	rm -r $(SOBJS) meshMain

realclean :
	rm -r $(SOBJS) meshMain

