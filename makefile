ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

include ${OCCA_DIR}/scripts/Makefile

# define variables
HDRDIR  = ./include

# set options for this machine
# specify which compilers to use for c and linking
CC	= mpicxx
LD	= mpicxx

# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I$(HDRDIR)  -g -DOCCA_VERSION_1_0 $(compilerFlags) $(flags)  $(paths)


# link flags to be used 
LDFLAGS	=  $(compilerFlags) -DOCCA_VERSION_1_0

# libraries to be linked in
LIBS	= -L$(OCCA_DIR)/lib $(links) -lm 

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c

# list of objects to be compiled
SOBJS    = \
src/meshHaloExchangeTri2D.o \
src/meshHaloSetupTri2D.o \
src/meshParallelConnectTri2D.o \
src/meshConnectTri2D.o \
src/meshParallelReaderTri2D.o \
src/meshReaderTri2D.o \
src/meshVTUTri2D.o \
src/meshMortonOrderingTri2D.o \
src/meshVolumeGeometricFactorsTri2D.o \
src/meshLoadReferenceNodesTri2D.o \
src/meshGradientTri2D.o \
src/readArray.o

all: meshMain

meshMain:$(SOBJS) src/meshMain.o 
	$(LD)  $(LDFLAGS) -o meshMainsrc/meshMain.o  $(SOBJS) $(LIBS)

meshHaloMain:$(SOBJS) src/meshHaloMainTri2D.o
	$(LD)  $(LDFLAGS) -o meshHaloMain src/meshHaloMainTri2D.o $(SOBJS) $(LIBS)

bandwidthTest:src/bandwidthTest.o
	$(LD) $(LDFLAGS) -o bandwidthTest src/bandwidthTest.o $(LIBS)

# what to do if user types "make clean"
clean :
	rm -r $(SOBJS) meshMain src/bandwidthTest.o bandwidthTest

realclean :
	rm -r $(SOBJS) meshMain src/bandwidthTest.o bandwidthTest

