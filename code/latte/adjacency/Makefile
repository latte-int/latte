CC = cc -fast
CDD_INCLUDE_DIR=/localapp/imosoft/sparc-sun-solaris2.7/include
CDD_LIB_DIR=/localapp/imosoft/sparc-sun-solaris2.7/lib

CFLAGS = -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE=\"cddlib\" -DVERSION=\"0.93\" -DHAVE_LIBGMP=1 -DSTDC_HEADERS=1  -I.  -I$(CDD_INCLUDE_DIR) -UGMPRATIONAL
LIBS = -L$(CDD_LIB_DIR) -lcdd -lgmp

all: ComputeAdjacency redcheck

ComputeAdjacency: ComputeAdjacency.o
	$(CC) $(CFLAGS)  -c ComputeAdjacency.c
	$(CC) $(CFLAGS)  -o ComputeAdjacency ComputeAdjacency.o $(LIBS)

redcheck: redcheck.o
	$(CC) $(CFLAGS) -c redcheck.c
	$(CC) $(CFLAGS) -o redcheck redcheck.o $(LIBS)

clean:
	rm -f *.o core redcheck ComputeAdjacency
