prefix = /h/ruriko/CDD/

ComputeAdjacency: ComputeAdjacency.o
	g++ -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE=\"cddlib\" -DVERSION=\"0.93\" -DHAVE_LIBGMP=1 -DSTDC_HEADERS=1  -I. -I. -I${prefix}/cddlib-093a/lib-src -I${prefix}/cddlib-093a/include -UGMPRATIONAL  -g -O2 -c `test -f 'ComputeAdjacency.c' || echo './'`ComputeAdjacency.c
	g++  -g -O2   -o ComputeAdjacency ComputeAdjacency.o -L${prefix}/cddlib-093a/lib-src -lcdd -lgmp

redcheck: redcheck.o
	g++ -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE=\"cddlib\" -DVERSION=\"0.93\" -DHAVE_LIBGMP=1 -DSTDC_HEADERS=1  -I. -I. -I${prefix}/cddlib-093a/lib-src-gmp -I${prefix}/cddlib-093a/include -UGMPRATIONAL  -g -O2 -c `test -f 'redcheck.c' || echo './'`redcheck.c
	g++  -g -O2   -o redcheck redcheck.o -L${prefix}/cddlib-093a/lib-src-gmp -lcddgmp -lgmp

clean:
	rm -f *.o core ComputeAdjacency
