CC = g++
CFLAGS = -I/localapp/imosoft/sparc-sun-solaris2.7/include

Libs = -lm -L/localapp/imosoft/sparc-sun-solaris2.7/lib -lntl

all: latte

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.phony: all clean

clean:
	@-echo "CLEANING..."
	@-rm -f *.o *~ core latte maximize minimize count ehrhart
	@-$(MAKE) -C barvinok clean
	@-$(MAKE) -C genFunction clean
	@-$(MAKE) -C vertices clean

objC = barvinok/dec.o barvinok/barvinok.o  barvinok/Cone.o barvinok/ConeDecom.o barvinok/Triangulation.o genFunction/piped.o genFunction/maple.o vertices/cdd.o preprocess.o rational.o cone.o print.o ramon.o timing.o ConeDeterminant.o dual.o PolyTree.o RudyResNTL.o ResSingle.o printsingle.cpp IntegralHull.o ReadingFile.o Grobner.o ConeInfo.o binarySearchIP.o CheckEmpty.o maximize.o

objD = barvinok/dec.o barvinok/barvinok.o  barvinok/Cone.o barvinok/ConeDecom.o barvinok/Triangulation.o genFunction/piped.o genFunction/maple.o vertices/cdd.o preprocess.o rational.o cone.o print.o ramon.o timing.o ConeDeterminant.o dual.o PolyTree.o RudyResNTL.o ResSingle.o printsingle.cpp IntegralHull.o ReadingFile.o Grobner.o ConeInfo.o binarySearchIP.o CheckEmpty.o  minimize.o

objE = barvinok/dec.o barvinok/barvinok.o  barvinok/Cone.o barvinok/ConeDecom.o barvinok/Triangulation.o genFunction/piped.o genFunction/maple.o vertices/cdd.o preprocess.o rational.o cone.o print.o ramon.o timing.o ConeDeterminant.o dual.o PolyTree.o RudyResNTL.o ResSingle.o printsingle.cpp IntegralHull.o ReadingFile.o Grobner.o ConeInfo.o binarySearchIP.o CheckEmpty.o count.o

objF = barvinok/dec.o barvinok/barvinok.o  barvinok/Cone.o barvinok/ConeDecom.o barvinok/Triangulation.o genFunction/piped.o genFunction/maple.o vertices/cdd.o preprocess.o rational.o cone.o print.o ramon.o timing.o ConeDeterminant.o dual.o PolyTree.o RudyResNTL.o ResSingle.o printsingle.cpp IntegralHull.o ReadingFile.o Grobner.o ConeInfo.o binarySearchIP.o CheckEmpty.o ehrhart.o


latte: maximize minimize count ehrhart redcheck

maximize: $(objC)
	$(CC) $(CFLAGS) $(LDFLAGS) $(objC) -o maximize $(Lib_Dirs) $(Libs)

minimize: $(objD)
	$(CC) $(CFLAGS) $(LDFLAGS) $(objD) -o minimize $(Lib_Dirs) $(Libs)

count: $(objE)
	$(CC) $(CFLAGS) $(LDFLAGS) $(objE) -o count $(Lib_Dirs) $(Libs)

ehrhart: $(objF)
	$(CC) $(CFLAGS) $(LDFLAGS) $(objF) -o ehrhart $(Lib_Dirs) $(Libs)

redcheck:
	$(MAKE) -C adjacency

clobber:
	rm -f *.o latte core

