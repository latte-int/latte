CC = gcc
CXX = g++
Libs = -lm -L/localapp/imosoft/sparc-sun-solaris2.7/alpha/ntl-5.3.1-gcc/lib -lntl -L/localapp/imosoft/sparc-sun-solaris2.7/lib -lgmp
#CC = cc
#CXX = CC
#Libs = -lm -L/localapp/imosoft/sparc-sun-solaris2.7/lib -lntl
CFLAGS = -I/localapp/imosoft/sparc-sun-solaris2.7/alpha/ntl-5.3.1-gcc/include -I/localapp/imosoft/sparc-sun-solaris2.7/include 


all: latte

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cc
	$(CXX) $(CFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

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
	$(CXX) $(CFLAGS) $(LDFLAGS) $(objC) -o maximize $(Lib_Dirs) $(Libs)

minimize: $(objD)
	$(CXX) $(CFLAGS) $(LDFLAGS) $(objD) -o minimize $(Lib_Dirs) $(Libs)

count: $(objE)
	$(CXX) $(CFLAGS) $(LDFLAGS) $(objE) -o count $(Lib_Dirs) $(Libs)

ehrhart: $(objF)
	$(CXX) $(CFLAGS) $(LDFLAGS) $(objF) -o ehrhart $(Lib_Dirs) $(Libs)

redcheck:
	$(MAKE) -C adjacency

clobber:
	rm -f *.o latte core

