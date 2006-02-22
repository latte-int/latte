bin_PROGRAMS = maximize minimize count ehrhart \
	test-todd

maximize_SOURCES = maximize.cpp
minimize_SOURCES = minimize.cpp
count_SOURCES = count.cpp
ehrhart_SOURCES = ehrhart.cpp

LDADD = liblatte.a
AM_CXXFLAGS := $(NTL_CFLAGS) $(AM_CXXFLAGS)
LIBS += $(NTL_LIBS)

noinst_LIBRARIES = liblatte.a
liblatte_a_SOURCES = \
	barvinok/dec.cpp barvinok/barvinok.cpp  barvinok/Cone.cpp \
	barvinok/ConeDecom.cpp barvinok/Triangulation.cpp \
	genFunction/piped.cpp genFunction/maple.cpp vertices/cdd.cpp \
	preprocess.cpp rational.cpp cone.cpp print.cpp ramon.cpp \
	timing.cpp ConeDeterminant.cpp dual.cpp PolyTree.cpp RudyResNTL.cpp \
	ResSingle.cpp printsingle.cpp IntegralHull.cpp ReadingFile.cpp \
	Grobner.cpp ConeInfo.cpp binarySearchIP.cpp CheckEmpty.cpp \
	banner.cpp \
	todd/todd.cpp todd/gmp_pow.cpp

BUILT_SOURCES = todd/todd.cpp
EXTRA_DIST = todd/todd.cpp
if HAVE_MAPLE
todd/todd.pre-cpp: todd/todd.maple
	(cd todd && $(MAPLEBIN) todd.maple )
.PRECIOUS: todd/todd.pre-cpp
todd/todd.cpp: todd/todd.pre-cpp todd/todd.sed
	sed -f todd/todd.sed $< > $@
endif

test_todd_SOURCES = todd/test-todd.cpp
test_todd_LDADD = $(LDADD) -lgmpxx

