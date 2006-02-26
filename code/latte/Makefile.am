bin_PROGRAMS = maximize minimize count ehrhart \
	test-todd simple-count

maximize_SOURCES = maximize.cpp
minimize_SOURCES = minimize.cpp
count_SOURCES = count.cpp
ehrhart_SOURCES = ehrhart.cpp
simple_count_SOURCES = simple-main.cpp

LDADD = liblatte.a -lgmpxx
AM_CXXFLAGS := $(NTL_CFLAGS) $(LIDIA_CFLAGS) $(AM_CXXFLAGS)
LIBS += $(NTL_LIBS) $(LIDIA_LIBS)

noinst_LIBRARIES = liblatte.a
liblatte_a_SOURCES = \
	barvinok/dec.cpp barvinok/barvinok.cpp  barvinok/Cone.cpp \
	barvinok/ConeDecom.cpp barvinok/Triangulation.cpp \
	genFunction/piped.cpp genFunction/maple.cpp vertices/cdd.cpp \
	preprocess.cpp rational.cpp cone.cpp print.cpp ramon.cpp \
	timing.cpp dual.cpp PolyTree.cpp RudyResNTL.cpp \
	ResSingle.cpp printsingle.cpp IntegralHull.cpp ReadingFile.cpp \
	Grobner.cpp ConeInfo.cpp binarySearchIP.cpp CheckEmpty.cpp \
	banner.cpp \
	todd/gmp_pow.cpp todd/todd-expansion.cpp \
	latte_gmp.cpp		\
	ExponentialSubst.cpp

noinst_HEADERS =							\
	CheckEmpty.h CheckRedundant.h 					\
	ConeInfo.h Grobner.h IntegralHull.h PolyTree.h ReadingFile.h	\
	RudyResNTL.h banner.h binarySearchIP.h cone.h dual.h flags.h	\
	myheader.h preprocess.h print.h ramon.h rational.h timing.h	\
	barvinok/Cone.h barvinok/ConeDecom.h				\
	barvinok/Triangulation.h					\
	barvinok/barvinok.h barvinok/dec.h genFunction/maple.h		\
	genFunction/piped.h vertices/cdd.h				\
	todd/gmp_pow.h todd/todd-expansion.h		\
	latte_gmp.h
	ExponentialSubst.cpp

test_todd_SOURCES = todd/test-todd.cpp


