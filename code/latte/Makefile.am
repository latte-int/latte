bin_PROGRAMS = maximize minimize count ehrhart

maximize_SOURCES = maximize.cpp
minimize_SOURCES = minimize.cpp
count_SOURCES = count.cpp
ehrhart_SOURCES = ehrhart.cpp

LDADD = liblatte.a
AM_CXXFLAGS := $(NTL_CFLAGS) $(LIDIA_CFLAGS) $(GMP_CFLAGS) $(AM_CXXFLAGS)
LIBS += $(NTL_LIBS) $(LIDIA_LIBS) $(GMP_LIBS)

noinst_LIBRARIES = liblatte.a
liblatte_a_SOURCES = \
	barvinok/dec.cpp barvinok/barvinok.cpp  barvinok/Cone.cpp \
	barvinok/ConeDecom.cpp barvinok/Triangulation.cpp \
	genFunction/piped.cpp genFunction/maple.cpp vertices/cdd.cpp \
	preprocess.cpp rational.cpp cone.cpp print.cpp ramon.cpp \
	timing.cpp dual.cpp PolyTree.cpp RudyResNTL.cpp \
	ResSingle.cpp Residue.cpp IntegralHull.cpp ReadingFile.cpp \
	Grobner.cpp ConeInfo.cpp binarySearchIP.cpp CheckEmpty.cpp \
	banner.cpp convert.cpp \
	todd/gmp_pow.cpp todd/todd-expansion.cpp \
	latte_gmp.cpp latte_random.cpp latte_system.cpp		\
	genFunction/IntCombEnum.cpp		\
	genFunction/matrix_ops.cpp		\
	ExponentialSubst.cpp 

noinst_HEADERS =							\
	CheckEmpty.h 		 					\
	ConeInfo.h Grobner.h IntegralHull.h PolyTree.h ReadingFile.h	\
	RudyResNTL.h banner.h binarySearchIP.h cone.h dual.h flags.h	\
	myheader.h preprocess.h print.h ramon.h rational.h timing.h	\
	Residue.h barvinok/Cone.h barvinok/ConeDecom.h			\
	barvinok/Triangulation.h					\
	barvinok/barvinok.h barvinok/dec.h genFunction/maple.h		\
	genFunction/piped.h vertices/cdd.h				\
	todd/gmp_pow.h todd/todd-expansion.h				\
	convert.h							\
	latte_gmp.h latte_ntl.h	latte_random.h latte_system.h		\
	genFunction/IntCombEnum.h     \
	genFunction/matrix_ops.h    \
	ExponentialSubst.h 

if HAVE_LIDIA
liblatte_a_SOURCES += 	genFunction/NTL_to_LiDIA.cpp
noinst_HEADERS += genFunction/NTL_to_LiDIA.h
endif

if HAVE_EXPERIMENTS
# Experimental stuff
bin_PROGRAMS += test-todd 

nodist_liblatte_a_SOURCES = \
	ExponentialApprox.cpp Irrational.cpp
nodist_noinst_HEADERS = \
	ExponentialApprox.h Irrational.h 

#nodist_simple_count_SOURCES = simple-main.cpp
nodist_test_todd_SOURCES = todd/test-todd.cpp

endif
