AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = maximize minimize count ehrhart latte2ine latte2ext

maximize_SOURCES = maximize.cpp
minimize_SOURCES = minimize.cpp
count_SOURCES = count.cpp
ehrhart_SOURCES = ehrhart.cpp
latte2ine_SOURCES = latte2ine.cpp
latte2ext_SOURCES = latte2ext.cpp

LDADD = liblatte.a
AM_CXXFLAGS = $(NTL_CFLAGS) $(LIDIA_CFLAGS) $(TOPCOM_CXXFLAGS) $(CDDLIB_CFLAGS) $(FORTYTWO_CXXFLAGS) $(GMP_CFLAGS)
LIBS += $(NTL_LIBS) $(LIDIA_LIBS) $(TOPCOM_LIBS) $(CDDLIB_LIBS) $(FORTYTWO_LIBS) $(GMP_LIBS)

noinst_LIBRARIES = liblatte.a
liblatte_a_SOURCES =							\
	barvinok/dec.cpp barvinok/barvinok.cpp  barvinok/Cone.cpp	\
	barvinok/Triangulation.cpp					\
	genFunction/piped.cpp genFunction/maple.cpp vertices/cdd.cpp	\
	preprocess.cpp rational.cpp cone.cpp print.cpp ramon.cpp	\
	timing.cpp dual.cpp PolyTree.cpp RudyResNTL.cpp			\
	ResSingle.cpp Residue.cpp IntegralHull.cpp ReadingFile.cpp	\
	Grobner.cpp ConeInfo.cpp binarySearchIP.cpp CheckEmpty.cpp	\
	banner.cpp convert.cpp						\
	todd/gmp_pow.cpp todd/todd-expansion.cpp			\
	latte_gmp.cpp latte_ntl.cpp latte_random.cpp latte_system.cpp	\
	genFunction/IntCombEnum.cpp					\
	genFunction/matrix_ops.cpp					\
	ExponentialSubst.cpp ExponentialEhrhart.cpp			\
	Irrational.cpp							\
	ReadLatteStyle.cpp						\
	triangulation/triangulate.cpp					\
	triangulation/RegularTriangulation.cpp				\
	triangulation/RegularTriangulationWithCdd.cpp

noinst_HEADERS =							\
	CheckEmpty.h							\
	ConeInfo.h Grobner.h IntegralHull.h PolyTree.h ReadingFile.h	\
	RudyResNTL.h banner.h binarySearchIP.h cone.h dual.h flags.h	\
	preprocess.h print.h ramon.h rational.h timing.h		\
	Residue.h barvinok/Cone.h					\
	barvinok/Triangulation.h					\
	barvinok/barvinok.h barvinok/dec.h genFunction/maple.h		\
	genFunction/piped.h vertices/cdd.h				\
	todd/gmp_pow.h todd/todd-expansion.h				\
	convert.h							\
	latte_gmp.h latte_ntl.h	latte_ntl_integer.h			\
	latte_random.h latte_system.h latte_cddlib.h			\
	genFunction/IntCombEnum.h					\
	genFunction/matrix_ops.h					\
	ExponentialSubst.h ExponentialEhrhart.h				\
	Irrational.h							\
	ReadLatteStyle.h Polyhedron.h					\
	triangulation/triangulate.h					\
	triangulation/RegularTriangulation.h				\
	triangulation/RegularTriangulationWithCdd.h

if HAVE_TOPCOM_LIB
liblatte_a_SOURCES += triangulation/TriangulationWithTOPCOM.cpp
noinst_HEADERS += triangulation/TriangulationWithTOPCOM.h
else
if HAVE_TOPCOM_BIN
liblatte_a_SOURCES += triangulation/TriangulationWithTOPCOM.cpp
noinst_HEADERS += triangulation/TriangulationWithTOPCOM.h
endif
endif

if HAVE_CDDLIB
liblatte_a_SOURCES += latte_cddlib.cpp
liblatte_a_SOURCES += triangulation/RegularTriangulationWithCddlib.cpp
noinst_HEADERS += triangulation/RegularTriangulationWithCddlib.h
endif

if HAVE_LIDIA
liblatte_a_SOURCES += genFunction/NTL_to_LiDIA.cpp
noinst_HEADERS += genFunction/NTL_to_LiDIA.h
endif

if HAVE_FORTYTWO_LIB
liblatte_a_SOURCES += latte_4ti2.cpp
liblatte_a_SOURCES += triangulation/RegularTriangulationWith4ti2.cpp
noinst_HEADERS += latte_4ti2.h
noinst_HEADERS += triangulation/RegularTriangulationWith4ti2.h

liblatte_a_SOURCES += DualizationWith4ti2.cpp
noinst_HEADERS += DualizationWith4ti2.h

liblatte_a_SOURCES += latte_4ti2_zsolve.cpp
noinst_HEADERS += latte_4ti2_zsolve.h
endif

if HAVE_EXPERIMENTS
# Experimental stuff

nodist_liblatte_a_SOURCES =				\
	ExponentialApprox.cpp				\
	TrivialSubst.cpp				\
	barvinok/SubspaceAvoidingDecomposition.cpp	\
	triangulation/RecursiveTriangulation.cpp	\
	triangulation/BoundaryTriangulation.cpp

nodist_noinst_HEADERS =					\
	ExponentialApprox.h				\
	TrivialSubst.h					\
	barvinok/SubspaceAvoidingDecomposition.h	\
	triangulation/RecursiveTriangulation.h		\
	triangulation/BoundaryTriangulation.h

bin_PROGRAMS += test-todd 
nodist_test_todd_SOURCES = todd/test-todd.cpp

bin_PROGRAMS += aggregate
nodist_aggregate_SOURCES = aggregation/aggregation.cpp

bin_PROGRAMS += normaliz
nodist_normaliz_SOURCES = normalize/normalize.cpp

endif

bin_PROGRAMS += triangulate
triangulate_SOURCES = triangulation/test-triangulation.cpp

SUBDIRS = adjacency

# Maple code
dist_pkgdata_DATA = simplify.add simplify2.add simplify3.add
