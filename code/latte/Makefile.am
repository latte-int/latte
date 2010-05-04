AUTOMAKE_OPTIONS = subdir-objects

SUBDIRS = adjacency . 
SUBDIRS += normalize

bin_PROGRAMS = latte-maximize latte-minimize count ehrhart latte2ine latte2ext

latte_maximize_SOURCES = maximize.cpp
latte_minimize_SOURCES = minimize.cpp
count_SOURCES = count.cpp
ehrhart_SOURCES = ehrhart.cpp
latte2ine_SOURCES = latte2ine.cpp
latte2ext_SOURCES = latte2ext.cpp

LDADD = liblatte.a
AM_CPPFLAGS = -I$(top_srcdir)/code
AM_CXXFLAGS = $(NTL_CFLAGS) $(LIDIA_CFLAGS) $(TOPCOM_CXXFLAGS) $(CDDLIB_CFLAGS) $(FORTYTWO_CXXFLAGS) $(GMP_CFLAGS)
LIBS += $(NTL_LIBS) $(LIDIA_LIBS) $(TOPCOM_LIBS) $(CDDLIB_LIBS) $(FORTYTWO_LIBS) $(GMP_LIBS)
AM_LDFLAGS = 

# For Gnulib relocatable-prog
AM_CPPFLAGS += -DINSTALLDIR=\"$(bindir)\"
if RELOCATABLE_VIA_LD
AM_LDFLAGS += `$(RELOCATABLE_LDFLAGS) $(bindir)`
endif
LDADD += ../gnulib/libgnu.a

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
	latte_relocatable.cpp						\
	genFunction/IntCombEnum.cpp					\
	genFunction/matrix_ops.cpp					\
	ExponentialSubst.cpp ExponentialEhrhart.cpp			\
	Irrational.cpp							\
	ReadLatteStyle.cpp						\
	cone_consumer.cpp						\
	triangulation/triangulate.cpp					\
	triangulation/RegularTriangulation.cpp				\
	triangulation/RegularTriangulationWithCdd.cpp			\
	ReadPolyhedron.cpp						\
	ProjectUp.cpp							\
	ReadSubcones.cpp

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
	latte_relocatable.h						\
	latte_random.h latte_system.h latte_cddlib.h			\
	genFunction/IntCombEnum.h					\
	genFunction/matrix_ops.h					\
	ExponentialSubst.h ExponentialEhrhart.h				\
	Irrational.h							\
	ReadLatteStyle.h Polyhedron.h					\
	cone_consumer.h							\
	triangulation/triangulate.h					\
	triangulation/RegularTriangulation.h				\
	triangulation/RegularTriangulationWithCdd.h			\
	ReadPolyhedron.h						\
	ProjectUp.h							\
	ReadSubcones.h

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
liblatte_a_SOURCES += VertexConesWith4ti2.cpp
noinst_HEADERS += latte_4ti2.h
noinst_HEADERS += triangulation/RegularTriangulationWith4ti2.h
noinst_HEADERS += VertexConesWith4ti2.h

liblatte_a_SOURCES += DualizationWith4ti2.cpp
noinst_HEADERS += DualizationWith4ti2.h

liblatte_a_SOURCES += latte_4ti2_zsolve.cpp
noinst_HEADERS += latte_4ti2_zsolve.h
endif

bin_PROGRAMS += polyhedron-to-cones
polyhedron_to_cones_SOURCES = polyhedron-to-cones.cpp

if HAVE_EXPERIMENTS
# Experimental stuff

nodist_liblatte_a_SOURCES =				\
	ExponentialApprox.cpp				\
	heap.cpp					\
	TrivialSubst.cpp				\
	barvinok/SubspaceAvoidingDecomposition.cpp	\
	triangulation/BoundaryTriangulation.cpp	

nodist_noinst_HEADERS =					\
	ExponentialApprox.h				\
	heap.h						\
	TrivialSubst.h					\
	barvinok/SubspaceAvoidingDecomposition.h	\
	triangulation/BoundaryTriangulation.h

bin_PROGRAMS += test-todd 
nodist_test_todd_SOURCES = todd/test-todd.cpp

bin_PROGRAMS += aggregate
nodist_aggregate_SOURCES = aggregation/aggregation.cpp

bin_PROGRAMS += triangulate-boundary
nodist_triangulate_boundary_SOURCES = triangulation/triangulate-boundary.cpp 

bin_PROGRAMS += complete-boundary-triangulation
nodist_complete_boundary_triangulation_SOURCES = triangulation/complete-boundary-triangulation.cpp 

if HAVE_CPLEX
AM_LDFLAGS += $(CPLEX_LDFLAGS)
AM_CXXFLAGS += $(CPLEX_CFLAGS)

nodist_liblatte_a_SOURCES += triangulation/SpecialSimplex.cpp 
nodist_noinst_HEADERS += triangulation/SpecialSimplex.h

bin_PROGRAMS += special-simplex
nodist_special_simplex_SOURCES = triangulation/special-simplex.cpp
special_simplex_LDFLAGS = $(CPLEX_LDFLAGS)
special_simplex_CXXFLAGS = $(AM_CXXFLAGS) $(CPLEX_CFLAGS)
endif

endif

bin_PROGRAMS += triangulate
triangulate_SOURCES = triangulation/test-triangulation.cpp

#bin_PROGRAMS += decompose
#decompose_SOURCES =     integration/decompose.cpp
#decompose_SOURCES +=    integration/PolyRep.cpp

#bin_PROGRAMS += integrate
#integrate_SOURCES = integration/decint.cpp
#integrate_SOURCES += integration/newIntegration.cpp
#integrate_SOURCES += integration/PolyRep.cpp
#integrate_SOURCES += integration/residue.cpp
#integrate_SOURCES += integration/multiply.h
#integrate_SOURCES +=     integration/timing.cpp

#bin_PROGRAMS += BurstIntegrate
#BurstIntegrate_SOURCES = integration/burstint.cpp
#BurstIntegrate_SOURCES += integration/newIntegration.cpp
#BurstIntegrate_SOURCES += integration/PolyRep.cpp
#BurstIntegrate_SOURCES += integration/PolyTrie.cpp
#BurstIntegrate_SOURCES += integration/residue.cpp
#BurstIntegrate_SOURCES += integration/multiply.h

#bin_PROGRAMS += BurstTest
#BurstTest_SOURCES = integration/burstTest.cpp
#BurstTest_SOURCES += integration/newIntegration.cpp
#BurstTest_SOURCES += integration/PolyRep.cpp
#BurstTest_SOURCES += integration/PolyTrie.cpp
#BurstTest_SOURCES += integration/residue.cpp
#BurstTest_SOURCES += integration/multiply.h
#BurstTest_SOURCES += integration/timing.cpp

bin_PROGRAMS += integrate_demo
integrate_demo_SOURCES =    integration/integrateDemo.cpp
integrate_demo_SOURCES +=   integration/PolyTrie.cpp
integrate_demo_SOURCES +=   integration/newIntegration.cpp
integrate_demo_SOURCES +=   integration/residue.cpp

check_PROGRAMS = ilio_test
ilio_test_SOURCES = genFunction/ilio.cpp

check_PROGRAMS += integrate_test
integrate_test_SOURCES =    integration/integrateTest.cpp
integrate_test_SOURCES +=   integration/PolyTrie.cpp
integrate_test_SOURCES +=   integration/newIntegration.cpp
integrate_test_SOURCES +=   integration/residue.cpp

check_PROGRAMS += multiply_test
multiply_test_SOURCES =    integration/PolyTrie.cpp
multiply_test_SOURCES +=   integration/multiplyTest.cpp

check_PROGRAMS += compare
compare_SOURCES = integration/versusMode.cpp
compare_SOURCES += integration/PolyTrie.cpp
compare_SOURCES += integration/PolyRep.cpp

TESTS =  integration/testMultiply.sh
TESTS += integration/testIntegration.sh
#TESTS += ilio_test

bench-table5: 
	integration/benchmark.sh 5 integration/table5.txt
bench-table6: 
	integration/benchmark.sh 6 integration/table6.txt
bench-table10: 
	integration/benchmark.sh 10 integration/table10.txt
bench-table11: 
	integration/benchmark.sh 11 integration/table11.txt
bench-table13: 
	integration/benchmark.sh 13 integration/table13.txt
bench-all:
	bench-table5 bench-table6 bench-table10 bench-table11 bench-table13
.PHONY: bench-all bench-table5 bench-table6 bench-table10 bench-table11 bench-table13

# Maple code
dist_pkgdata_DATA = simplify.add simplify2.add simplify3.add
