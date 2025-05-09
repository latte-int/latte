AUTOMAKE_OPTIONS = subdir-objects

SUBDIRS = .
SUBDIRS += normalize

TESTS =

bin_PROGRAMS = latte-maximize latte-minimize count ehrhart latte2ine latte2ext

latte_maximize_SOURCES = maximize.cpp
latte_minimize_SOURCES = minimize.cpp
count_SOURCES = countDriver.cpp
ehrhart_SOURCES = ehrhart.cpp
latte2ine_SOURCES = latte2ine.cpp
latte2ext_SOURCES = latte2ext.cpp


LDADD = liblatte.la
AM_CPPFLAGS = -I$(top_srcdir)/code
AM_CXXFLAGS = $(NTL_CFLAGS) $(LIDIA_CFLAGS) $(TOPCOM_CXXFLAGS) $(CDDLIB_CFLAGS) $(FORTYTWO_CXXFLAGS) $(GMP_CFLAGS)
LIBS += $(NTL_LIBS) $(LIDIA_LIBS) $(TOPCOM_LIBS) $(CDDLIB_LIBS) $(FORTYTWO_LIBS) $(GMP_LIBS)
AM_LDFLAGS = 

# For Gnulib relocatable-prog
AM_CPPFLAGS += -DINSTALLDIR=\"$(bindir)\"
if RELOCATABLE_VIA_LD
AM_LDFLAGS += `$(RELOCATABLE_LDFLAGS) $(bindir)`
endif
LDADD += ../gnulib/libgnu.la

lib_LTLIBRARIES = liblatte.la
# The following line is necessary for the Mac OS X enable-shared build.
# Otherwise, get "dyld: lazy symbol binding failed: Symbol not found: _set_program_name_and_installdir"
# when running the binary.  --mkoeppe
liblatte_la_LIBADD = ../gnulib/libgnu.la
liblatte_la_SOURCES =							\
	barvinok/dec.cpp barvinok/barvinok.cpp  barvinok/ComputeOmega.cpp	\
	count.cpp \
	LattException.cpp \
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
	count.h 								\
	LattException.h						\
	ConeInfo.h Grobner.h IntegralHull.h PolyTree.h ReadingFile.h	\
	RudyResNTL.h banner.h binarySearchIP.h cone.h dual.h flags.h	\
	preprocess.h print.h ramon.h rational.h timing.h		\
	Residue.h barvinok/ComputeOmega.h					\
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
	ReadSubcones.h \
	integration/GeneralMonomialSum.h 


if HAVE_TOPCOM_LIB
liblatte_la_SOURCES += triangulation/TriangulationWithTOPCOM.cpp
noinst_HEADERS += triangulation/TriangulationWithTOPCOM.h
else
if HAVE_TOPCOM_BIN
liblatte_la_SOURCES += triangulation/TriangulationWithTOPCOM.cpp
noinst_HEADERS += triangulation/TriangulationWithTOPCOM.h
endif
endif

if HAVE_CDDLIB
liblatte_la_SOURCES += latte_cddlib.cpp
liblatte_la_SOURCES += triangulation/RegularTriangulationWithCddlib.cpp
noinst_HEADERS += triangulation/RegularTriangulationWithCddlib.h
endif

if HAVE_LIDIA
liblatte_la_SOURCES += genFunction/NTL_to_LiDIA.cpp
noinst_HEADERS += genFunction/NTL_to_LiDIA.h
endif

if HAVE_FORTYTWO_LIB
liblatte_la_SOURCES += latte_4ti2.cpp
liblatte_la_SOURCES += triangulation/RegularTriangulationWith4ti2.cpp
liblatte_la_SOURCES += VertexConesWith4ti2.cpp
noinst_HEADERS += latte_4ti2.h
noinst_HEADERS += triangulation/RegularTriangulationWith4ti2.h
noinst_HEADERS += VertexConesWith4ti2.h

liblatte_la_SOURCES += DualizationWith4ti2.cpp
noinst_HEADERS += DualizationWith4ti2.h

liblatte_la_SOURCES += latte_4ti2_zsolve.cpp
noinst_HEADERS += latte_4ti2_zsolve.h
endif

bin_PROGRAMS += polyhedron-to-cones
polyhedron_to_cones_SOURCES = polyhedron-to-cones.cpp

if HAVE_EXPERIMENTS
#Experimental stuff
nodist_liblatte_la_SOURCES =				\
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

nodist_liblatte_la_SOURCES += triangulation/SpecialSimplex.cpp 
nodist_noinst_HEADERS += triangulation/SpecialSimplex.h
bin_PROGRAMS += special-simplex
nodist_special_simplex_SOURCES = triangulation/special-simplex.cpp
special_simplex_LDFLAGS = $(CPLEX_LDFLAGS)
special_simplex_CXXFLAGS = $(AM_CXXFLAGS) $(CPLEX_CFLAGS)
endif
endif #end HAVE_EXPERIMENTS

bin_PROGRAMS += triangulate
triangulate_SOURCES = triangulation/test-triangulation.cpp

check_PROGRAMS =

if HAVE_LIDIA
check_PROGRAMS += ilio_test
ilio_test_SOURCES = genFunction/ilio.cpp
TESTS += ilio_test
endif

# Maple code
dist_pkgdata_DATA = simplify.add simplify2.add simplify3.add

###### VOLUME COMPUTATION AND RATIONAL VALUATIONS ##########


#main valuation headers and source
noinst_HEADERS += valuation/PolytopeValuation.h
noinst_HEADERS += valuation/valuation.h
noinst_HEADERS += valuation/Perturbation.h

liblatte_la_SOURCES += valuation/PolytopeValuation.cpp
liblatte_la_SOURCES += valuation/valuation.cpp
liblatte_la_SOURCES += valuation/Perturbation.cpp

#Polytope building headers and source.
noinst_HEADERS += buildPolytopes/BuildPolytope.h
noinst_HEADERS += buildPolytopes/BuildGraphPolytope.h
noinst_HEADERS += buildPolytopes/BuildHypersimplexEdgePolytope.h
noinst_HEADERS += buildPolytopes/BuildRandomPolytope.h
noinst_HEADERS += buildPolytopes/GraphMaker.h
noinst_HEADERS += buildPolytopes/BuildRandomPolynomials.h

liblatte_la_SOURCES  += buildPolytopes/BuildPolytope.cpp
liblatte_la_SOURCES  += buildPolytopes/BuildGraphPolytope.cpp
liblatte_la_SOURCES  += buildPolytopes/BuildHypersimplexEdgePolytope.cpp
liblatte_la_SOURCES  += buildPolytopes/BuildRandomPolytope.cpp
liblatte_la_SOURCES  += buildPolytopes/GraphMaker.cpp
liblatte_la_SOURCES  += buildPolytopes/BuildRandomPolynomials.cpp

#top knapsack headers and sources. 
noinst_HEADERS += top-knapsack/PeriodicFunction.h
noinst_HEADERS += top-knapsack/TopKnapsack.h
liblatte_la_SOURCES += top-knapsack/TopKnapsack.cpp
liblatte_la_SOURCES += top-knapsack/PeriodicFunction.cpp

#interpolation. to delete this.
noinst_HEADERS += interpolation/PolynomialInterpolation.h



#main program to find volumes and integrals of polytopes.
bin_PROGRAMS += integrate
integrate_SOURCES = valuation/Driver.cpp

#main program to find the top k terms of the Ehrhart polynomial for knapsacks
bin_PROGRAMS += top-ehrhart-knapsack
top_ehrhart_knapsack_SOURCES = top-knapsack/Driver.cpp
top_ehrhart_knapsack_LDADD = $(LDADD)

#convert a cdd H-rep to latte H-rep (by clearing the fractions)
bin_PROGRAMS += ConvertCDDineToLatte
ConvertCDDineToLatte_SOURCES = convertCDD_ine_to_latte.cpp
ConvertCDDineToLatte_LDADD = $(LDADD)

#convert a cdd V-rep to latte V-rep
bin_PROGRAMS += ConvertCDDextToLatte
ConvertCDDextToLatte_SOURCES = convertCDD_ext_to_latte.cpp
ConvertCDDextToLatte_LDADD = $(LDADD)


#Run ./valuation/test/integrateHyperrectangleTest.sh 
#  to start a maple script that will test integration of polynomials over rectangles in many dimensions.
check_PROGRAMS += test-hyperrectangle-integration
test_hyperrectangle_integration_SOURCES = valuation/test/testIntegrationHyperrectanglesDriver.cpp 

#A collection of examples of how to use latte as a lib.
check_PROGRAMS += test-library-examples
test_library_examples_SOURCES = valuation/testLibraryExamples.cpp
TESTS += test-library-examples

### test-library-examples bleeds lots of temporary files and result files into the current directory.
### Delete them.
DISTCLEANFILES = 								\
	4ti2_triangulation_output Check_emp.lp Check_emp.lps Check_emp.out 	\
	tri.ext tri.ine tri.iad tri.ead tri.icd tri.ecd				\
	latte_stats lifted_cone_for_4ti2_triangulation numOfLatticePoints 	\
	numOfUnimodularCones test-library-examples.hrep totalTime

#hack...have_Programmer_extra was giving errors
#if HAVE_PROGRAMMER_EXTRA
#if HAVE_DATABASE

#we should delete this, it only makes sure the BuidPolytope base class works.~Brandon
## bin_PROGRAMS += testBuildPolytope
## testBuildPolytope_SOURCES = valuation/test/testBuildPolytopes.cpp
## testBuildPolytope_SOURCES += buildPolytopes/BuildPolytope.cpp



bin_PROGRAMS += ehrhart3
ehrhart3_SOURCES = ehrhart3.cpp
ehrhart3_SOURCES += interpolation/PolynomialInterpolation.cpp


#counts the number of linear forms in a polynomisal.
bin_PROGRAMS += count-linear-forms-from-polynomial
count_linear_forms_from_polynomial_SOURCES = valuation/countLinearFormsFromPolynomialDriver.cpp

#we should remove this at some point. 
## bin_PROGRAMS += testPolynomialInterpolation
## testPolynomialInterpolation_SOURCES = interpolation/PolynomialInterpolation.cpp
## testPolynomialInterpolation_SOURCES += interpolation/testDriver.cpp

#we should remove this at some point.
## bin_PROGRAMS += my-test-ehrhart
## my_test_ehrhart_SOURCES = buildPolytopes/Driver.cpp
## my_test_ehrhart_SOURCES += buildPolytopes/BuildRandomPolytope.cpp
## my_test_ehrhart_SOURCES += buildPolytopes/BuildHypersimplexEdgePolytope.cpp
## my_test_ehrhart_SOURCES += buildPolytopes/GraphMaker.cpp
## my_test_ehrhart_SOURCES += buildPolytopes/BuildGraphPolytope.cpp

#converts a really big fraction (say 70000 digits) to decimal.
#bin_PROGRAMS += convert-rational-to-decimal
#convert_rational_to_decimal_SOURCES = valuation/test/convertRationalToDecimalDriver.cpp
#convert_rational_to_decimal_LDADD = $(LDADD)

#check_PROGRAMS += test-rationalNTL
#test_rationalNTL_SOURCES = testRationalNTLDriver.cpp 
#test_rationalNTL_SOURCES += rational.cpp
#test_rationalNTL_LDADD = $(LDADD)


#A maple script can find the volume of many simplices at once using the PolytopeValuation class.
check_PROGRAMS += test-volume-for-SL
test_volume_for_SL_SOURCES = valuation/test/testVolumeForSLDriver.cpp

#endif #end of  HAVE_PROGRAMMER_EXTRA



if HAVE_DATABASE

#This is used to make many polytopes for the new research paper.
bin_PROGRAMS += build-database-polytopes
build_database_polytopes_SOURCES = valuation/test/buildDatabaseRandomPolytopesDriver.cpp
build_database_polytopes_SOURCES += sqlite/IntegrationDB.cpp
build_database_polytopes_SOURCES += sqlite/SqliteDB.cpp
build_database_polytopes_LDFLAGS = -lsqlite3

#This is used to test integration classes from the database
bin_PROGRAMS += integrate-database-polytopes
integrate_database_polytopes_SOURCES = valuation/test/runIntegrationTestFromDatabaseDriver.cpp
integrate_database_polytopes_SOURCES += sqlite/IntegrationDB.cpp
integrate_database_polytopes_SOURCES += sqlite/SqliteDB.cpp
integrate_database_polytopes_LDFLAGS = -lsqlite3

#This is used to test volume classes from the database
bin_PROGRAMS += find-volume-database-polytopes
find_volume_database_polytopes_SOURCES = valuation/test/runVolumeTestFromDatabaseDriver.cpp
find_volume_database_polytopes_SOURCES += sqlite/VolumeDB.cpp
find_volume_database_polytopes_SOURCES += sqlite/SqliteDB.cpp
find_volume_database_polytopes_LDFLAGS = -lsqlite3

#This is used to test integration classes from the database
bin_PROGRAMS += fix-dual-database-polytopes
fix_dual_database_polytopes_SOURCES = valuation/test/fixDualPolytopesInDatabaseDriver.cpp
fix_dual_database_polytopes_SOURCES += sqlite/IntegrationDB.cpp
fix_dual_database_polytopes_SOURCES += sqlite/SqliteDB.cpp
fix_dual_database_polytopes_LDFLAGS = -lsqlite3

#This is used update the lawrence valuations.
bin_PROGRAMS += fix-lawrence-valuation
fix_lawrence_valuation_SOURCES = valuation/test/fixLawrenceValuationDriver.cpp
fix_lawrence_valuation_SOURCES += sqlite/IntegrationDB.cpp
fix_lawrence_valuation_SOURCES += sqlite/SqliteDB.cpp
fix_lawrence_valuation_LDFLAGS = -lsqlite3

#This is used to test integration classes from the database
bin_PROGRAMS += print-database-time-table 
print_database_time_table_SOURCES = valuation/test/printDatabaseTimeTablesDriver.cpp
print_database_time_table_SOURCES += sqlite/IntegrationDB.cpp
print_database_time_table_SOURCES += sqlite/VolumeDB.cpp
print_database_time_table_SOURCES += sqlite/SqliteDB.cpp
print_database_time_table_LDFLAGS = -lsqlite3

endif #HAVE_DATABASE


################ INTEGRATION ####################

liblatte_la_SOURCES +=				\
	integration/PolyRep.cpp			\
	integration/residue.cpp			\
	integration/newIntegration.cpp		\
	integration/PolyTrie.cpp		\
	top-ehrhart/TopEhrhart.cpp		

noinst_HEADERS +=				\
	integration/multiply.h			\
	integration/iterators.h			\
	integration/iterators.hpp		\
	integration/consumers.h			\
	integration/PolyRep.h			\
	integration/PolyTrie.h			\
	integration/newIntegration.h		\
	integration/burstTrie.h			\
	integration/burstTrie.hpp		\
	integration/blockReps.h			\
	integration/residue.h			\
	top-ehrhart/TopEhrhart.h

#if HAVE_PROGRAMMER_EXTRA
if HAVE_DATABASE
bin_PROGRAMS += BlockIntegrate
BlockIntegrate_SOURCES = integration/blockInt.cpp

bin_PROGRAMS += BurstIntegrate
BurstIntegrate_SOURCES = integration/burstInt.cpp

#bin_PROGRAMS += BurstBench
#BurstBench_SOURCES = integration/burstBench.cpp

#bin_PROGRAMS += BlockBench
#BlockBench_SOURCES = integration/blockBench.cpp

bin_PROGRAMS += integrate_demo
integrate_demo_SOURCES =    integration/integrateDemo.cpp

check_PROGRAMS += compare
compare_SOURCES = integration/versusMode.cpp

endif #end of HAVE_PROGRAMMER_EXTRA

check_PROGRAMS += integrate_test
integrate_test_SOURCES =    integration/integrateTest.cpp
TESTS += integration/testIntegration.sh

check_PROGRAMS += multiply_test
multiply_test_SOURCES =   integration/multiplyTest.cpp
TESTS +=  integration/testMultiply.sh

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

