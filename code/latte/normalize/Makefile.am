AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS =
noinst_LIBRARIES =
nodist_noinst_HEADERS =

noinst_LIBRARIES += libnormaliz.a
nodist_libnormaliz_a_SOURCES = ReductionTest.cpp normalize.cpp
nodist_noinst_HEADERS += ReductionTest.h
AM_CPPFLAGS = -I$(srcdir)/..
if HAVE_CPLEX
nodist_libnormaliz_a_SOURCES += ReductionTestWithCPLEX.cpp 
nodist_noinst_HEADERS += ReductionTestWithCPLEX.h
endif

LDADD = libnormaliz.a ../liblatte.a
AM_CXXFLAGS = $(NTL_CFLAGS) $(LIDIA_CFLAGS) $(TOPCOM_CXXFLAGS) $(CDDLIB_CFLAGS) $(FORTYTWO_CXXFLAGS) $(GMP_CFLAGS)
LIBS += $(NTL_LIBS) $(LIDIA_LIBS) $(TOPCOM_LIBS) $(CDDLIB_LIBS) $(FORTYTWO_LIBS) $(GMP_LIBS)

if HAVE_CPLEX
AM_LDFLAGS = $(CPLEX_LDFLAGS)
AM_CXXFLAGS += $(CPLEX_CFLAGS)
endif

bin_PROGRAMS += normaliz
nodist_normaliz_SOURCES = normalize_main.cpp

if HAVE_CPLEX
bin_PROGRAMS += check-representation
nodist_check_representation_SOURCES = check-representation.cpp
endif

bin_PROGRAMS += normaliz_wrapper
nodist_normaliz_wrapper_SOURCES = normaliz_wrapper.cpp
normaliz_wrapper_LDADD = $(LDADD) -l4ti2util
