AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS =
noinst_LIBRARIES =
nodist_noinst_HEADERS =

noinst_LIBRARIES += libnormaliz.a
libnormaliz_a_SOURCES = ReductionTest.cpp
nodist_noinst_HEADERS += ReductionTest.h
AM_CPPFLAGS = -I$(srcdir)/..
if HAVE_CPLEX
libnormaliz_a_SOURCES += ReductionTestWithCPLEX.cpp
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
nodist_normaliz_SOURCES = normalize.cpp

# if HAVE_CPLEX
# bin_PROGRAMS += check-representation
# nodist_check_representation_SOURCES = check-representation.cpp
# endif

