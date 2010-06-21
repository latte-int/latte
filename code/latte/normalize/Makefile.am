AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS =
noinst_LIBRARIES =
noinst_HEADERS =

if HAVE_FORTYTWO_LIB
noinst_LIBRARIES += libnormalize.a
libnormalize_a_SOURCES = ReductionTest.cpp normalize.cpp
noinst_HEADERS += ReductionTest.h
AM_CPPFLAGS = -I$(srcdir)/..
AM_LDFLAGS = 
LDADD = libnormalize.a ../liblatte.a

# For Gnulib relocatable-prog
AM_CPPFLAGS += -DINSTALLDIR=\"$(bindir)\"
if RELOCATABLE_VIA_LD
AM_LDFLAGS += `$(RELOCATABLE_LDFLAGS) $(bindir)`
endif
LDADD += ../../gnulib/libgnu.a

if HAVE_CPLEX
libnormalize_a_SOURCES += ReductionTestWithCPLEX.cpp 
noinst_HEADERS += ReductionTestWithCPLEX.h
endif

AM_CXXFLAGS = $(NTL_CFLAGS) $(LIDIA_CFLAGS) $(TOPCOM_CXXFLAGS) $(CDDLIB_CFLAGS) $(FORTYTWO_CXXFLAGS) $(GMP_CFLAGS)
LIBS += $(NTL_LIBS) $(LIDIA_LIBS) $(TOPCOM_LIBS) $(CDDLIB_LIBS) $(FORTYTWO_LIBS) $(GMP_LIBS)

if HAVE_CPLEX
AM_LDFLAGS += $(CPLEX_LDFLAGS)
AM_CXXFLAGS += $(CPLEX_CFLAGS)
endif

bin_PROGRAMS += hilbert-from-rays
hilbert_from_rays_SOURCES = normalize_main.cpp

if HAVE_CPLEX
bin_PROGRAMS += check-representation
check_representation_SOURCES = check-representation.cpp
endif

bin_PROGRAMS += hilbert-from-rays-symm
hilbert_from_rays_symm_SOURCES = normaliz_wrapper.cpp
hilbert_from_rays_symm_LDADD = $(LDADD) -l4ti2util

endif
