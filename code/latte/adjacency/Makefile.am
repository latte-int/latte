if HAVE_CDDLIB
libexec_PROGRAMS = ComputeAdjacency
ComputeAdjacency_SOURCES = ComputeAdjacency.cpp
endif 

AM_CXXFLAGS = -DGMPRATIONAL $(CDDLIB_CFLAGS) $(GMP_CFLAGS)
LIBS += $(CDDLIB_LIBS) $(GMP_LIBS)
