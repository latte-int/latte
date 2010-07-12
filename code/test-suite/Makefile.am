#noinst_PROGRAMS = set-max-cpu-time

MAXRUNTIME=10

TESTS_ENVIRONMENT=env MAXRUNTIME=$(MAXRUNTIME) $(SHELL)

BUILT_SOURCES = test.pl

TESTS =

TESTS +=					\
	test-homog.sh				\
	test-default.sh				\
	test-irr-exp-maxdet1000.sh		\
	test-allprimal-exp-maxdet1000.sh

## Tests with defaults of older versions
TESTS += 					\
	test-cdd.sh				\
	test-cdd-cddlib.sh

if HAVE_FORTYTWO_LIB
TESTS +=					\
	test-4ti2.sh				\
	test-irr-exp-maxdet1000-4ti2.sh		\
	test-allprimal-exp-maxdet1000-4ti2.sh
endif

EXTRA_DIST = $(TESTS)


