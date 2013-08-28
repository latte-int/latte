#noinst_PROGRAMS = set-max-cpu-time

MAXRUNTIME=10

TESTS_ENVIRONMENT=env MAXRUNTIME=$(MAXRUNTIME)

BUILT_SOURCES = test.pl

TESTS =

TESTS +=					\
	test-homog.sh				\
	test-default.sh				\
	test-irr-exp-maxdet1000.sh		\
	test-allprimal-exp-maxdet1000.sh

#### Using 4ti2 is actually the new default, no need to test it separately.
## if HAVE_FORTYTWO_LIB
## TESTS +=					\
## 	test-4ti2.sh				\
## 	test-irr-exp-maxdet1000-4ti2.sh		\
## 	test-allprimal-exp-maxdet1000-4ti2.sh
## endif

## Tests with defaults of older versions
TESTS += 					\
	test-cdd.sh				\
	test-cdd-cddlib.sh

EXTRA_DIST = $(TESTS)

clean-local:
	-rm -rf results

installcheck-local:
	$(TESTS_ENVIRONMENT) USE_INSTALLED_BINARY=t $(srcdir)/test-default.sh
	$(TESTS_ENVIRONMENT) USE_INSTALLED_BINARY=t $(srcdir)/test-cdd.sh
