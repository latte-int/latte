nobase_dist_pkgdata_DATA = 											\
	Conebyconeapproximations_08_11_2010.mpl RealBarvinok-mars-exemples-2014-03-10.mpl m-knapsack.mpl	\
	README

# not installed:
EXTRA_DIST = testMapleLatteKnapsack.mpl

### For automake's old "serial test harness"
##TESTS_ENVIRONMENT = ./test-maple.sh -D CHECK_EXAMPLES=true -I $(srcdir)

LOG_COMPILER = ./test-maple.sh -D CHECK_EXAMPLES=true -I $(srcdir)

TESTS =

if HAVE_MAPLE
TESTS += Conebyconeapproximations_08_11_2010.mpl RealBarvinok-mars-exemples-2014-03-10.mpl 
TESTS += m-knapsack.mpl testMapleLatteKnapsack.mpl
endif


CLEANFILES = knapsackEquation.latte knapsackEquation.latte.topKehrhart.mpl

