#!/usr/bin/perl -w

use File::Basename;

# count always prints the number of lattice points in the file
# numOfLatticePoints. In order to check the output of the test,
# we must open the file and verify the number of lattice points
# stored there.
$OUTPUT_FILE_NAME = "numOfLatticePoints";

# add file names and lattice points to an array
# We obtained the results marked "non-authoritative" ourselves, 
# using various flavors of the algorithm; we hope they are correct.
# The other ones are obtained from other sources; we believe 
# they are correct. 
my @files_nums = (
   # FILENAME::CORRECT-ANSWER::TIME-LIMIT
   "example1::3",
   "example2::2",
   "example3::2",
   "24_cell::33",
   "cuww1::1",
   "magic4x4::8",
   "magic5x5::::86400",
   "hickerson/hickerson-7::8",
   "hickerson/hickerson-8::22",
   "hickerson/hickerson-9::9",
   "hickerson/hickerson-10::24",
   "hickerson/hickerson-11::12",
   "hickerson/hickerson-12::38",
   "hickerson/hickerson-13::14::1800",
   "hickerson/hickerson-14::32::86400",
   "hickerson/hickerson-15::20::1800",
   "hickerson/hickerson-16::54::86400",
   "hickerson/hickerson-17::18::86400",
   "hickerson/hickerson-18::44::86400",
   "hickerson/hickerson-19::20::86400",
   "hickerson/hickerson-20::74::86400",
   "hickerson/hickerson-24::96::86400",
   "hickerson/hickerson-28::92::86400",
   "hickerson/hickerson-32::122::86400",
   "hickerson/hickerson-36::138::86400",
   "mcallister/HivePolytopeImpossible::86400",
    "mcallister/HivePolytopeMinutes1::14438636", #not-authoritative
    "mcallister/HivePolytopeMinutes2::4890291",  #not-authoritative
    "mcallister/HivePolytopeMinutes3::1990152",  #not-authoritative
    "mcallister/HivePolytopeMinutes4::0",        #not-authoritative
    "mcallister/HivePolytopeMinutes5::2129924",	 #not-authoritative
    "mcallister/HivePolytopeSeconds::5231",      #not-authoritative
    "yoshida/24_cell_latte::33",
    "yoshida/3x3x3_semi_10.equ::1273125",
    "yoshida/3x3x3_semi_11.equ::2467302",
    "yoshida/3x3x3_semi_12.equ::4547458",
    "yoshida/3x3x3_semi_13.equ::8027223",
    "yoshida/3x3x3_semi_14.equ::13648170",
    "yoshida/3x3x3_semi_15.equ::22454470",
    "yoshida/3x3x3_semi_16.equ::35884827",
    "yoshida/3x3x3_semi_17.equ::55883718",
    "yoshida/3x3x3_semi_18.equ::85034962",
    "yoshida/3x3x3_semi_1.equ::12",
    "yoshida/3x3x3_semi_2.equ::132",
    "yoshida/3x3x3_semi_3.equ::847",
    "yoshida/3x3x3_semi_4.equ::3921",
    "yoshida/3x3x3_semi_5.equ::14286",
    "yoshida/3x3x3_semi_6.equ::43687",
    "yoshida/3x3x3_semi_7.equ::116757",
    "yoshida/3x3x3_semi_8.equ::280656",
    "yoshida/3x3x3_semi_9.equ::619219",
    "yoshida/4flow_1.equ::34441480172695101274",
    "yoshida/4flow_2.equ::28493245103068590026",
    "yoshida/4flow_3.equ::91608082255943644656",
    "yoshida/4x4Table1::1225914276768514",
    "yoshida/4x4Table10::63313191414342827754566531364533378588986467",
    "yoshida/4x4Table11::209",
    "yoshida/4x4Table2::993810896945891",
    "yoshida/4x4Table3::25387360604030",
    "yoshida/4x4Table4::13571026063401838164668296635065899923152079",
    "yoshida/4x4Table5::646911395459296645200004000804003243371154862",
    "yoshida/4x4Table6::319720249690111437887229255487847845310463475",
    "yoshida/4x4Table7::322773560821008856417270275950599107061263625",
    "yoshida/4x4Table8::6977523720740024241056075121611021139576919",
    "yoshida/4x4Table9::861316343280649049593236132155039190682027614",
    "yoshida/4x5_1::316052820930116909459822049052149787748004963058022997262397",
    "yoshida/4x5_2::23196436596128897574829611531938753",
    "yoshida/4x5_3::23196436596128897574829611531938753",
    "yoshida/5flow_1.equ::6817997013081449330251623043931489475270",
    "yoshida/5flow_2.equ::277145720781272784955528774814729345461",
    "yoshida/5flow_3.equ::710305971948234346520365668331191134724",
    "yoshida/aardallenstra.equ::0",
    "yoshida/cube::16",
    "yoshida/cube2::8",
    "yoshida/digraph4_10.equ::179777378508547",
    "yoshida/digraph4_1.equ::223",
    "yoshida/digraph4_2.equ::330",
    "yoshida/digraph4_3.equ::3002",
    "yoshida/digraph4_4.equ::785528058",
    "yoshida/digraph4_5.equ::20673947895",
    "yoshida/digraph4_6.equ::14100406254",
    "yoshida/digraph4_7.equ::1906669380",
    "yoshida/digraph4_8.equ::19470466783680",
    "yoshida/digraph4_9.equ::106036300535520",
    "yoshida/digraph5_10.equ::65348330279808617817420057",
    "yoshida/digraph5_1.equ::14805",
    "yoshida/digraph5_2.equ::6950747024",
    "yoshida/digraph5_3.equ::222850218035543",
    "yoshida/digraph5_4.equ::563408416219655157542748",
    "yoshida/digraph5_5.equ::1108629405144880240444547243",
    "yoshida/digraph5_6.equ::3997121684242603301444265332",
    "yoshida/digraph5_7.equ::160949617742851302259767600",
    "yoshida/digraph5_8.equ::15711217216898158096466094",
    "yoshida/digraph5_9.equ::102815492358112722152328",
    "yoshida/hyp_simp_4_1.equ::4",
    "yoshida/hyp_simp_4_2.equ::6",
    "yoshida/hyp_simp_4_3.equ::4",
    "yoshida/hyp_simp_5_1.equ::5",
    "yoshida/hyp_simp_5_2.equ::10",
    "yoshida/hyp_simp_5_3.equ::10",
    "yoshida/hyp_simp_5_4.equ::5",
    "yoshida/hyp_simp_6_1.equ::6",
    "yoshida/hyp_simp_6_2.equ::15",
    "yoshida/hyp_simp_6_3.equ::20",
    "yoshida/hyp_simp_6_4.equ::15",
    "yoshida/hyp_simp_6_5.equ::6",
    "yoshida/hyp_simp_7_1.equ::7",
    "yoshida/hyp_simp_7_2.equ::21",
    "yoshida/hyp_simp_7_3.equ::35",
    "yoshida/hyp_simp_7_4.equ::35",
    "yoshida/hyp_simp_7_5.equ::21",
    "yoshida/hyp_simp_7_6.equ::7",
    "yoshida/knapsack1.equ::42",
    "yoshida/knapsack2.equ::92378",
    "yoshida/knapsackbaby1.equ::7",
    "yoshida/knapsackbaby2.equ::2",
    "yoshida/mountExample1::35353",
    "yoshida/mountExample2::3187528",
    "yoshida/mountExample3::97080796::360",
    "yoshida/mountExample4::1326849651::360",
    "yoshida/mountExample5::::86400",
    "yoshida/test4x4_1::665711555567792389878908993624629379187969880179721169068827951",
    "yoshida/test4x4_2::63292704423941655080293971395348848807454253204720526472462015",
    "yoshida/test4x4_3::43075357146173570492117291685601604830544643769252831337342557",
    "yoshida/tru_cube_latte::0",
    "yoshida/tru_simplex_latte::0",
    "yoshida/3x3x4_1.equ",
### Following are identical to hickerson/*
##     "yoshida/HD::::86400",
##     "yoshida/HD1",
##     "yoshida/HD10::::86400",
##     "yoshida/HD11::::86400",
##     "yoshida/HD12::::86400",
##     "yoshida/HD13::::86400",
##     "yoshida/HD14::::86400",
##     "yoshida/HD15::::86400",
##     "yoshida/HD16::::86400",
##     "yoshida/HD17::::86400",
##     "yoshida/HD18::::86400",
##     "yoshida/HD2",
##     "yoshida/HD3",
##     "yoshida/HD4",
##     "yoshida/HD5",
##     "yoshida/HD6",
##     "yoshida/HD7::::1800",
##     "yoshida/HD8::::86400",
##     "yoshida/HD9::::1800",
    "yoshida/cube_test1",
    "yoshida/cube_test1.ine",
    "yoshida/cuww1_1.equ",
    "yoshida/cuww2_1.equ",
    "yoshida/cuww3_1.equ",
    "yoshida/cuww4_1.equ",
    "yoshida/cuww5_1.equ",
    "yoshida/dean1",
    "yoshida/dean2::::86400",
    "yoshida/dean3::::86400",
    "yoshida/prob10_1.equ::::86400",
    "yoshida/prob1_1.equ::::720",
    "yoshida/prob2_1.equ::::720",
    "yoshida/prob3_1.equ::::720",
    "yoshida/prob4_1.equ::::720",
    "yoshida/prob5_1.equ::::720",
    "yoshida/prob6_1.equ::::1800",
    "yoshida/prob7_1.equ::::1800"
);

#$MAXRUNTIME = 1800;
$MAXRUNTIME = 60;

chop($LATTEDIR = `cd \`dirname $0\`; cd ../..; pwd`);
$PARAMETERS = $ARGV[0];
$COMMAND = "ulimit -t $MAXRUNTIME; ulimit -c 0; env LD_LIBRARY_PATH=/localapp/imosoft/sparc-sun-solaris2.7/alpha/gmp-4.1.4-gcc33/lib:\$LD_LIBRARY_PATH $LATTEDIR/code/latte/count $PARAMETERS";

$EXAMPLESDIR = "$LATTEDIR/EXAMPLES";
chop($DATE = `/bin/date +%Y-%m-%d`);
chop($HOSTNAME = `hostname`);
$LOGDIR = "log-$DATE-count $PARAMETERS-pid$$\@$HOSTNAME";

print "Logging to '$LOGDIR/log'\n";
mkdir($LOGDIR);
chdir($LOGDIR);

open SUMMARY, ">summary" or die;
print SUMMARY "Running with time limit ", $MAXRUNTIME, "\n";

foreach $file_num (@files_nums)
{
   # data is of the form filename::num_lattice_points
   @tmp = split('::', $file_num);
   $file_name = $tmp[0];
   $num_lattice_pts = $tmp[1];
   $time_limit = $tmp[2];

   unlink($OUTPUT_FILE_NAME);
   print "$file_name: \t";
   print SUMMARY "$file_name: \t";
   if ($time_limit && ($time_limit > $MAXRUNTIME)) {
       print "Skipped.";
       #print SUMMARY "Skipped.";
   }
   else {
       $ret_val = system( "$COMMAND $EXAMPLESDIR/$file_name >> log 2>&1 ");
       if ($ret_val != 0) {
	   print "ERROR STATUS $ret_val";
	   print SUMMARY "ERROR STATUS $ret_val";
       } else {
	   open IN, "<$OUTPUT_FILE_NAME";
	   if (!IN) {
	       print "Can't open $OUTPUT_FILE_NAME";
	       print SUMMARY "Can't open $OUTPUT_FILE_NAME";
	   }
	   else {
	       while ($line = <IN>) {
		   chomp $line;
		   # line should contain the number of lattice points
		   #print "Found numOfLatticePoints = $line\n";
		   print "Result: $line  ";
		   print SUMMARY "Result: $line  ";
		   if ($num_lattice_pts) {
		       if ($line != $num_lattice_pts) {
			   print "WRONG";
			   print SUMMARY "WRONG";
		       } else {
			   print "GOOD";
			   print SUMMARY "GOOD";
			   if (open IN, "<totalTime") {
			       $time = <IN>;
			       chomp $time;
			       print " Time: $time sec";
			       print SUMMARY " Time: $time sec";
			   }
		       }
		   }
		   else {
		       if (open IN, "<totalTime") {
			   $time = <IN>;
			   chomp $time;
			   print " Time: $time sec";
			   print SUMMARY " Time: $time sec";
		       }
		   }
	       }
	       close(IN);
	   }
       }
   }
   print "\n";
   print SUMMARY "\n";
}
