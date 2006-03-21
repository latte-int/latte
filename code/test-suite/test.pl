#!/usr/bin/perl -w

# count always prints the number of lattice points in the file
# numOfLatticePoints. In order to check the output of the test,
# we must open the file and verify the number of lattice points
# stored there.
$OUTPUT_FILE_NAME = "numOfLatticePoints";

# add file names and lattice points to an array
my @files_nums = ("example1::3",
   "example2::2",
   "example3::2",
   "24_cell::33",
   "hickerson/hickerson-7::8",
   "hickerson/hickerson-8::22",
   "hickerson/hickerson-9::9",
   "hickerson/hickerson-10::24",
   "hickerson/hickerson-11::12",
   "hickerson/hickerson-12::38",
   "hickerson/hickerson-13::14",
   "hickerson/hickerson-14::32",
   "hickerson/hickerson-15::20",
   "hickerson/hickerson-16::54",
   "hickerson/hickerson-17::18",
   "hickerson/hickerson-18::44",
   "hickerson/hickerson-19::20",
   "hickerson/hickerson-20::74",
   "hickerson/hickerson-24::96",
   "hickerson/hickerson-28::92",
   "hickerson/hickerson-32::122",
   "hickerson/hickerson-36::138");

foreach $file_num (@files_nums)
{
   # data is of the form filename::num_lattice_points
   @tmp = split('::', $file_num);
   $file_name = $tmp[0];
   $num_lattice_pts = $tmp[1];

   $ret_val = system("./../code/latte/count --maxdet=10 --exp $file_name");
   if ($ret_val != 0) {
      print "count aborted with error on $file_name\n";
   } else {
      open IN, "<$OUTPUT_FILE_NAME" or die "Can't open $OUTPUT_FILE_NAME : $!\n";
      while ($line = <IN>) {
         chomp $line;
         # line should contain the number of lattice points
         print "Found numOfLatticePoints = $line\n";
         if ($line != $num_lattice_pts) {
            print "count found incorrect number of lattice points on $file_name\n";
         } else {
            print "count sucessfully tested example: $file_name\n";
         }
      }
      close(IN);
   }
}
