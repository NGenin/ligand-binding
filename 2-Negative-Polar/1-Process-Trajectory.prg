use File::Slurp;
use strict;
use autodie;
use warnings qw(all);
use File::Copy qw(copy);
use Math::Complex;
use Cwd;

# clear the screen
system 'clear';
print "\n";

#Variable Declaration
my $v;
my @abs_filename;
my $prmtop;
my $dir = getcwd;
my $abs_filename_prmtop;
my $nb_traj;
my $ch_arg=$ARGV[1];
my $sel;
my $y;
my $z;
my $range_frames;
my $div;
my $divide_trigger=1;
my $start_1;
my $last_1;
my $start_2;
my $last_2;
my $start_3;
my $last_3;
my $sel_ctp;
my $sel_MG;
my $type_MG;

#Index of Asparagine and Glutamine residues belonging to Surface Exposed protein walls
#Custom-change to your specific systems
my @sel_asp1 = qw(152 163 165 1440 2986 2982 2947 2950 2939 953 932 928 3340 3250 3223 3220 3266 3226 3228 3263);

my @sel_asp2 = qw(543 2484 1240 3391 1271 1273 1169 1180 3366 3411 3437 1850 1851 1853 1344 2983 3202 1352 3093 1959);

my @sel_asp3 = qw(1272 885 805 1 1057 1428 1773 804 2073 2830 1761 1394 531 361 2400 2524 2523 2559 2548 2182);

my @sel_glu1 = qw(709 3364 2948 3057 3063 942 946 3319 3243 3245 3114 1355 1340 3393 3390 150 149 161 223 214);

my @sel_glu2 = qw(237 230 71 3201 1064 927 1372 1370 884 674 1763 2153 2493 1848 1447 80 1432 1425 1424 312);

my @sel_glu3 = qw(1280 1206 3050 1363 1766 860 2151 2843 2842 1774 1733 1687 1665 2557 2555 2138 893 650);


#Concatenate Indexes in array, and divide into six groups to avoid memory crashing
#Custom-change to your specific systems
my @sel_Z=(\@sel_asp1,\@sel_asp2,\@sel_asp3,\@sel_glu1,\@sel_glu2,\@sel_glu3);
my @sel_cut= qw(20 20 20 20 20 18);

#Simulation frame lengths
#Custom-change to your specific systems
my @sel_frames= qw(511 299 488 554 494 124 110 115 108 202 200 26 26 36 57 56 32 31 215 63 63 32 142 39 77 54 134 115 151 154 249 248 205 148 105 378 324 149 109 109 109 111 117 111 111 249 266 509 202 476 848 445 783 422 204 201 61 344 62 72 472 112 187 118 134 123 124 237 202 200 200 200 500 504 502 737 730);



########################################################################################################
################################################  FILE LOOP ############################################
########################################################################################################

my $datestring_start = localtime();
my $outfile="output-global-polar.txt";
#OUTPUT
open (FILE, "> $outfile") || die;

#Protein trajectory files
#Custom-change to your specific systems
for (my $fold_num=1; $fold_num<78; $fold_num++) {
print "SIM ID=\n";

#251.5ns
if ($fold_num==1) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-1/stripped.nc";
$type_MG=1;
print "aMD-new3-1\n";
}
#75.75ns
if ($fold_num==2) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new6-3/stripped.nc";
$type_MG=1;
print "aMD-new6-3\n";
}
#122ns
if ($fold_num==3) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new6-4/stripped.nc";
$type_MG=1;
print "aMD-new6-4\n";
}
#138.5ns
if ($fold_num==4) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new6-5/stripped.nc";
$type_MG=1;
print "aMD-new6-5\n";
}
#123.5ns
if ($fold_num==5) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new6-6/stripped.nc";
$type_MG=1;
print "aMD-new6-6\n";
}
#31ns
if ($fold_num==6) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new7-1/stripped.nc";
$type_MG=1;
print "aMD-new7-1\n";
}
#27.5ns
if ($fold_num==7) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new7-2/stripped.nc";
$type_MG=1;
print "aMD-new7-2\n";
}
#28.5ns
if ($fold_num==8) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new7-3/stripped.nc";
$type_MG=1;
print "aMD-new7-3\n";
}
#27ns
if ($fold_num==9) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new7-4/stripped.nc";
$type_MG=1;
print "aMD-new7-4\n";
}
#101 ns
if ($fold_num==10) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-close-0/stripped.nc";
$type_MG=1;
print "TL-close-0\n";
}
#100 ns
if ($fold_num==11) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-close-1/stripped.nc";
$type_MG=1;
print "TL-close-1\n";
}
#13ns
if ($fold_num==12) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-1/stripped.nc";
$type_MG=1;
print "TL-open1-1\n";
}
#13ns
if ($fold_num==13) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-2/stripped.nc";
$type_MG=1;
print "TL-open1-2\n";
}
#18ns
if ($fold_num==14) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-4/stripped.nc";
$type_MG=1;
print "TL-open1-4\n";
}
#28.5ns
if ($fold_num==15) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-5/stripped.nc";
$type_MG=1;
print "TL-open1-5\n";
}
#28ns
if ($fold_num==16) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-6/stripped.nc";
$type_MG=1;
print "TL-open1-6\n";
}
#16ns
if ($fold_num==17) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-1/stripped.nc";
$type_MG=1;
print "TL-open2-1\n";
}
#15.5ns
if ($fold_num==18) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-2/stripped.nc";
$type_MG=1;
print "TL-open2-2\n";
}
#107.5ns 
if ($fold_num==19) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-4/stripped.nc";
$type_MG=1;
print "TL-open2-4\n";
}
#31.5ns
if ($fold_num==20) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-5/stripped.nc";
$type_MG=1;
print "TL-open2-5\n";
}
#31.5ns
if ($fold_num==21) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-6/stripped.nc";
$type_MG=1;
print "TL-open2-6\n";
}
#16ns
if ($fold_num==22) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-1/stripped.nc";
$type_MG=1;
print "TL-open3-1\n";
}
#71ns
if ($fold_num==23) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-2/stripped.nc";
$type_MG=1;
print "TL-open3-2\n";
}
#19.5ns
if ($fold_num==24) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-3/stripped.nc";
$type_MG=1;
print "TL-open3-3\n";
}
#38.5ns
if ($fold_num==25) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-5/stripped.nc";
$type_MG=1;
print "TL-open3-5\n";
}
#27ns
if ($fold_num==26) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-6/stripped.nc";
$type_MG=1;
print "TL-open3-6\n";
}
#67 ns
if ($fold_num==27) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-1/stripped.nc";
$type_MG=2;
print "iso-CH3-1\n";
}
#57.5 ns
if ($fold_num==28) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-2/stripped.nc";
$type_MG=2;
print "iso-CH3-2\n";
}
#75.5 ns
if ($fold_num==29) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-3/stripped.nc";
$type_MG=2;
print "iso-CH3-3\n";
}
#77 ns
if ($fold_num==30) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-4/stripped.nc";
$type_MG=2;
print "iso-CH3-4\n";
}
#124.5 ns
if ($fold_num==31) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-5/stripped.nc";
$type_MG=2;
print "iso-CH3-5\n";
}
#124 ns
if ($fold_num==32) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-6/stripped.nc";
$type_MG=2;
print "iso-CH3-6\n";
}
#102.5 ns
if ($fold_num==33) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-1/stripped.nc";
$type_MG=2;
print "trans-CH3-1\n";
}
#74 ns
if ($fold_num==34) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-2/stripped.nc";
$type_MG=2;
print "trans-CH3-2\n";
}
#52.5 ns
if ($fold_num==35) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-3/stripped.nc";
$type_MG=2;
print "trans-CH3-3\n";
}
#189 ns
if ($fold_num==36) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-5/stripped.nc";
$type_MG=2;
print "trans-CH3-5\n";
}
#162 ns
if ($fold_num==37) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-6/stripped.nc";
$type_MG=2;
print "trans-CH3-6\n";
}
#74.5 ns
if ($fold_num==38) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-3-2/stripped.nc";
$type_MG=2;
print "trans-CH3-3-2\n";
}
#54.5 ns
if ($fold_num==39) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-1/stripped.nc";
$type_MG=2;
print "trans-CH3-m4-1\n";
}
#54.5 ns
if ($fold_num==40) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-2/stripped.nc";
$type_MG=2;
print "trans-CH3-m4-2\n";
}
#54.5 ns
if ($fold_num==41) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-3/stripped.nc";
$type_MG=2;
print "trans-CH3-m4-3\n";
}
#55.5 ns
if ($fold_num==42) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-5/stripped.nc";
$type_MG=2;
print "trans-CH3-m4-5\n";
}
#58.5 ns
if ($fold_num==43) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-6/stripped.nc";
$type_MG=2;
print "trans-CH3-m4-6\n";
}
#55.5 ns
if ($fold_num==44) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-7/stripped.nc";
$type_MG=2;
print "trans-CH3-m4-7\n";
}
#22.2 ns
if ($fold_num==45) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new-1/stripped.nc";
$type_MG=2;
print "trans-new-1\n";
}
#49.8 ns
if ($fold_num==46) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new-2/stripped.nc";
$type_MG=2;
print "trans-new-2\n";
}
#53.5 ns
if ($fold_num==47) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-1/stripped.nc";
$type_MG=2;
print "trans-new5-1\n";
}
#101.8 ns
if ($fold_num==48) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-2/stripped.nc";
$type_MG=2;
print "trans-new5-2\n";
}
#40.4 ns
if ($fold_num==49) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-3/stripped.nc";
$type_MG=2;
print "trans-new5-3\n";
}
#95.2 ns
if ($fold_num==50) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-5/stripped.nc";
$type_MG=2;
print "trans-new5-5\n";
}
#169.6 ns
if ($fold_num==51) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-6/stripped.nc";
$type_MG=2;
print "trans-new5-6\n";
}
#89 ns
if ($fold_num==52) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new7-3/stripped.nc";
$type_MG=2;
print "trans-new7-3\n";
}
#156.6 ns
if ($fold_num==53) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new7-5/stripped.nc";
$type_MG=2;
print "trans-new7-5\n";
}
#75 ns
if ($fold_num==54) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1/stripped.nc";
$type_MG=2;
print "trans-new8-1\n";
}
#40.8 ns
if ($fold_num==55) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1-3/stripped.nc";
$type_MG=2;
print "trans-new8-1-3\n";
}
#40.2 ns
if ($fold_num==56) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1-4/stripped.nc";
$type_MG=2;
print "trans-new8-1-4\n";
}
#12.2 ns
if ($fold_num==57) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1-5/stripped.nc";
$type_MG=2;
print "trans-new8-1-5\n";
}
#68.8 ns
if ($fold_num==58) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-3/stripped.nc";
print "trans-new8-3\n";
}
#12.4 ns
if ($fold_num==59) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1-7/stripped.nc";
$type_MG=2;
print "trans-new8-1-7\n";
}
#14.2 ns
if ($fold_num==60) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-2/stripped.nc";
$type_MG=2;
print "trans-new8-2\n";
}
#94.2 ns
if ($fold_num==61) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-6/stripped.nc";
$type_MG=2;
print "trans-new8-6\n";
}
#56ns
if ($fold_num==62) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/stripped.nc";
$type_MG=3;
print "aMD-new-4-alt-rep3\n";
}
#93.5ns
if ($fold_num==63) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep5/stripped.nc";
$type_MG=3;
print "aMD-new-4-alt-rep5\n";
}
#59ns
if ($fold_num==64) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep1/stripped.nc";
$type_MG=3;
print "aMD-new-4-alt3-rep1\n";
}
#67ns
if ($fold_num==65) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep2/stripped.nc";
$type_MG=3;
print "aMD-new-4-alt3-rep2\n";
}
#61.5ns
if ($fold_num==66) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep3/stripped.nc";
$type_MG=3;
print "aMD-new-4-alt3-rep3\n";
}
#62ns
if ($fold_num==67) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep4/stripped.nc";
$type_MG=3;
print "aMD-new-4-alt3-rep4\n";
}
#118.5ns
if ($fold_num==68) {
$prmtop= "out15-parmed2.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep5/stripped.nc";
$type_MG=3;
print "aMD-new-4-alt3-rep5\n";
}
#101ns
if ($fold_num==69) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep6/stripped.nc";
$type_MG=3;
print "aMD-new-4-alt3-rep6\n";
}
#100ns
if ($fold_num==70) {
$prmtop= "out15-parmed.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/syner+3-5/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/syner+4-6/stripped.nc";
$type_MG=1;
print "syner+4-6\n";
}
#100ns
if ($fold_num==71) {
$prmtop= "out15-parmed.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/syner+3-5/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/syner+3-5/stripped.nc";
$type_MG=1;
print "syner+3-5\n";
}
#100ns
if ($fold_num==72) {
$prmtop= "out15-parmed.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/syner+3-5/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/syner+2-new/stripped.nc";
$type_MG=1;
print "syner+2-new\n";
}
#250ns
if ($fold_num==73) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-2/stripped.nc";
$type_MG=1;
print "aMD-new3-2\n";
}
#252ns
if ($fold_num==74) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-3/stripped.nc";
$type_MG=1;
print "aMD-new3-3\n";
}
#251ns
if ($fold_num==75) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-4/stripped.nc";
$type_MG=1;
print "aMD-new3-4\n";
}
#368.5ns
if ($fold_num==76) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-5/stripped.nc";
$type_MG=1;
print "aMD-new3-5\n";
}
#365ns
if ($fold_num==77) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-6/stripped.nc";
$type_MG=1;
print "aMD-new3-6\n";
}

#Residue Index of the non lingand-bound Mg2+ atoms to be excluded from the calculations
if ($type_MG==1){
$sel_MG='4308 4309 4297';
}
if ($type_MG==2){
$sel_MG='4307 4308 4298';
}
if ($type_MG==3){
$sel_MG='4306 4307 4297';
}

$sel_ctp='resname ctp';

#Divide Protein Trajectories into several parts in order to avoid memory crashing
$range_frames=$sel_frames[$fold_num-1];
if ($range_frames>400) {
$div=int($range_frames/2);
$divide_trigger=2;
$start_1=0;
$last_1=$div;
$start_2=$div+1;
$last_2=$range_frames-1;
}
if ($fold_num==77) {
$divide_trigger=3;
$start_1=0;
$last_1=250;
$start_2=251;
$last_2=500;
$start_3=501;
$last_3=729;
}


for (my $divide=0; $divide<$divide_trigger; $divide++) {

for (my $cut=1; $cut<7; $cut++) {
$y=$sel_cut[$cut-1];

#Deal with special case
if (($fold_num==77) and ($divide_trigger>1)){
$sel_ctp='resname ctp and not resid 4325';
$sel_MG=$sel_MG . ' 4324';
}
#END Deal with special case

################################################################################################## 
#############################################  CUT   #############################################
##################################################################################################

my $outfile="scr-stat.vmd";
#OUTPUT
open (FILE2, "> $outfile") || die;
print (FILE2 "set outFile out-stat.dat\n");
print (FILE2 "set out [open \$outFile w]\n");
#print dummy value to prevent void file error (Died)
print (FILE2 "puts \$out \"0 0\"\n");
#TRAJECTORY
print (FILE2 "mol new $abs_filename_prmtop type parm7 first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n");
$v=$abs_filename[0];
if ($divide_trigger==1){
print (FILE2 "mol addfile $v type netcdf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n");
}
if ($divide_trigger==2){
if ($divide==0){
print (FILE2 "mol addfile $v type netcdf first $start_1 last $last_1 step 1 filebonds 1 autobonds 1 waitfor all\n");
}
if ($divide==1){
print (FILE2 "mol addfile $v type netcdf first $start_2 last $last_2 step 1 filebonds 1 autobonds 1 waitfor all\n");
}
}
if ($divide_trigger==3){
if ($divide==0){
print (FILE2 "mol addfile $v type netcdf first $start_1 last $last_1 step 1 filebonds 1 autobonds 1 waitfor all\n");
}
if ($divide==1){
print (FILE2 "mol addfile $v type netcdf first $start_2 last $last_2 step 1 filebonds 1 autobonds 1 waitfor all\n");
}
if ($divide==2){
print (FILE2 "mol addfile $v type netcdf first $start_3 last $last_3 step 1 filebonds 1 autobonds 1 waitfor all\n");
}
}
#ALIGN
print (FILE2 "set startFrame 0\n");
#Number of frames
print (FILE2 "set nFrames [molinfo top get numframes]\n");
if ($cut==1){
print (FILE2 "puts \$out \"5000 \$nFrames\"\n");
}
#use frame 0 for the reference
print (FILE2 "set reference [atomselect top \"protein\" frame 0]\n");
print (FILE2 "for {set frame 1} {\$frame < \$nFrames} {incr frame} {\n");
#get the correct frame
print (FILE2 "set compare [atomselect top \"protein\" frame \$frame]\n");
#compute the transformation
print (FILE2 "set trans_mat [measure fit \$compare \$reference]\n");
#move sel
print (FILE2 "set move_sel [atomselect top \"all\" frame \$frame]\n");
#do the alignment
print (FILE2 "\$move_sel move \$trans_mat\n");
print (FILE2 "}\n");
#PRECISION
#print (FILE2 "set tcl_precision 12\n");


#STAT
#####################################  FRAME LOOP #########################################
print (FILE2 "for {set frame 0} {\$frame < \$nFrames} {incr frame} {\n");
for (my $x=0; $x<$y; $x++) {
$sel=$sel_Z[$cut-1][$x];
print (FILE2 "set sel $sel\n");
#Deal with Asparagines
if ($cut<4){
print (FILE2 "set sel1 [atomselect top \"(($sel_ctp) and (within 3 of resid $sel)) or ((resname MG and not resid $sel_MG) and (within 4.5 of (resid $sel and name OD1 OD2)))\" frame \$frame]\n");
}
#Deal with Glutamines
if ($cut>=4){
print (FILE2 "set sel1 [atomselect top \"(($sel_ctp) and (within 3 of resid $sel)) or ((resname MG and not resid $sel_MG) and (within 4.5 of (resid $sel and name OE1 OE2)))\" frame \$frame]\n");
}
#Sort
print (FILE2 "set sel2 [\$sel1 get resid]\n");
print (FILE2 "set sel3 [lsort -unique \$sel2]\n");
print (FILE2 "set sel4 [llength \$sel3]\n");
print (FILE2 "if {\$sel4 != \"0\"} {\n");
print (FILE2 "puts \$out \"\$sel\"\n");
print (FILE2 "}\n");
}
print (FILE2 "}\n");
print (FILE2 "\n");
#CLOSE OUTPUT
print (FILE2 "close \$out\n");
#exit VMD
print (FILE2 "exit\n");
close (FILE2);
my $cmd = "vmd -dispdev text -nt -e scr-stat.vmd";
system($cmd);
#unlink "scr-stat.vmd";


#PRINT COMBINED INFO
my @pdb_input = read_file("out-stat.dat") or die;
my $array_size=scalar @pdb_input;
my $line;
for (my $count=0; $count<$array_size; $count++) {
        $line= @pdb_input[$count]; 
	print (FILE"$line");
}

}
################################################################################################## 
############################################# END  CUT loop ######################################
##################################################################################################

}
#end divide
$divide_trigger=1;




print (FILE"\n");



}
########################################################################################################
################################################ END FILE LOOP #########################################
########################################################################################################
close (FILE);


#Compute Statistics:

#GLOBAL STAT 
my $outfile="output-global-polar-final.txt";
open (FILE3, "> $outfile") || die;
my @pdb_input = read_file("output-global-polar.txt") or die;
my $array_size=scalar @pdb_input;
my @line_handle;
my $resid;
my $zone;
my $stat_sum;
my $stat_frames_sum=0;
my $stat_frames;

print (FILE3"STAT GLOBAL\n");

for (my $count=0; $count<$array_size; $count++) {
@line_handle = split ( /\s+/, @pdb_input[$count] );
$zone= @line_handle[0];
if ($zone==5000) {
$stat_frames= @line_handle[1];
$stat_frames_sum=$stat_frames_sum+$stat_frames;
}
}
print (FILE3"nb_frames $stat_frames_sum\n");
print (FILE3"\n");


for (my $cut=1; $cut<7; $cut++) {
$y=$sel_cut[$cut-1];

for (my $x=0; $x<$y; $x++) {
$stat_sum=0;
$sel=$sel_Z[$cut-1][$x];

for (my $count=0; $count<$array_size; $count++) {
@line_handle = split ( /\s+/, @pdb_input[$count] );
$resid= @line_handle[0];
if ($resid==$sel){
$stat_sum=$stat_sum+1;
}
}

$stat_sum=($stat_sum/$stat_frames_sum)*100;
$stat_sum=sprintf "%.2f", $stat_sum;
print (FILE3"$sel $stat_sum\n");
}

print (FILE3"\n");
}

close (FILE3);

my $datestring_end = localtime();
print "START = $datestring_start\n";
print "END = $datestring_end\n";

exit;



