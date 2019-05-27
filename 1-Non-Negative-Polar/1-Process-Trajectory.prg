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


#Index of residues belonging to Surface Exposed protein walls
#Custom-change to your specific systems

my @sel1 = qw(3067 2988 3037 3010 2917 1283 1253 1287 151 2639 25 2638 1190 1219 1212 3470 950 628 3249 3335 526 458 461 832 3385 1765 1764 1767 1270 1260 239 229 2618 1442 1419 1077 863 1403 1396 859 3026 1706 1300 1335 515 2443 760 1667 1975 2195 3363);

my @sel2 = qw(1193 1230 1238 2992 2989 2956 225 3099 2575 1462 76 155 104 216 3046 1357 3458 3459 3234 3247 3222 3268 907 947 961 894 1319 1379 728 1675 793 791 1098 487 1423 2991 103 2569 14 1275 1133 913 911 3362 539 583 659 757 813);

my @sel3 = qw(2993 3039 115 2984 160 231 19 84 5 73 1441 1456 3196 1019 3087 1074 3140 1345 3094 3062 3112 962 902 614 720 1303 3447 3479 718 3229 608 601 3441);

my @sel4 = qw(1965 790 1113 1125 3028 1770 2187 1089 1191 1172 1170 1200 2964 2969 2996 2619 2580 1417 128 3438 1128 1955 2440 482 524 519 453 1082 1978 360 2436 1982 3130);

my @sel5 = qw(2588 166 121 108 3076 6 2615 1454 1464 1465 1418 3124 3034 2977 2974 2953 1289 1234 3427 3472 716 1322 1302 1149 1327 1135 1338 3060 900 957 3333 3334 3264 3251 1110 1095 455 724 1849 1968 1989 1962 838 751 736 3195 138 3075 3081 1030 1020 1011);

my @sel6 = qw(1060 2949 2943 1246 1217 1248 1202 1844 943 897 3337 649 1134 1152 1707 1118 1115 855 1081 1086 3047 3052 1207 1841 975 970 3324 921 3240 626 2495 533 2551 477 475 1369 1402 1391 109 2926 1284 1157 1840 1276 1373 1389 1915 1392 3253 631 643 3043);

my @sel7 = qw(2641 55 27 243 2585 3072 226 3035 3005 1266 2976 2978 2998 1195 1186 1176 3417 3424 3434 3475 707 1346 3065 3066 959 956 1014 1012 1001 1004 762 620 619 1103 1376 1108 843 525 969 968 864 892);

my @sel8 = qw(1384 1987 1298 1132 1343 3262 2494 3021 2944 2990 12 1445 877 870 862 1078 1410 135 110 3125 3107 3022 1354 1365 1332 1112 784 3457 3061 3095 621 3237 3255 3224 529 1990 1274 3106 1359 1177 1184 1247 1205 1160 1526 1267 1337 761 775 704 732 2154 1674);

my @sel9 = qw(960 3120 3110 3101 117 144 3000 3042 3015 874 56 240 3137 9 120 1121 3396 1198 1237 612 735 705 3227 3252 630 627 594 646 998 939 3221 3430 3431 1279 1292 1362 1356 1353 839 1956 1377 851);

my @sel10 = qw(789 3420 770 528 3020 3003 2613 2594 491 2491 215 219 1449 2562 3089 1252 1250 1226 1222 1326 1148 3109 1034 1013 1063 1295 1311 708 3311 3257 591 678 697 693 1719 1405 1119 1312 1313 2439 772 2591 32 680 2504 2074 3055 856 1071 358 1723 516 366);

#Concatenate Indexes in array, and divide into six groups to avoid memory crashing
#Custom-change to your specific systems
my @sel_Z=(\@sel1,\@sel2,\@sel3,\@sel4,\@sel5,\@sel6,\@sel7,\@sel8,\@sel9,\@sel10);
my @sel_cut= qw(51 49 33 33 52 52 42 53 42 53);

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
print "aMD-new3-1\n";
}
#75.75ns
if ($fold_num==2) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new6-3/stripped.nc";
print "aMD-new6-3\n";
}
#122ns
if ($fold_num==3) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new6-4/stripped.nc";
print "aMD-new6-4\n";
}
#138.5ns
if ($fold_num==4) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new6-5/stripped.nc";
print "aMD-new6-5\n";
}
#123.5ns
if ($fold_num==5) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new6-6/stripped.nc";
print "aMD-new6-6\n";
}
#31ns
if ($fold_num==6) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new7-1/stripped.nc";
print "aMD-new7-1\n";
}
#27.5ns
if ($fold_num==7) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new7-2/stripped.nc";
print "aMD-new7-2\n";
}
#28.5ns
if ($fold_num==8) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new7-3/stripped.nc";
print "aMD-new7-3\n";
}
#27ns
if ($fold_num==9) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new7-4/stripped.nc";
print "aMD-new7-4\n";
}
#101 ns
if ($fold_num==10) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-close-0/stripped.nc";
print "TL-close-0\n";
}
#100 ns
if ($fold_num==11) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-close-1/stripped.nc";
print "TL-close-1\n";
}
#13ns
if ($fold_num==12) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-1/stripped.nc";
print "TL-open1-1\n";
}
#13ns
if ($fold_num==13) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-2/stripped.nc";
print "TL-open1-2\n";
}
#18ns
if ($fold_num==14) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-4/stripped.nc";
print "TL-open1-4\n";
}
#28.5ns
if ($fold_num==15) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-5/stripped.nc";
print "TL-open1-5\n";
}
#28ns
if ($fold_num==16) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-6/stripped.nc";
print "TL-open1-6\n";
}
#16ns
if ($fold_num==17) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-1/stripped.nc";
print "TL-open2-1\n";
}
#15.5ns
if ($fold_num==18) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-2/stripped.nc";
print "TL-open2-2\n";
}
#107.5ns 
if ($fold_num==19) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-4/stripped.nc";
print "TL-open2-4\n";
}
#31.5ns
if ($fold_num==20) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-5/stripped.nc";
print "TL-open2-5\n";
}
#31.5ns
if ($fold_num==21) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-6/stripped.nc";
print "TL-open2-6\n";
}
#16ns
if ($fold_num==22) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-1/stripped.nc";
print "TL-open3-1\n";
}
#71ns
if ($fold_num==23) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-2/stripped.nc";
print "TL-open3-2\n";
}
#19.5ns
if ($fold_num==24) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-3/stripped.nc";
print "TL-open3-3\n";
}
#38.5ns
if ($fold_num==25) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-5/stripped.nc";
print "TL-open3-5\n";
}
#27ns
if ($fold_num==26) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-6/stripped.nc";
print "TL-open3-6\n";
}
#67 ns
if ($fold_num==27) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-1/stripped.nc";
print "iso-CH3-1\n";
}
#57.5 ns
if ($fold_num==28) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-2/stripped.nc";
print "iso-CH3-2\n";
}
#75.5 ns
if ($fold_num==29) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-3/stripped.nc";
print "iso-CH3-3\n";
}
#77 ns
if ($fold_num==30) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-4/stripped.nc";
print "iso-CH3-4\n";
}
#124.5 ns
if ($fold_num==31) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-5/stripped.nc";
print "iso-CH3-5\n";
}
#124 ns
if ($fold_num==32) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-6/stripped.nc";
print "iso-CH3-6\n";
}
#102.5 ns
if ($fold_num==33) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-1/stripped.nc";
print "trans-CH3-1\n";
}
#74 ns
if ($fold_num==34) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-2/stripped.nc";
print "trans-CH3-2\n";
}
#52.5 ns
if ($fold_num==35) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-3/stripped.nc";
print "trans-CH3-3\n";
}
#189 ns
if ($fold_num==36) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-5/stripped.nc";
print "trans-CH3-5\n";
}
#162 ns
if ($fold_num==37) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-6/stripped.nc";
print "trans-CH3-6\n";
}
#74.5 ns
if ($fold_num==38) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-3-2/stripped.nc";
print "trans-CH3-3-2\n";
}
#54.5 ns
if ($fold_num==39) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-1/stripped.nc";
print "trans-CH3-m4-1\n";
}
#54.5 ns
if ($fold_num==40) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-2/stripped.nc";
print "trans-CH3-m4-2\n";
}
#54.5 ns
if ($fold_num==41) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-3/stripped.nc";
print "trans-CH3-m4-3\n";
}
#55.5 ns
if ($fold_num==42) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-5/stripped.nc";
print "trans-CH3-m4-5\n";
}
#58.5 ns
if ($fold_num==43) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-6/stripped.nc";
print "trans-CH3-m4-6\n";
}
#55.5 ns
if ($fold_num==44) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-7/stripped.nc";
print "trans-CH3-m4-7\n";
}
#22.2 ns
if ($fold_num==45) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new-1/stripped.nc";
print "trans-new-1\n";
}
#49.8 ns
if ($fold_num==46) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new-2/stripped.nc";
print "trans-new-2\n";
}
#53.5 ns
if ($fold_num==47) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-1/stripped.nc";
print "trans-new5-1\n";
}
#101.8 ns
if ($fold_num==48) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-2/stripped.nc";
print "trans-new5-2\n";
}
#40.4 ns
if ($fold_num==49) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-3/stripped.nc";
print "trans-new5-3\n";
}
#95.2 ns
if ($fold_num==50) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-5/stripped.nc";
print "trans-new5-5\n";
}
#169.6 ns
if ($fold_num==51) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-6/stripped.nc";
print "trans-new5-6\n";
}
#89 ns
if ($fold_num==52) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new7-3/stripped.nc";
print "trans-new7-3\n";
}
#156.6 ns
if ($fold_num==53) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new7-5/stripped.nc";
print "trans-new7-5\n";
}
#75 ns
if ($fold_num==54) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1/stripped.nc";
print "trans-new8-1\n";
}
#40.8 ns
if ($fold_num==55) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1-3/stripped.nc";
print "trans-new8-1-3\n";
}
#40.2 ns
if ($fold_num==56) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1-4/stripped.nc";
print "trans-new8-1-4\n";
}
#12.2 ns
if ($fold_num==57) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1-5/stripped.nc";
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
print "trans-new8-1-7\n";
}
#14.2 ns
if ($fold_num==60) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-2/stripped.nc";
print "trans-new8-2\n";
}
#94.2 ns
if ($fold_num==61) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-6/stripped.nc";
print "trans-new8-6\n";
}
#56ns
if ($fold_num==62) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/stripped.nc";
print "aMD-new-4-alt-rep3\n";
}
#93.5ns
if ($fold_num==63) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep5/stripped.nc";
print "aMD-new-4-alt-rep5\n";
}
#59ns
if ($fold_num==64) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep1/stripped.nc";
print "aMD-new-4-alt3-rep1\n";
}
#67ns
if ($fold_num==65) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep2/stripped.nc";
print "aMD-new-4-alt3-rep2\n";
}
#61.5ns
if ($fold_num==66) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep3/stripped.nc";
print "aMD-new-4-alt3-rep3\n";
}
#62ns
if ($fold_num==67) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep4/stripped.nc";
print "aMD-new-4-alt3-rep4\n";
}
#118.5ns
if ($fold_num==68) {
$prmtop= "out15-parmed2.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep5/stripped.nc";
print "aMD-new-4-alt3-rep5\n";
}
#101ns
if ($fold_num==69) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep6/stripped.nc";
print "aMD-new-4-alt3-rep6\n";
}
#100ns
if ($fold_num==70) {
$prmtop= "out15-parmed.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/syner+3-5/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/syner+4-6/stripped.nc";
print "syner+4-6\n";
}
#100ns
if ($fold_num==71) {
$prmtop= "out15-parmed.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/syner+3-5/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/syner+3-5/stripped.nc";
print "syner+3-5\n";
}
#100ns
if ($fold_num==72) {
$prmtop= "out15-parmed.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/syner+3-5/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/syner+2-new/stripped.nc";
print "syner+2-new\n";
}
#250ns
if ($fold_num==73) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-2/stripped.nc";
print "aMD-new3-2\n";
}
#252ns
if ($fold_num==74) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-3/stripped.nc";
print "aMD-new3-3\n";
}
#251ns
if ($fold_num==75) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-4/stripped.nc";
print "aMD-new3-4\n";
}
#368.5ns
if ($fold_num==76) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-5/stripped.nc";
print "aMD-new3-5\n";
}
#365ns
if ($fold_num==77) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-6/stripped.nc";
print "aMD-new3-6\n";
}

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


for (my $divide=0; $divide<$divide_trigger; $divide++) {

for (my $cut=1; $cut<11; $cut++) {
$y=$sel_cut[$cut-1];

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
print (FILE2 "set sel1 [atomselect top \"(resname ctp) and (within 3 of resid $sel)\" frame \$frame]\n");
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
############################################# END  CUT loop #########################################
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

#LABEL:


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


for (my $cut=1; $cut<11; $cut++) {
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



