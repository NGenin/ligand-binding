use File::Slurp;
use strict;
use autodie;
use warnings qw(all);
use File::Copy qw(copy);
use Math::Complex;
use Cwd;
use List::MoreUtils qw( minmax );
use List::MoreUtils qw( uniq );


# clear the screen
system 'clear';
print "\n";

#Variable Declaration
my $v;
my @abs_filename;
my $prmtop;
my $abs_filename_prmtop;
my $sel_ctp;
my $sel_res;
my $range_frames;
my $div;
my $divide_trigger=1;
my $start_1;
my $last_1;
my $start_2;
my $last_2;
my $start_3;
my $last_3;
my $start_4;
my $last_4;
my $start_5;
my $last_5;
my $start_6;
my $last_6;
my $cat;

#Index of residues belonging to Surface Exposed protein walls
#Custom-change to your specific systems
my $sel1_res = '749 688 699 740 734 698 742 695 691 750 701 692 748 687 702 745 686 694 738 689 746 731 737 690 706 741 703 710';
my $sel2_res = '595 3232 541 3310 538 2507 2506 3233 596 618 3235 754 3236 3315 3312 3314 3239 753 752 537 3230 3309 3336 3316 3327 593 542 3323 592 682 3225 758 616 3267 3248 759 3317 617 755 3318 3322 3343 3338 622 683 597 610 611 609 613 602 600 605 606 603 599 604';
my $sel3_res = '910 909 922 1043 901 905 919 958 912 951 903 1044 914 904 963 954 896 916 915 945 967 917 908 923 895 1048 955 1047 918 949 941 948 926 931 944 952 899 965 966';
my $sel4_res = '1371 1347 1367 1323 3059 1324 1310 1308 1309 1307 1305 1306 1320 1315 1316 1314 1318';
my $sel5_res = '1371 699 740 734 742 701 702 1347 710 738 1323 731 737 706 741 703 1324 1310 1308 1309 1307 1305 1306 1320 1315 1316 1314 1318';
my $sel6_res = '3454 3484 3444 3487 3455 3480 3450 3453 3483 3442 3436 3445 3481 3473 3485 3451 3426 3486 3449 3476 3429 3419 3448 3422 3474 3423 3446 3435 3428 3421 3477 3433 3425 714 713 712 3471 3465 3462 3467 3463 3468 3460 3461';
my $sel7_res = '1210 1209 1211 1244 1245 1241 1214 1215 1179 1213 1242 1208 1146 1182';
my $sel8_res = '1144 1197 1224 1221 1259 1225 1341 1255 1258 1218 1220 1151 1290 1254 1269 1256 1262 1263 1233 1265 1232 1147 1261 1227 1268 1257 1194 1192 1264 1342 1188 1235';
my $sel9_res = '3001 3004 3008 3009 2971 2918 2968 3012 3011 2994 2995 2963 2920 2921 2925 2919 2972 3002 3007 2970 2975 2957 2966 2965 3006 2962 2959 2973 2954 2955 3013 2960 2952 2924 2922';
my $sel10_res = '3103 130 123 3088 3126 1399 127 125 3102 1414 124 116 3040 146 118 3036 119 142 3041 131 3017 113 3105 111 114 3016 2985';
my $sel11_res = '7 10 3082 11 1443 3194 8 3192 3078 3084 3080 2617 875 1463 2640 3200 2571 1461 3197 3085 3069 865 3199 1453 13 866 1439 876 3149 873 867 3083 2614 871 1448 3074 3073 1438 3068 3070';
my $sel12_res = '24 2608 232 21 2637 2607 2589 233 2609 2605 2610 20 2598 2612 238 18 23 2587 2611 2603 2604 235 2597 217 29 2600 26 218 2601 221 2602 28 2596 2599 222';
my $sel13_res = '634 639 635 523 642 633 769 632 768 774 638 1380 1383 1378 1374 1382 766';
my $sel14_res = '726 730 725 722 729 811 723 733 1328';
my $sel15_res = '1159 1153';
my $sel16_res = '1336 2997 1286 1285 1715 3024 1716 1282 1709 1712 1294 1281 3023 1154 1717 3049 1293 1288 1278 1721 3027';
my $sel17_res = '1401 1412 3090 134 1413 3032 3031 133 1358 129 137 3045 3044 3097 140';
my $sel18_res = '1022 1005 1009 1023 3108 1015 1002 1017 1010 1026 3122 3071 1006 1008 1033 1018 1029 1021 1027';

#Index of Ligands
my $sel_ctp_1='4311 4313 4315 4317 4319 4321 4323 4325 4327 4329 4331 4333 4335 4337 4339 4341 4343 4345 4347 4349 4351 4353 4355 4357 4359 4361 4363 4365 4367 4369 4371 4373 4375 4377 4379 4381 4383 4385 4387 4389';
my $sel_ctp_2='4311 4313 4315 4317 4319 4321 4323 4327 4329 4331 4333 4335 4337 4339 4341 4343 4345 4347 4349 4351 4353 4355 4357 4359 4361 4363 4365 4367 4369 4371 4373 4375 4377 4379 4381 4383 4385 4387 4389';
my $sel_ctp_3='4311 4313 4315 4317 4319 4321 4323 4327 4329 4331 4333 4335 4337 4339 4341 4343 4345 4347 4349 4351 4353 4355 4357 4359 4361 4363 4365 4367 4369 4371 4373 4375 4379 4381 4383 4385 4387 4389';
my $sel_ctp_4='4313 4315 4317 4319 4321 4323 4325 4327 4329 4331 4333 4335 4337 4339 4341 4343 4345 4347 4349 4351 4353 4355 4357 4359 4361 4363 4365 4367 4369 4371 4373 4375 4377 4379 4381 4383 4385 4387 4389';
my $sel_ctp_5='4312 4314 4316 4318 4320 4322 4324 4326 4328 4330 4332 4334 4336 4338 4340 4342 4344 4346 4348 4350 4352 4354 4356 4358 4360 4362 4364 4366 4368 4370 4372 4374 4376 4378 4380 4382 4384 4386 4388';
my $sel_ctp_6='4311 4313 4315 4317 4319 4321 4323 4325 4327 4329 4331 4333 4335 4337 4339 4341 4343 4345 4347 4349 4351 4353 4355 4357 4359 4361 4363 4365 4367 4369 4371 4373 4375 4377 4379 4381 4383 4385 4387';
my $sel_ctp_7='4317 4319 4321 4323 4325 4327 4329 4331 4333 4335 4337 4339 4341 4343 4345 4347 4349 4351 4353 4355 4357 4359 4361 4363 4365 4367 4369 4371 4373 4375 4377 4379 4381 4383 4385 4387 4389';
my $sel_ctp_8='4315 4317 4319 4321 4323 4325 4327 4329 4331 4333 4335 4337 4339 4341 4343 4345 4347 4349 4351 4353 4355 4357 4359 4361 4363 4365 4367 4369 4371 4373 4375 4377 4379 4381 4383 4385 4387 4389';
my $sel_ctp_9='4313 4315 4317 4319 4321 4323 4325 4327 4329 4331 4333 4335 4337 4339 4341 4343 4345 4347 4349 4351 4353 4355 4357 4359 4361 4363 4365 4367 4369 4371 4373 4375 4377 4379 4381 4383 4385 4387 4389';

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
$cat=1;
print "aMD-new3-1\n";
}
#75.75ns
if ($fold_num==2) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new6-3/stripped.nc";
$cat=2;
print "aMD-new6-3\n";
}
#122ns
if ($fold_num==3) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new6-4/stripped.nc";
$cat=2;
print "aMD-new6-4\n";
}
#138.5ns
if ($fold_num==4) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new6-5/stripped.nc";
$cat=2;
print "aMD-new6-5\n";
}
#123.5ns
if ($fold_num==5) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new6-6/stripped.nc";
$cat=2;
print "aMD-new6-6\n";
}
#31ns
if ($fold_num==6) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new7-1/stripped.nc";
$cat=3;
print "aMD-new7-1\n";
}
#27.5ns
if ($fold_num==7) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new7-2/stripped.nc";
$cat=3;
print "aMD-new7-2\n";
}
#28.5ns
if ($fold_num==8) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new7-3/stripped.nc";
$cat=3;
print "aMD-new7-3\n";
}
#27ns
if ($fold_num==9) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new7-4/stripped.nc";
$cat=3;
print "aMD-new7-4\n";
}
#101 ns
if ($fold_num==10) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-close-0/stripped.nc";
$cat=2;
print "TL-close-0\n";
}
#100 ns
if ($fold_num==11) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-close-1/stripped.nc";
$cat=2;
print "TL-close-1\n";
}
#13ns
if ($fold_num==12) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-1/stripped.nc";
$cat=4;
print "TL-open1-1\n";
}
#13ns
if ($fold_num==13) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-2/stripped.nc";
$cat=4;
print "TL-open1-2\n";
}
#18ns
if ($fold_num==14) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-4/stripped.nc";
$cat=4;
print "TL-open1-4\n";
}
#28.5ns
if ($fold_num==15) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-5/stripped.nc";
$cat=4;
print "TL-open1-5\n";
}
#28ns
if ($fold_num==16) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-6/stripped.nc";
$cat=4;
print "TL-open1-6\n";
}
#16ns
if ($fold_num==17) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-1/stripped.nc";
$cat=4;
print "TL-open2-1\n";
}
#15.5ns
if ($fold_num==18) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-2/stripped.nc";
$cat=4;
print "TL-open2-2\n";
}
#107.5ns 
if ($fold_num==19) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-4/stripped.nc";
$cat=4;
print "TL-open2-4\n";
}
#31.5ns
if ($fold_num==20) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-5/stripped.nc";
$cat=4;
print "TL-open2-5\n";
}
#31.5ns
if ($fold_num==21) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-6/stripped.nc";
$cat=4;
print "TL-open2-6\n";
}
#16ns
if ($fold_num==22) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-1/stripped.nc";
$cat=4;
print "TL-open3-1\n";
}
#71ns
if ($fold_num==23) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-2/stripped.nc";
$cat=4;
print "TL-open3-2\n";
}
#19.5ns
if ($fold_num==24) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-3/stripped.nc";
$cat=4;
print "TL-open3-3\n";
}
#38.5ns
if ($fold_num==25) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-5/stripped.nc";
$cat=4;
print "TL-open3-5\n";
}
#27ns
if ($fold_num==26) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-6/stripped.nc";
$cat=4;
print "TL-open3-6\n";
}
#67 ns
if ($fold_num==27) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-1/stripped.nc";
$cat=5;
print "iso-CH3-1\n";
}
#57.5 ns
if ($fold_num==28) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-2/stripped.nc";
$cat=5;
print "iso-CH3-2\n";
}
#75.5 ns
if ($fold_num==29) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-3/stripped.nc";
$cat=5;
print "iso-CH3-3\n";
}
#77 ns
if ($fold_num==30) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-4/stripped.nc";
$cat=5;
print "iso-CH3-4\n";
}
#124.5 ns
if ($fold_num==31) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-5/stripped.nc";
$cat=5;
print "iso-CH3-5\n";
}
#124 ns
if ($fold_num==32) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-6/stripped.nc";
$cat=5;
print "iso-CH3-6\n";
}
#102.5 ns
if ($fold_num==33) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-1/stripped.nc";
$cat=5;
print "trans-CH3-1\n";
}
#74 ns
if ($fold_num==34) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-2/stripped.nc";
$cat=5;
print "trans-CH3-2\n";
}
#52.5 ns
if ($fold_num==35) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-3/stripped.nc";
$cat=5;
print "trans-CH3-3\n";
}
#189 ns
if ($fold_num==36) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-5/stripped.nc";
$cat=5;
print "trans-CH3-5\n";
}
#162 ns
if ($fold_num==37) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-6/stripped.nc";
$cat=5;
print "trans-CH3-6\n";
}
#74.5 ns
if ($fold_num==38) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-3-2/stripped.nc";
$cat=5;
print "trans-CH3-3-2\n";
}
#54.5 ns
if ($fold_num==39) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-1/stripped.nc";
$cat=5;
print "trans-CH3-m4-1\n";
}
#54.5 ns
if ($fold_num==40) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-2/stripped.nc";
$cat=5;
print "trans-CH3-m4-2\n";
}
#54.5 ns
if ($fold_num==41) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-3/stripped.nc";
$cat=5;
print "trans-CH3-m4-3\n";
}
#55.5 ns
if ($fold_num==42) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-5/stripped.nc";
$cat=5;
print "trans-CH3-m4-5\n";
}
#58.5 ns
if ($fold_num==43) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-6/stripped.nc";
$cat=5;
print "trans-CH3-m4-6\n";
}
#55.5 ns
if ($fold_num==44) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-7/stripped.nc";
$cat=5;
print "trans-CH3-m4-7\n";
}
#22.2 ns
if ($fold_num==45) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new-1/stripped.nc";
$cat=5;
print "trans-new-1\n";
}
#49.8 ns
if ($fold_num==46) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new-2/stripped.nc";
$cat=5;
print "trans-new-2\n";
}
#53.5 ns
if ($fold_num==47) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-1/stripped.nc";
$cat=5;
print "trans-new5-1\n";
}
#101.8 ns
if ($fold_num==48) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-2/stripped.nc";
$cat=5;
print "trans-new5-2\n";
}
#40.4 ns
if ($fold_num==49) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-3/stripped.nc";
$cat=5;
print "trans-new5-3\n";
}
#95.2 ns
if ($fold_num==50) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-5/stripped.nc";
$cat=5;
print "trans-new5-5\n";
}
#169.6 ns
if ($fold_num==51) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-6/stripped.nc";
$cat=5;
print "trans-new5-6\n";
}
#89 ns
if ($fold_num==52) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new7-3/stripped.nc";
$cat=5;
print "trans-new7-3\n";
}
#156.6 ns
if ($fold_num==53) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new7-5/stripped.nc";
$cat=5;
print "trans-new7-5\n";
}
#75 ns
if ($fold_num==54) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1/stripped.nc";
$cat=5;
print "trans-new8-1\n";
}
#40.8 ns
if ($fold_num==55) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1-3/stripped.nc";
$cat=5;
print "trans-new8-1-3\n";
}
#40.2 ns
if ($fold_num==56) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1-4/stripped.nc";
$cat=5;
print "trans-new8-1-4\n";
}
#12.2 ns
if ($fold_num==57) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1-5/stripped.nc";
$cat=5;
print "trans-new8-1-5\n";
}
#68.8 ns
if ($fold_num==58) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-3/stripped.nc";
$cat=5;
print "trans-new8-3\n";
}
#12.4 ns
if ($fold_num==59) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1-7/stripped.nc";
$cat=5;
print "trans-new8-1-7\n";
}
#14.2 ns
if ($fold_num==60) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-2/stripped.nc";
$cat=5;
print "trans-new8-2\n";
}
#94.2 ns
if ($fold_num==61) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-6/stripped.nc";
$cat=5;
print "trans-new8-6\n";
}
#56ns
if ($fold_num==62) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/stripped.nc";
$cat=6;
print "aMD-new-4-alt-rep3\n";
}
#93.5ns
if ($fold_num==63) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep5/stripped.nc";
$cat=6;
print "aMD-new-4-alt-rep5\n";
}
#59ns
if ($fold_num==64) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep1/stripped.nc";
$cat=6;
print "aMD-new-4-alt3-rep1\n";
}
#67ns
if ($fold_num==65) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep2/stripped.nc";
$cat=6;
print "aMD-new-4-alt3-rep2\n";
}
#61.5ns
if ($fold_num==66) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep3/stripped.nc";
$cat=6;
print "aMD-new-4-alt3-rep3\n";
}
#62ns
if ($fold_num==67) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep4/stripped.nc";
$cat=6;
print "aMD-new-4-alt3-rep4\n";
}
#118.5ns
if ($fold_num==68) {
$prmtop= "out15-parmed2.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep5/stripped.nc";
$cat=6;
print "aMD-new-4-alt3-rep5\n";
}
#101ns
if ($fold_num==69) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep6/stripped.nc";
$cat=6;
print "aMD-new-4-alt3-rep6\n";
}
#100ns
if ($fold_num==70) {
$prmtop= "out15-parmed.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/syner+3-5/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/syner+4-6/stripped.nc";
$cat=7;
print "syner+4-6\n";
}
#100ns
if ($fold_num==71) {
$prmtop= "out15-parmed.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/syner+3-5/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/syner+3-5/stripped.nc";
$cat=8;
print "syner+3-5\n";
}
#100ns
if ($fold_num==72) {
$prmtop= "out15-parmed.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/syner+3-5/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/syner+2-new/stripped.nc";
$cat=9;
print "syner+2-new\n";
}
#250ns
if ($fold_num==73) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-2/stripped.nc";
$cat=1;
print "aMD-new3-2\n";
}
#252ns
if ($fold_num==74) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-3/stripped.nc";
$cat=1;
print "aMD-new3-3\n";
}
#251ns
if ($fold_num==75) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-4/stripped.nc";
$cat=1;
print "aMD-new3-4\n";
}
#368.5ns
if ($fold_num==76) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-5/stripped.nc";
$cat=1;
print "aMD-new3-5\n";
}
#365ns
if ($fold_num==77) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-6/stripped.nc";
$cat=1;
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
if ($fold_num==77) {
$divide_trigger=3;
$start_1=0;
$last_1=250;
$start_2=251;
$last_2=500;
$start_3=501;
$last_3=729;
}

#Assign Ligand Index depending on Simulation Trajectory
if ($cat==1){
$sel_ctp=$sel_ctp_1;
}
if ($cat==2){
$sel_ctp=$sel_ctp_2;
}
if ($cat==3){
$sel_ctp=$sel_ctp_3;
}
if ($cat==4){
$sel_ctp=$sel_ctp_4;
}
if ($cat==5){
$sel_ctp=$sel_ctp_5;
}
if ($cat==6){
$sel_ctp=$sel_ctp_6;
}
if ($cat==7){
$sel_ctp=$sel_ctp_7;
}
if ($cat==8){
$sel_ctp=$sel_ctp_8;
}
if ($cat==9){
$sel_ctp=$sel_ctp_9;
}


for (my $divide=0; $divide<$divide_trigger; $divide++) {

###################################################################################################### 
#############################################  ZONE LOOP #############################################
######################################################################################################
for (my $zone=1; $zone<19; $zone++) {

#Deal with special case
if (($fold_num==77) and ($divide_trigger>1)){
$sel_ctp=$sel_ctp_2;
}
#END Deal with special case

#ASSIGN selection
if ($zone==1){
$sel_res=$sel1_res;
}
if ($zone==2){
$sel_res=$sel2_res;
}
if ($zone==3){
$sel_res=$sel3_res;
}
if ($zone==4){
$sel_res=$sel4_res;
}
if ($zone==5){
$sel_res=$sel5_res;
}
if ($zone==6){
$sel_res=$sel6_res;
}
if ($zone==7){
$sel_res=$sel7_res;
}
if ($zone==8){
$sel_res=$sel8_res;
}
if ($zone==9){
$sel_res=$sel9_res;
}
if ($zone==10){
$sel_res=$sel10_res;
}
if ($zone==11){
$sel_res=$sel11_res;
}
if ($zone==12){
$sel_res=$sel12_res;
}
if ($zone==13){
$sel_res=$sel13_res;
}
if ($zone==14){
$sel_res=$sel14_res;
}
if ($zone==15){
$sel_res=$sel15_res;
}
if ($zone==16){
$sel_res=$sel16_res;
}
if ($zone==17){
$sel_res=$sel17_res;
}
if ($zone==18){
$sel_res=$sel18_res;
}
#END ASSIGN selection


my $outfile="scr-stat.vmd";
#OUTPUT
open (FILE2, "> $outfile") || die;
#normal file
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
if ($zone==1){
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
print (FILE2 "set sel1 [atomselect top \"(resid $sel_ctp) and (within 3 of resid $sel_res)\" frame \$frame]\n");
print (FILE2 "set sel2 [\$sel1 get resid]\n");
print (FILE2 "set sel3 [lsort -unique \$sel2]\n");
print (FILE2 "set sel4 [llength \$sel3]\n");
#set dummy value for sel11 and sel 12
print (FILE2 "set sel11 \"\"\n");
print (FILE2 "set sel12 0\n");
#Add NTPs bound at channel NTP but not at channel directly and not already bound at another zone
print (FILE2 "set sel5 [atomselect top \"((resid $sel_ctp) and (within 3 of protein))\" frame \$frame]\n");
print (FILE2 "set sel6 [\$sel5 get resid]\n");
print (FILE2 "set sel7 [lsort -unique \$sel6]\n");
print (FILE2 "set sel8 [llength \$sel7]\n");
print (FILE2 "if {\$sel4 > 0 && \$sel8 > 0} {\n");
print (FILE2 "set sel9 [atomselect top \"((resid $sel_ctp) and (within 4 of resid \$sel3)) and not resid \$sel3 \$sel7\" frame \$frame]\n");
print (FILE2 "set sel10 [\$sel9 get resid]\n");
print (FILE2 "set sel11 [lsort -unique \$sel10]\n");
print (FILE2 "set sel12 [llength \$sel11]\n");
print (FILE2 "}\n");
#total:
print (FILE2 "set nb [expr \$sel4 + \$sel12]\n");
print (FILE2 "puts \$out \"$zone \$nb \$sel3 \$sel11\"\n");
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
###################################################################################################### 
############################################# END  ZONE loop #########################################
######################################################################################################




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
my $outfile="output-global-polar-final2.txt";
open (FILE3, "> $outfile") || die;
my @pdb_input = read_file("output-global-polar.txt") or die;
my $array_size=scalar @pdb_input;
my @line_handle;
my $stat_sum;
my $stat_frames_sum=0;
my $stat_frames;
my $buffer;
my $sel1;
my $sel2;
my ($min,$max);
my $sel_next;
my $sel_CH_frame;
my $nb_frames_step;
my @line_handle_next;
my $nb_resid;
my $pos;
my $resid_total_size=0;

print (FILE3"STAT GLOBAL\n");

#NB FRAMES
for (my $count=0; $count<$array_size; $count++) {
@line_handle = split ( /\s+/, @pdb_input[$count] );
$sel1= @line_handle[0];
if ($sel1==5000) {
$stat_frames= @line_handle[1];
$stat_frames_sum=$stat_frames_sum+$stat_frames;
}
}
print (FILE3"nb_frames $stat_frames_sum\n");
print (FILE3"\n");

#ZONE STATS
for (my $zone=1; $zone<19; $zone++) {
$stat_sum=0;
my @array_max;

for (my $count=0; $count<$array_size; $count++) {
@line_handle = split ( /\s+/, @pdb_input[$count] );
$sel1= @line_handle[0];
$sel2= @line_handle[1];
if ($sel1==$zone){
$stat_sum=$stat_sum+$sel2;
push (@array_max, ($sel2));
}
}

($min,$max) = minmax(@array_max);
$stat_sum=($stat_sum/$stat_frames_sum);
$stat_sum=sprintf "%.3f", $stat_sum;
print (FILE3"avg zone $zone stat $stat_sum\n");
print (FILE3"max zone $zone max $max\n");
}
print (FILE3"\n");


#CHANNEL STATS
#Here the protein binding sites are further processed 
#into 5 protein channels: CH2, CH3A, CH3B, CH3C and CH3D
for (my $ch=1; $ch<6; $ch++) {
$stat_sum=0;
$nb_frames_step=0;
my @array_max;

for (my $count=0; $count<$array_size; $count++) {
@line_handle = split ( /\s+/, @pdb_input[$count] );
$sel1= @line_handle[0];
$sel2= @line_handle[1];

#get number of frames for the current traj
if ($sel1==5000){
$nb_frames_step=$sel2+1;
}

#CH2
if ($ch==1){
my @resid;
my $resid;
$resid_total_size=0;
if ($sel1==1){
#ZONE 1 contribution
$nb_resid=$sel2;
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle[$pos];
push (@resid, $resid);
}
#ADD ZONE 2 contribution
@line_handle_next = split ( /\s+/, @pdb_input[$count+$nb_frames_step] );
$nb_resid= @line_handle_next[1];
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle_next[$pos];
push (@resid, $resid);
}
#ADD ZONE 3 contribution
@line_handle_next = split ( /\s+/, @pdb_input[$count+(2*$nb_frames_step)] );
$nb_resid= @line_handle_next[1];
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle_next[$pos];
push (@resid, $resid);
}
#ADD ZONE 4 contribution
@line_handle_next = split ( /\s+/, @pdb_input[$count+(3*$nb_frames_step)] );
$nb_resid= @line_handle_next[1];
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle_next[$pos];
push (@resid, $resid);
}
#ADD ZONE 13 contribution
@line_handle_next = split ( /\s+/, @pdb_input[$count+(12*$nb_frames_step)] );
$nb_resid= @line_handle_next[1];
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle_next[$pos];
push (@resid, $resid);
}
#TOTAL contribution
my @resid_total= uniq(@resid);
$resid_total_size=scalar @resid_total;
##TOTAL frame
$stat_sum=$stat_sum+$resid_total_size;
push (@array_max, ($resid_total_size));
}
}

#CH3A
if ($ch==2){
my @resid;
my $resid;
$resid_total_size=0;
if ($sel1==5){
#ZONE 5 contribution
$nb_resid=$sel2;
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle[$pos];
push (@resid, $resid);
}
#ADD ZONE 14 contribution
@line_handle_next = split ( /\s+/, @pdb_input[$count+(9*$nb_frames_step)] );
$nb_resid= @line_handle_next[1];
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle_next[$pos];
push (@resid, $resid);
}
#TOTAL contribution
my @resid_total= uniq(@resid);
$resid_total_size=scalar @resid_total;
##TOTAL frame
$stat_sum=$stat_sum+$resid_total_size;
push (@array_max, ($resid_total_size));
}
}

#CH3B
if ($ch==3){
my @resid;
my $resid;
$resid_total_size=0;
if ($sel1==6){
#ZONE 6 contribution
$nb_resid=$sel2;
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle[$pos];
push (@resid, $resid);
}
#ADD ZONE 7 contribution
@line_handle_next = split ( /\s+/, @pdb_input[$count+$nb_frames_step] );
$nb_resid= @line_handle_next[1];
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle_next[$pos];
push (@resid, $resid);
}
#ADD ZONE 15 contribution
@line_handle_next = split ( /\s+/, @pdb_input[$count+(9*$nb_frames_step)] );
$nb_resid= @line_handle_next[1];
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle_next[$pos];
push (@resid, $resid);
}
#TOTAL contribution
my @resid_total= uniq(@resid);
$resid_total_size=scalar @resid_total;
##TOTAL frame
$stat_sum=$stat_sum+$resid_total_size;
push (@array_max, ($resid_total_size));
}
}

#CH3C
if ($ch==4){
my @resid;
my $resid;
$resid_total_size=0;
if ($sel1==8){
#ZONE 8 contribution
$nb_resid=$sel2;
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle[$pos];
push (@resid, $resid);
}
#ADD ZONE 9 contribution
@line_handle_next = split ( /\s+/, @pdb_input[$count+$nb_frames_step] );
$nb_resid= @line_handle_next[1];
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle_next[$pos];
push (@resid, $resid);
}
#ADD ZONE 16 contribution
@line_handle_next = split ( /\s+/, @pdb_input[$count+(8*$nb_frames_step)] );
$nb_resid= @line_handle_next[1];
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle_next[$pos];
push (@resid, $resid);
}
#TOTAL contribution
my @resid_total= uniq(@resid);
$resid_total_size=scalar @resid_total;
##TOTAL frame
$stat_sum=$stat_sum+$resid_total_size;
push (@array_max, ($resid_total_size));
}
}

#CH3D
if ($ch==5){
my @resid;
my $resid;
$resid_total_size=0;
if ($sel1==10){
#ZONE 10 contribution
$nb_resid=$sel2;
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle[$pos];
push (@resid, $resid);
}
#ADD ZONE 11 contribution
@line_handle_next = split ( /\s+/, @pdb_input[$count+$nb_frames_step] );
$nb_resid= @line_handle_next[1];
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle_next[$pos];
push (@resid, $resid);
}
#ADD ZONE 12 contribution
@line_handle_next = split ( /\s+/, @pdb_input[$count+(2*$nb_frames_step)] );
$nb_resid= @line_handle_next[1];
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle_next[$pos];
push (@resid, $resid);
}
#ADD ZONE 17 contribution
@line_handle_next = split ( /\s+/, @pdb_input[$count+(7*$nb_frames_step)] );
$nb_resid= @line_handle_next[1];
$pos=1;
for (my $scan_line=0; $scan_line<$nb_resid; $scan_line++){
$pos++;
$resid = @line_handle_next[$pos];
push (@resid, $resid);
}
#TOTAL contribution
my @resid_total= uniq(@resid);
$resid_total_size=scalar @resid_total;
##TOTAL frame
$stat_sum=$stat_sum+$resid_total_size;
push (@array_max, ($resid_total_size));
}
}


}
#END LINE LOOP


($min,$max) = minmax(@array_max);
$stat_sum=($stat_sum/$stat_frames_sum);
$stat_sum=sprintf "%.3f", $stat_sum;
print (FILE3"avg CH $ch stat $stat_sum\n");
print (FILE3"max CH $ch max $max\n");
}
#END CHANNEL LOOP


print (FILE3"\n");


close (FILE3);

my $datestring_end = localtime();
print "START = $datestring_start\n";
print "END = $datestring_end\n";

exit;



