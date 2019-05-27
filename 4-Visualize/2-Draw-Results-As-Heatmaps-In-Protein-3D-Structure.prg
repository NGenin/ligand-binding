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

my $trigger_heatmap;
my $trigger_res;
my $trigger_ctp;

#Visualization parameters:
#trigger_heatmap=1 triggers the display of the binding sites coloration for binding propensity scores
$trigger_heatmap=0;
#trigger_res=1 triggers the display of the protein binding sites divided into macro-regions or zones, without heatmap coloration for the binding propensities
$trigger_res=1;
#trigger_ctp=1 triggers the display of the ligands
$trigger_ctp=0;



#Variable Declaration and Argument Extraction
#Arg[0] is the protein trajectory ID
my $v;
my @abs_filename;
my $prmtop;
my $dir = getcwd;
my $abs_filename_prmtop;
my $nb_traj;
my $fold_num=$ARGV[0];
my $sel;
my $sel_resid;
my $size;
my $occupancy;
my $stat;
my $colorid;

#List of protein ligand-binding amino acids, and their respective ligand-binding scores.
#Divided into 20 macro-regions or zones (@sel_ 1 to 19 and ext).
#In addition specific amino acids belonging to flexible domains or junction regions, are assigned to arrays @sel_flex and @sel_junc.

my @sel1_res = qw(749 688 699 734 698 742 695 691 750 701 692 748 687 702 745 686 694 710 738 689 751 735 746 731 690 706 741 705 703 766 767 769 770 768 740 737 717 715 719);
my @sel1_stat = qw(7.80 7.56 5.13 4.37 4.32 3.41 2.89 2.86 2.60 2.26 2.25 2.11 2.04 2.03 1.91 1.76 1.38 1.37 1.34 1.16 1.15 1.11 1.06 1.02 0.83 0.69 0.68 0.67 0.61 32.25 4.03 3.69 3.52 1.79 4.48 0.97 0.42 1.83 0.27);

my @sel2_res = qw(595 3232 541 3310 538 2507 677 2506 2490 3233 596 3234 3235 594 754 2487 3236 756 2492 3315 3312 3314 3313 3239 753 3231 752 537 3230 3309 3336 3316 679 3327 593 542 3323 592 682 3225 758 620 628 630 3267 3248 676 759 3317 755 3247 3318 3322 678 3343 3338 622 683 3334 765 762 764 763 597 618 610 611 614 612 609 613 615 602 608 616 600 605 617 606 601 603 599 604);
my @sel2_stat = qw(37.00 31.06 30.53 29.68 28.32 24.12 23.41 21.35 16.67 15.65 15.00 14.73 12.67 12.41 11.99 10.19 10.06 9.93 9.62 9.61 8.69 7.75 6.83 5.75 5.44 5.38 5.29 5.27 5.05 4.85 4.84 4.63 4.59 4.35 4.28 4.20 4.12 4.12 4.01 3.84 3.67 3.17 3.10 2.82 2.82 2.81 2.70 2.48 2.17 2.12 2.10 2.09 2.03 2.00 1.96 1.94 1.92 1.86 1.86 2.08 0.98 0.67 0.46 14.46 14.22 9.07 6.31 5.83 5.38 5.02 3.47 3.26 3.19 3.08 2.92 2.52 2.40 2.16 2.03 1.34 1.22 0.93 0.91);

my @sel3_res = qw(910 909 922 1043 901 905 919 958 912 951 903 900 1044 914 904 963 954 906 962 896 916 915 902 945 897 959 967 917 908 961 923 895 1048 964 955 1047 918 949 941 921 913 960 948 926 931 944 952 907 899 965 975 966 947);
my @sel3_stat = qw(32.33 31.85 31.80 27.13 26.66 23.62 22.41 18.10 15.36 15.30 11.89 11.89 11.55 11.39 10.81 10.72 10.09 9.46 8.22 7.76 7.50 7.22 6.83 6.77 6.10 5.54 5.03 5.00 4.45 4.41 4.40 4.26 3.89 3.69 3.26 3.18 3.14 3.09 2.84 2.57 2.40 2.35 2.25 2.25 2.00 1.96 1.87 1.82 1.73 1.71 1.48 1.10 1.09);

my @sel4_res = qw(1371 1347 1367 1323 3059 1346 1324 1310 1309 1308 1311 1312 1307 1305 1306 1320 1313 1315 1319 1316 1314 1322 1318 1317 1137 1140 1141 1136 1138 1135);
my @sel4_stat = qw(5.88 1.72 1.41 1.17 0.98 0.59 0.38 8.10 4.09 4.26 3.40 3.36 0.95 0.45 0.41 4.41 2.61 2.11 1.64 1.55 1.48 1.36 1.34 1.14 5.76 5.61 1.97 1.30 1.02 0.68);

my @sel5_res = qw(1371 699 734 742 701 702 1347 710 738 1323 735 731 706 741 705 703 1324 1310 1309 1308 1311 1312 1307 1305 1306 1320 1313 1315 1319 1316 1314 1322 1318 1317 740 737 715 717 719 1137 1140 1141 1136 1138 1135);
my @sel5_stat = qw(5.88 5.13 4.37 3.41 2.26 2.03 1.72 1.37 1.34 1.17 1.11 1.02 0.69 0.68 0.67 0.61 0.38 8.10 4.09 4.26 3.40 3.36 0.95 0.45 0.41 4.41 2.61 2.11 1.64 1.55 1.48 1.36 1.34 1.14 4.48 0.97 1.83 0.42 0.27 5.76 5.61 1.97 1.30 1.02 0.68);

my @sel6_res = qw(3454 3484 3444 3487 3455 3480 3450 3453 3483 3442 3436 3445 3481 3473 3485 3451 3426 3486 3449 3482 3479 3452 3476 3430 3429 3419 3475 3447 3448 3422 3474 3423 3446 3435 3428 3421 3477 3433 3425 714 713 3431 3427 712 3424 3420 3471 3465 3466 3462 3459 3467 3463 3468 3460 3461 715 717 719);
my @sel6_stat = qw(16.14 15.95 14.29 12.42 10.35 6.76 6.54 5.48 5.19 3.63 3.37 3.34 3.14 3.09 3.03 2.96 2.63 2.58 2.42 2.27 2.21 2.04 1.72 1.65 1.53 1.46 1.45 1.37 1.22 1.15 1.12 1.04 0.94 0.81 0.69 0.67 0.65 0.60 0.48 1.83 1.48 1.00 0.96 0.57 0.54 0.49 6.49 4.91 1.34 1.10 0.81 0.77 0.56 0.55 0.40 0.29 1.83 0.42 0.27);

my @sel7_res = qw(1210 1209 1211 1244 1245 1241 1214 1215 1179 1213 1242 1208 1146 1182 1216 1143 1189 1187 1186 1185 1188);
my @sel7_stat = qw(3.88 2.71 1.61 1.60 1.57 0.99 0.89 0.78 0.69 0.67 0.58 0.53 0.49 0.47 2.49 0.75 0.58 0.54 0.49 0.27 0.25);

my @sel8_res = qw(1144 1197 1224 1221 1259 1225 1341 1255 1258 1193 1218 1220 1151 1290 1254 1269 1253 1230 1256 1262 1263 1233 1229 1265 1232 1147 1291 1267 1191 1217 1261 1260 1231 1289 1227 1268 1190 1257 1266 1194 1192 1264 1342 1235 1216 1137 1140 1141 1136 1138 1135 1143 1189 1187 1186 1185 1188);
my @sel8_stat = qw(8.06 5.58 4.12 3.69 3.64 3.66 2.38 2.34 1.81 1.73 1.62 1.48 1.38 1.33 1.21 1.15 1.09 1.03 0.93 0.91 0.90 0.85 0.84 0.83 0.78 0.72 0.71 0.61 0.60 0.60 0.59 0.58 0.58 0.55 0.53 0.51 0.48 0.45 0.45 0.43 0.38 0.31 0.28 0.21 2.49 5.76 5.61 1.97 1.30 1.02 0.68 0.75 0.58 0.54 0.49 0.27 0.25);

my @sel9_res = qw(2917 3001 3004 3008 3009 3005 2971 2918 2968 3012 3011 2994 2995 2963 2969 2920 2921 2925 2919 2972 3002 3007 2970 2964 2967 2975 2957 2977 2966 2965 2989 3006 2962 2959 2973 2976 2954 2955 3000 3013 2960 2991 2999 3052 2953 2952 2961 2924 2922 2993);
my @sel9_stat = qw(8.76 7.16 5.95 5.46 5.08 4.69 4.26 4.17 3.49 2.02 2.02 1.60 1.46 1.33 1.28 1.16 1.10 1.06 1.04 1.01 0.95 0.94 0.93 0.89 0.89 0.83 0.67 0.65 0.62 0.61 0.51 0.48 0.41 0.41 0.40 0.40 0.39 0.36 0.34 0.33 0.32 0.32 0.31 0.29 0.28 0.27 0.25 0.24 0.24 0.24);

my @sel10_res = qw(3103 130 123 3088 3126 1399 127 125 3102 1414 124 116 3040 126 146 122 118 3036 119 3039 142 3042 3104 225 3035 115 120 3038 121 3124 3041 131 3017 3037 113 3105 2984 111 114 110 3016 2985 147 148);
my @sel10_stat = qw(20.51 18.24 15.54 14.90 12.17 8.87 8.22 6.55 5.52 4.76 4.19 3.41 3.41 2.94 2.79 2.77 2.38 2.34 2.14 2.14 2.10 2.09 2.07 1.95 1.88 1.59 1.51 1.51 1.48 1.48 1.42 1.33 1.33 1.18 1.11 1.02 0.98 0.97 0.88 0.73 0.69 0.60 0.59 0.58);

my @sel11_res = qw(7 10 2639 3082 2641 11 1443 3194 8 3192 3078 3084 3080 2617 875 1463 3198 2640 1456 3200 1441 2571 1461 874 3197 3196 3085 3069 865 1454 3199 1453 13 1460 1462 866 1439 876 3149 3195 3077 873 867 3083 1464 3081 1442 2614 871 1448 3074 3073 3076 1438 3068 3070);
my @sel11_stat = qw(36.63 21.79 21.76 15.61 14.76 12.58 11.64 11.55 10.83 10.56 10.35 8.87 8.81 8.32 7.16 6.87 6.64 5.86 5.68 5.61 5.18 5.06 4.92 4.90 4.58 4.54 3.97 3.79 3.75 3.66 3.54 3.33 3.25 3.08 2.76 2.74 2.70 2.65 2.53 2.29 2.29 2.15 2.05 2.04 2.03 1.85 1.54 1.41 1.21 1.13 1.20 1.03 1.00 0.97 0.76 0.57);

my @sel12_res = qw(24 2608 232 21 2588 2637 2607 2589 233 2606 2609 2605 2610 20 2598 2612 238 18 25 23 2591 2587 2611 2603 231 2593 2604 235 2597 217 29 2600 26 27 218 2601 221 19 2590 2602 28 2596 2599 222);
my @sel12_stat = qw(74.60 73.68 52.07 48.85 44.76 43.22 41.42 39.37 30.00 26.36 23.55 23.27 21.13 19.92 13.88 13.56 13.54 11.28 9.44 6.96 6.15 5.60 5.22 4.53 4.39 4.15 4.11 3.77 3.67 3.67 3.43 3.22 2.93 2.83 2.68 2.42 2.38 2.20 2.09 1.72 1.56 1.52 1.51 1.46);

my @sel13_res = qw(634 639 641 640 635 523 642 894 633 632 638 774 527 771 530 529 525 1380 528 1383 1378 1381 1374 524 526 1382 766 767 769 770 765 768 762 764 763 597 618 610 611 614 612 609 613 615 602 608 616 600 605 617 606 601 603 599 604 1093 1106 1097 1099 1098 1096 1092 1094 1100 1091 1095 1104 1103 1107 1109 1088 1102 1110 1108 740 737 1379 1375);
my @sel13_stat = qw(17.73 14.43 12.99 12.96 5.86 5.27 5.24 4.85 4.37 1.79 1.60 1.60 1.40 1.28 1.20 1.20 1.11 1.04 0.96 0.76 0.84 0.68 0.63 0.50 0.38 0.27 32.25 4.03 3.69 3.52 2.08 1.79 0.98 0.67 0.46 14.46 14.22 9.07 6.31 5.83 5.38 5.02 3.47 3.26 3.19 3.08 2.92 2.52 2.40 2.16 2.03 1.34 1.22 0.93 0.91 19.11 12.86 12.26 11.98 7.00 6.41 7.10 5.70 4.65 3.11 2.87 1.86 1.17 0.99 0.91 0.98 0.67 0.54 0.50 4.48 0.97 2.15 0.84);

my @sel14_res = qw(730 729 811 1300 1328 782 733 779 781 1093 1106 1097 1099 1098 1096 1092 1094 1100 1091 1095 1104 1103 1107 1109 1088 1102 1110 1108 787 785 783 786 788 791 790 789 1310 1309 1308 1311 1312 1307 1305 1306 726 722 725 723 1379 1375);
my @sel14_stat = qw(2.29 1.56 1.52 1.37 0.71 0.67 0.61 0.61 0.35 19.11 12.86 12.26 11.98 7.00 6.41 7.10 5.70 4.65 3.11 2.87 1.86 1.17 0.99 0.91 0.98 0.67 0.54 0.50 7.47 3.75 2.54 2.00 1.61 1.09 0.53 0.46 8.10 4.09 4.26 3.40 3.36 0.95 0.45 0.41 3.90 2.02 2.12 1.42 2.15 0.84);

my @sel15_res = qw(1159 1153 3471 3465 3466 3462 3459 3467 3463 3468 3460 3461 1158 1297 1299 1296 1330 1298 1331 726 722 725 723);
my @sel15_stat = qw(3.00 0.24 6.49 4.91 1.34 1.10 0.81 0.77 0.56 0.55 0.40 0.29 2.99 7.73 5.23 2.08 2.02 1.21 0.58 3.90 2.02 2.12 1.42);

my @sel16_res = qw(1336 2997 1286 1285 2996 1715 3024 1716 1282 1709 1712 1294 1713 1281 1292 3023 1711 1154 1717 3049 1710 1293 1288 1287 1332 1278 3025 3047 1335 1721 3027 1158 1296 1331);
my @sel16_stat = qw(2.35 2.18 3.26 2.22 1.81 1.79 1.77 1.72 1.53 1.46 1.46 1.46 1.38 1.16 1.04 0.96 0.95 0.88 0.65 0.56 0.52 0.46 0.46 0.41 0.38 0.33 0.33 0.32 0.21 0.20 0.16 2.99 2.08 0.58);

my @sel17_res = qw(1401 1412 3090 134 1413 3091 3092 3032 3031 133 1358 129 137 1357 3045 3044 3097 1403 3034 1354 140 1123 1406 1407 1405 1408);
my @sel17_stat = qw(6.22 4.88 4.59 4.34 3.39 1.94 1.82 1.54 1.25 0.93 0.71 0.70 0.68 0.68 0.60 0.59 0.52 0.46 0.36 0.29 0.25 4.01 1.05 0.64 0.38 0.36);

my @sel18_res = qw(1126 1963 1124 1116 1966 1961 1967 839 1678 1117 1964 1125 1113 1120 1122 1962 1128 1681 1129 1127 847 1679 1970 1680 1682 1969 1114 854 1351 1684 1677 1123 1093 1106 1097 1099 1098 1096 1092 1094 1100 1091 1095 1104 1103 1107 1109 1088 1102 1110 1108 787 785 783 786 788 791 790 789 1406 1407 1405 1408 844 840 842 843 836 846 1297 1299 1296 1330 1298 1331);
my @sel18_stat = qw(18.21 13.32 10.32 10.04 10.00 9.03 8.84 7.50 7.12 4.71 4.50 4.26 4.04 3.23 2.97 2.12 2.09 2.02 1.65 1.48 1.20 1.00 0.75 0.69 0.58 0.42 0.38 0.34 0.33 0.26 0.18 4.01 19.11 12.86 12.26 11.98 7.00 6.41 7.10 5.70 4.65 3.11 2.87 1.86 1.17 0.99 0.91 0.98 0.67 0.54 0.50 7.47 3.75 2.54 2.00 1.61 1.09 0.53 0.46 1.05 0.64 0.38 0.36 22.29 9.24 8.71 3.66 2.81 1.72 7.73 5.23 2.08 2.02 1.21 0.58);

my @sel19_res = qw(486 2188 2442 2409 488 451 845 453 490 841 1983 484 485 1984 2191 2187 2185 838 2258 455 1985 489 2441 487 2443 837 2190 2186 482 2259 491 1093 1106 1097 1099 1098 1096 1092 1094 1100 1091 1095 1104 1103 1107 1109 1088 1102 1110 1108 766 844 840 842 843 836 846);
my @sel19_stat = qw(51.31 50.62 40.79 38.61 28.26 22.87 22.75 18.81 18.07 15.97 14.35 13.87 12.30 11.68 10.15 7.18 6.50 6.38 4.12 2.65 2.61 2.18 1.92 1.80 1.13 1.06 1.01 0.76 0.63 0.58 0.58 19.11 12.86 12.26 11.98 7.00 6.41 7.10 5.70 4.65 3.11 2.87 1.86 1.17 0.99 0.91 0.98 0.67 0.54 0.50 32.25 22.29 9.24 8.71 3.66 2.81 1.72);

my @sel_ext_res = qw(1022 1005 1019 1009 1023 3108 1015 1002 1017 3120 1010 1026 3122 3071 1006 1008 1016 3067 998 1033 1020 1018 1012 3075 1029 1030 1021 1027);
my @sel_ext_stat = qw(9.72 5.73 5.03 4.78 4.08 3.62 3.39 3.09 2.70 2.70 1.73 1.59 1.47 1.37 1.36 1.29 1.22 1.10 1.04 1.02 0.96 0.90 0.87 0.86 0.81 0.73 0.54 0.47);

my @sel_flex_res = qw(597 618 610 611 614 612 609 613 615 602 608 616 600 605 617 606 601 603 599 604 1310 1309 1308 1311 1312 1307 1305 1306 1320 1313 1315 1319 1316 1314 1322 1318 1317 3471 3465 3466 3462 3459 3467 3463 3468 3460 3461 1093 1106 1097 1099 1098 1096 1092 1094 1100 1091 1095 1104 1103 1107 1109 1088 1102 1110 1108 1406 1407 1405 1408 844 840 842 843 836 846 787 785 783 786 788 791 790 789 767 769 770 765 768 762 764 763 766 1123);
my @sel_flex_stat = qw(14.46 14.22 9.07 6.31 5.83 5.38 5.02 3.47 3.26 3.19 3.08 2.92 2.52 2.40 2.16 2.03 1.34 1.22 0.93 0.91 8.10 4.09 4.26 3.40 3.36 0.95 0.45 0.41 4.41 2.61 2.11 1.64 1.55 1.48 1.36 1.34 1.14 6.49 4.91 1.34 1.10 0.81 0.77 0.56 0.55 0.40 0.29 19.11 12.86 12.26 11.98 7.00 6.41 7.10 5.70 4.65 3.11 2.87 1.86 1.17 0.99 0.91 0.98 0.67 0.54 0.50 1.05 0.64 0.38 0.36 22.29 9.24 8.71 3.66 2.81 1.72 7.47 3.75 2.54 2.00 1.61 1.09 0.53 0.46 4.03 3.69 3.52 2.08 1.79 0.98 0.67 0.46 32.25 4.01);

my @sel_junc_res = qw(715 717 719 1137 1140 1141 1136 1138 1135 1143 1216 1189 1187 1186 1185 1188 1158 1296 1297 1299 1330 1298 1331 740 737 1379 1375 726 722 725 723);
my @sel_junc_stat = qw(1.83 0.42 0.27 5.76 5.61 1.97 1.30 1.02 0.68 0.75 2.49 0.58 0.54 0.49 0.27 0.25 2.99 2.08 7.73 5.23 2.02 1.21 0.58 4.48 0.97 2.15 0.84 3.90 2.02 2.12 1.42);


my @sel_res=(\@sel1_res,\@sel2_res,\@sel3_res,\@sel4_res,\@sel5_res,\@sel6_res,\@sel7_res,\@sel8_res,\@sel9_res,\@sel10_res,\@sel11_res,\@sel12_res,\@sel13_res,\@sel14_res,\@sel15_res,\@sel16_res,\@sel17_res,\@sel18_res,\@sel19_res,\@sel_ext_res,\@sel_flex_res,\@sel_junc_res);

my @sel_stat=(\@sel1_stat,\@sel2_stat,\@sel3_stat,\@sel4_stat,\@sel5_stat,\@sel6_stat,\@sel7_stat,\@sel8_stat,\@sel9_stat,\@sel10_stat,\@sel11_stat,\@sel12_stat,\@sel13_stat,\@sel14_stat,\@sel15_stat,\@sel16_stat,\@sel17_stat,\@sel18_stat,\@sel19_stat,\@sel_ext_stat,\@sel_flex_stat,\@sel_junc_stat);



my @sel_colorid= qw(23 26 9 14 7 10 23 10 26 9 26 10 5 15 9 5 15 12 6 6 4 3);

#List of Protein Trajectories
#A
#250ns
if ($fold_num==5) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-1/stripped.nc";
$nb_traj=1;
print "aMD-new3-1\n";
}
#250ns
if ($fold_num==6) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-2/stripped.nc";
$nb_traj=1;
print "aMD-new3-2\n";
}
#252ns
if ($fold_num==7) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-3/stripped.nc";
$nb_traj=1;
print "aMD-new3-3\n";
}
#251ns
if ($fold_num==8) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-4/stripped.nc";
$nb_traj=1;
print "aMD-new3-4\n";
}
#368.5ns
if ($fold_num==9) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-5/stripped.nc";
$nb_traj=1;
print "aMD-new3-5\n";
}
#365ns
if ($fold_num==10) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new3-6/stripped.nc";
$nb_traj=1;
print "aMD-new3-6\n";
}

#B
#75.75ns
if ($fold_num==17) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new6-3/stripped.nc";
$nb_traj=1;
print "aMD-new6-3\n";
}
#122ns
if ($fold_num==18) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new6-4/stripped.nc";
$nb_traj=1;
print "aMD-new6-4\n";
}
#138.5ns
if ($fold_num==19) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new6-5/stripped.nc";
$nb_traj=1;
print "aMD-new6-5\n";
}
#123.5ns
if ($fold_num==20) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new6-6/stripped.nc";
$nb_traj=1;
print "aMD-new6-6\n";
}

#C
#28.5ns
if ($fold_num==141) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new7-1/stripped.nc";
$nb_traj=1;
print "aMD-new7-1\n";
}
#27.5ns
if ($fold_num==142) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new7-2/stripped.nc";
$nb_traj=1;
print "aMD-new7-2\n";
}
#28.5ns
if ($fold_num==143) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new7-3/stripped.nc";
$nb_traj=1;
print "aMD-new7-3\n";
}
#27ns
if ($fold_num==144) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new7-4/stripped.nc";
$nb_traj=1;
print "aMD-new7-4\n";
}

#C'
#101 ns
if ($fold_num==21) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-close-0/stripped.nc";
$nb_traj=1;
print "TL-close-0\n";
}
#100 ns
if ($fold_num==22) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-close-1/stripped.nc";
$nb_traj=1;
print "TL-close-1\n";
}

#D
#13ns
if ($fold_num==49) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-1/stripped.nc";
$nb_traj=1;
print "TL-open1-1\n";
}
#13ns
if ($fold_num==50) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-2/stripped.nc";
$nb_traj=1;
print "TL-open1-2\n";
}
#18ns
if ($fold_num==51) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-4/stripped.nc";
$nb_traj=1;
print "TL-open1-4\n";
}
#28.5ns
if ($fold_num==52) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-5/stripped.nc";
$nb_traj=1;
print "TL-open1-5\n";
}
#28ns
if ($fold_num==53) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open1-6/stripped.nc";
$nb_traj=1;
print "TL-open1-6\n";
}
#16ns
if ($fold_num==54) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-1/stripped.nc";
$nb_traj=1;
print "TL-open2-1\n";
}
#15.5ns
if ($fold_num==55) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-2/stripped.nc";
$nb_traj=1;
print "TL-open2-2\n";
}
#107.5ns 
if ($fold_num==56) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-4/stripped.nc";
$nb_traj=1;
print "TL-open2-4\n";
}
#31.5ns
if ($fold_num==57) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-5/stripped.nc";
$nb_traj=1;
print "TL-open2-5\n";
}
#31.5ns
if ($fold_num==58) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open2-6/stripped.nc";
$nb_traj=1;
print "TL-open2-6\n";
}
#16ns
if ($fold_num==59) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-1/stripped.nc";
$nb_traj=1;
print "TL-open3-1\n";
}
#71ns
if ($fold_num==60) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-2/stripped.nc";
$nb_traj=1;
print "TL-open3-2\n";
}
#19.5ns
if ($fold_num==61) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-3/stripped.nc";
$nb_traj=1;
print "TL-open3-3\n";
}
#38.5ns
if ($fold_num==62) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-5/stripped.nc";
$nb_traj=1;
print "TL-open3-5\n";
}
#27ns
if ($fold_num==63) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/TL-open3-6/stripped.nc";
$nb_traj=1;
print "TL-open3-6\n";
}

#E
#67 ns
if ($fold_num==78) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-1/stripped.nc";
$nb_traj=1;
print "iso-CH3-1\n";
}
#57.5 ns
if ($fold_num==79) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-2/stripped.nc";
$nb_traj=1;
print "iso-CH3-2\n";
}
#75.5 ns
if ($fold_num==80) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-3/stripped.nc";
$nb_traj=1;
print "iso-CH3-3\n";
}
#77 ns
if ($fold_num==81) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-4/stripped.nc";
$nb_traj=1;
print "iso-CH3-4\n";
}
#124.5 ns
if ($fold_num==82) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-5/stripped.nc";
$nb_traj=1;
print "iso-CH3-5\n";
}
#124 ns
if ($fold_num==83) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/iso-CH3-6/stripped.nc";
$nb_traj=1;
print "iso-CH3-6\n";
}

#F
#102.5 ns
if ($fold_num==95) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-1/stripped.nc";
$nb_traj=1;
print "trans-CH3-1\n";
}
#74 ns
if ($fold_num==96) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-2/stripped.nc";
$nb_traj=1;
print "trans-CH3-2\n";
}
#52.5 ns
if ($fold_num==97) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-3/stripped.nc";
$nb_traj=1;
print "trans-CH3-3\n";
} 
#189 ns
if ($fold_num==99) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-5/stripped.nc";
$nb_traj=1;
print "trans-CH3-5\n";
}
#162 ns
if ($fold_num==100) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-6/stripped.nc";
$nb_traj=1;
print "trans-CH3-6\n";
}
#74.5 ns
if ($fold_num==101) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-3-2/stripped.nc";
$nb_traj=1;
print "trans-CH3-3-2\n";
}

#F'
#54.5 ns
if ($fold_num==151) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-1/stripped.nc";
$nb_traj=1;
print "trans-CH3-m4-1\n";
}
#54.5 ns
if ($fold_num==152) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-2/stripped.nc";
$nb_traj=1;
print "trans-CH3-m4-2\n";
}
#54.5 ns
if ($fold_num==153) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-3/stripped.nc";
$nb_traj=1;
print "trans-CH3-m4-3\n";
}
#55.5 ns
if ($fold_num==155) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-5/stripped.nc";
$nb_traj=1;
print "trans-CH3-m4-5\n";
}
#58.5 ns
if ($fold_num==156) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-6/stripped.nc";
$nb_traj=1;
print "trans-CH3-m4-6\n";
}
#55.5 ns
if ($fold_num==157) {
$prmtop= "strip.prmtop";  
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-CH3-m4-7/stripped.nc";
$nb_traj=1;
print "trans-CH3-m4-7\n";
}


#G
#22.2 ns
if ($fold_num==103) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new-1/stripped.nc";
$nb_traj=1;
print "trans-new-1\n";
}
#49.8 ns
if ($fold_num==104) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new-2/stripped.nc";
$nb_traj=1;
print "trans-new-2\n";
}


#H
#53.5 ns
if ($fold_num==109) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-1/stripped.nc";
$nb_traj=1;
print "trans-new5-1\n";
}
#101.8 ns
if ($fold_num==110) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-2/stripped.nc";
$nb_traj=1;
print "trans-new5-2\n";
}
#40.4 ns
if ($fold_num==111) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-3/stripped.nc";
$nb_traj=1;
print "trans-new5-3\n";
}
#95.2 ns
if ($fold_num==113) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-5/stripped.nc";
$nb_traj=1;
print "trans-new5-5\n";
}
#169.6 ns 
if ($fold_num==114) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new5-6/stripped.nc";
$nb_traj=1;
print "trans-new5-6\n";
}

#I
#89 ns
if ($fold_num==116) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new7-3/stripped.nc";
$nb_traj=1;
print "trans-new7-3\n";
}
#156.6 ns
if ($fold_num==117) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new7-5/stripped.nc";
$nb_traj=1;
print "trans-new7-5\n";
}
#75 ns
if ($fold_num==118) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1/stripped.nc";
$nb_traj=1;
print "trans-new8-1\n";
}

#J
#40.8 ns
if ($fold_num==120) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1-3/stripped.nc";
$nb_traj=1;
print "trans-new8-1-3\n";
}
#40.2 ns
if ($fold_num==121) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1-4/stripped.nc";
$nb_traj=1;
print "trans-new8-1-4\n";
}
#12.2 ns
if ($fold_num==122) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1-5/stripped.nc";
$nb_traj=1;
print "trans-new8-1-5\n";
}
#68.8 ns
if ($fold_num==123) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-3/stripped.nc";
$nb_traj=1;
print "trans-new8-3\n";
}
#12.4 ns
if ($fold_num==124) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-1-7/stripped.nc";
$nb_traj=1;
print "trans-new8-1-7\n";
}
#14.2 ns
if ($fold_num==126) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-2/stripped.nc";
$nb_traj=1;
print "trans-new8-2\n";
}
#94.2 ns
if ($fold_num==127) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/iso-CH3-1/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/trans-new8-6/stripped.nc";
$nb_traj=1;
print "trans-new8-6\n";
}


#K
#56ns
if ($fold_num==64) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/stripped.nc";
$nb_traj=1;
print "aMD-new-4-alt-rep3\n";
}
#93.5ns
if ($fold_num==65) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep5/stripped.nc";
$nb_traj=1;
print "aMD-new-4-alt-rep5\n";
}

#L
#59ns
if ($fold_num==72) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep1/stripped.nc";
$nb_traj=1;
print "aMD-new-4-alt3-rep1\n";
}
#67ns
if ($fold_num==73) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep2/stripped.nc";
$nb_traj=1;
print "aMD-new-4-alt3-rep2\n";
}
#61.5ns
if ($fold_num==74) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep3/stripped.nc";
$nb_traj=1;
print "aMD-new-4-alt3-rep3\n";
}
#62ns
if ($fold_num==75) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep4/stripped.nc";
$nb_traj=1;
print "aMD-new-4-alt3-rep4\n";
}
#118.5ns 
if ($fold_num==76) {
$prmtop= "out15-parmed2.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep5/stripped.nc";
$nb_traj=1;
print "aMD-new-4-alt3-rep5\n";
}
#101ns
if ($fold_num==77) {
$prmtop= "strip.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/aMD-new-4-alt-rep3/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/aMD-new-4-alt3-rep6/stripped.nc";
$nb_traj=1;
print "aMD-new-4-alt3-rep6\n";
}

#M
#100ns 
if ($fold_num==137) {
$prmtop= "out15-parmed.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/syner+3-5/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/syner+4-6/stripped.nc";
$nb_traj=1;
print "syner+4-6\n";
}
#100ns
if ($fold_num==136) {
$prmtop= "out15-parmed.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/syner+3-5/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/syner+3-5/stripped.nc";
$nb_traj=1;
print "syner+3-5\n";
}
#100ns
if ($fold_num==135) {
$prmtop= "out15-parmed.prmtop"; 
$abs_filename_prmtop = "/home/ng/Desktop/DATA-sims/syner+3-5/strip.prmtop";
$abs_filename[0]= "/home/ng/Desktop/DATA-sims/syner+2-new/stripped.nc";
$nb_traj=1;
print "syner+2-new\n";
}


my $outfile="visu.vmd";
open (FILE2, "> $outfile") || die;
print (FILE2 "#!/usr/local/bin/vmd\n");
print (FILE2 "# VMD script written by save_state \$Revision: 1.47 \$\n");
print (FILE2 "# VMD version: 1.9.2\n");
print (FILE2 "set viewplist {}\n");
print (FILE2 "set fixedlist {}\n");
print (FILE2 "proc vmdrestoremymaterials {} {\n");
print (FILE2 "  set mlist { Opaque Transparent BrushedMetal Diffuse Ghost Glass1 Glass2 Glass3 Glossy HardPlastic MetallicPastel Steel Translucent Edgy EdgyShiny EdgyGlass Goodsell AOShiny AOChalky AOEdgy BlownGlass GlassBubble RTChrome }\n");
print (FILE2 "  set mymlist [material list]\n");
print (FILE2 "  foreach mat \$mlist {\n");
print (FILE2 "    if { [lsearch \$mymlist \$mat] == -1 } { \n");
print (FILE2 "      material add \$mat\n");
print (FILE2 "    }\n");
print (FILE2 "  }\n");
print (FILE2 "  material change ambient Opaque 0.000000\n");
print (FILE2 "  material change diffuse Opaque 0.650000\n");
print (FILE2 "  material change specular Opaque 0.500000\n");
print (FILE2 "  material change shininess Opaque 0.534020\n");
print (FILE2 "  material change mirror Opaque 0.000000\n");
print (FILE2 "  material change opacity Opaque 1.000000\n");
print (FILE2 "  material change outline Opaque 0.000000\n");
print (FILE2 "  material change outlinewidth Opaque 0.000000\n");
print (FILE2 "  material change transmode Opaque 0.000000\n");
print (FILE2 "  material change ambient Transparent 0.000000\n");
print (FILE2 "  material change diffuse Transparent 0.650000\n");
print (FILE2 "  material change specular Transparent 0.500000\n");
print (FILE2 "  material change shininess Transparent 0.534020\n");
print (FILE2 "  material change mirror Transparent 0.000000\n");
print (FILE2 "  material change opacity Transparent 0.300000\n");
print (FILE2 "  material change outline Transparent 0.000000\n");
print (FILE2 "  material change outlinewidth Transparent 0.000000\n");
print (FILE2 "  material change transmode Transparent 0.000000\n");
print (FILE2 "  material change ambient BrushedMetal 0.080000\n");
print (FILE2 "  material change diffuse BrushedMetal 0.390000\n");
print (FILE2 "  material change specular BrushedMetal 0.340000\n");
print (FILE2 "  material change shininess BrushedMetal 0.150000\n");
print (FILE2 "  material change mirror BrushedMetal 0.000000\n");
print (FILE2 "  material change opacity BrushedMetal 1.000000\n");
print (FILE2 "  material change outline BrushedMetal 0.000000\n");
print (FILE2 "  material change outlinewidth BrushedMetal 0.000000\n");
print (FILE2 "  material change transmode BrushedMetal 0.000000\n");
print (FILE2 "  material change ambient Diffuse 0.000000\n");
print (FILE2 "  material change diffuse Diffuse 0.620000\n");
print (FILE2 "  material change specular Diffuse 0.000000\n");
print (FILE2 "  material change shininess Diffuse 0.530000\n");
print (FILE2 "  material change mirror Diffuse 0.000000\n");
print (FILE2 "  material change opacity Diffuse 1.000000\n");
print (FILE2 "  material change outline Diffuse 0.000000\n");
print (FILE2 "  material change outlinewidth Diffuse 0.000000\n");
print (FILE2 "  material change transmode Diffuse 0.000000\n");
print (FILE2 "  material change ambient Ghost 0.000000\n");
print (FILE2 "  material change diffuse Ghost 0.000000\n");
print (FILE2 "  material change specular Ghost 1.000000\n");
print (FILE2 "  material change shininess Ghost 0.230000\n");
print (FILE2 "  material change mirror Ghost 0.000000\n");
print (FILE2 "  material change opacity Ghost 0.100000\n");
print (FILE2 "  material change outline Ghost 0.000000\n");
print (FILE2 "  material change outlinewidth Ghost 0.000000\n");
print (FILE2 "  material change transmode Ghost 0.000000\n");
print (FILE2 "  material change ambient Glass1 0.000000\n");
print (FILE2 "  material change diffuse Glass1 0.500000\n");
print (FILE2 "  material change specular Glass1 0.650000\n");
print (FILE2 "  material change shininess Glass1 0.530000\n");
print (FILE2 "  material change mirror Glass1 0.000000\n");
print (FILE2 "  material change opacity Glass1 0.150000\n");
print (FILE2 "  material change outline Glass1 0.000000\n");
print (FILE2 "  material change outlinewidth Glass1 0.000000\n");
print (FILE2 "  material change transmode Glass1 0.000000\n");
print (FILE2 "  material change ambient Glass2 0.520000\n");
print (FILE2 "  material change diffuse Glass2 0.760000\n");
print (FILE2 "  material change specular Glass2 0.220000\n");
print (FILE2 "  material change shininess Glass2 0.590000\n");
print (FILE2 "  material change mirror Glass2 0.000000\n");
print (FILE2 "  material change opacity Glass2 0.680000\n");
print (FILE2 "  material change outline Glass2 0.000000\n");
print (FILE2 "  material change outlinewidth Glass2 0.000000\n");
print (FILE2 "  material change transmode Glass2 0.000000\n");
print (FILE2 "  material change ambient Glass3 0.150000\n");
print (FILE2 "  material change diffuse Glass3 0.250000\n");
print (FILE2 "  material change specular Glass3 0.750000\n");
print (FILE2 "  material change shininess Glass3 0.800000\n");
print (FILE2 "  material change mirror Glass3 0.000000\n");
print (FILE2 "  material change opacity Glass3 0.500000\n");
print (FILE2 "  material change outline Glass3 0.000000\n");
print (FILE2 "  material change outlinewidth Glass3 0.000000\n");
print (FILE2 "  material change transmode Glass3 0.000000\n");
print (FILE2 "  material change ambient Glossy 0.000000\n");
print (FILE2 "  material change diffuse Glossy 0.650000\n");
print (FILE2 "  material change specular Glossy 1.000000\n");
print (FILE2 "  material change shininess Glossy 0.880000\n");
print (FILE2 "  material change mirror Glossy 0.000000\n");
print (FILE2 "  material change opacity Glossy 1.000000\n");
print (FILE2 "  material change outline Glossy 0.000000\n");
print (FILE2 "  material change outlinewidth Glossy 0.000000\n");
print (FILE2 "  material change transmode Glossy 0.000000\n");
print (FILE2 "  material change ambient HardPlastic 0.000000\n");
print (FILE2 "  material change diffuse HardPlastic 0.560000\n");
print (FILE2 "  material change specular HardPlastic 0.280000\n");
print (FILE2 "  material change shininess HardPlastic 0.690000\n");
print (FILE2 "  material change mirror HardPlastic 0.000000\n");
print (FILE2 "  material change opacity HardPlastic 1.000000\n");
print (FILE2 "  material change outline HardPlastic 0.000000\n");
print (FILE2 "  material change outlinewidth HardPlastic 0.000000\n");
print (FILE2 "  material change transmode HardPlastic 0.000000\n");
print (FILE2 "  material change ambient MetallicPastel 0.000000\n");
print (FILE2 "  material change diffuse MetallicPastel 0.260000\n");
print (FILE2 "  material change specular MetallicPastel 0.550000\n");
print (FILE2 "  material change shininess MetallicPastel 0.190000\n");
print (FILE2 "  material change mirror MetallicPastel 0.000000\n");
print (FILE2 "  material change opacity MetallicPastel 1.000000\n");
print (FILE2 "  material change outline MetallicPastel 0.000000\n");
print (FILE2 "  material change outlinewidth MetallicPastel 0.000000\n");
print (FILE2 "  material change transmode MetallicPastel 0.000000\n");
print (FILE2 "  material change ambient Steel 0.250000\n");
print (FILE2 "  material change diffuse Steel 0.000000\n");
print (FILE2 "  material change specular Steel 0.380000\n");
print (FILE2 "  material change shininess Steel 0.320000\n");
print (FILE2 "  material change mirror Steel 0.000000\n");
print (FILE2 "  material change opacity Steel 1.000000\n");
print (FILE2 "  material change outline Steel 0.000000\n");
print (FILE2 "  material change outlinewidth Steel 0.000000\n");
print (FILE2 "  material change transmode Steel 0.000000\n");
print (FILE2 "  material change ambient Translucent 0.000000\n");
print (FILE2 "  material change diffuse Translucent 0.700000\n");
print (FILE2 "  material change specular Translucent 0.600000\n");
print (FILE2 "  material change shininess Translucent 0.300000\n");
print (FILE2 "  material change mirror Translucent 0.000000\n");
print (FILE2 "  material change opacity Translucent 0.800000\n");
print (FILE2 "  material change outline Translucent 0.000000\n");
print (FILE2 "  material change outlinewidth Translucent 0.000000\n");
print (FILE2 "  material change transmode Translucent 0.000000\n");
print (FILE2 "  material change ambient Edgy 0.000000\n");
print (FILE2 "  material change diffuse Edgy 0.660000\n");
print (FILE2 "  material change specular Edgy 0.000000\n");
print (FILE2 "  material change shininess Edgy 0.750000\n");
print (FILE2 "  material change mirror Edgy 0.000000\n");
print (FILE2 "  material change opacity Edgy 1.000000\n");
print (FILE2 "  material change outline Edgy 0.620000\n");
print (FILE2 "  material change outlinewidth Edgy 0.940000\n");
print (FILE2 "  material change transmode Edgy 0.000000\n");
print (FILE2 "  material change ambient EdgyShiny 0.000000\n");
print (FILE2 "  material change diffuse EdgyShiny 0.660000\n");
print (FILE2 "  material change specular EdgyShiny 0.960000\n");
print (FILE2 "  material change shininess EdgyShiny 0.750000\n");
print (FILE2 "  material change mirror EdgyShiny 0.000000\n");
print (FILE2 "  material change opacity EdgyShiny 1.000000\n");
print (FILE2 "  material change outline EdgyShiny 0.760000\n");
print (FILE2 "  material change outlinewidth EdgyShiny 0.940000\n");
print (FILE2 "  material change transmode EdgyShiny 0.000000\n");
print (FILE2 "  material change ambient EdgyGlass 0.000000\n");
print (FILE2 "  material change diffuse EdgyGlass 0.660000\n");
print (FILE2 "  material change specular EdgyGlass 0.500000\n");
print (FILE2 "  material change shininess EdgyGlass 0.750000\n");
print (FILE2 "  material change mirror EdgyGlass 0.000000\n");
print (FILE2 "  material change opacity EdgyGlass 0.620000\n");
print (FILE2 "  material change outline EdgyGlass 0.620000\n");
print (FILE2 "  material change outlinewidth EdgyGlass 0.940000\n");
print (FILE2 "  material change transmode EdgyGlass 0.000000\n");
print (FILE2 "  material change ambient Goodsell 0.520000\n");
print (FILE2 "  material change diffuse Goodsell 1.000000\n");
print (FILE2 "  material change specular Goodsell 0.000000\n");
print (FILE2 "  material change shininess Goodsell 0.000000\n");
print (FILE2 "  material change mirror Goodsell 0.000000\n");
print (FILE2 "  material change opacity Goodsell 1.000000\n");
print (FILE2 "  material change outline Goodsell 4.000000\n");
print (FILE2 "  material change outlinewidth Goodsell 0.900000\n");
print (FILE2 "  material change transmode Goodsell 0.000000\n");
print (FILE2 "  material change ambient AOShiny 0.000000\n");
print (FILE2 "  material change diffuse AOShiny 0.850000\n");
print (FILE2 "  material change specular AOShiny 0.200000\n");
print (FILE2 "  material change shininess AOShiny 0.530000\n");
print (FILE2 "  material change mirror AOShiny 0.000000\n");
print (FILE2 "  material change opacity AOShiny 1.000000\n");
print (FILE2 "  material change outline AOShiny 0.000000\n");
print (FILE2 "  material change outlinewidth AOShiny 0.000000\n");
print (FILE2 "  material change transmode AOShiny 0.000000\n");
print (FILE2 "  material change ambient AOChalky 0.000000\n");
print (FILE2 "  material change diffuse AOChalky 0.850000\n");
print (FILE2 "  material change specular AOChalky 0.000000\n");
print (FILE2 "  material change shininess AOChalky 0.530000\n");
print (FILE2 "  material change mirror AOChalky 0.000000\n");
print (FILE2 "  material change opacity AOChalky 1.000000\n");
print (FILE2 "  material change outline AOChalky 0.000000\n");
print (FILE2 "  material change outlinewidth AOChalky 0.000000\n");
print (FILE2 "  material change transmode AOChalky 0.000000\n");
print (FILE2 "  material change ambient AOEdgy 0.000000\n");
print (FILE2 "  material change diffuse AOEdgy 0.900000\n");
print (FILE2 "  material change specular AOEdgy 0.200000\n");
print (FILE2 "  material change shininess AOEdgy 0.530000\n");
print (FILE2 "  material change mirror AOEdgy 0.000000\n");
print (FILE2 "  material change opacity AOEdgy 1.000000\n");
print (FILE2 "  material change outline AOEdgy 0.620000\n");
print (FILE2 "  material change outlinewidth AOEdgy 0.930000\n");
print (FILE2 "  material change transmode AOEdgy 0.000000\n");
print (FILE2 "  material change ambient BlownGlass 0.040000\n");
print (FILE2 "  material change diffuse BlownGlass 0.340000\n");
print (FILE2 "  material change specular BlownGlass 1.000000\n");
print (FILE2 "  material change shininess BlownGlass 1.000000\n");
print (FILE2 "  material change mirror BlownGlass 0.000000\n");
print (FILE2 "  material change opacity BlownGlass 0.100000\n");
print (FILE2 "  material change outline BlownGlass 0.000000\n");
print (FILE2 "  material change outlinewidth BlownGlass 0.000000\n");
print (FILE2 "  material change transmode BlownGlass 1.000000\n");
print (FILE2 "  material change ambient GlassBubble 0.250000\n");
print (FILE2 "  material change diffuse GlassBubble 0.340000\n");
print (FILE2 "  material change specular GlassBubble 1.000000\n");
print (FILE2 "  material change shininess GlassBubble 1.000000\n");
print (FILE2 "  material change mirror GlassBubble 0.000000\n");
print (FILE2 "  material change opacity GlassBubble 0.040000\n");
print (FILE2 "  material change outline GlassBubble 0.000000\n");
print (FILE2 "  material change outlinewidth GlassBubble 0.000000\n");
print (FILE2 "  material change transmode GlassBubble 1.000000\n");
print (FILE2 "  material change ambient RTChrome 0.000000\n");
print (FILE2 "  material change diffuse RTChrome 0.650000\n");
print (FILE2 "  material change specular RTChrome 0.500000\n");
print (FILE2 "  material change shininess RTChrome 0.530000\n");
print (FILE2 "  material change mirror RTChrome 0.700000\n");
print (FILE2 "  material change opacity RTChrome 1.000000\n");
print (FILE2 "  material change outline RTChrome 0.000000\n");
print (FILE2 "  material change outlinewidth RTChrome 0.000000\n");
print (FILE2 "  material change transmode RTChrome 0.000000\n");
print (FILE2 "}\n");
print (FILE2 "vmdrestoremymaterials\n");
# Atom selection macros
print (FILE2 "atomselect macro at {resname ADE A THY T\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro acidic {resname ASP GLU\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro cyclic {resname HIS PHE PRO TRP TYR\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro acyclic {protein and not cyclic\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro aliphatic {resname ALA GLY ILE LEU VAL\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro alpha {protein and name CA\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro amino {protein\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro aromatic {resname HIS PHE TRP TYR\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro basic {resname ARG HIS LYS HSP\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro bonded {numbonds > 0\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro buried {resname ALA LEU VAL ILE PHE CYS MET TRP\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro cg {resname CYT C GUA G\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro charged {basic or acidic\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro hetero {not (protein or nucleic)\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro hydrophobic {resname ALA LEU VAL ILE PRO PHE MET TRP\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro small {resname ALA GLY SER\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro medium {resname VAL THR ASP ASN PRO CYS ASX PCA HYP\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro large {protein and not (small or medium)\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro neutral {resname VAL PHE GLN TYR HIS CYS MET TRP ASX GLX PCA HYP\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro polar {protein and not hydrophobic\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro purine {resname ADE A GUA G\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro pyrimidine {resname CYT C THY T URA U\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro surface {protein and not buried\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro lipid {resname DLPE DMPC DPPC GPC LPPC PALM PC PGCL POPC POPE\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro lipids {lipid\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro ion {resname AL BA CA CAL CD CES CLA CL CO CS CU CU1 CUA HG IN IOD K MG MN3 MO3 MO4 MO5 MO6 NA NAW OC7 PB POT PT RB SOD TB TL WO4 YB ZN ZN1 ZN2\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro ions {ion\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro sugar {resname AGLC\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro solvent {not (protein or sugar or nucleic or lipid)\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro carbon {name \"C.*\" and not ion\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro hydrogen {name \"[0-9]?H.*\"\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro nitrogen {name \"N.*\"\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro oxygen {name \"O.*\"\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro sulfur {name \"S.*\" and not ion\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro noh {not hydrogen\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro heme {resname HEM HEME\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro conformationall {altloc \"\"\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro conformationA {altloc \"\" or altloc \"A\"\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro conformationB {altloc \"\" or altloc \"B\"\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro conformationC {altloc \"\" or altloc \"C\"\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro conformationD {altloc \"\" or altloc \"D\"\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro conformationE {altloc \"\" or altloc \"E\"\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro conformationF {altloc \"\" or altloc \"F\"\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro drude {type DRUD or type LP\n");
print (FILE2 "}\n");
print (FILE2 "atomselect macro unparametrized beta<1\n");
print (FILE2 "atomselect macro addedmolefacture {occupancy 0.8}\n");
# Display settings
print (FILE2 "display eyesep       0.065000\n");
print (FILE2 "display focallength  2.000000\n");
print (FILE2 "display height       6.000000\n");
print (FILE2 "display distance     -2.000000\n");
print (FILE2 "display projection   Perspective\n");
print (FILE2 "display nearclip set 0.500000\n");
print (FILE2 "display farclip  set 10.000000\n");
print (FILE2 "display depthcue   on\n");
print (FILE2 "display cuestart   0.500000\n");
print (FILE2 "display cueend     10.000000\n");
print (FILE2 "display cuestart   0.500000\n");
print (FILE2 "display cueend     10.000000\n");
print (FILE2 "display cuedensity 0.320000\n");
print (FILE2 "display cuemode    Exp2\n");
print (FILE2 "display shadows off\n");
print (FILE2 "display ambientocclusion off\n");
print (FILE2 "display aoambient 0.800000\n");
print (FILE2 "display aodirect 0.300000\n");
print (FILE2 "display dof off\n");
print (FILE2 "display dof_fnumber 64.000000\n");
print (FILE2 "display dof_focaldist 0.700000\n");
print (FILE2 "axes location off\n");
print (FILE2 "light 0 on\n");
print (FILE2 "light 1 on\n");
print (FILE2 "light 2 on\n");
print (FILE2 "light 3 on\n");
print (FILE2 "mol new $abs_filename_prmtop type parm7 first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n");
for (my $x=0; $x<$nb_traj; $x++) {
$v=$abs_filename[$x];
print (FILE2 "mol addfile $v type netcdf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n");
}
print (FILE2 "mol delrep 0 top\n");

#Convert stat to occupancy
if ($trigger_heatmap==1){
for (my $zone=0; $zone<22; $zone++) {
$size = @{ $sel_stat[$zone] };
for (my $x=0; $x<$size; $x++) {
$sel_resid=$sel_res[$zone][$x];
$stat=$sel_stat[$zone][$x];
$occupancy=$stat/100;
$occupancy= sprintf "%.3f", $occupancy;
print (FILE2 "set sel_resid [atomselect top \"resid $sel_resid\"]\n");
print (FILE2 "\$sel_resid set occupancy $occupancy\n");
}
}
}

#visu var
my $n=0;


#protein
print (FILE2 "mol representation Lines 6.000000\n");
print (FILE2 "mol color ColorID 16\n");
print (FILE2 "mol selection {protein and not resid 3718 to 4177}\n");
print (FILE2 "mol material Opaque\n");
print (FILE2 "mol addrep top\n");
print (FILE2 "mol selupdate $n top 0\n");
print (FILE2 "mol colupdate $n top 0\n");
print (FILE2 "mol scaleminmax top $n 0.000000 0.000000\n");
print (FILE2 "mol smoothrep top $n 0\n");
print (FILE2 "mol drawframes top $n {now}\n");
print (FILE2 "mol clipplane center 0 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  0 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 0 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 0 $n top {0}\n");
print (FILE2 "mol clipplane center 1 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  1 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 1 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 1 $n top {0}\n");
print (FILE2 "mol clipplane center 2 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  2 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 2 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 2 $n top {0}\n");
print (FILE2 "mol clipplane center 3 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  3 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 3 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 3 $n top {0}\n");
print (FILE2 "mol clipplane center 4 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  4 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 4 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 4 $n top {0}\n");
print (FILE2 "mol clipplane center 5 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  5 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 5 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 5 $n top {0}\n");

#
$n++;
print (FILE2 "mol representation QuickSurf 1.000000 0.500000 1.000000 1.000000\n");
print (FILE2 "mol color ColorID 2\n");
print (FILE2 "mol selection {nucleic}\n");
print (FILE2 "mol material Opaque\n");
print (FILE2 "mol addrep top\n");
print (FILE2 "mol selupdate $n top 0\n");
print (FILE2 "mol colupdate $n top 0\n");
print (FILE2 "mol scaleminmax top $n 0.000000 0.000000\n");
print (FILE2 "mol smoothrep top $n 0\n");
print (FILE2 "mol drawframes top $n {now}\n");
print (FILE2 "mol clipplane center 0 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  0 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 0 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 0 $n top {0}\n");
print (FILE2 "mol clipplane center 1 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  1 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 1 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 1 $n top {0}\n");
print (FILE2 "mol clipplane center 2 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  2 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 2 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 2 $n top {0}\n");
print (FILE2 "mol clipplane center 3 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  3 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 3 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 3 $n top {0}\n");
print (FILE2 "mol clipplane center 4 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  4 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 4 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 4 $n top {0}\n");
print (FILE2 "mol clipplane center 5 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  5 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 5 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 5 $n top {0}\n");

#
$n++;
print (FILE2 "mol representation QuickSurf 1.000000 0.500000 1.000000 1.000000\n");
print (FILE2 "mol color ColorID 6\n");
print (FILE2 "mol selection {resid 3718 to 4177}\n");
print (FILE2 "mol material Opaque\n");
print (FILE2 "mol addrep top\n");
print (FILE2 "mol selupdate $n top 0\n");
print (FILE2 "mol colupdate $n top 0\n");
print (FILE2 "mol scaleminmax top $n 0.000000 0.000000\n");
print (FILE2 "mol smoothrep top $n 0\n");
print (FILE2 "mol drawframes top $n {now}\n");
print (FILE2 "mol clipplane center 0 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  0 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 0 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 0 $n top {0}\n");
print (FILE2 "mol clipplane center 1 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  1 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 1 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 1 $n top {0}\n");
print (FILE2 "mol clipplane center 2 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  2 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 2 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 2 $n top {0}\n");
print (FILE2 "mol clipplane center 3 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  3 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 3 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 3 $n top {0}\n");
print (FILE2 "mol clipplane center 4 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  4 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 4 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 4 $n top {0}\n");
print (FILE2 "mol clipplane center 5 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  5 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 5 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 5 $n top {0}\n");

if ($trigger_heatmap==1){
for (my $zone=0; $zone<22; $zone++) {
my $sel_incr;
$size = @{ $sel_stat[$zone] };
for (my $x=0; $x<$size; $x++) {
$sel_resid=$sel_res[$zone][$x];
if ($x==0){
$sel_incr=$sel_resid;
}
if ($x!=0){
$sel_incr=$sel_incr . " " . $sel_resid;
}
}
#
$n++;
print (FILE2 "mol representation VDW 1.000000 12.000000\n");
print (FILE2 "mol color Occupancy\n");
print (FILE2 "mol selection {resid $sel_incr}\n");
print (FILE2 "mol material AOEdgy\n");
print (FILE2 "mol addrep top\n");
print (FILE2 "mol selupdate $n top 0\n");
print (FILE2 "mol colupdate $n top 0\n");
print (FILE2 "mol scaleminmax top $n 0.000000 0.746000\n");
print (FILE2 "mol smoothrep top $n 0\n");
print (FILE2 "mol drawframes top $n {now}\n");
print (FILE2 "mol clipplane center 0 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  0 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 0 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 0 $n top {0}\n");
print (FILE2 "mol clipplane center 1 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  1 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 1 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 1 $n top {0}\n");
print (FILE2 "mol clipplane center 2 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  2 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 2 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 2 $n top {0}\n");
print (FILE2 "mol clipplane center 3 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  3 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 3 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 3 $n top {0}\n");
print (FILE2 "mol clipplane center 4 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  4 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 4 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 4 $n top {0}\n");
print (FILE2 "mol clipplane center 5 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  5 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 5 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 5 $n top {0}\n");
}
}


if ($trigger_res==1){
for (my $zone=0; $zone<22; $zone++) {
my $sel_incr;
$colorid=$sel_colorid[$zone];
$size = @{ $sel_res[$zone] };
for (my $x=0; $x<$size; $x++) {
$sel_resid=$sel_res[$zone][$x];
if ($x==0){
$sel_incr=$sel_resid;
}
if ($x!=0){
$sel_incr=$sel_incr . " " . $sel_resid;
}
}
#
$n++;
if ($zone!=3){
print (FILE2 "mol representation VDW 1.000000 12.000000\n");
}
if ($zone==3){
print (FILE2 "mol representation VDW 1.100000 12.000000\n");
}
print (FILE2 "mol color ColorID $colorid\n");
print (FILE2 "mol selection {resid $sel_incr}\n");
if (($zone!=16) and ($zone!=17) and ($zone!=18) and ($zone!=20) and ($zone!=21)){
print (FILE2 "mol material AOEdgy\n");
}
if ($zone==16){
print (FILE2 "mol material Opaque\n");
}
if (($zone==17) or ($zone==18) or ($zone==20) or ($zone==21)){
print (FILE2 "mol material RTChrome\n");
}
print (FILE2 "mol addrep top\n");
print (FILE2 "mol selupdate $n top 0\n");
print (FILE2 "mol colupdate $n top 0\n");
print (FILE2 "mol scaleminmax top $n 0.000000 0.746000\n");
print (FILE2 "mol smoothrep top $n 0\n");
print (FILE2 "mol drawframes top $n {now}\n");
print (FILE2 "mol clipplane center 0 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  0 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 0 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 0 $n top {0}\n");
print (FILE2 "mol clipplane center 1 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  1 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 1 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 1 $n top {0}\n");
print (FILE2 "mol clipplane center 2 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  2 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 2 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 2 $n top {0}\n");
print (FILE2 "mol clipplane center 3 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  3 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 3 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 3 $n top {0}\n");
print (FILE2 "mol clipplane center 4 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  4 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 4 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 4 $n top {0}\n");
print (FILE2 "mol clipplane center 5 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  5 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 5 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 5 $n top {0}\n");
if ($zone==4){
print (FILE2 "mol showrep top $n 0\n");
}
}
}


if ($trigger_ctp==1){
#
$n++;
print (FILE2 "mol representation VDW 1.000000 12.000000\n");
print (FILE2 "mol color Name\n");
print (FILE2 "mol selection {resname ctp}\n");
print (FILE2 "mol material RTChrome\n");
print (FILE2 "mol addrep top\n");
print (FILE2 "mol selupdate $n top 0\n");
print (FILE2 "mol colupdate $n top 0\n");
print (FILE2 "mol scaleminmax top $n 0.000000 0.746000\n");
print (FILE2 "mol smoothrep top $n 0\n");
print (FILE2 "mol drawframes top $n {now}\n");
print (FILE2 "mol clipplane center 0 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  0 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 0 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 0 $n top {0}\n");
print (FILE2 "mol clipplane center 1 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  1 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 1 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 1 $n top {0}\n");
print (FILE2 "mol clipplane center 2 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  2 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 2 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 2 $n top {0}\n");
print (FILE2 "mol clipplane center 3 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  3 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 3 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 3 $n top {0}\n");
print (FILE2 "mol clipplane center 4 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  4 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 4 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 4 $n top {0}\n");
print (FILE2 "mol clipplane center 5 $n top {0.0 0.0 0.0}\n");
print (FILE2 "mol clipplane color  5 $n top {0.5 0.5 0.5 }\n");
print (FILE2 "mol clipplane normal 5 $n top {0.0 0.0 1.0}\n");
print (FILE2 "mol clipplane status 5 $n top {0}\n");
}


print (FILE2 "mol rename top $prmtop\n");
#Set visualization angle and zoom
print (FILE2 "set viewpoints([molinfo top]) {{{1 0 0 -89.398} {0 1 0 -104.763} {0 0 1 -85.4039} {0 0 0 1}} {{0.512352 -0.85857 0.0189647 0} {0.223002 0.111588 -0.968469 0} {0.829328 0.500381 0.248643 0} {0 0 0 1}} {{0.0474999 0 0 0} {0 0.0474999 0 0} {0 0 0.0474999 0} {0 0 0 1}} {{1 0 0 -0.13} {0 1 0 -0.56} {0 0 1 0} {0 0 0 1}}}\n");
#
print (FILE2 "lappend viewplist [molinfo top]\n");
print (FILE2 "set topmol [molinfo top]\n");
# done with molecule 0
print (FILE2 "foreach v \$viewplist {\n");
print (FILE2 "  molinfo \$v set {center_matrix rotate_matrix scale_matrix global_matrix} \$viewpoints(\$v)\n");
print (FILE2 "}\n");
print (FILE2 "foreach v \$fixedlist {\n");
print (FILE2 "  molinfo \$v set fixed 1\n");
print (FILE2 "}\n");
print (FILE2 "unset viewplist\n");
print (FILE2 "unset fixedlist\n");
print (FILE2 "mol top \$topmol\n");
print (FILE2 "unset topmol\n");
print (FILE2 "proc vmdrestoremycolors {} {\n");
print (FILE2 "color scale colors RWB {1.0 0.0 0.0} {1.0 1.0 1.0} {0.0 0.0 1.0}\n");
print (FILE2 "color scale colors BWR {0.0 0.0 1.0} {1.0 1.0 1.0} {1.0 0.0 0.0}\n");
print (FILE2 "color scale colors RGryB {1.0 0.0 0.0} {0.5 0.5 0.5} {0.0 0.0 1.0}\n");
print (FILE2 "color scale colors BGryR {0.0 0.0 1.0} {0.5 0.5 0.5} {1.0 0.0 0.0}\n");
print (FILE2 "color scale colors RGB {1.0 0.0 0.0} {0.0 1.0 0.0} {0.0 0.0 1.0}\n");
print (FILE2 "color scale colors BGR {0.0 0.0 1.0} {0.0 1.0 0.0} {1.0 0.0 0.0}\n");
print (FILE2 "color scale colors RWG {1.0 0.0 0.0} {1.0 1.0 1.0} {0.0 1.0 0.0}\n");
print (FILE2 "color scale colors GWR {0.0 1.0 0.0} {1.0 1.0 1.0} {1.0 0.0 0.0}\n");
print (FILE2 "color scale colors GWB {0.0 1.0 0.0} {1.0 1.0 1.0} {0.0 0.0 1.0}\n");
print (FILE2 "color scale colors BWG {0.0 0.0 1.0} {1.0 1.0 1.0} {0.0 1.0 0.0}\n");
print (FILE2 "color scale colors BlkW {0.0 0.0 0.0} {0.5 0.5 0.5} {1.0 1.0 1.0}\n");
print (FILE2 "color scale colors WBlk {1.0 1.0 1.0} {0.5 0.5 0.5} {0.0 0.0 0.0}\n");
#scale method:
print (FILE2 "  color scale method RGB\n");
#midpoint:
print (FILE2 "  color scale midpoint 0.09\n");
#
print (FILE2 "  set colorcmds {\n");
#white background:
print (FILE2 "    {color Display {Background} white}\n");
#
print (FILE2 "    {color Display {BackgroundTop} black}\n");
print (FILE2 "    {color Display {BackgroundBot} blue2}\n");
print (FILE2 "    {color Display {FPS} white}\n");
print (FILE2 "    {color Name {LPA} green}\n");
print (FILE2 "    {color Name {LPB} green}\n");
print (FILE2 "    {color Name {M} pink}\n");
print (FILE2 "    {color Name {K} cyan}\n");
print (FILE2 "    {color Name {E} purple}\n");
print (FILE2 "    {color Type {LP} green}\n");
print (FILE2 "    {color Type {DRUD} pink}\n");
print (FILE2 "    {color Type {M} pink}\n");
print (FILE2 "    {color Type {K} cyan}\n");
print (FILE2 "    {color Type {E} purple}\n");
print (FILE2 "    {color Element {X} cyan}\n");
print (FILE2 "    {color Element {Ac} ochre}\n");
print (FILE2 "    {color Element {Ag} ochre}\n");
print (FILE2 "    {color Element {Al} ochre}\n");
print (FILE2 "    {color Element {Am} ochre}\n");
print (FILE2 "    {color Element {Ar} ochre}\n");
print (FILE2 "    {color Element {As} ochre}\n");
print (FILE2 "    {color Element {At} ochre}\n");
print (FILE2 "    {color Element {Au} ochre}\n");
print (FILE2 "    {color Element {B} ochre}\n");
print (FILE2 "    {color Element {Ba} ochre}\n");
print (FILE2 "    {color Element {Be} ochre}\n");
print (FILE2 "    {color Element {Bh} ochre}\n");
print (FILE2 "    {color Element {Bi} ochre}\n");
print (FILE2 "    {color Element {Bk} ochre}\n");
print (FILE2 "    {color Element {Br} ochre}\n");
print (FILE2 "    {color Element {Ca} ochre}\n");
print (FILE2 "    {color Element {Cd} ochre}\n");
print (FILE2 "    {color Element {Ce} ochre}\n");
print (FILE2 "    {color Element {Cf} ochre}\n");
print (FILE2 "    {color Element {Cl} ochre}\n");
print (FILE2 "    {color Element {Cm} ochre}\n");
print (FILE2 "    {color Element {Co} ochre}\n");
print (FILE2 "    {color Element {Cr} ochre}\n");
print (FILE2 "    {color Element {Cs} ochre}\n");
print (FILE2 "    {color Element {Cu} ochre}\n");
print (FILE2 "    {color Element {Db} ochre}\n");
print (FILE2 "    {color Element {Ds} ochre}\n");
print (FILE2 "    {color Element {Dy} ochre}\n");
print (FILE2 "    {color Element {Er} ochre}\n");
print (FILE2 "    {color Element {Es} ochre}\n");
print (FILE2 "    {color Element {Eu} ochre}\n");
print (FILE2 "    {color Element {F} ochre}\n");
print (FILE2 "    {color Element {Fe} ochre}\n");
print (FILE2 "    {color Element {Fm} ochre}\n");
print (FILE2 "    {color Element {Fr} ochre}\n");
print (FILE2 "    {color Element {Ga} ochre}\n");
print (FILE2 "    {color Element {Gd} ochre}\n");
print (FILE2 "    {color Element {Ge} ochre}\n");
print (FILE2 "    {color Element {He} ochre}\n");
print (FILE2 "    {color Element {Hf} ochre}\n");
print (FILE2 "    {color Element {Hg} ochre}\n");
print (FILE2 "    {color Element {Ho} ochre}\n");
print (FILE2 "    {color Element {Hs} ochre}\n");
print (FILE2 "    {color Element {I} ochre}\n");
print (FILE2 "    {color Element {In} ochre}\n");
print (FILE2 "    {color Element {Ir} ochre}\n");
print (FILE2 "    {color Element {K} ochre}\n");
print (FILE2 "    {color Element {Kr} ochre}\n");
print (FILE2 "    {color Element {La} ochre}\n");
print (FILE2 "    {color Element {Li} ochre}\n");
print (FILE2 "    {color Element {Lr} ochre}\n");
print (FILE2 "    {color Element {Lu} ochre}\n");
print (FILE2 "    {color Element {Md} ochre}\n");
print (FILE2 "    {color Element {Mg} ochre}\n");
print (FILE2 "    {color Element {Mn} ochre}\n");
print (FILE2 "    {color Element {Mo} ochre}\n");
print (FILE2 "    {color Element {Mt} ochre}\n");
print (FILE2 "    {color Element {Na} ochre}\n");
print (FILE2 "    {color Element {Nb} ochre}\n");
print (FILE2 "    {color Element {Nd} ochre}\n");
print (FILE2 "    {color Element {Ne} ochre}\n");
print (FILE2 "    {color Element {Ni} ochre}\n");
print (FILE2 "    {color Element {No} ochre}\n");
print (FILE2 "    {color Element {Np} ochre}\n");
print (FILE2 "    {color Element {Os} ochre}\n");
print (FILE2 "    {color Element {Pa} ochre}\n");
print (FILE2 "    {color Element {Pb} ochre}\n");
print (FILE2 "    {color Element {Pd} ochre}\n");
print (FILE2 "    {color Element {Pm} ochre}\n");
print (FILE2 "    {color Element {Po} ochre}\n");
print (FILE2 "    {color Element {Pr} ochre}\n");
print (FILE2 "    {color Element {Pt} ochre}\n");
print (FILE2 "    {color Element {Pu} ochre}\n");
print (FILE2 "    {color Element {Ra} ochre}\n");
print (FILE2 "    {color Element {Rb} ochre}\n");
print (FILE2 "    {color Element {Re} ochre}\n");
print (FILE2 "    {color Element {Rf} ochre}\n");
print (FILE2 "    {color Element {Rg} ochre}\n");
print (FILE2 "    {color Element {Rh} ochre}\n");
print (FILE2 "    {color Element {Rn} ochre}\n");
print (FILE2 "    {color Element {Ru} ochre}\n");
print (FILE2 "    {color Element {Sb} ochre}\n");
print (FILE2 "    {color Element {Sc} ochre}\n");
print (FILE2 "    {color Element {Se} ochre}\n");
print (FILE2 "    {color Element {Sg} ochre}\n");
print (FILE2 "    {color Element {Si} ochre}\n");
print (FILE2 "    {color Element {Sm} ochre}\n");
print (FILE2 "    {color Element {Sn} ochre}\n");
print (FILE2 "    {color Element {Sr} ochre}\n");
print (FILE2 "    {color Element {Ta} ochre}\n");
print (FILE2 "    {color Element {Tb} ochre}\n");
print (FILE2 "    {color Element {Tc} ochre}\n");
print (FILE2 "    {color Element {Te} ochre}\n");
print (FILE2 "    {color Element {Th} ochre}\n");
print (FILE2 "    {color Element {Ti} ochre}\n");
print (FILE2 "    {color Element {Tl} ochre}\n");
print (FILE2 "    {color Element {Tm} ochre}\n");
print (FILE2 "    {color Element {U} ochre}\n");
print (FILE2 "    {color Element {V} ochre}\n");
print (FILE2 "    {color Element {W} ochre}\n");
print (FILE2 "    {color Element {Xe} ochre}\n");
print (FILE2 "    {color Element {Y} ochre}\n");
print (FILE2 "    {color Element {Yb} ochre}\n");
print (FILE2 "    {color Element {Zr} ochre}\n");
print (FILE2 "    {color Resname {HIE} silver}\n");
print (FILE2 "    {color Resname {U5} green}\n");
print (FILE2 "    {color Resname {A} white}\n");
print (FILE2 "    {color Resname {U} pink}\n");
print (FILE2 "    {color Resname {C} cyan}\n");
print (FILE2 "    {color Resname {G} purple}\n");
print (FILE2 "    {color Resname {A3} lime}\n");
print (FILE2 "    {color Resname {DG5} mauve}\n");
print (FILE2 "    {color Resname {DA} ochre}\n");
print (FILE2 "    {color Resname {DG} iceblue}\n");
print (FILE2 "    {color Resname {DT} black}\n");
print (FILE2 "    {color Resname {DC} yellow2}\n");
print (FILE2 "    {color Resname {DG3} yellow3}\n");
print (FILE2 "    {color Resname {DC5} green2}\n");
print (FILE2 "    {color Resname {DC3} green3}\n");
print (FILE2 "    {color Resname {MG} cyan2}\n");
print (FILE2 "    {color Resname {CA} cyan3}\n");
print (FILE2 "    {color Resname {SUL} blue2}\n");
print (FILE2 "    {color Resname {Na+} blue3}\n");
print (FILE2 "    {color Resname {ZK} violet}\n");
print (FILE2 "    {color Resname {ZHE} violet2}\n");
print (FILE2 "    {color Resname {ZR} magenta}\n");
print (FILE2 "    {color Resname {ZD} magenta2}\n");
print (FILE2 "    {color Resname {ZE} red2}\n");
print (FILE2 "    {color Resname {K+} red3}\n");
print (FILE2 "    {color Resname {Cl-} orange2}\n");
print (FILE2 "    {color Resname {ctp} blue}\n");
print (FILE2 "    {color Chain {X} blue}\n");
print (FILE2 "    {color Segname {} blue}\n");
print (FILE2 "    {color Conformation {all} blue}\n");
print (FILE2 "    {color Molecule {0} blue}\n");
print (FILE2 "    {color Molecule {$prmtop} blue}\n");
print (FILE2 "   {color Structure {3_10_Helix} blue}\n");
print (FILE2 "    {color Surface {Grasp} gray}\n");
print (FILE2 "    {color Labels {Springs} orange}\n");
print (FILE2 "    {color Stage {Even} gray}\n");
print (FILE2 "    {color Stage {Odd} silver}\n");
print (FILE2 "  }\n");
print (FILE2 "  foreach colcmd \$colorcmds {\n");
print (FILE2 "    set val [catch {eval \$colcmd}]\n");
print (FILE2 "  }\n");
print (FILE2 "  color change rgb 0 0.0 0.0 1.0\n");
print (FILE2 "  color change rgb 2 0.3499999940395355 0.3499999940395355 0.3499999940395355\n");
print (FILE2 "  color change rgb 3 1.0 0.5 0.0\n");
print (FILE2 "  color change rgb 4 1.0 1.0 0.0\n");
print (FILE2 "  color change rgb 5 0.5 0.5 0.20000000298023224\n");
print (FILE2 "  color change rgb 6 0.6000000238418579 0.6000000238418579 0.6000000238418579\n");
print (FILE2 "  color change rgb 7 0.0 1.0 0.0\n");
print (FILE2 "  color change rgb 9 1.0 0.6000000238418579 0.6000000238418579\n");
print (FILE2 "  color change rgb 11 0.6499999761581421 0.0 0.6499999761581421\n");
print (FILE2 "  color change rgb 12 0.5 0.8999999761581421 0.4000000059604645\n");
print (FILE2 "  color change rgb 13 0.8999999761581421 0.4000000059604645 0.699999988079071\n");
print (FILE2 "  color change rgb 14 0.5 0.30000001192092896 0.0\n");
print (FILE2 "  color change rgb 15 0.5 0.5 0.75\n");
print (FILE2 "  color change rgb 17 0.8799999952316284 0.9700000286102295 0.019999999552965164\n");
print (FILE2 "  color change rgb 18 0.550000011920929 0.8999999761581421 0.019999999552965164\n");
print (FILE2 "  color change rgb 19 0.0 0.8999999761581421 0.03999999910593033\n");
print (FILE2 "  color change rgb 20 0.0 0.8999999761581421 0.5\n");
print (FILE2 "  color change rgb 21 0.0 0.8799999952316284 1.0\n");
print (FILE2 "  color change rgb 22 0.0 0.7599999904632568 1.0\n");
print (FILE2 "  color change rgb 23 0.019999999552965164 0.3799999952316284 0.6700000166893005\n");
print (FILE2 "  color change rgb 24 0.009999999776482582 0.03999999910593033 0.9300000071525574\n");
print (FILE2 "  color change rgb 25 0.27000001072883606 0.0 0.9800000190734863\n");
print (FILE2 "  color change rgb 26 0.44999998807907104 0.0 0.8999999761581421\n");
print (FILE2 "  color change rgb 27 0.8999999761581421 0.0 0.8999999761581421\n");
print (FILE2 "  color change rgb 28 1.0 0.0 0.6600000262260437\n");
print (FILE2 "  color change rgb 29 0.9800000190734863 0.0 0.23000000417232513\n");
print (FILE2 "  color change rgb 30 0.8100000023841858 0.0 0.0\n");
print (FILE2 "  color change rgb 31 0.8899999856948853 0.3499999940395355 0.0\n");
print (FILE2 "  color change rgb 32 0.9599999785423279 0.7200000286102295 0.0\n");
print (FILE2 "}\n");
print (FILE2 "vmdrestoremycolors\n");
print (FILE2 "label textsize 1.0\n");

#ALIGN TRAJECTORY
print (FILE2 "set startFrame 0\n");
#Number of frames
print (FILE2 "set nFrames [molinfo top get numframes]\n");
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

close (FILE2);

#LAUNCH VMD
my $cmd = "vmd -e visu.vmd >vmd-out.txt";
system($cmd);
unlink "visu.vmd";

exit;




