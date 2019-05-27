use File::Slurp;
#use strict;
use autodie;
use warnings qw(all);
use File::Copy qw(copy);
use Math::Complex;
use Cwd;

# clear the screen
system 'clear';
print "\n";

#This script associate the letter amino acid code to each amino acid index being processed for ligand-binding
#and sorts by descending order their propensity score to bind the ligand.

#The string 'ZONE' identifies the beginning of the amino acid list in output-global-polar-final.txt
#The character X is used as a trigger to separate specific protein zones
#In which case, 'X' is to be optionally written by the user in the output-global-polar-final.txt
#between each zone (in a new line between the zones)
#The output file of step 1 is to be saved as output-global-polar-final2.txt
#Otherwise, just enter 'X' after the last amino acid line in output-global-polar-final2.txt

#Variable Declaration
my @pdb_input = read_file("output-global-polar-final2.txt") or die;
my $array_size=scalar @pdb_input;
my @line_handle;
my $resid;
my $stat_sum;
my $stat;
my $sel;
my $trigger=0;
my $zone=1;
my @table;
my $resname;
my $resid_check;
my $letter_code;
my $done=0;
my $jump=0;
my $sel_table=0;
my @res;
my @table_amended_1;

#EXTRACT INFO
for (my $count=0; $count<$array_size; $count++) {
@line_handle = split ( /\s+/, @pdb_input[$count] );
$resid= @line_handle[0];
$stat= @line_handle[1];

if ($resid eq 'X') {
$trigger = 0;
push(@table,('X X'));
push(@res,('9999'));
}

if ($trigger == 1) {
$sel=$resid . " " . $stat;
push(@table,($sel));
push(@res,($resid));
}
if ($resid eq 'ZONE') {
$trigger=1;
} 

}


#ASSOCIATE RESID TO RESNAME
my @pdb_input_ref = read_file("Example-Protein-Coordinates.pdb") or die;
my $array_size_ref=scalar @pdb_input_ref;
$size=scalar @res;
for (my $x=0; $x<$size; $x++) {
$resid_check=$res[$x];
$sel_table=$table[$x];
$done=0;

if ($resid_check==9999){
push(@table_amended_1,('X X'));
$jump=1;
}

if ($jump==0){
for (my $count=1; $count<$array_size_ref; $count++) {
$resid=substr($pdb_input_ref[$count], 22, 4);
$resid=$resid+0;
if (($resid == $resid_check) and ($done==0)){
$resname=substr($pdb_input_ref[$count], 17, 3);
if ($resname eq 'ALA'){
$letter_code = 'A';
}
elsif ($resname eq 'ASN'){
$letter_code = 'N';
}
elsif ($resname eq 'ASP'){
$letter_code = 'D';
}
elsif ($resname eq 'ARG'){
$letter_code = 'R';
}
elsif ($resname eq 'CYS'){
$letter_code = 'C';
}
elsif ($resname eq 'GLN'){
$letter_code = 'Q';
}
elsif ($resname eq 'GLU'){
$letter_code = 'E';
}
elsif ($resname eq 'GLY'){
$letter_code = 'G';
}
elsif ($resname eq 'HIE'){
$letter_code = 'H';
}
elsif ($resname eq 'HIP'){
$letter_code = 'HP';
}
elsif ($resname eq 'ILE'){
$letter_code = 'I';
}
elsif ($resname eq 'LEU'){
$letter_code = 'L';
}
elsif ($resname eq 'LYS'){
$letter_code = 'K';
}
elsif ($resname eq 'MET'){
$letter_code = 'M';
}
elsif ($resname eq 'PHE'){
$letter_code = 'F';
}
elsif ($resname eq 'PRO'){
$letter_code = 'P';
}
elsif ($resname eq 'SER'){
$letter_code = 'S';
}
elsif ($resname eq 'THR'){
$letter_code = 'T';
}
elsif ($resname eq 'TRP'){
$letter_code = 'W';
}
elsif ($resname eq 'TYR'){
$letter_code = 'Y';
}
elsif ($resname eq 'VAL'){
$letter_code = 'V';
}
else {
print "ERROR RESNAME NOT FOUND line_pdb **$count** **$resid** **$resid_check** rank=**$x**\n";
exit;
}
$done=1;
}

}
#end line loop
$sel=$letter_code . $sel_table;
push(@table_amended_1,($sel));
}
#end jump
$jump=0;

}


#SORT by STAT
my @zone;
my $rank=1;
my $prev_zone_rank=0;
my $size2;
print "\n";
print "SORT\n";
print "\n";
$size=scalar @table_amended_1;
for (my $y=0; $y<$size; $y++) {
$resid=$res[$y];

if($resid == 9999){
my @zone;

for (my $z=$prev_zone_rank; $z<$y; $z++) {
$prev_zone_rank=$y+1;
push(@zone,($table_amended_1[$z]));
}
my @new_zone = sort { (split(' ',$b))[1] <=> (split(' ',$a))[1]} @zone;
$size2=scalar @new_zone;
print "ZONE $rank\n";
for (my $w=0; $w<$size2; $w++) {
print "@new_zone[$w]\n";
}

print "\n";
$rank++;
}


}



exit;






