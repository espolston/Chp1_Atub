#!perl -w

use strict;
use warnings;

my $input_fam = shift(@ARGV) or die; #incomplete fam file produced by 1_merged.sbatch 
my $pheno_info = shift(@ARGV) or die; #sex code and phenotype value info 
my $output = shift(@ARGV) or die; #output fam file
#my $dat = shift(@ARGV) or die; #output dat file for within
unlink(qq{$output});

#perl fam_edit.pl /scratch/midway3/espolston/binary_commongarden.fam /scratch/midway3/espolston/Common_Garden_Seq_Stats-moremeta.tsv /scratch/midway3/espolston/commongarden_final.fam 
#Need to take 1_merged.sbatch .fam files output and edit them 
#col 1: family code (FID) (P#_Ag or cg_ID)
#col 2: indiv code (IID)
#col 3: pair/cluster 
#col 4: zeros
#col 5: sex code ('1' = male, '2' = female, '0' = unknown) from google sheets files
#col 6: phenotype value ('1' = control/natural, '2' = case/agricultural)


#Leaving col 1-4 unchanged

open (A, "<", $input_fam) 
	or die "Could not open input " . $input_fam;
	
open (B, "<", $pheno_info) 
	or die "Could not open input " . $pheno_info;

#Open pheno info file and match individual ids to the sex and environment code	
my %indiv = ();
my %dats = ();

my $iid_col = 0;
my $sex_col = 3;
my $env_col = 4;
my $pair_col = 2;

while(my $ph_line = <B>){
	chomp $ph_line;
	#Dont actually need to skip first line here!
	last unless $ph_line =~ /\S/;
	#print "line 42 " . $ph_line . "\n";
	my @b = split(' ', $ph_line);
	
	if($b[0] ne "Sample"){
		my $IID = $b[$iid_col];
		my $envir_code = $b[$env_col];
		my $pair = $b[$pair_col];
		my $sex_code = $b[$sex_col];
		
		#print "line 52 (iid, sex, pair, envircode)" . $IID . "\t" . $sex_code . "\t" . $pair . "\t" . $envir_code ."\n";
		
		#removed envir to envir_code and sex to sex code section because new pheno data not Ag and Nat
		
		$indiv{$IID} = $sex_code . " " . $envir_code . " " . $pair;
		#$dats{$IID} = $IID . " " . $IID . " " . $pair;  #FIDs in the first column, IIDs in the second column, and cluster names in the third column
	
	}#else{
	#	$iid_col = grep{/Sample/} @b;
	#	$sex_col = grep{/Sex/} @b;
	#	$env_col = grep{/Env/} @b;
	#	$pair_col = grep{/Pair/} @b;
	#}
}
close B; 

open (C, ">", $output)
	or die "Could not open output " . $output;
#open(D, ">", $dat)
#	or die "Could not open dat output " . $dat;

#Open generic fam file and update sex and environment code based on pheno info file
while(my $line = <A>){
	chomp $line;
	#print "line 87 " . $line . "\n";
	my @a = split(/\s+/, $line);
	my $col_1 = $a[0];
	my $col_2 = $a[1];
	my $col_3 = $a[2];
	my $col_4 = $a[3];
	my $col_5 = $a[4];
	my $col_6 = $a[5];
	
	
	if($col_5 != 0 || $col_6 != -9){
		print "Col 5 or 6 overwriting previous update";
		die "Fam file line: " . $line;
	}elsif($col_1 =~ /P/){ #doing same thing but for drought
		my @drought_set = split(/\_/, $col_1);
		$col_1 = $drought_set[0] . "_" . $drought_set[1];
		$col_2 = $drought_set[2];
		$col_5 = 0;
		my $envir = $drought_set[1];
		my @p = split('P', $drought_set[0]);
		$col_3 = $p[1];
		my $envir_code = 0;
		
		if($envir eq "Ag" || $envir eq "Agricultural"){
			$envir_code = 2;
		}elsif($envir eq "Nat" || $envir eq "Natural"){
			$envir_code = 1;
		}elsif($envir eq "NA"){
			$envir_code = -9;
		}else{
			print "Environmental codes not in correct format";
			die "fam file line: " . $line;
		}
		
		$col_6 = $envir_code;
		
		#if(!exists($dats{$col_1})){
		#	$dats{$col_1} = $col_1 . " " . $col_1. " " . $pair;
		#	print D $col_1 . "\t" . $col_1 . "\t" . $pair . "\n";
		#}else{
		#	die "IID repeated for dat file w drought" . $col_1;
		#}
		
	}else{
		if(exists($indiv{$col_1})){
			my @vals = split(/\s+/, $indiv{$col_1});
			$col_5 = $vals[0];
			$col_6 = $vals[1];
			$col_3 = $vals[2];
			#if(exists($dats{$col_1})){
			#	my @vals_dat = split(/\s+/,  $dats{$col_1});
			#	print D $vals_dat[0] . "\t" . $vals_dat[1] . "\t" . $vals_dat[2] . "\n";
			#}else{
			#	die "IID not found for dat file" . $col_1;
			#}
		}else{
			die "IID not found for fam file" . $col_1;
		}
	}
	
	print C $col_1 . "\t" . $col_2 . "\t" . $col_3 . "\t" . $col_4 . "\t" . $col_5 . "\t" . $col_6 . "\n"; 
}
close A;
close C; 
#close D;
