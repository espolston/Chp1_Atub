#!perl -w

use strict;
use warnings;

my $sample_file = shift(@ARGV) or die; #indivs from vcf
my $input_fam = shift(@ARGV) or die; #fam file produced by fam_edit_merged.pl
my $output = shift(@ARGV) or die; #output fam file
unlink(qq{$output});

#perl pixy_pop.pl /scratch/midway3/espolston/sample_ids.txt  /scratch/midway3/espolston/vcf/dcgm.fam /scratch/midway3/espolston/dcg_popfile.txt
#col 1: family code (FID)
#col 2: indiv code (IID)
#col 3: pair/cluster 
#col 4: zeros
#col 5: sex code ('1' = male, '2' = female, '0' = unknown) from google sheets files
#col 6: phenotype value ('1' = control/natural, '2' = case/agricultural)

open (A, "<", $input_fam) 
	or die "Could not open input " . $input_fam;

#Open match individual ids from vcf to the environment code	in fam file
my %indiv = ();
my %dats = ();

my $fid_col = 0;
my $iid_col = 1;
my $env_col = 5;

#get indiv and enviros from fam file
while(my $ph_line = <A>){
	chomp $ph_line;
	last unless $ph_line =~ /\S/;
	#print "Fam line:" . $ph_line . "\n";
	my @a = split(' ', $ph_line);

	my $IID = $a[$iid_col];
	my $FID = $a[$fid_col];
	my $envir = $a[$env_col];
	my $col_1 = "Blank";
	my $envir_code = "Blank";
	
	if($FID =~ /P/){
		$col_1 = $FID . "_" . $IID . "_T";
		#print "COL 1 for drought samples worked" . $col_1;
	}else{
		$col_1 = $FID;
	}
		
	if($envir == 2){
		$envir_code = "AG";
	}elsif($envir == 1){
		$envir_code = "NAT";
	}else{
		print "Environmental codes not in correct format";
		die "Fam file line error " . $ph_line;
	}
	
	print "(indiviudal, envircode)" . $col_1 . "\t" . $envir_code ."\n";

	$indiv{$col_1} = $col_1 . " " . $envir_code;
	
}
close A; 

open (B, "<", $sample_file) 
	or die "Could not open input " . $sample_file;

open (C, ">", $output)
	or die "Could not open output " . $output;

#Open vcf header for order
while(my $line = <B>){
	chomp $line;
	my $individual = "Blank";
	my $environment = "Blank";
	
	if(exists($indiv{$line})){
		my @vals = split(/\s+/, $indiv{$line});
		$individual = $vals[0];
		$environment = $vals[1];
	}else{
		die "order vcf and output error";
	}
	
	print C $individual . "\t" . $environment . "\n"; 
	print "adding to output " . $individual . "\t" . $environment . "\n";
}
close B;
close C; 
