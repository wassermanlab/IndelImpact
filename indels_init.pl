#!/bin/env perl

# This is a simple wrapper script to run to initialize indel alternative score computation from a VCF file. 

# Usage:
# i => VCF file with indels
# b => BED-formatted file containing TFBSs (from MANTA)
# m => PFM/PWM matrix file (e.g. from JASPAR 2014)
# o => Final output. 

use strict;
use warnings;

use Getopt::Long;

my $input;
my $bed_file;
my $matrices;
my $out_file;

GetOptions(
	'i=s' => \$input,
	'b=s' => \$bed_file,
	'm=s' => \$matrices,
	'o=s' => \$out_file
);

my $temp_vcf = "vcf_temp_indelrunner.vcf";

print "Overlapping with TFBSs...\n";
system("perl vcf_overlap_bed.pl -v $input -b $bed_file -o $temp_vcf");
print "Calculating alternative scores and alt/ref ratios...\n";
system("perl indels_main.pl -s human -i $temp_vcf -m $matrices -o OUT > $out_file");
print "Removing temporary files...\n";
system("rm $temp_vcf");
print "Done.\n";
exit;
