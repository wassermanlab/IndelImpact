#!/bin/env perl

# Script to overlap VCF input with BED input. Used for computing overlaps between indel (or SNV) entries and TFBSs in corresponding BED format. 

# Usage:
# v => VCF input file
# b => BED-formatted file to overlap with
# o => output file

use strict;
use warnings;

use Getopt::Long;

my $vcf_file;
my $bed_file;
my $out_file;

GetOptions(
	'-v=s' => \$vcf_file,
	'-b=s' => \$bed_file,
	'-o=s' => \$out_file
);

my $temp_bed_insertions = "temp_insertions.bed";
my $temp_bed_noninsert = "temp_noninsertions.bed";
my $temp_out_file = "temp_bedtools_overlap.bed";

## Convert VCF to BED and split into insertions and noninsertions (typically deletions) file. 
open(TEMPINS, ">$temp_bed_insertions");
open(TEMPNONINS, ">$temp_bed_noninsert");
open(VCF, $vcf_file);
while (my $line = <VCF>) {
	chomp $line;
	next unless $line;
	next if ($line =~ /^#/);
	my @info = split "\t", $line;
	(my $chrom = $info[0]) =~ s/^chr//;
	my $pos = $info[1];
	my $ref_str = $info[3];
	my $alt_str = $info[4];
	## pos-1 because BED is zero-indexed. 
	if (length $alt_str > length $ref_str) {
		print TEMPINS $chrom,"\t",$pos-1,"\t",$pos+1,"\t$ref_str\t$alt_str\n"; ## other fields to be added
	} else {
		print TEMPNONINS $chrom,"\t",$pos-1,"\t",$pos+(length $ref_str),"\t$ref_str\t$alt_str\n";
	}
}
close(VCF);

## Opening here for the sheer purpose of clearing the file if it already exists. 
open(TEMPOUT, ">$temp_out_file");
## Require overlap of 100% for insertions to ensure that it actually overlaps the insertion
`bedtools intersect -a $temp_bed_insertions -b $bed_file -wa -wb -f 1 > $temp_out_file`;
`bedtools intersect -a $temp_bed_noninsert -b $bed_file -wa -wb >> $temp_out_file`;
close(TEMPOUT);

## Convert back to VCF format and output back to out_file. 
open(TEMPOUTREAD, $temp_out_file);
open(OUT, ">$out_file");
while (my $line = <TEMPOUTREAD>) {
	chomp $line;
	next unless $line;
	my @info = split "\t", $line;
	my $chrom = "chr".$info[0];
	my $start_pos = $info[1]+1;
	my $ref_str = $info[3];
	my $alt_str = $info[4];
	my $tfbs_start = $info[6]+1;
	my $tfbs_end = $info[7];
	my $tf_id = $info[8];
	my $scores = $info[9];
	my $strand = $info[10];
	my @score_arr = split ";", $scores;
	print OUT $chrom,"\t",$start_pos,"\t.\t",$ref_str,"\t",$alt_str,"\t",$tfbs_start,"\t",$tfbs_end,"\t",$tf_id,"\t",$score_arr[0],"\t",$score_arr[1]*0.01,"\t",$strand,"\n";
	## STUB: Use the structure of temp_out_file to find which columns to output to the VCF to provide all of the properties needed 
	## for the TFBS, etc.
}
close(OUT);
close(TEMPOUTREAD);

# Delete temporary files.
`rm $temp_out_file`;
`rm $temp_bed_insertions`;
`rm $temp_bed_noninsert`;
exit;
