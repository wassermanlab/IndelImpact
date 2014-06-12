#!/bin/env perl

=head1 NAME

indels_main.pl

=head1 SYNOPSIS

  indels_main
            -s species
            -i in_file
            [-f in_format]
            ([-m matrix_file] [-pwm])
            | (([-db jaspar_db]
                [-c collection] [-t tax_group(s)] [-ic min_ic])
                | ([-id matrix_ids] | [-n matrix_names]))
            [-th threshold]
            [-d score_diff]
            [-o out_file]

=head1 OPTIONS

  Options may be abbreviation where unique. Where an option may take
  multiple values, the values may be specified as a comma separated string
  or as multiple options switches or some combination thereof, e.g. the
  multiple values for the tax groups (-t) option may be specified as any of
  the following:
      -t "vertebrates, insects"
      -t vertebrates -t insects
      -t vertebrates -t "insects, nematodes"
 
  -s species        = Name of species.
  -i in_file        = Input variant file to analyze. This may be a VCF or
                      a PolyPhen annotation file.
  -f in_format      = Format of input variation file. Either VCF or
                      TSD:
                        default = VCF
  -m matrix_file    = Optional file name of the TFBS profile matrices with
                      which to search. This options overrides the JASPAR DB
                      options below.
  -pwm              = If specified, and the -m option is provided, then the
                      associated file contains PWMs rather than PFMs.
  -db               = Optional JASPAR database name.
                          default = 'JASPAR_2010'
  -id matrix_ids    = Optional JASPAR ID(s) of specific TFBS profile(s) to
                      use. May specify multiple values (see above).
                          e.g. -id 'MA0046'
  -n matrix_names   = Optional JASPAR name(s) of specific TFBS profile(s) to
                      use. May specify multiple values (see above).
                          e.g. -n 'HNF1A'
  -c collection     = Optional JASPAR collection of TFBS profiles with which
                      to search the sequence(s)
                          default = 'CORE'
  -t tax_groups     = Optional taxonomic supergroups(s) of TFBS profiles to
                      use. May specify multiple values (see above).
                          e.g. -t "vertebrates,insects,nematodes"
                          default = 'vertebrates'
  -ic min_ic        = Optional minimum information content (specificity of
                      TFBS profile matrices to use
  -th threshold     = Optional minimum score of putative TFBS hit to report.
                      Specify as an absolute score, i.e. 14.1 or as a
                      relative score, i.e. '75%'.
                          default = '80%'
  -d score_diff     = Optional minimum score difference of TFBS affected by
                      variations.
                          default = '0%'
  -o out_dir        = Optional directory for output files, current repertory 
  					  by default.
  -p				= Optional Calculation method Using p-values instead of 
  					  relative scores.	

=head1 DESCRIPTION

For each variation in the input VCF file, compute which TFBSs are affected
by the variation and output the results. TFBS profile(s) may be specified
either by a matrix file (PFM or PWM), or read from a JASPAR database. If read
from a database a specific profile may be specified by either a matrix ID or
matrix name. Otherwise a set of matrices may be specified using a combination
of the collection name, minimum infomation content and/or tax groups.

=head1 ALGORITHM

For each variation fead from the input VCF file, fetch the sequence surrounding
the variation from Ensembl. For each of the TFBS profiles scan each of the
sequence alleles and score the predicted binding site matrix at all positions
which overlap the variation. Determine the highest scoring binding site for each
of the alleles. If one of the allele binding sites scores above the TFBS score
threshold and the difference between the scores for the different alleles is
greater than or equal to the specified score difference, report these TFBSs
in an output file.

=cut

BEGIN {
  push @INC, ('/raid2/local/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi', '/raid2/local/lib/perl5/site_perl/5.8.8', '/raid2/local/lib/perl5', '/raid2/local/lib/perl5/site_perl', '/raid2/local/lib64/perl5/vendor_perl/5.8.8/x86_64-linux-thread-multi', '/raid2/local/lib/perl5/vendor_perl/5.8.8', '/raid2/local/lib/perl5/vendor_perl', '/raid2/local/lib64/perl5/5.8.8/x86_64-linux-thread-multi', '/raid2/local/lib/perl5/5.8.8', '/ssd1/apps/oPOSSUM3/cgi-bin/', '/usr/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi/', '/home/azhang/lib/', '/raid2/local/src/ensembl-64/ensembl/modules/','/raid2/local/bin');
}

use warnings;
use strict;
use threads;
use lib '/usr/local/src/ensembl-current/ensembl/modules';

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Registry;
use TFBS::Matrix::PFM;
use TFBS::DB::JASPAR5;
use Bio::SeqIO;
use Statistics::R;
use POSIX;


use constant VAR_FILE_FORMAT    => 'VCF';

use constant COLLECTION         => "CORE";
use constant TAX_GROUPS         => ("vertebrates");

use constant MINIMAL_TFBS_THRESHOLD => "10%";
use constant MINIMAL_TFBS_THRESHOLD_NUM => "0.10";
use constant TFBS_THRESHOLD         => "85%";
use constant TFBS_SCORE_DIFF        => "0%";

use constant JASPAR_DB_NAME     => "JASPAR_2014";
use constant JASPAR_DB_HOST     => "vm5.cmmt.ubc.ca";
use constant JASPAR_DB_USER     => "jaspar_r";
use constant JASPAR_DB_PASS     => "";

use constant ENSEMBL_DB_HOST    => "vm2.cmmt.ubc.ca";
use constant ENSEMBL_DB_USER    => "ensembl_r";
use constant ENSEMBL_DB_PASS    => "";

use Parallel::Iterator qw/iterate_as_array/;

my $species;
my $in_file;
my $in_format;
my $matrix_file;
my $jaspar_db;
my @matrix_ids;
my @matrix_names;
my $collection;
my @tax_groups;
my $min_ic;
my $threshold_str;
my $pwm_flag;
my $min_score_diff_str;
my $out_dir;
my $p_val_method;
GetOptions(
	   's=s'  => \$species,
	   'i=s'  => \$in_file,
	   'f=s'  => \$in_format,
	   'm=s'  => \$matrix_file,
	   'pwm'  => \$pwm_flag,
	   'db=s' => \$jaspar_db,
	   'id=s' => \@matrix_ids,
	   'n=s'  => \@matrix_names,
	   'c=s'  => \$collection,
	   't=s'  => \@tax_groups,
	   'ic=i' => \$min_ic,
	   'th=s' => \$threshold_str,
	   'd=s'  => \$min_score_diff_str,
	   'o=s'  => \$out_dir,
	   'p'    => \$p_val_method,
	  );

unless ($out_dir){
  $out_dir = ".";	
}
unless (-e $out_dir){
  mkdir "$out_dir";	
}
# 
# Three directories used to store/pass information between threads. Do not use these directories 
# for other purposes. These directories may all be made children of another temporary directory, 
# which can be specified by the user (has yet to be implemented). Consider using a database instead. 
# SQLite is a favourable option, for its speed compared to MySQL. 
#
my $out_dir_scoring = "tempScores";
my $temp_matrix_dir = "tempMatrix";
my $temp_scores_dir = "allScoresDir";
unless (-e $out_dir_scoring) {
  mkdir "$out_dir_scoring";
}
unless (-e $temp_matrix_dir) {
  mkdir "$temp_matrix_dir";
}
unless (-e $temp_scores_dir) {
  mkdir "$temp_scores_dir";
}
unless ($in_file) {
  pod2usage(
	    -verbose => 1,
	    -msg     => "No input variation file specified\n"
	   );
}

unless ($in_format) {
  $in_format = VAR_FILE_FORMAT;
}

$in_format = uc $in_format;

unless ($in_format eq 'VCF' || $in_format eq 'PPH') {
  pod2usage(
	    -verbose => 1,
	    -msg     => "Input file format should be one of VCF or PPH\n"
	   );
}

unless ($species) {
  pod2usage(
	    -verbose => 1,
	    -msg     => "No species name specified\n"
	   );
}

if ($matrix_file) {
  if (@matrix_names || @matrix_ids || $collection || $min_ic || @tax_groups) {
    pod2usage(
	      -verbose => 1,
	      -msg     =>
	      "Please specify EITHER a matrix file name (-m) OR one"
              . " or more of matrix_ids (-id), matrix_names (-n),"
              . " collection (-c), min_ic (-ic), tax_groups (-t)\n"
	     );
  }
} else {
  $jaspar_db  = JASPAR_DB_NAME unless $jaspar_db;

  if (@matrix_names || @matrix_ids) {
    if (@matrix_names && @matrix_ids) {
      pod2usage(
                -verbose => 1,
                -msg     =>
		"Please specify only one of matrix_ids (-id) or"
		. " matrix_names (-n)\n"
	       );
    }

    if ($collection || $min_ic || @tax_groups) {
      pod2usage(
                -verbose => 1,
                -msg     =>
		"Please specify EITHER specific matrices (-id/-n) OR a"
		. " set matrices with some combination of collection (-c),"
		. " min_ic (-ic) and/or tax_groups (-t)\n"
	       );
    }

    if (@matrix_ids) {
      @matrix_ids = split(/\s*,\s*/, join(',', @matrix_ids));
    }

    if (@matrix_names) {
      @matrix_names = split(/\s*,\s*/, join(',', @matrix_names));
    }
  } else {
    $collection = COLLECTION unless $collection;

    @tax_groups = TAX_GROUPS unless @tax_groups;

    if (@tax_groups) {
      @tax_groups = split(/\s*,\s*/, join(',', @tax_groups));
    }
  }
}

$threshold_str      = TFBS_THRESHOLD unless $threshold_str;
$min_score_diff_str = TFBS_SCORE_DIFF unless $min_score_diff_str;

my $threshold;
my $threshold_type = 'absolute';
if ($threshold_str =~ /(\S+)%/) {
  $threshold = $1 / 100;
  $threshold_type = 'relative';
} else {
  $threshold = $threshold_str;
}

my $min_score_diff;
my $score_diff_type = 'absolute';
if ($min_score_diff_str =~ /(\S+)%/) {
  $min_score_diff = $1 / 100;
  $score_diff_type = 'relative';
} else {
  $min_score_diff = $min_score_diff_str;
}

#print "\nFetching matrix set...\n";

my $matrix_set = get_matrix_set();
unless ($matrix_set && $matrix_set->size > 0) {
  die "Error fetching matrix set\n";
}

my $max_profile_width = max_profile_width($matrix_set);

#print "\nChecking VCF file...\n";

if ($in_format eq 'VCF') {
  check_VCF($in_file);
} else {
  die "PPH search option not implemented yet!\n";
}

#print "Searching for TFBSs affected by variations...\n";

if ($in_format eq 'VCF') {
  vcf_search_tfbss($in_file, $matrix_set);
} else {
  #pph_search_tfbss($in_file, $matrix_set, $out_file);
  die "PPH search option not implemented yet!\n";
}

#print "Done.\n\n";

exit;
# 
# Load the ensembl registry. 
#
sub load_ensembl_registry
  {
    Bio::EnsEMBL::Registry->load_registry_from_db(
						  -host    => ENSEMBL_DB_HOST,
						  -user    => ENSEMBL_DB_USER,
						  -pass    => ENSEMBL_DB_PASS
						 );
  }

#
# Process VCF file.
#
# For each line in the file, extract the corresponding reference sequence with
# enough flanking nucleotides to search all TFBS matrices for each position
# where the matrix overlaps the variant nucleotides. Scan both the reference
# sequence and the variant sequence. In cases where the matrix scores above
# threshold for one sequence and below threshold for the other and the score
# difference is >= than the minimum score difference, write this information
# to the output file.
# 
# Keep in mind that the min_score_diff, as well as chromosome M positions, are not
# supported yet. Chromosome M positions induce errors upon sequence 
# extraction. 
#
sub vcf_search_tfbss
  {
    my ($file, $matrix_set) = @_;
        
    load_ensembl_registry();
    my $slice_adp = Bio::EnsEMBL::Registry->get_adaptor(
							$species, 'core', 'slice'
						       );
    unless ($slice_adp) {
      die "Error getting Ensembl SliceAdaptor\n";
    }

    open(FH, $file) || die "Error opening input VCF file $file\n";


    my %tab_scores;
    my %tab_info;

    my $miter = $matrix_set->Iterator();
    # 
    # Precomputation before iterating through the file. Necessary to 
    # prevent unnecessary conditional checks within the critical portion, 
    # as well as to supply parameters for all_scores_two. 
    #
    my ($matrixfilestr, $matrixIDstr, $startsitestr, $endsitestr) = precompute_parameters($miter);
    #
    # Splits into array now, to avoid doing this computation later. Note that this should not, 
    # and cannot, be done without making a copy of the array later, as scalar values are preserved
    # for each worker using the Storable module. 
    #
    my @startsitestrarr = split ",", $startsitestr;
    my @endsitestrarr = split ",", $endsitestr;
    # 
    # Iterate over all the lines in the file, producing the reference and alternative
    # scores. Meanwhile, batch + submit jobs to the cluster compute nodes. Note that 
    # batching is necessary due to the fixed time interval after which the cluster can 
    # dequeue jobs off of the wait queue. 
    my $line_num = 0;
    ## Print header.
    print "Chr\tPosition\tRef\tAlt\tID\tRefMatrix\trTFBSstart\trTFBSend\trTFBSscore\trTFBSrscore\trTFBSstrand\t";
    print "AltMatrix\taTFBSscore\taTFBSrscore\tARratio\taTFBSstrand\taTFBSstart\taTFBSend\n";
    while (my $line = <FH>) {
	++$line_num;
	if ($line_num % 100 == 0) {print $line_num,"\n";}
	chomp $line;
	next unless $line;
	next if $line =~ /^#/;

	## Load attributes.
	my @elem = split /\t/, $line;
	my $chrom    = $elem[0];
	my $position = $elem[1];
	my $ref  = uc $elem[3];
	my $alt_str  = uc $elem[4]; ## only one alternative allele
	my $tfbs_start = $elem[5];
	my $tfbs_end = $elem[6];
	my $tf_id = $elem[7];
	my $ref_abs_score = $elem[8];
	my $ref_rel_score = $elem[9];
	my $ref_strand = $elem[10];
	my $ref_len = length $ref;
	if ($chrom =~ /^chr(\S+)/) {
	    $chrom = $1;
	}
	my $alt = $alt_str;
	my $rel_position = 2*$max_profile_width + length($alt);
	my $seq_start = $position - $rel_position + 1;
	my $seq_end   = $position + $rel_position - 1;
	my $slice = $slice_adp->fetch_by_region("chromosome", $chrom, $seq_start, $seq_end );
	unless ($slice) {
	    print STDERR
		"ERROR line $line_num: chr$chrom:$position $ref/$alt_str",
		" - error fetching sequence slice",
		" chr$chrom:$seq_start-$seq_end\n";
	    next;
	}
	my $ref_seq = $slice->seq;
	my $strand = 1;
	my $ref_from_seq = substr($ref_seq, $position - $seq_start, $ref_len);
	if ($ref_from_seq ne $ref) {
	    my $rc_ref = reverse_complement($ref);
	    unless ($ref_from_seq eq $rc_ref) {
		print STDERR
		    "ERROR line $line_num: chr$chrom:$position",
		    " $ref/$alt_str - reference allele mis-match:",
		    " $ref -> $ref_from_seq\n";
		next;
	    } else {
		$strand = -1;
	    }
	}

	my $alt_len = length $alt;

	## Produce alternative sequence. May be refactored out.
        my $alt_seq = $ref_seq;
        if ($strand == -1) {
                substr($alt_seq, $rel_position - 1 - $ref_len + 1, $ref_len, reverse_complement($alt));
        } else {
                substr($alt_seq, $rel_position - 1, $ref_len, $alt);
        }
	
	## Call C++ code to compute scores with the alternative sequence across all relevant positions. 
	my $alt_scorefile = score_caller(\@startsitestrarr, \@endsitestrarr, $rel_position, $alt_len, $matrixfilestr, $alt_seq, $matrixIDstr, $chrom, $position, $temp_scores_dir, $alt);
        my %lookuphash_alt = all_scores_hash($alt_scorefile);
	my $miter = $matrix_set->Iterator();
	while (my $matrix = $miter->next) {
	    my $tf_name = $matrix->ID;
	    unless ($tf_name eq $tf_id) {
		next;
	    }

	    my $matrixpath = "$temp_matrix_dir/"."$tf_name".".txt";
	    my $ref_start_site;

	    ## In theory, should never be on the negative strand. Just in case.  
	    if ($strand == -1) {
		$ref_start_site = $rel_position - $matrix->length + 2 - $ref_len;
	    } else {
		$ref_start_site = $rel_position - $matrix->length + 1;
	    }

	    #
	    # Retreive results of alternative score calculations. 
	    #
	    my $alt_top_score;
	    my $alt_top_position;
	    my $alt_top_strand;
	    my $alt_top_rel_score;
	    my @alt_seq_array = @{$lookuphash_alt{$matrix->ID}};
	    my $alt_end_pos =  $rel_position + $alt_len + $matrix->length - 2;
	    ($alt_top_score, $alt_top_rel_score, $alt_top_strand, $alt_top_position) = top_score_search(\@alt_seq_array, $ref_start_site, $alt_end_pos, $rel_position, $matrix, $position);
	    
	    ## Print results. Can be refactored. 
	    print $chrom,"\t",$position,"\t",$ref,"\t",$alt,"\t.\t";
	    print $tf_id,"\t",$tfbs_start,"\t",$tfbs_end,"\t",$ref_abs_score,"\t",$ref_rel_score,"\t",$ref_strand,"\t";
	    print $matrix->ID,"\t",$alt_top_score,"\t",$alt_top_rel_score,"\t",$alt_top_rel_score/$ref_rel_score,"\t",($alt_top_strand==1)?"+":"-","\t",$alt_top_position,"\t",$alt_top_position+$matrix->length-1,"\n";
	    my $output = "$out_dir_scoring/"."out"."$tf_name"."$chrom"."$position"."$alt".".txt";
#		print STDOUT $output,"\n";
	    my $erroutput = "tempError/"."out"."$tf_name"."$chrom"."$position"."$alt".".txt";
	    
	}
	
    }
  
close(FH);
    
  } 

## Outdated procedure. 
sub output_printing{
  my $tf_names = shift;
  my $file_tab = shift;
  my $tab_scores = shift;
  my $header='';
  print STDOUT "Creating files output...\n";
  foreach my $tf_name (keys(%$tf_names)) {
    my $file_name = "$tf_name"."_"."$in_file";
    my $path_name = "$out_dir"."/"."$tf_name"."/"."$file_name";
    if (!open(OUT,">$path_name")) {
      next;
    }
    open(OUT,">$path_name") or die ("$0 cannot open $path_name");
    foreach my $line_num (1..(scalar @$file_tab)) {
      foreach my $alt (keys(%{${${$tab_scores}{$tf_name}}{$line_num}})) {
	my $line = $$file_tab[$line_num-1];
	my @lineray = split /\t/, $line;
	$lineray[4] = $alt;
	$line = join "\t", @lineray;
	chomp $line;
	if ($line =~ /^#CHROM/) {
	  $line = "$line"."\tREF_SCORE\tREF_REL_SCORE\tALT_SCORE\tALT_REL_SCORE\tREF_TOP_SITE\tREF_TOP_STRAND\tALT_TOP_SITE\tALT_TOP_STRAND\tALT/REF_LOG";
	  $header = $line;
	} elsif ($$tab_scores{$tf_name}{$line_num}{$alt}) {
	  my $liste = $$tab_scores{$tf_name}{$line_num}{$alt};
	  $line = "$line"."\t"."$$liste[0]"."%\t"."$$liste[1]"."%\t".
	    "$$liste[2]"."\t"."$$liste[3]"."\t"."$$liste[4]"."\t"."$$liste[5]"."\t".
	      "$$liste[6]"."\t"."$$liste[7]"."\t"."$$liste[8]";
	} elsif (($line !~ /^#/) && ($line)) {
	  $line = "$line"."\t\.\t\.\t\.";	
	}
	print OUT $line,"\n";
      }
    }
    close(OUT);
    my @rank_title = ("ALT/REF_LOG","Z_ALT_LOG","PART_DIST_LOG","TOT_DIST_LOG","EMPIRIC_P-VALUE_BELOW_REF","EMPIRIC_P-VALUE");
    $file_name = "RANKING"."_"."$tf_name"."_"."$in_file";
    $path_name = "$out_dir"."/"."$tf_name"."/"."$file_name";
    open(RANK,">$path_name") or die ("$0 cannot open $path_name");
    for (my $i = 0; $i < 1; $i++) {
      my %rank_tab;
      foreach my $line_num (keys(%{$$tab_scores{$tf_name}})) {
	foreach my $alt (keys(%{$$tab_scores{$tf_name}{$line_num}})) {
	  if ($$tab_scores{$tf_name}{$line_num}{$alt}) {
	    my $items = $$tab_scores{$tf_name}{$line_num}{$alt};
	    $rank_tab{$$items[$i+2]}{$line_num}{$alt} = 1;
	  }
	}
      }
      print RANK $rank_title[$i],"\n";
      print RANK $header,"\n";
      foreach my $value (sort(numeric_sort keys(%rank_tab))) {
	foreach my $item (sort(keys%{$rank_tab{$value}})) {
	  foreach my $alt (keys%{$rank_tab{$value}{$item}}) {
	    my $line = $$file_tab[$item-1];
	    my $liste = $$tab_scores{$tf_name}{$item}{$alt};	
	    my @lineray = split("\t", $line);
	    $lineray[4] = $alt;
	    my $line1 = join "\t", @lineray;
	    chomp $line1;
	    $line = "$line1"."\t"."$$liste[0]"."%\t"."$$liste[1]"."%\t".
	      "$$liste[2]";
	    print RANK $line,"\n";
	  }
	} 	
      }
      print RANK "\n\n";
    }
    close(RANK);
  }   
}
				

# Read through VCF file and check whether it passes overall criteria for
# analysus. Currently the only check is whether any two variations are within
# the maximum matrix width nucleotides from one another. If so output a warning
# message, as the script does not currently take into account cases where more
# than one variant positions may affect a single TFBS.
#
# NOTE: This check ASSUMES the VCF file is sorted by chromosomal position!!!
#
sub check_VCF
  {
    my ($file) = @_;

    my $ok = 1;

    open(FH, $file) || die "Error opening input VCF file $file\n";

    my $last_chrom = ' ';
    my $last_pos = -1;
    my $line_num = 0;
    while (my $line = <FH>) {
      $line_num++;

      chomp $line;

      # skip blank lines
      next unless $line;

      # skip comment lines
      next if $line =~ /^#/;

      my @elem = split /\t/, $line;

      my $chrom    = $elem[0];
      my $position = $elem[1];

      if ($chrom eq $last_chrom) {
	if ($position - $last_pos + 1 <= $max_profile_width) {
#	  print
#	    "\n\nWARNING: VCF contains at least one case where",
#	      " multiple variations may affect a single TFBS.",
#		" This script does not currently take into",
#		  " consideration such cases, i.e. only the effect of",
#                   " individual variations on a single TFBS are reported.",
#		      "\n\n";

	  $ok = 0;

	  last;
	}
      }

      $last_chrom = $chrom;
      $last_pos = $position;
    }

    close(FH);

    return $ok;
  }
# 
# Retrieve the set of matrices, either from the JASPAR database if no matrix
# file is specified, or from the relevant matrix file. 
#
sub get_matrix_set
  {
    my $matrix_set;

    if ($matrix_file) {
      if ($pwm_flag) {
	$matrix_set = read_PWMs($matrix_file);
	die "Error reading PWMs from $matrix_file\n" unless $matrix_set;
      } else {
	$matrix_set = read_PFMs($matrix_file);
	die "Error reading PFMs from $matrix_file\n" unless $matrix_set;
      }
    } else {
      my $db = TFBS::DB::JASPAR5->connect(
					  "dbi:mysql:" . $jaspar_db . ":" . JASPAR_DB_HOST,
					  JASPAR_DB_USER, JASPAR_DB_PASS
					 );

      die "Error connecting to JASPAR database $jaspar_db\n" unless $db;
      my %matrix_args = (
			 '-matrixtype'   => 'PWM'
			);

      $matrix_args{-ID}          = \@matrix_ids if @matrix_ids;
      $matrix_args{-name}        = \@matrix_names if @matrix_names;
      $matrix_args{-collection}  = $collection if $collection;
      $matrix_args{-tax_group}   = \@tax_groups if @tax_groups;
      $matrix_args{-min_ic}      = $min_ic if $min_ic;

      $matrix_set = $db->get_MatrixSet(%matrix_args);

      unless ($matrix_set && $matrix_set->size() > 0) {
	my $msg = "Error fetching JASPAR matrices from $jaspar_db with"
	  . " the following criteria:\n";

	foreach my $key (keys %matrix_args) {
	  $msg .= sprintf "\t$key\t=> %s\n", $matrix_args{$key};
	}

	die "$msg\n";
      }
    }
    
    return $matrix_set;
  }
# 
# Reads in matrices from the specified matrix file. Matrices are in position frequency
# matrix form. 
#
sub read_PFMs
  {
    my ($file) = @_;

    open(FH, $file) || die "Error opening PFM file $file\n";

    my $matrix_set = TFBS::MatrixSet->new();

    my $name          = '';
    my $matrix_string = '';
    my $line_count    = 0;
    while (my $line = <FH>) {
      chomp $line;
      next unless $line;
      if ($line =~ /^>\s*(\S+)/) {
	$name = $1;
      } else {
	if ($line =~ /^\s*[ACGT]\s*\[\s*(.*)\s*\]/) {
	  # line of the form: A [ # # # ... # ]
	  $matrix_string .= "$1\n";
	} elsif ($line =~ /^\s*\d+/) {
	  # line of the form: # # # ... #
	  $matrix_string .= "$line\n";
	} else {
	  next;
	}
	$line_count++;

	if ($line_count == 4) {
	  my $pfm = TFBS::Matrix::PFM->new(
					   -matrixstring => $matrix_string,
					   -name         => $name,
					   -ID           => $name
					  );

	  $matrix_set->add_Matrix($pfm->to_PWM);

	  $line_count    = 0;
	  $name          = '';
	  $matrix_string = '';
	}
      }
    }
    close(FH);

    return $matrix_set;
  }
# 
# Reads in position weight matrices from the specified matrix file. 
#
sub read_PWMs
  {
    my ($file) = @_;

    open(FH, $file) || die "Error opening PWM file $file\n";

    my $matrix_set = TFBS::MatrixSet->new();

    my $name          = '';
    my $matrix_string = '';
    my $line_count    = 0;
    while (my $line = <FH>) {
      chomp $line;
      next if !$line;
      if ($line =~ /^>\s*(\S+)/) {
	$name = $1;
      } else {
	if ($line =~ /^\s*[ACGT]\s*\[\s*(.*)\s*\]/) {
	  # line of the form: A [ # # # ... # ]
	  $matrix_string .= "$1\n";
	} elsif ($line =~ /^\s*\d+/) {
	  # line of the form: # # # ... #
	  $matrix_string .= "$line\n";
	} else {
	  next;
	}
	$line_count++;

	if ($line_count == 4) {
	  my $pwm = TFBS::Matrix::PWM->new(
					   -matrixstring => $matrix_string,
					   -name         => $name,
					   -ID           => $name
					  );

	  $matrix_set->add_Matrix($pwm);

	  $line_count    = 0;
	  $name          = '';
	  $matrix_string = '';
	}
      }
    }
    close(FH);

    return $matrix_set;
  }
# 
# Computes the maximum matrix width (# of columns) among all of the matrices
# in the matrix set. 
#
sub max_profile_width
  {
    my ($matrix_set) = @_;

    my $max_width = 0;

    my $iter = $matrix_set->Iterator;
    while (my $matrix = $iter->next) {
      if ($matrix->length > $max_width) {
	$max_width = $matrix->length;
      }
    }

    return $max_width;
  }
#
# (VERY SLOW) Score calculation function. Give it sequence informations and a pwm matrix, 
# it returns an item top_site with all needed score information.	
#
sub top_site
  {
    my ($matrix, $seq, $search_start, $search_end) = @_;
    my $site_set = $matrix->search_seq( # this right here is the trouble
				       -seqstring  => $seq,
				       -threshold  => MINIMAL_TFBS_THRESHOLD,
				       -subpart    =>  {
							-start  => $search_start,
							-end    => $search_end
						       }
				      );
    unless ($site_set && $site_set->size > 0) {
      return undef;
    }
    my $top_site;
    my $top_score = -999;
    my $iter = $site_set->Iterator;
    while (my $site = $iter->next) {
      if ($site->score > $top_score) {
	$top_score = $site->score;
	$top_site = $site;
      }
    }
    return $top_site;
  }
# 
# Computes the reverse, and the complement, of a particular DNA sequence. Should 
# be self explanatory. 
#
sub reverse_complement
  {
    my $seq = shift;

    $seq =~ tr/acgtACGT/tgcaTGCA/;
 
    return reverse $seq;
  }

#	 
# Really simple functions to sort numeric items.
# (They are actually not used in the script but 
# i keep them for checking data if necessary)
#
sub par_num { return $a <=> $b }
# 
# Keep in mind that this sort will not appreciate NA's in the 
# output. Does not affect the output values themselves, only the
# ordering of them. NA's arise from p-value calculations in which
# the reference/alternative is the lowest/highest score. 
#
sub numeric_sort {
  if ($a < $b) {
    return -1;
  } elsif ($a == $b) {
    return 0;
  } else {
    return 1;
  }
}
# 
# Prototype procedure, to do indel computations quickly, in C++. 
# Not implemented yet, because of certain issues...
#
sub top_scores_two {
  my ($matrixfiles, $seq, $startlist, $endlist, $namestr, $chrom, $pos) = @_;
  `./littleprogramindels $matrixfiles $seq $startlist $endlist $namestr $chrom $pos $temp_scores_dir`;
  return $temp_scores_dir."/"."$chrom"."$pos".".txt";
}
# 
# Returns a path to a list of arrays for each matrix. Each array stores 
# the computation for the score at every position along the sequence,
# with the given matrix. 
#
sub all_scores_two {
  my ($matrixfiles, $seq, $startlist, $endlist, $namestr, $chrom, $pos, $temp_scores_dir, $allele) = @_;
  #print STDOUT $endlist,"\n\n";
  `./littleprogram $matrixfiles $seq $startlist $endlist $namestr $chrom $pos $temp_scores_dir $allele`;
  return $temp_scores_dir."/"."$chrom"."$pos"."$allele".".txt";
}
# 
# Retrieves the appropriate entry from the file created by the procedure above. 
# Note that this is desperately slow, due to grep having to iterate through. Do not use.
#
sub all_scores_retrieve {
  my ($scorefile, $matrix_name) = @_;
  my $line = `grep $matrix_name $scorefile`;
  #print STDOUT "LINE: ", $line,"\n";
  my @retarr = split ",", $line;
  shift(@retarr);
  return @retarr;
}
#
# Improvement to the above procedure -- reads the file returned by all_scores_two, storing
# the entries into a (matrix, *array) hash. 
#
sub all_scores_hash {
  my ($scorefile) = @_;
  my %ret_hash;
  open(scoreFH, $scorefile);
  while (my $line = <scoreFH>) {
    next unless $line;
    chomp $line;
    my @array = split ",", $line;
    $ret_hash{shift(@array)} = \@array;
  }
  close(scoreFH);
  return %ret_hash;
}
#
# Perl computation of the all_scores procedure. Once again, do not use unless absolutely necessary. 
# Does not work for indels. 
#
sub all_scores
  {
    my ($matrix, $seq, $search_start, $search_end) = @_;
    my $base;	
    ## note that limit is EXCLUSIVE, NOT INCLUSIVE, as $position was originally 1-indexed from database
    my $mat = $matrix->matrix;
    my $fscore;
    my $bscore;
    my @seq_array;
    my @sequence = split "", $seq;
    for ($base = $search_start; $base <= ($search_end - $matrix->length + 1); ++$base) {
      my $pos;
      my $nuc;
      $fscore = 0;
      $bscore = 0;
      for ($pos = 0; $pos < $matrix->length; ++$pos) {
	$nuc = $sequence[$base+$pos-1];
	$nuc =~ tr/ACGT/0123/;
	$fscore += $mat->[$nuc][$pos];
	$bscore += $mat->[3-$nuc][$matrix->length - 1 - $pos];
      } 
      $seq_array[2*$base] = $fscore;
      $seq_array[2*$base+1] = $bscore;
    }
	
    return @seq_array;
  }
#
# Searches for the top score among a specified array, indices, and 
# matrices. No score-swapping in this version. Since this version only searches
# alternative sequences, it ONLY provides the top hit. 
#
sub top_score_search {
  my ($arrayref, $start, $end, $index, $matrix, $position) = @_;
  my $matrix_entry;
  my $j;
  my $matrixref = $matrix->matrix;
  my $comp_matrix_entry;
  my $entryFor;
  my $entryRev;
  my $best_position;
  my $best_strand;
  my $best_score = -1024;
  my $best_rel_score;
#  my @good_poss;
#  my @good_strand;
#  my @good_score;
#  my @good_rel_score;
#  my $thres = 0.85*($matrix->max_score - $matrix->min_score) + $matrix->min_score;
  for ($j = $start; $j <= ($end - $matrix->length + 1); ++$j) {
      $matrix_entry = $index - $j;
      $comp_matrix_entry = $matrix->length - 1 - $matrix_entry;
      $entryFor = @$arrayref[2*$j];
      $entryRev = @$arrayref[2*$j+1];
      if (($entryFor > $entryRev) && ($entryFor > $best_score)) {
	  $best_score = (sprintf "%.3f", $entryFor);
	  $best_rel_score = (sprintf '%.3f', ($entryFor-$matrix->min_score)/($matrix->max_score-$matrix->min_score));
          $best_strand = 1;
          $best_position = $j-$start+$position-$matrix->length+1;
      } elsif ($entryRev > $best_score) {
	  $best_score = (sprintf "%.3f", $entryRev);
          $best_rel_score = (sprintf '%.3f', ($entryRev-$matrix->min_score)/($matrix->max_score-$matrix->min_score));
	  $best_strand = -1;
          $best_position = $j-$start+$position-$matrix->length+1;
      }
#      if ($entryFor > $thres) {
#	  push(@good_poss, $j-$start+$position-$matrix->length+1);
#	  push(@good_score, (sprintf '%.3f', $entryFor));
#	  push(@good_strand, 1);
#	  push(@good_rel_score, (sprintf '%.3f', ($entryFor-$matrix->min_score)/($matrix->max_score-$matrix->min_score)));
 #     }
  #    if ($entryRev > $thres) {
#	  push(@good_poss, $j-$start+$position-$matrix->length+1);
#	  push(@good_score, (sprintf '%.3f', $entryRev));
#	  push(@good_strand, -1);
#	  push(@good_rel_score, (sprintf '%.3f', ($entryRev-$matrix->min_score)/($matrix->max_score-$matrix->min_score)));
 #     }
  }
	
  return ($best_score, $best_rel_score, $best_strand, $best_position);

}
# 
# Computation for the top score for a given position, matrix, start and end. Is faster than
# the slow top_site computation. Note however that top_site returns the coordinates of the 
# top site (in a top_site object wrapper) -- this merely returns the score. Simple O(m^2) 
# implementation. 
# Works for indels. This should be the only time when this is used. 
#
sub top_score 
  {
    my ($matrix, $seq, $search_start, $search_end) = @_;
    my $base;	
    ## note that limit is EXCLUSIVE, NOT INCLUSIVE, as $position was originally 1-indexed from database
    my $limit = $search_end - $matrix->length; 
    my $mat = $matrix->matrix;
    my $top_score_ret = -1024;
    my $top_site_ret = 0;
    my $top_strand_ret = 1;
    my $fscore;
    my $bscore;
    ## search_start is counted from a 1-indexed perspective (as $position was), and hence, 
    ## the -1 is necessary. 
    for ($base = $search_start; $base <= $limit+1; ++$base) {
      my $pos;
      my $nuc;
      $fscore = 0;
      $bscore = 0;
      for ($pos = 0; $pos < $matrix->length; ++$pos) {
	$nuc = substr($seq, $base+$pos-1, 1);
	if ($nuc eq "N") {
				# fasta format N equals any base
	  $fscore += ($mat->[0][$pos] + $mat->[1][$pos] + $mat->[2][$pos] + $mat->[3][$pos])/4;
	  $bscore += ($mat->[0][$matrix->length - 1 - $pos] + $mat->[1][$matrix->length - 1 - $pos] + $mat->[2][$matrix->length - 1 - $pos] + $mat->[3][$matrix->length - 1 - $pos])/4;
	} else {
	  $nuc =~ tr/ACGT/0123/;
	  $fscore += $mat->[$nuc][$pos];
	  $bscore += $mat->[3-$nuc][$matrix->length - 1 - $pos];
	}
      } 
      if ($fscore > $bscore && $fscore > $top_score_ret) {
	$top_score_ret = $fscore;
	$top_site_ret = $base - $search_start;
	$top_strand_ret = 1;
      } elsif ($bscore > $top_score_ret) {
	$top_score_ret = $bscore;
	$top_site_ret = $base-$search_start;
	$top_strand_ret = -1;
      }
    }	
    $top_score_ret = sprintf "%.3f", $top_score_ret;
    my @retarr = ($top_score_ret, $top_site_ret, $top_strand_ret);
    return \@retarr;
  }

sub score_caller
{
    my ($startsite_ref, $endsite_ref, $rel_pos, $len, $matrixfilestr, $ref_seq, $matrixIDstr, $chrom, $position, $temp_scores_dir, $allele) = @_;
    
    my @startsitestrcopy = @{$startsite_ref};
    my @endsitestrcopy = @{$endsite_ref};
    foreach my $startsite (@startsitestrcopy) {
	$startsite += $rel_pos - $len;
    }
    foreach my $endsite (@endsitestrcopy) {
	$endsite += $rel_pos + $len;
    }
    my $startsitestr2 = join ",", @startsitestrcopy;
    my $endsitestr2 = join ",", @endsitestrcopy;
    # 
    # Retrieve a hash of (matrix,array) pairs corresponding to the calculation
    # of all_scores for particular matrices. 
    #
    ## TO BE DONE: Instead, just retrieve entries from MANTA to fill the properties here. 
    return all_scores_two($matrixfilestr, $ref_seq, $startsitestr2, $endsitestr2, $matrixIDstr, $chrom, $position, $temp_scores_dir, $allele);
}

sub precompute_parameters
{
    my ($miter) = @_;
    my %tf_names;
    my $matrixfilestr = "";
    my $matrixIDstr = "";
    my $startsitestr = "";
    my $endsitestr = "";
    while (my $matrix = $miter->next) {
	my $tf_name = $matrix->ID;
	$tf_names{$tf_name}=1;
	unless (-e "$out_dir/$tf_name"){
	    mkdir "$out_dir/$tf_name" or die $!; ## Make output directory for each matrix, to store plots, etc.
	} 
	unless (-e "$out_dir/$tf_name/Rel_score"){
	    mkdir "$out_dir/$tf_name/Rel_score" or die $!;	
	}
	my $raw = $matrix->rawprint;
	my $matrixpath = "$temp_matrix_dir/"."$tf_name".".txt";
	unless (-e $matrixpath) {
	    open(OUTMATRIX, ">".$matrixpath);
	    print OUTMATRIX ">"."$tf_name\n"."$raw";
	    close(OUTMATRIX);
	}
	
	#
	# Computes the reference start and end sites based SOLELY on matrix
	# length --> this means that these are not the actual start or end 
	# sites. The values of these will be updated when processing the file. 
	#
	my $ref_start_site = 2 - $matrix->length;
	my $ref_end_site = $matrix->length - 2;
	#
	# Computes information to pass on to the threads, for execution of the 
	# all_scores_two procedure. Possibly a better way to implement this. 
	# Do NOT use arrays, as they are passed by reference to the individual 
	# threads, and the threads each retrieve copies of the original array 
	# upon each iteration. OPTIMIZE ME. 
	#
	if ($matrixfilestr eq "") {
	    $matrixfilestr = $matrixpath;
	    $matrixIDstr = $tf_name;
	    $startsitestr = "$ref_start_site";
	    $endsitestr = "$ref_end_site";
	} else {
	    $matrixfilestr .= ",".$matrixpath;
	    $matrixIDstr .= ",".$tf_name;
	    $startsitestr .= ",$ref_start_site";
	    $endsitestr .= ",$ref_end_site";
	}
    }
    return ($matrixfilestr, $matrixIDstr, $startsitestr, $endsitestr);
}
