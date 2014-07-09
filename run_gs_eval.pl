#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use FAlite;
use DataBrowser;
use Getopt::Std;
use vars qw($opt_h $opt_5 $opt_m $opt_c $opt_d $opt_i $opt_a $opt_s $opt_3 $opt_b $opt_D $opt_A $opt_U $opt_E);
getopts('h:5:m:c:d:i:a:s:3:b:D:A:U:E:');

#==========================#
# Genesmith Wrapper Script #
#==========================#
my $UP      = 0;     # Order Upstream
my $START   = 0;     # Order Start, includes canonical ATG at the end, length = 3 + order
my $EXON    = 0;     # Order CDS, always 3 states
my $DON     = 0;     # Order Donor, starts with canonical GT, length = 2 + order
my $INTRON  = 0;     # Order Intron, always 3 states
my $ACCEP   = 0;     # Order Acceptor, ends with canonical AG, length = 2 + order
my $STOP    = 0;     # Order Stop, starts with stop codon, length = 3 + order
my $DOWN    = 0;     # Order Downstream
my $BRANCH  = "NO";  # Order for the branch states created by <hmmgen_branch.pl>
my $L_DON   = 2;     # Quantity of Donor States
my $L_ACCEP = 2;     # Quantity of Acceptor States 
my $L_UP    = 500;   # Length of Upstream region parsed for training
my $L_DOWN  = 500;   # Length of Downstream region parsed for training


die "
usage: $0 [options] <GFF> <FASTA>

universal parameters:
  -5 <order>   upstream state info                        Default = $UP
  -m <order>   start codon state info                     Default = $START
  -c <order>   coding state info                          Default = $EXON
  -d <order>   donor site state info                      Default = $DON
  -i <order>   intron body state info                     Default = $INTRON
  -a <order>   acceptor site state info                   Default = $ACCEP
  -s <order>   stop codon state info                      Default = $STOP
  -3 <order>   downstream state info                      Default = $DOWN
  -b <string>  OPTIONAL branch states, enter 'Y' for use  Default = $BRANCH
  -D <length>  Donor Site Length                          Default = $L_DON
  -A <length>  Acceptor Site Length                       Default = $L_ACCEP
  -U <length>  upstream training                          Default = $L_UP
  -E <length>  downtream training                         Default = $L_DOWN
  -h           help (format and details)
" unless @ARGV == 2;

$UP	     = $opt_5 if $opt_5;
$START   = $opt_m if $opt_m;
$EXON    = $opt_c if $opt_c;
$DON     = $opt_d if $opt_d;
$INTRON  = $opt_i if $opt_i;
$ACCEP   = $opt_a if $opt_a;
$STOP    = $opt_s if $opt_s;
$DOWN    = $opt_3 if $opt_3;
$BRANCH  = $opt_b if $opt_b;
$L_DON   = $opt_D if $opt_D;
$L_ACCEP = $opt_A if $opt_A;
$L_UP    = $opt_U if $opt_U;
$L_DOWN  = $opt_E if $opt_E;

my ($GFF, $FASTA) = @ARGV;

# Sanity Check - options
if ($START   !~ /^\d+$/ or
    $DON     !~ /^\d+$/ or 
    $ACCEP   !~ /^\d+$/ or
    $STOP    !~ /^\d+$/ or
    $UP      !~ /^\d+$/ or
    $EXON    !~ /^\d+$/ or
    $INTRON  !~ /^\d+$/ or
    $DOWN    !~ /^\d+$/ or
    ($BRANCH !~ /^NO$/  and $BRANCH !~ /^Y$/) or
    $L_DON   !~ /^\d+$/ or
    $L_ACCEP !~ /^\d+$/ or    
    $L_UP    !~ /^\d+$/ or
    $L_DOWN  !~ /^\d+$/   ) {
    die "Invalid input [options]\n";
}
# OPTIONAL states
if ($BRANCH !~ /^Y$/ and $BRANCH !~ /^NO$/) {die "Invalid BRANCH input [options]\n";}

# Get Training set Filename without File handle
my ($taxa)  = $GFF =~ /\/?(\w+\.*\w+\d*)\.gff$/;


#----------------------------#
# Get Test and Training Sets #
#----------------------------#
# print ">>> Creating TEST and TRAINING sets\n";
`test_train_sets.pl $GFF $FASTA`;

my $sets = 1;  # Stores Total number of Sets
foreach my $fh (glob("./$taxa\_*")) {
	if ($fh =~ /test/) {
		my ($set) = $fh =~ /_test(\d+)\.\w+/;
		$sets = $set + 1 if $set >= $sets;
	}
}
# print "\t$sets sets\n\n";


#-------------#
# Create HMMs #
#-------------#
my $PATH   = "./HMM/";
# print ">>> Generating HMMs\n";
for (my $i=0; $i < $sets; $i++) {
	my $file = "$taxa\_train$i";
	my $em_path = $PATH . $file . "/";
	
	my $gff_fh = "./$taxa\_train$i\.gff";
	my $fa_fh  = "./$taxa\_train$i\.fa";
	if (!-d $em_path) {
		my $cmd = "run_hmmgen.pl";
		$cmd   .= " -b $BRANCH" if $BRANCH =~ /^Y$/;
		$cmd   .= " -D $L_DON -A $L_ACCEP -U $L_UP -E $L_DOWN";
		$cmd   .= " $gff_fh $fa_fh";
		`$cmd`;
# 		print "\t$em_path\tDIRECTORY CREATED\n";
	} else {
# 		print "\t$em_path\tDIRECTORY EXISTS\n";
	}

	my $cmd = "hmm_assmbl.pl";
	$cmd   .= " -D $L_DON -A $L_ACCEP -U $L_UP -E $L_DOWN";
	$cmd   .= " -5 $UP -m $START -c $EXON -d $DON -i $INTRON -a $ACCEP -s $STOP";
	$cmd   .= " -b $BRANCH" if $BRANCH =~ /^Y$/;
	$cmd   .= " $gff_fh $fa_fh";
	`$cmd`;
# 	print "\t$cmd\n\n";
}


#-------------------#
# Running Genesmith #
#-------------------#
# print "\n>>> Running Genesmith\n";
my $output_fh = "$taxa\_pred.gff";
my $exp_fa    = "$taxa\_exp.fa";
my $exp_gff   = "$taxa\_exp.gff";
`rm $output_fh` if -e $output_fh;
`rm $exp_fa`    if -e $exp_fa;
`rm $exp_gff`   if -e $exp_gff;
`touch $output_fh`;
`touch $exp_fa`;
`touch $exp_gff`;

foreach my $fh (glob("./$taxa\_*\.hmm")) {
	my ($set, $st_quant) = $fh =~ /\.\/$taxa\_train(\d+)_(\d+)\.hmm/;
	my $test_fa  = "$taxa\_test$set\.fa";
	my $test_gff = "$taxa\_test$set\.gff";
	`genesmith $fh $test_fa >> $output_fh`;
	`cat $test_fa  >> $exp_fa`;
	`cat $test_gff >> $exp_gff`;
# 	print "\tgenesmith $fh\ $test_fa\ >> $output_fh\n";
}


#----------------------------#
# Evaluate Gene Predications #
#----------------------------#
# print "\n>>> Evaluate Predictions\n";
my $results = `evaluator.pl $exp_fa $exp_gff $output_fh`;
chomp($results);

my @eval_stats = split("\n", $results);
my $nuc_counts = $eval_stats[0];
my $nuc_stats  = $eval_stats[1];
my $cds_counts = $eval_stats[2];
my $kog_counts = $eval_stats[3];

print $taxa,       "\t";
print $nuc_stats,  "\t\n";
print $cds_counts, "\t";
print $kog_counts, "\t\n";


### Remove Extra Files
