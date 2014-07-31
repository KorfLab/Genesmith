#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use FAlite;
use DataBrowser;
use Getopt::Std;
use vars qw($opt_h $opt_P $opt_S);
getopts('h:P:S:');

#========================================#
# Genesmith Prediction using HMMER score #
#========================================#

# Conditions for Optimal S.cerevisiae Gene Model
my $UP      = 0;     # Order Upstream
my $START   = 0;     # Order Start, includes canonical ATG at the end
my $EXON    = 5;     # Order CDS, always 3 states
my $DON     = 1;     # Order Donor, starts with canonical GT
my $INTRON  = 5;     # Order Intron, always 3 states
my $ACCEP   = 1;     # Order Acceptor, ends with canonical AG
my $STOP    = 0;     # Order Stop, starts with stop codon
my $DOWN    = 0;     # Order Downstream
my $L_DON   = 4;     # Quantity of Donor States
my $L_ACCEP = 3;     # Quantity of Acceptor States 
my $L_UP    = 500;   # Length of Upstream region parsed for training
my $L_DOWN  = 200;   # Length of Downstream region parsed for training

# Default Coefficients for HMMER and SW-alignment scores
my $HMMER_COEFF = 1.0;
my $SW_COEFF    = 1.0;

die "
usage: $0 [options] <GFF> <FASTA>

parameters:
   -P <digit>  fractional coefficient for HMMER          Default = $HMMER_COEFF
   -S <digit>  fractional coefficient for SW-alignment   Default = $SW_COEFF
   -h          help (format and details)
" unless @ARGV == 2;

$HMMER_COEFF = $opt_P if $opt_P;
$SW_COEFF    = $opt_S if $opt_S;

my ($GFF, $FASTA) = @ARGV;

# Sanity Check
if ($HMMER_COEFF !~ /\d?.?\d+/ or $SW_COEFF !~ /\d?.?\d+/) {
    die "Invalid input [options]\n";
}

# Get Training set Filename without File handle
my ($taxa)  = $GFF =~ /\/?(\w+\.*\w+\d*)\.gff$/;


#----------------------------#
# Get Test and Training Sets #
#----------------------------#
# print ">>> Creating TEST and TRAINING sets\n";
#`test_train_sets.pl $GFF $FASTA`;

# my $sets = 1;  # Stores Total number of Sets
# foreach my $fh (glob("./$taxa\_*")) {
# 	if ($fh =~ /test/) {
# 		my ($set) = $fh =~ /_test(\d+)\.\w+/;
# 		$sets = $set + 1 if $set >= $sets;
# 	}
# }
# print "\t$sets sets\n\n";
my $sets = 5;

#-------------#
# Create HMMs #
#-------------#
my $PATH   = "./HMM/";
print ">>> Generating HMMs\n";
for (my $i=0; $i < $sets; $i++) {
	my $file = "$taxa\_train$i";
	my $em_path = $PATH . $file . "/";
	
	my $gff_fh = "./$taxa\_train$i\.gff";
	my $fa_fh  = "./$taxa\_train$i\.fa";
	if (!-d $em_path) {
		my $cmd = "run_hmmgen.pl";
		$cmd   .= " -D $L_DON -A $L_ACCEP -U $L_UP -E $L_DOWN";
		$cmd   .= " $gff_fh $fa_fh";
		`$cmd`;
		print "\t$em_path\tDIRECTORY CREATED\n";
	} else {
		print "\t$em_path\tDIRECTORY EXISTS\n";
	}

	my $cmd = "hmm_assmbl.pl";
	$cmd   .= " -D $L_DON -A $L_ACCEP -U $L_UP -E $L_DOWN";
	$cmd   .= " -5 $UP -m $START -c $EXON -d $DON -i $INTRON -a $ACCEP -s $STOP";
	$cmd   .= " $gff_fh $fa_fh";
	`$cmd`;
	print "\t$cmd\n\n";
}


#-------------------#
# Running Genesmith #
#-------------------#
print "\n>>> Running Genesmith\n";
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
	
	# Create temp FASTA file with one ID/seq
	open(IN, "<$test_fa") or die "Error reading $test_fa\n";
	my $fasta = new FAlite(\*IN);
	while (my $entry = $fasta->nextEntry) {
		my ($fa_id)   = $entry->def =~ /^>(\S+)$/;
		my $profile = "./KOG/$fa_id\/$fa_id\.hmm";
		my $protein = "./KOG/$fa_id\/$taxa\.prot";
		open(ONE, ">one_id.fa") or die "Error writing into one_id.fa\n";
		print ONE $entry;
		close ONE;
		`genesmith -P $HMMER_COEFF -S $SW_COEFF -p $profile -s $protein $fh ./one_id.fa >> $output_fh`;
	}
	`cat $test_fa  >> $exp_fa`;
	`cat $test_gff >> $exp_gff`;
	print "\tset:  ", $set, "\n";
}

#----------------------------#
# Evaluate Gene Predications #
#----------------------------#
print "\n>>> Evaluate Predictions\n";
my $results = `evaluator.pl $exp_fa $exp_gff $output_fh`;
chomp($results);

my @eval_stats = split("\n", $results);
my $nuc_counts = $eval_stats[0];
my $nuc_stats  = $eval_stats[1];
my $cds_counts = $eval_stats[2];
my $kog_counts = $eval_stats[3];

my $opt_info = "> -5 $UP -m $START -c $EXON -d $DON -i $INTRON -a $ACCEP -s $STOP -3 $DOWN";
print $opt_info,   "\n";
print $taxa,       "\t";
print $nuc_stats,  "\t\n";
print $cds_counts, "\t";
print $kog_counts, "\t\n";


### Remove Model Files
my $model_fh = "$taxa\_train*.hmm";
`rm $model_fh`;

