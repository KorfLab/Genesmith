#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use FAlite;
use DataBrowser;

#==========================#
# Genesmith Wrapper Script #
#==========================#
die "$0 <GFF> <FASTA>\n" unless @ARGV == 2;
my ($GFF, $FASTA) = @ARGV;
my ($taxa)  = $GFF =~ /\/?(\w+\.*\w+\d*)\.gff$/;


#----------------------------#
# Get Test and Training Sets #
#----------------------------#
print ">>> Creating TEST and TRAINING sets\n";
#`test_train_sets.pl $GFF $FASTA`;

my $sets = 1;  # Stores Total number of Sets
foreach my $fh (glob("./$taxa\_*")) {
	if ($fh =~ /test/) {
		my ($set) = $fh =~ /_test(\d+)\.\w+/;
		$sets = $set + 1 if $set >= $sets;
	}
}
print "\t$sets sets\n\n";


#-------------#
# Create HMMs #
#-------------#
my $CDS    = 5;  # CDS order
my $INTRON = 5;  # Intron order
my $STOP   = 5;  # Stop codon order
my $PATH   = "./HMM/";
print ">>> Generating HMMs\n";
for (my $i=0; $i < $sets; $i++) {
	my $file = "$taxa\_train$i";
	my $em_path = $PATH . $file . "/";
	
	my $gff_fh = "./$taxa\_train$i\.gff";
	my $fa_fh  = "./$taxa\_train$i\.fa";
	if (!-d $em_path) {
		`run_hmmgen.pl $gff_fh $fa_fh`;
		print "\t$em_path\tDIRECTORY CREATED\n";
	} else {
		print "\t$em_path\tDIRECTORY EXISTS\n";
	}

	my $cmd = "hmm_assmbl.pl";
	$cmd   .= " -c $CDS -i $INTRON -s $STOP";
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
	`genesmith $fh $test_fa >> $output_fh`;
	`cat $test_fa  >> $exp_fa`;
	`cat $test_gff >> $exp_gff`;
	print "\tgenesmith $fh\ $test_fa\ >> $output_fh\n";
}

### Evaluate Gene Predicitions

### Remove Extra Files

