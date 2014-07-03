#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use FAlite;
use HMMstar;
use DataBrowser;
use Getopt::Std;
use vars qw($opt_h $opt_s);
getopts('h:s:');


#=====================================#
# Purpose                             #
# - Split inputted FASTA/GFF into     #
#   combinations of test and training #
#   sets                              #
#=====================================#
my $SETS = 5;

die "
usage: $0 [options] <GFF> <FASTA>

options:
  -s <int/str>  # sets for cross validation   Default = $SETS
                enter <all> if you want each
                test set to equal 1
  -h            help (view usage statement)
" unless @ARGV == 2;
my ($GFF, $FASTA) = @ARGV;
$SETS = $opt_s if $opt_s;


# Of the Validated KOGs create test and training sets (4 combos)
my %seqs;   #hash of all FASTA sequences
my ($taxa) = $GFF =~ /\/?(\w+\.*\w+\d*)\.gff$/;
my @val_kogs;

open (IN, "<$FASTA") or die "Error reading FASTA file\n";
my $fasta = new FAlite(\*IN);
while(my $entry = $fasta->nextEntry) {
	my ($id)    = $entry->def =~ /^>+(\S+)/;
	$seqs{$id}  = uc($entry->seq);
	push(@val_kogs, $id);
}
close IN;

print "#Valid KOGs: ", scalar(@val_kogs), "\n";

$SETS = scalar(@val_kogs) if $SETS eq 'all';
my $test_quant = int(scalar(@val_kogs)/$SETS);
test_train_seqs($GFF, \%seqs, \@val_kogs, $test_quant, $taxa, $SETS);


#===================================================================#
# SUBROUTINES														#
#===================================================================#


# Create Test and Training Sets of FASTA seqs
sub test_train_seqs{
	my ($GFF, $seqs, $val_kogs, $test_quant, $taxa, $sets) = @_;
	my $quant;
	
	# FASTA files
	for(my $i=0; $i < $sets; $i++) {
		print "SET: $i\t";
		# Determine Start Position
		my $pos = $test_quant * $i;
		
		# Determine number of iterations
		if ($i == $sets - 1) {$quant = scalar(@$val_kogs) - $pos;}
		else                 {$quant = $test_quant;}
		print $pos, "\t", $quant, "\n";
		
		open(TEST,     ">$taxa\_test$i\.fa" )  or die "Error writing into TEST\n";
		open(TRAIN,    ">$taxa\_train$i\.fa")  or die "Error writing into TRAIN\n";
		open(TESTGFF,  ">$taxa\_test$i\.gff" ) or die "Error writing into TESTGFF\n";
		open(TRAINGFF, ">$taxa\_train$i\.gff") or die "Error writing into TRAINGFF\n";

		
		for(my $f=0; $f < scalar(@$val_kogs); $f++) {
			my $id  = $val_kogs->[$f];
			my $seq = $seqs->{$id};
			my $ann = `grep "$id" $GFF`;
			
			if ($f >= $pos and $f < ($pos + $quant)) {
				print TEST ">$id\n", $seq, "\n";
				print TESTGFF $ann;
			} else {
				print TRAIN ">$id\n", $seq, "\n";
				print TRAINGFF $ann;
			}
		}
		close TEST;
		close TRAIN;
		close TESTGFF;
		close TRAINGFF;
		print "COMPLETE\n\n";	
	}
}
