#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use STATE;
use HMMstar;
use DataBrowser;
use Getopt::Std;
use vars qw($opt_h $opt_b);
getopts('hb');

#========================================================#
# PURPOSE                                                #
# - Generate all possible state emissions orders 0-5     #
# - Individual state emissions stored in "HMM" directory #
#========================================================#

# DEFAULT settings [options]
my $MAX_ORDER  = 5;      # The Highest order of emissions to be generated
my $L_DON      = 9;      # Quantity of Donor States [2-9]
#my $L_ACCEP    = 30;     # Quantity of Acceptor States [10-30]
my @l_accep    = (2, 5, 10, 15, 20, 25, 30);
my $INTERGENIC = 1000;   # Length of Upstream/Downstream regions parsed for training [100-1000]

die "
usage: $0 [options] <GFF> <FASTA>

universal parameters:
  -b   option to train Branch State emissions
  -h   help
" unless @ARGV == 2;
my ($GFF, $FASTA) = @ARGV;

#------------------------------------------------------------------------------#
# OPTIONAL STATES: branch states 

if ($opt_b) {
	for (my $i=0; $i <= $MAX_ORDER; $i++) {
		for (my $d=2; $d < $L_DON+1;$d++) {
			for (my $a=0; $a < @l_accep; $a++) {
				`hmmgen_branch.pl -i $i -b $i -D $d -A $l_accep[$a] $GFF $FASTA`;
			}
		}
	}
}
#------------------------------------------------------------------------------#

# Generate Gene model states from order 0 - 5
for (my $i=0; $i <= $MAX_ORDER; $i++) {
	for(my $e=100; $e < $INTERGENIC+1; $e+=100) {
		`hmmgen_gene.pl -5 $i -m $i -c $i -d $i -i $i -a $i -s $i -3 $i -D 2 -A 2 -U $e -E $e $GFF $FASTA`;
	}
}

for (my $i=0; $i <= $MAX_ORDER; $i++) {
	for(my $d=2; $d < $L_DON+1; $d++) {
		for(my $a=0; $a < @l_accep; $a++) {
			next if $d == 2 and $a == 2;
			`hmmgen_gene.pl -5 $i -m $i -c $i -d $i -i $i -a $i -s $i -3 $i -D $d -A $l_accep[$a] $GFF $FASTA`;
		}
	}
}

