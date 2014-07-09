#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use STATE;
use HMMstar;
use DataBrowser;
use Getopt::Std;
use vars qw($opt_b $opt_D $opt_A $opt_U $opt_E);
getopts('h:b:D:A:U:E:');

#========================================================#
# PURPOSE                                                #
# - Generate all possible state emissions orders 0-5     #
# - Individual state emissions stored in "HMM" directory #
#========================================================#

# DEFAULT settings [options]
my $MAX_ORDER = 7;     # The Highest order of emissions to be generated
my $BRANCH    = "NO";  # Order for the branch states created by <hmmgen_branch.pl>
my $L_DON     = 2;     # Quantity of Donor States
my $L_ACCEP   = 2;     # Quantity of Acceptor States 
my $L_UP      = 500;   # Length of Upstream region parsed for training
my $L_DOWN    = 500;   # Length of Downstream region parsed for training


die "
usage: $0 [options] <GFF> <FASTA>

universal parameters:
  -b <string>  OPTIONAL branch states, enter 'Y' for use  Default = $BRANCH
  -D <length>  Donor Site Length                          Default = $L_DON
  -A <length>  Acceptor Site Length                       Default = $L_ACCEP
  -U <length>  upstream training                          Default = $L_UP
  -E <length>  downtream training                         Default = $L_DOWN
  -h           help (format and details)
" unless @ARGV == 2;

$BRANCH  = $opt_b if $opt_b;
$L_DON   = $opt_D if $opt_D;
$L_ACCEP = $opt_A if $opt_A;
$L_UP    = $opt_U if $opt_U;
$L_DOWN  = $opt_E if $opt_E;

my ($GFF, $FASTA) = @ARGV;

# Sanity Check - options
if ($BRANCH  !~ /^\w+$/ or
    $L_DON   !~ /^\d+$/ or
    $L_ACCEP !~ /^\d+$/ or    
    $L_UP    !~ /^\d+$/ or
    $L_DOWN  !~ /^\d+$/   ) {
    die "Invalid input [options]\n";
}

#------------------------------------------------------------------------------#
# OPTIONAL STATES: branch states 

if ($BRANCH =~ /^Y$/) {
	for (my $i=0; $i <= $MAX_ORDER; $i++) {
		my $cmd = "hmmgen_branch.pl ";
		$cmd   .= "-i $i -b $i -D $L_DON -A $L_ACCEP ";
		$cmd   .= "$GFF $FASTA";
		`$cmd`;
	}
}
#------------------------------------------------------------------------------#

# Generate Gene model states from order 0 - 5
for (my $i=0; $i <= $MAX_ORDER; $i++) {
	my $cmd = "hmmgen_gene.pl ";
	$cmd   .= "-5 $i -m $i -c $i -d $i -i $i -a $i -s $i -3 $i ";
	$cmd   .= "-D $L_DON -A $L_ACCEP -U $L_UP -E $L_DOWN ";
	$cmd   .= "$GFF $FASTA";
	`$cmd`;
} 

