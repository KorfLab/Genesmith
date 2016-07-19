#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use Getopt::Std;
use vars qw($opt_h $opt_1 $opt_S $opt_5 $opt_m $opt_c $opt_d $opt_i $opt_a $opt_s $opt_3 $opt_B $opt_b $opt_D $opt_A $opt_U $opt_E);
getopts('h1S5:m:c:d:i:a:s:3:Bb:D:A:U:E:');


#========================================#
# Genesmith Prediction using HMMER score #
#========================================#
my $TAXA        = "C.elegans";
my $UP          = 0;     # Order Upstream
my $START       = 0;     # Order Start, includes canonical ATG at the end, length = 3 + order
my $EXON        = 0;     # Order CDS, always 3 states
my $DON         = 0;     # Order Donor, starts with canonical GT, length = 2 + order
my $INTRON      = 0;     # Order Intron, always 3 states
my $ACCEP       = 0;     # Order Acceptor, ends with canonical AG, length = 2 + order
my $STOP        = 0;     # Order Stop, starts with stop codon, length = 3 + order
my $DOWN        = 0;     # Order Downstream
my $BRANCH      = 0;     # Order for the branch states created by <hmmgen_branch.pl>
my $L_DON       = 5;     # Quantity of Donor States
my $L_ACCEP     = 10;    # Quantity of Acceptor States 
my $L_UP        = 1000;  # Length of Upstream region parsed for training
my $L_DOWN      = 1000;  # Length of Downstream region parsed for training
my $HMMER_COEFF = 1.0;
my $SW_COEFF    = 1.0;


die "
usage: $0 [options] <DIR_PATH_test_train_sets> <wb_gene_summary_list.txt>

universal parameters:
  -1           Convert for Standard to Basic HMM with 1 CDS state
  -S           Option to exclude start and stop states from basic model
  -B           OPTION to include branch states
  -b <order>   branch motif state info                    Default = $BRANCH
  -5 <order>   upstream state info                        Default = $UP
  -m <order>   start codon state info                     Default = $START
  -c <order>   coding state info                          Default = $EXON
  -d <order>   donor site state info                      Default = $DON
  -i <order>   intron body state info                     Default = $INTRON
  -a <order>   acceptor site state info                   Default = $ACCEP
  -s <order>   stop codon state info                      Default = $STOP
  -3 <order>   downstream state info                      Default = $DOWN
  -D <length>  Donor Site Length                          Default = $L_DON
  -A <length>  Acceptor Site Length                       Default = $L_ACCEP
  -U <length>  upstream training                          Default = $L_UP
  -E <length>  downtream training                         Default = $L_DOWN
  -h           help (format and details)
" unless @ARGV == 2;

$UP	         = $opt_5 if $opt_5;
$START       = $opt_m if $opt_m;
$EXON        = $opt_c if $opt_c;
$DON         = $opt_d if $opt_d;
$INTRON      = $opt_i if $opt_i;
$ACCEP       = $opt_a if $opt_a;
$STOP        = $opt_s if $opt_s;
$DOWN        = $opt_3 if $opt_3;
$BRANCH      = $opt_b if $opt_b;
$L_DON       = $opt_D if $opt_D;
$L_ACCEP     = $opt_A if $opt_A;
$L_UP        = $opt_U if $opt_U;
$L_DOWN      = $opt_E if $opt_E;

my ($DIR, $LIST) = @ARGV;

# Sanity Check - options
if ($START       !~ /^\d+$/    or
    $DON         !~ /^\d+$/    or 
    $ACCEP       !~ /^\d+$/    or
    $STOP        !~ /^\d+$/    or
    $UP          !~ /^\d+$/    or
    $EXON        !~ /^\d+$/    or
    $INTRON      !~ /^\d+$/    or
    $DOWN        !~ /^\d+$/    or
    $BRANCH      !~ /^\d+$/    or
    $L_DON       !~ /^\d+$/    or
    $L_ACCEP     !~ /^\d+$/    or    
    $L_UP        !~ /^\d+$/    or
    $L_DOWN      !~ /^\d+$/    or
    $opt_S and !$opt_1           ) {
    die "Invalid input [options]\n";
}
$DIR .= "/" if $DIR !~ /\/$/;

my @l_don      = (2, 4, 6, 9);            # set of donor state lengths to run
my @l_accep    = (2, 5, 10, 15, 20, 30);  # set of acceptor state lengths to run

if ($opt_S) {
	for(my $D=0; $D < @l_don; $D++) {
		for(my $A=0; $A < @l_accep; $A++) {
			my $cmd = "run_wb_eval.pl ";
			$cmd   .= "-1 " if $opt_1;
			$cmd   .= "-S " if $opt_S;
			$cmd   .= "-U $L_UP -E $L_DOWN -D $l_don[$D] -A $l_accep[$A] ";    # Constant baseline parameters chosen for select species
			$cmd   .= "-5 $UP -c $EXON -d $DON -i $INTRON -a $ACCEP -3 $DOWN ";
			$cmd   .= "$DIR $LIST";
			print $cmd, "\n";
		}
	}
} else {
	for(my $D=0; $D < @l_don; $D++) {
		for(my $A=0; $A < @l_accep; $A++) {
			my $cmd = "run_wb_eval.pl ";
			$cmd   .= "-1 " if $opt_1;
			$cmd   .= "-U $L_UP -E $L_DOWN -D $l_don[$D] -A $l_accep[$A] ";    # Constant baseline parameters chosen for select species
			$cmd   .= "-5 $UP -c $EXON -d $DON -i $INTRON -a $ACCEP -3 $DOWN ";
			$cmd   .= "-m $START -s $STOP ";
			$cmd   .= "$DIR $LIST";
			print $cmd, "\n";
		}
	}
}

