#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use FAlite;
use DataBrowser;
use Getopt::Std;
use vars qw($opt_h $opt_1 $opt_S $opt_t $opt_5 $opt_m $opt_c $opt_d $opt_i $opt_a $opt_s $opt_3 $opt_B $opt_b $opt_D $opt_A $opt_U $opt_E $opt_p $opt_P $opt_w $opt_W);
getopts('h1St:5:m:c:d:i:a:s:3:Bb:D:A:U:E:pP:wW:');

#========================================#
# Genesmith Prediction using HMMER score #
#========================================#
my $UP          = 0;     # Order Upstream
my $START       = 0;     # Order Start, includes canonical ATG at the end, length = 3 + order
my $EXON        = 0;     # Order CDS, always 3 states
my $DON         = 0;     # Order Donor, starts with canonical GT, length = 2 + order
my $INTRON      = 0;     # Order Intron, always 3 states
my $ACCEP       = 0;     # Order Acceptor, ends with canonical AG, length = 2 + order
my $STOP        = 0;     # Order Stop, starts with stop codon, length = 3 + order
my $DOWN        = 0;     # Order Downstream
my $BRANCH      = 0;  # Order for the branch states created by <hmmgen_branch.pl>
my $L_DON       = 5;     # Quantity of Donor States
my $L_ACCEP     = 10;     # Quantity of Acceptor States 
my $L_UP        = 500;   # Length of Upstream region parsed for training
my $L_DOWN      = 500;   # Length of Downstream region parsed for training
my $HMMER_COEFF = 1.0;
my $SW_COEFF    = 1.0;
my $TRANS       = 1.0;


die "
usage: $0 [options] <GFF> <FASTA>

universal parameters:
  -1           Convert for Standard to Basic HMM with 1 CDS state
  -S           Option to exclude start and stop states from basic model
  -p           OPTION to include HMMER score
  -w           OPTION to include SW-alignment score
  -P <digit>   fractional coefficient for HMMER           Default = $HMMER_COEFF
  -W <digit>   fractional coefficient for SW-alignment    Default = $SW_COEFF
  -t <digit>   Weight for CDS transitions                 Default = $TRANS
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
$HMMER_COEFF = $opt_P if $opt_P;
$SW_COEFF    = $opt_W if $opt_W;
$TRANS       = $opt_t if $opt_t;

my ($GFF, $FASTA) = @ARGV;

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
    $HMMER_COEFF !~ /\d?.?\d+/ or
    $SW_COEFF    !~ /\d?.?\d+/ or
#     $TRANS       !~ /\d?.?\d+/ or
    $opt_S and !$opt_1           ) {
    die "Invalid input [options]\n";
}

# Get Training set Filename without File handle
my ($taxa)  = $GFF =~ /\/?(\w+\.*\w+\d*)\.gff$/;


#-------------#
# Create HMMs #
#-------------#
my $optinfo = "D$L_DON\-A$L_ACCEP\-U$L_UP\-E$L_DOWN\-5$UP\-m$START\-c$EXON\-d$DON\-i$INTRON\-a$ACCEP\-s$STOP\-3$DOWN";
$optinfo   .= "\-1"  if $opt_1 and !$opt_S;
$optinfo   .= "\-1S" if $opt_1 and $opt_S;
$optinfo   .= "\-B-b$BRANCH"      if $opt_B; 
$optinfo   .= "\-p-P$HMMER_COEFF" if $opt_p;
$optinfo   .= "\-w-W$SW_COEFF"    if $opt_w;

my $sets    = 5;           # Number of sets (5 sets for 5-fold cross validation)
my $PATH    = "./HMM/";
my $DIR     = "../Data/";  # Directory path where test and training sets are

# print ">>> Generating HMMs\n";
for (my $i=0; $i < $sets; $i++) {
	my $file = "$taxa\_train$i";
	my $em_path = $PATH . $file . "/";
	
	my $gff_fh = "$DIR$taxa\_train$i\.gff";
	my $fa_fh  = "$DIR$taxa\_train$i\.fa";
	if (!-d $em_path) {
		my $cmd = "run_hmmgen.pl";
		$cmd   .= " -b " if $opt_B;
		$cmd   .= " $gff_fh $fa_fh";
		`$cmd`;
# 		print "\t$em_path\tDIRECTORY CREATED\n";
	} else {
# 		print "\t$em_path\tDIRECTORY EXISTS\n";
	}

	my $hmm_fh = "$taxa\_$optinfo\__$i\.hmm";	
	my $cmd    = "hmm_assemble.pl";
	$cmd      .= " -1"        if $opt_1 and !$opt_S;
	$cmd      .= " -1S"       if $opt_1 and $opt_S;
	$cmd      .= " -t $TRANS" if $opt_t;
	$cmd      .= " -D $L_DON -A $L_ACCEP -U $L_UP -E $L_DOWN";
	$cmd      .= " -5 $UP -m $START -c $EXON -d $DON -i $INTRON -a $ACCEP -s $STOP -3 $DOWN";
	$cmd      .= " -B -b $BRANCH" if $opt_B;
	$cmd      .= " $gff_fh $fa_fh";
	`$cmd > $hmm_fh`;
# 	print "\t$cmd\n\n";
}


#-------------------#
# Running Genesmith #
#-------------------#
# print "\n>>> Running Genesmith\n";
my $output_fh = "$taxa\_$optinfo\_pred.gff";
my $exp_fa    = "$taxa\_$optinfo\_exp.fa";
my $exp_gff   = "$taxa\_$optinfo\_exp.gff";
# `rm $output_fh` if -e $output_fh;
# `rm $exp_fa`    if -e $exp_fa;
# `rm $exp_gff`   if -e $exp_gff;
`touch $output_fh`;
`touch $exp_fa`;
`touch $exp_gff`;

for (my $i=0; $i < $sets; $i++) {
	my $fh = "$taxa\_$optinfo\__$i\.hmm";
	my $test_fa  = "$DIR$taxa\_test$i\.fa";
	my $test_gff = "$DIR$taxa\_test$i\.gff";
	my @kogs     = `grep ">" $test_fa`;
	foreach my $kog (@kogs) {
		chomp($kog);
		my ($fa_id) = $kog =~ /^>(\S+)/;
		my $profile = "./KOG/$fa_id\/$fa_id\.hmm";
		my $protein = "./KOG/$fa_id\/$taxa\.prot";
		my $dna     = "./KOG/$fa_id\/$taxa\.dna";
		my $cmd = "genesmith ";
		$cmd   .= "-P $HMMER_COEFF -p $profile " if $opt_p;
		$cmd   .= "-S $SW_COEFF    -s $protein " if $opt_w;
		$cmd   .= "$fh $dna >> $output_fh";
		`$cmd`;
	}
	`cat $test_fa  >> $exp_fa`;
	`cat $test_gff >> $exp_gff`;
# 	print "\tset:  ", $i, "\n";
}


#----------------------------#
# Evaluate Gene Predications #
#----------------------------#
# print "\n>>> Evaluate Predictions\n";
my $results = `eval_quick_compare.pl $exp_gff $output_fh`;
chomp($results);

my @eval_stats = split("\n", $results);
my $header     = ">$taxa\-$optinfo";
print $header,   "\n";
foreach my $stat (@eval_stats) {
	print $stat, "\n";
}

### Remove Model and Prediction Files
my $model_fh = "$taxa\_$optinfo\__*\.hmm";
`rm $model_fh`;
`rm $output_fh`;
`rm $exp_fa`;
`rm $exp_gff`;
