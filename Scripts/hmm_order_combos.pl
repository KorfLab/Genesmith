#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use Getopt::Std;
use vars qw($opt_h $opt_1 $opt_S);
getopts('h1S');

my $SPECIES = "C.elegans";

#=================================================================#
# Purpose:                                                        #
#  - create every HMM combination (state order, length, quantity) #
#  - run an evaluation on every combination                       #
#=================================================================#
die "
$0 [options] <DIR_PATH_test_train_sets> <wb_gene_summary_list.txt>

universal parameters:
  -1           Convert for Standard to Basic HMM with 1 CDS state
  -S           Option to exclude start and stop states from basic model
  -h           help
" unless @ARGV == 2;
my ($DIR, $LIST) = @ARGV;
$DIR .= "/" if $DIR !~ /\/$/;

my $INTERGENIC = 1000; # max = 1000, min = 500
my $L_DON      = 9;    # max = 9, min = 5
my $L_ACCEP    = 30;   # max = 30, min = 10, supports states in increments of 5
# my @l_don      = (2, 4, 6);           # set of donor state lengths to run
# my @l_accep    = (2, 5, 10, 15, 20);  # set of acceptor state lengths to run

if ($opt_S) {
	for(my $I=1000; $I < $INTERGENIC+1; $I+=1) {
		for(my $D=9; $D < $L_DON+1; $D++) {
			for(my $A=30; $A < $L_ACCEP+1; $A++) {
				for (my $u=0; $u < 2; $u++) {
					for (my $c=3; $c < 6; $c++) {
						for (my $d=0; $d < 2; $d++) {
							for (my $i=3; $i < 6; $i++) {
								for (my $a=0; $a < 2; $a++) {
									for (my $e=0; $e < 2; $e++) {
										my $cmd = "run_wb_eval.pl ";
										$cmd   .= "-1 " if $opt_1;
										$cmd   .= "-S " if $opt_S;
										$cmd   .= "-U $I -E $I -D $D -A $A ";    # Constant baseline parameters chosen for select species
										$cmd   .= "-5 $u -c $c -d $d -i $i -a $a -3 $e ";
										$cmd   .= "$DIR $LIST";
										print $cmd, "\n";
									}
								}
							}
						}
					}
				}
			}
		}
	}
} else {
	for(my $I=1000; $I < $INTERGENIC+1; $I+=1) {
		for(my $D=9; $D < $L_DON+1; $D++) {
			for(my $A=30; $A < $L_ACCEP+1; $A++) {
				for (my $u=0; $u < 2; $u++) {
					for (my $m=0; $m < 2; $m++) {
						for (my $c=3; $c < 6; $c++) {
							for (my $d=0; $d < 2; $d++) {
								for (my $i=3; $i < 6; $i++) {
									for (my $a=0; $a < 2; $a++) {
										for (my $s=0; $s < 2; $s++) {
											for (my $e=0; $e < 2; $e++) {
												my $cmd = "run_wb_eval.pl ";
												$cmd   .= "-1 " if $opt_1;
												$cmd   .= "-U $I -E $I -D $D -A $A ";    # Constant baseline parameters chosen for select species
												$cmd   .= "-5 $u -c $c -d $d -i $i -a $a -3 $e ";
												$cmd   .= "-m $m -s $s ";
												$cmd   .= "$DIR $LIST";
												print $cmd, "\n";
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


__END__
# Example run on Saetta
time run_wb_eval.pl -D 6 -A 15 -c 3 -i 3 -m 1 C.elegans ~/scratch/Basic/

real    9m36.512s
user    9m34.521s

