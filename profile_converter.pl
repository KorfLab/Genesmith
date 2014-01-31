#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';

#=========================================================#
# PURPOSE:                                                #
# - Takes all HMMER profiles stored in the <hmm_profiles> #
#   and converts them to HMMER 3 format                   #
#      * <hmm_profiles> and all its profiles must be      #
#        added manually                                   #
# - The reformatted profiles are then stored in a new     #
#   directory <hmmer3_profiles>                           #
#=========================================================#

my $hmmer2path = "/Users/rdandekar/Work/KorfProgramProjects/Genesmith/hmm_profiles/*.hmm";
my $hmmer3path = "/Users/rdandekar/Work/KorfProgramProjects/Genesmith/hmmer3_profiles/";
`mkdir $hmmer3path`;

foreach my $fh (glob($hmmer2path)) {
	my ($kog)   = $fh =~ /\/(KOG\d+\.hmm)$/;
	my $profile = $hmmer3path . $kog;
	`hmmconvert --outfmt 3/b $fh > $profile`;
}
