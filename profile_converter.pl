#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';

#=====================================================================#
# PURPOSE:                                                            #
# - Only converts KOGs Profiles                                       #
# - Takes all HMMER profiles stored in the HMMER 2 profiles Directory #
#   and converts them to HMMER 3 format                               #
#      * Directories for HMMER 2 & 3 profiles must be                 #
#        added manually                                               #
# - The reformatted profiles are then stored in a new                 #
#   directory <HMMER 3 directory PATH>                                #
#=====================================================================#

die "usage: profile_converter.pl <HMMER 2 profile PATH> <HMMER 3 profile PATH>\n" unless @ARGV == 2;
my ($hmmer2path, $hmmer3path) = @ARGV;

$hmmer2path .= "*.hmm";
foreach my $fh (glob($hmmer2path)) {
	my ($kog)   = $fh =~ /\/(KOG\d+\.hmm)$/;
	my $profile = $hmmer3path . $kog;
	`hmmconvert --outfmt 3/b $fh > $profile`;
}
