#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use DataBrowser;


# Convert all HMMER2 profiles to HMMER3 format within KOG directory in Genesmith
my $MAINPATH = "/Users/rdandekar/Work/KorfProgramProjects/Genesmith/KOG/";
my %profiles;

foreach my $fh (glob($MAINPATH . "*")) {
	my ($kog_id) = $fh =~ /\/Genesmith\/KOG\/(\S+)/;
	$profiles{$kog_id} = $fh . "/$kog_id\.hmm" if $kog_id ne "warnings";
}

foreach my $id (keys %profiles) {
	my $new_profile = $id . ".hmm";
	my ($path) = $profiles{$id} =~ /(\S+\/$id\/)$new_profile/; 
	`hmmconvert --outfmt 3/b $profiles{$id} > $new_profile`;
	`mv $new_profile $profiles{$id}`;
}