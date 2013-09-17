#!/usr/bin/perl
use strict; use warnings;
use Getopt::Std;

die "usage: evaluator.pl <ann> <prediction>\n" unless @ARGV == 2;

my $ann = read_file($ARGV[0]);
my $obs = read_file($ARGV[1]);

my %match;
for (my $i = 0; $i < @$ann; $i++) {
	$match{$ann->[$i]}{$obs->[$i]}++;
}

use DataBrowser; browse(\%match);

foreach my $s1 (sort keys %match) {
	my $found = 0;
	my $missed = 0;
	foreach my $s2 (keys %{$match{$s1}}) {
		if ($s1 eq $s2) {$found += $match{$s1}{$s2}}
		else            {$missed += $match{$s1}{$s2}}
	}
	printf "%s\t%.3f\n", $s1, $found / ($found + $missed);
}

sub read_file {
	my ($filename) = @_;
	
	my @ann;
	open(IN, $filename) or die;
	while (<IN>) {
		next if /^#/;
		next unless /\S/;
		my @f = split;
		my ($beg, $end, $state);
		if (@f == 3) {
			($beg, $end, $state) = @f;
		} elsif (@f == 8) {
			($beg, $end, $state) = ($f[3], $f[4], $f[2]);
		} else {
			die "malformatted annotation";
		}
		
		for (my $i = $beg - 1; $i < $end; $i++) {
			$ann[$i] = $state;
		}
	}
	close IN;
		
	return \@ann;
}

