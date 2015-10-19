#!/usr/bin/perl
use strict;
use warnings;
use fatal;


die "convert_ann_to_gff.pl <ZFF>" unless @ARGV == 1;
my ($zff_file) = @ARGV;

my %Feature = (
	'Exon'  => 'CDS',
	'Einit' => 'CDS',
	'Eterm' => 'CDS',
	'Esngl' => 'CDS',
);
my $seq;
my %H;

open(ZFF, "<$zff_file") or die "Error reading GFF\n";
while (<ZFF>) {
	if (/^>(\S+)/) {
		$seq = $1;
	} else {
		my @f = split;
		if (@f == 4) {
			my $strand = $f[1] < $f[2] ? '+' : '-';
			if ($strand eq '-') {($f[1], $f[2]) = ($f[2], $f[1])}
			print join("\t", $seq, 'snap', $Feature{$f[0]}, $f[1], $f[2], '.',
				$strand, $f[3]), "\n";
		} elsif (@f == 9) {
			print join("\t", $seq, 'verified', $Feature{$f[0]}, $f[1], $f[2], '.',
				$f[3], $f[8]), "\n";
		} else {die "input does not appear to be ZFF"}
	}
}
