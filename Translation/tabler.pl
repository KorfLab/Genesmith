#!/usr/bin/perl
use strict;
use warnings;

my %d;
while (<>) {
	my ($field, $val) = $_ =~ /(\S+)\s+= (\S+)/;
	$d{$field} = $val;
}

my %gc;
for (my $i = 0; $i < 64; $i++) {
	my $b1 = substr($d{Base1}, $i, 1);
	my $b2 = substr($d{Base2}, $i, 1);
	my $b3 = substr($d{Base3}, $i, 1);
	$gc{"$b1$b2$b3"} = substr($d{AAs}, $i, 1);
}

foreach my $codon (sort keys %gc) {
	print "'$gc{$codon}',";
}
print "\n";
