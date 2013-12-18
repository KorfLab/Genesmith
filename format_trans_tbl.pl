#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use DataBrowser;

# Reformats Translation Table into two columns (AA and its matching Codon) 
die "usage: $0 <trans_table>\n" unless @ARGV == 1;
my ($TBL) = @ARGV;

open (IN,  "<$TBL") or die "Error reading input\n";
open (OUT, ">codon_tbl.txt") or die "Error writing output\n";

my $aa_seq;
my $starts;
my $base1;
my $base2;
my $base3;
while (<IN>) {
	chomp;
	my ($name, $seq) = $_ =~ /^(\S+)\s+=\s+(\S+)/;
	$aa_seq = $seq if $name =~ /AAs/;
	$starts = $seq if $name =~ /Starts/;
	$base1  = $seq if $name =~ /Base1/;
	$base2  = $seq if $name =~ /Base2/;
	$base3  = $seq if $name =~ /Base3/;
}

my %codon_key;
for (my $i=0; $i < length($aa_seq); $i++) {
	my $aa     = substr($aa_seq, $i, 1);
	my $status = substr($starts, $i, 1);
	my $first  = substr($base1, $i, 1);
	my $second = substr($base2, $i, 1);
	my $third  = substr($base3, $i, 1);
	my $codon  = $first . $second . $third;
	
	print OUT $aa, "\t", $codon, "\n";
}

close IN;
close OUT;