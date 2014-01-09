#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use DataBrowser;

die "usage: format_gff.pl <GFF>\n" unless @ARGV == 1;
my ($GFF) = @ARGV;
my ($FH)  = $GFF =~ /(\w+\d+_\w+)_\w+.gff/;

my $prev_feat = "none";
my $start_id;
my $end_id;
my $cds_start;
my $cds_end;

open (IN, "<$GFF")           or die "Error reading <$GFF>\n";
open (OUT, ">$FH\_pred.gff") or die "Error writing into OUT\n";
while (<IN>) {
	next unless /\S/;
	my ($id, $source, $feat, $start, $end, $score, $strand, $attrib) = split;
	
	if ($prev_feat =~ /^GU0/ and $feat =~ /^start0/) {
		$start_id  = $id;
		$cds_start = $start;
	}
	if ($prev_feat =~ /^cds/ and $feat =~ /^D/) {
		$cds_end = $end;
		print OUT $id,        "\t",
			      $source,    "\t",
			      "CDS",      "\t",
		          $cds_start, "\t",
		          $cds_end,   "\t",
		          $score,     "\t",
		          $strand,    "\t",
		          $id,        "\n";
	}
	if ($prev_feat =~ /^A/ and $feat =~ /^cds/) {
		$cds_start = $start;
	}
	if ($prev_feat =~ /^stop1/ and $feat =~ /^stop2/) {
		$end_id  = $id;
		$cds_end = $end;
		if ($start_id eq $end_id) {
			print OUT $id,        "\t",
			          $source,    "\t",
			          "CDS",      "\t",
			          $cds_start, "\t",
			          $cds_end,   "\t",
			          $score,     "\t",
			          $strand,    "\t",
			          $id,        "\n";
		}
	}
	$prev_feat = $feat;
}
close IN;
close OUT;

