#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';

#==========================================================#
# Purpose: Convert GFF file to GTF file                    #
#    - for Eval software                                   #
#    - all genes in GFF must be on the positive (+) strand #
#    - Use GFF format with 8 columns (no phase) not 9      #
#    - All genes used have only one mRNA transcript        #
#==========================================================#
die "$0 <GFF> <wb_gene_summary_list.txt>\n" unless @ARGV == 2;
my ($GFF, $LIST) = @ARGV;

# Step 1: Get CDS counts
my %data;
open(my $IN, "<", $GFF) or die;
while (my $line = <$IN>) {
	my ($gene, $so, $fea, $beg, $end, $dot, $strand, $text) = split(/\t/, $line);
	$data{$gene}++;
}
close $IN;


# Step 2: Get Transcript info
my %transcripts;
open(my $DB, "<", $LIST) or die;
while (my $line = <$DB>) {
	chomp($line);
	my @info = split(/\t/, $line);
	my $gene = $info[0];
	my ($mrna, $cds_count) = split("=", $info[6]);
	
	if (defined $data{$gene}) {
		$transcripts{$gene} = $mrna;
	} else {
		next;
	}
}
close $DB;


# Step 3: Convert to GTF
my $prev  = "NA";
my $count = 1;
open(my $IN2, "<", $GFF) or die;
while (my $line = <$IN2>) {
	chomp($line);
	my ($gene, $so, $fea, $beg, $end, $dot, $strand, $text) = split(/\t/, $line);
	my $cds_count = $data{$gene};
	my $mrna      = $transcripts{$gene};
	
	if ($prev eq $gene) {
		$count++;
		$prev = $gene;
		if ($count == $cds_count) {
			# stop codon
			my $ex_end   = $end - 3;
			my $stop_beg = $end - 2;
			my $new_feat = "stop_codon";
			$new_feat    = "start_codon" if $strand eq "-";
			print "$gene\t$so\t$fea\t$beg\t$ex_end\t.\t$strand\t.\t", "gene_id \"$gene\"; transcript_id \"$mrna\";\n";
			print "$gene\t$so\t$new_feat\t$stop_beg\t$end\t.\t$strand\t.\t", "gene_id \"$gene\"; transcript_id \"$mrna\";\n";
		} else {
			print "$gene\t$so\t$fea\t$beg\t$end\t.\t$strand\t.\t", "gene_id \"$gene\"; transcript_id \"$mrna\";\n";
		}
	} else {
		# start codon
		$count = 1;
		$prev = $gene;
		my $start_end = $beg + 2;
		my $new_feat  = "start_codon";
		$new_feat     = "stop_codon" if $strand eq "-";
		print "$gene\t$so\t$new_feat\t$beg\t$start_end\t.\t$strand\t.\t", "gene_id \"$gene\"; transcript_id \"$mrna\";\n";
		print "$gene\t$so\t$fea\t$beg\t$end\t.\t$strand\t.\t", "gene_id \"$gene\"; transcript_id \"$mrna\";\n";
	}
	
}
close $IN2;

