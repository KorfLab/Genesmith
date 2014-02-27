#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use HMMstar;
use DataBrowser;
use FAlite;

#=============================================#
# Purpose:                                    #
#   Create output files containing            #
#   -  KOG gene lengths (genomic and protein) #
#   -  CDS/Intron lengths                     #
#=============================================#
die "usage: gff_lengths.pl <GFF> <dnaFASTA> <aaFASTA>\n" unless @ARGV == 3;
my ($GFF, $FASTA, $AA) = @ARGV;
my ($TAXA) = $GFF =~ /(\w\.\w)\w+/;
$TAXA =~ s/\.//;

### Get CDS and Intron Lengths
my $genome = new HMMstar::Genome($FASTA, $GFF);
my %gene_struc;

foreach my $contig ($genome->contigs) {
	foreach my $cds ($contig->coding_sequences) {
		my $ex_count = 1;		# KOG specific exon count
		my $in_count = 1;		# KOG specific intron count
		
		# Exon
		foreach my $exon ($cds->exons) {
			$gene_struc{$cds->name}{Exons}{$ex_count} = length($exon->sequence);
			$gene_struc{$cds->name}{tot_exons}			 = $ex_count;
			$ex_count++;
		}
		# Intron
		foreach my $intron ($cds->introns) {
			# some KOGs only have 1 CDS and no Introns
			if (defined $intron->sequence) {
				$gene_struc{$cds->name}{Introns}{$in_count} = length($intron->sequence);
				$gene_struc{$cds->name}{tot_introns}		 = $in_count;
				$in_count++;
			}
		}		
	}
}
#browse(\%gene_struc);


### Get Gene Lengths
my %genomic_len = gene_lengths($FASTA);
my %protein_len = gene_lengths($AA);
#browse(\%protein_len);


### Print Outputs Files
open(KOG,    ">$TAXA\_gene_feat.txt") or die "Error writing into KOG\n";
open(CDS,    ">$TAXA\_cds_len.txt")   or die "Error writing into CDS\n";
open(INTRON, ">$TAXA\_in_len.txt")    or die "Error writing into INTON\n";

my $ex_quant;
my $in_quant;
my $bp_len;
my $aa_len;
foreach my $kog (keys %genomic_len) {
	$ex_quant   = $gene_struc{$kog}{tot_exons};
	$in_quant   = $gene_struc{$kog}{tot_introns};
	$in_quant   = 0 if(!defined $in_quant);
	$bp_len     = $genomic_len{$kog};
	$aa_len     = $protein_len{$kog};
	# KOG lengths
	print KOG $TAXA,     "\t", $bp_len,   "\t", $aa_len, "\t", 
	          $ex_quant, "\t", $in_quant, "\n";
	          
	# CDS lengths
	foreach my $cds (sort {$a<=>$b} keys %{$gene_struc{$kog}{Exons}}) {
		my $cds_len = $gene_struc{$kog}{Exons}{$cds};
		print CDS $TAXA, "\t", $cds, "\t", $cds_len, "\n";
	}
	# Intron lengths
	if ($in_quant == 0) {
		#print INTRON $TAXA, "\t0\t0\n";
	} else {
		foreach my $intron (sort {$a<=>$b} keys %{$gene_struc{$kog}{Introns}}) {
			my $in_len = $gene_struc{$kog}{Introns}{$intron};
			print INTRON $TAXA, "\t", $intron, "\t", $in_len, "\n";
		}
	}
}

close KOG;
close CDS;
close INTRON;

#=============#
# SUBROUTINES #
#=============#

# Creates hash all Gene Lengths for a given FASTA file
sub gene_lengths{
	my ($FH) = @_;
	my %gene_lengths;
	
	open (SEQ, "<$FH") or die "Error reading $FH\n";
	my $fasta = new FAlite(\*SEQ);
	
	while (my $entry = $fasta->nextEntry) {
		my ($id) = $entry->def =~ /^>(\S+)$/;
		$gene_lengths{$id} = length($entry->seq);
	}
	return %gene_lengths;
}
