#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use FAlite;
use HMMstar;
use DataBrowser;


#=====================================#
# Purpose                             #
# - Split inputted FASTA/GFF into     #
#   combinations of test and training #
#   sets                              #
#=====================================#
die "usage: $0 <GFF> <FASTA>\n" unless @ARGV == 2;
my ($GFF, $FASTA) = @ARGV;

# Extract Gene Sequence Structure
my $genome = new HMMstar::Genome($FASTA, $GFF);
my @type   = qw(Upstream Exons Introns Downstream);
my %gene_struc;

# Create hash of Gene Structure (Exon, Introns, etc.) Seqs
foreach my $contig ($genome->contigs) {
	foreach my $cds ($contig->coding_sequences) {
		my $offset = 500;		# seq size used to train upstream/downstream regions
		my $ex_count = 1;		# KOG specific exon count
		my $in_count = 1;		# KOG specific intron count
		
		# Upstream
		my $upstream = $cds->start_site->sequence($offset,0); chop($upstream);
		$gene_struc{$cds->name}{$type[0]} = lc($upstream);
		
		# Exon
		foreach my $exon ($cds->exons) {
			$gene_struc{$cds->name}{$type[1]}{$ex_count} = lc($exon->sequence);
			$gene_struc{$cds->name}{tot_exons}			 = $ex_count;
			$ex_count++;
		}
		# Intron
		foreach my $intron ($cds->introns) {
			# some KOGs only have 1 CDS and no Introns
			if (defined $intron->sequence) {
				$gene_struc{$cds->name}{$type[2]}{$in_count} = lc($intron->sequence);
				$gene_struc{$cds->name}{tot_introns}		 = $in_count;
				$in_count++;
			}
		}
		
		# Downstream
		my $downstream = $cds->stop_site->sequence(0,$offset + 2);
		$downstream = substr($downstream, 3, length($downstream)-3);
		$gene_struc{$cds->name}{$type[3]} = lc($downstream);
	}
}



#-------------------------------------#
# Validate                            #
# - Remove KOGs with Alt Splice Sites #
#   form the Training Set             #
# - Create Test and Training sets     #
#-------------------------------------#
my $alt_splice = validate(\%gene_struc, \@type);
foreach my $id (keys %$alt_splice) {delete $gene_struc{$id};}
#browse($alt_splice);
#browse(\%gene_struc);


# Of the Validated KOGs create test and training sets (4 combos)
my %seqs;   #hash of all FASTA sequences
my $taxa;

open (IN, "<$FASTA") or die "Error reading FASTA file\n";
my $fasta = new FAlite(\*IN);
while(my $entry = $fasta->nextEntry) {
	my ($id)    = $entry->def =~ /^>(\S+)$/;
	($taxa)     = $id =~ /^(\w\w)\w+/;
	$seqs{$id}  = uc($entry->seq);
}
close IN;

my @val_kogs;
foreach my $id (keys %gene_struc) {push(@val_kogs, $id);}
my $sets       = 4;
my $test_quant = int(scalar(@val_kogs)/$sets);

test_train_seqs($GFF, \%seqs, \@val_kogs, $test_quant, $taxa, $sets);


#===================================================================#
# SUBROUTINES														#
#===================================================================#

# Validates Gene Structures and stores genes that do not follow basic model
sub validate{
	my ($gene_struc, $type) = @_;
	my %alt_splice;

	foreach my $id (keys %$gene_struc) {
		my $ex_max = $gene_struc{$id}{tot_exons};    # Total exons
		my $in_max = $gene_struc{$id}{tot_introns};  # Total introns
		for (my $i=1; $i < $ex_max + 1; $i++) {
			my $ex_seq = $gene_struc{$id}{Exons}{$i};
			if ($i == 1 and defined $in_max) {
				my $start = substr($ex_seq, 0, 3);
				if ($start =~ /atg/) {
					my $in_seq  = $gene_struc{$id}{Introns}{$i};
					my $don     = substr($in_seq, 0, 2);
					my $accep   = substr($in_seq, -2, 2);
					if ($don !~ /gt/) {
						$alt_splice{$id}{don}{$i} = $don;
					}
					if ($accep !~ /ag/) {
						$alt_splice{$id}{accep}{$i} = $accep;
					}
					if ($start eq $ex_seq) { 
						$alt_splice{$id}{start_for_ex_1} = $start;
					}
				} else {
					$alt_splice{$id}{first_ex} = $start;
				}
			} elsif (!defined $in_max) {
				my $start = substr($ex_seq, 0, 3);
				my $stop  = substr($ex_seq, -3, 3);
				if ($start !~ /atg/ and $stop !~ /taa|tga|tag/) {
					$alt_splice{$id}{one_cds} = "$start..$stop";
				}
			} elsif ($i == $ex_max) {
				my $stop = substr($ex_seq, -3, 3);
				if ($stop !~ /taa|tga|tag/) {
					$alt_splice{$id}{last_ex} = $stop;
				}
			} else {
				my $in_seq  = $gene_struc{$id}{Introns}{$i};
				my $don     = substr($in_seq, 0, 2);
				my $accep   = substr($in_seq, -2, 2);
				if ($don !~ /gt/) {
					$alt_splice{$id}{don}{$i} = $don;
				}
				if ($accep !~ /ag/) {
					$alt_splice{$id}{accep}{$i} = $accep;
				}
			}
		}
	}
	return \%alt_splice;
}

# Create Test and Training Sets of FASTA seqs
sub test_train_seqs{
	my ($GFF, $seqs, $val_kogs, $test_quant, $taxa, $sets) = @_;
	my $quant;
	
	# FASTA files
	for(my $i=0; $i < $sets; $i++) {
		print "SET: $i\t";
		# Determine Start Position
		my $pos = $test_quant * $i;
		
		# Determine number of iterations
		if ($i == $sets - 1) {$quant = scalar(@$val_kogs) - $pos;}
		else                 {$quant = $test_quant;}
		print $pos, "\t", $quant, "\n";
		
		open(TEST,     ">$taxa\_test$i\.fa" )  or die "Error writing into TEST\n";
		open(TRAIN,    ">$taxa\_train$i\.fa")  or die "Error writing into TRAIN\n";
		open(TESTGFF,  ">$taxa\_test$i\.gff" ) or die "Error writing into TESTGFF\n";
		open(TRAINGFF, ">$taxa\_train$i\.gff") or die "Error writing into TRAINGFF\n";

		
		for(my $f=0; $f < scalar(@$val_kogs); $f++) {
			my $id  = $val_kogs->[$f];
			my $seq = $seqs->{$id};
			my $ann = `grep "$id" $GFF`;
			
			if ($f >= $pos and $f <= ($pos + $quant)) {
				print TEST ">$id\n", $seq, "\n";
				print TESTGFF $ann;
			} else {
				print TRAIN ">$id\n", $seq, "\n";
				print TRAINGFF $ann;
			}
		}
		close TEST;
		close TRAIN;
		close TESTGFF;
		close TRAINGFF;
		print "COMPLETE\n\n";	
	}
}
