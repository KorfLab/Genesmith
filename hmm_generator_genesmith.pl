#!/usr/bin/perl
use strict;
use warnings;
use fatal;
use HMMstar;
use Getopt::Std;
our ($opt_h, $opt_v);
getopts('hv');

# Global Variables
my $DATE = `date`; chomp($DATE);

die "
usage: $0 [options] <gff> <fasta>
options:
  -h  help
  -v  verbose
" unless @ARGV == 2;

my ($GFF, $DNA) = @ARGV;

#-----------------------------------#
#	Generate State Emission Counts	#
#-----------------------------------#
my $genome = new HMMstar::Genome($DNA, $GFF);
my %states;				# HMM hash
my $state_count = 0;

my %downst_em;
my %upst_em;
my %start_em;
my %stop_em;
my %donor_em;
my %acceptor_em;
my %intron_em;
my %exon_em;


foreach my $contig ($genome->contigs) {
	foreach my $cds ($contig->coding_sequences) {
		my $offset = 500;
		
		# Upstream inergenic region
		my $upstream = $cds->start_site->sequence($offset,0); chop($upstream);
		get_emission_counts($upstream, length($upstream), 1, \%upst_em);
		get_state_info('Genomic_S', length($upstream), 1, \%states);
		
		# Start site
		my $start = $cds->start_site->sequence(0,2);
		get_emission_counts($start, length($start), 0, \%start_em);
		get_state_info('Start', length($start), 0, \%states);
		
		
		#Coding 2nd order
		foreach my $exon ($cds->exons) {
			my $order = 2;
			my $body = substr($exon->sequence, 0, $exon->length);
			get_emission_counts($body, length($body), 2, \%exon_em);
			get_state_info('Coding', length($body), 2, \%states);
		}
		
		#Introns
		foreach my $intron ($cds->introns) {
			# donor
			my $d = substr($intron->sequence, 0, 2);
			get_emission_counts($d, length($d), 0, \%donor_em);
			get_state_info('Don', length($d), 0, \%states);
			
			# intron body as 1st order
			my $body = substr($intron->sequence, 2, $intron->length -4);
			get_emission_counts($body, length($body), 1, \%intron_em);
			get_state_info('Intron', length($body), 1, \%states);
			
			# acceptor
			my $a = substr($intron->sequence, -2, 2);
			get_emission_counts($a, length($a), 0, \%acceptor_em);
			get_state_info('Accep', length($a), 0, \%states);
		}
		
		
		# Stop site
		my $stop = $cds->stop_site->sequence(0,2);
		get_emission_counts($stop, length($stop), 0, \%stop_em);
		get_state_info('Stop', length($stop), 0, \%states);
		
		# Downstream inergenic region
		my $downstream = $cds->stop_site->sequence(0,$offset + 2);
		$downstream = substr($downstream, 3, length($downstream)-3);
		get_emission_counts($downstream, length($downstream), 1, \%downst_em);
		get_state_info('Genomic_E', length($downstream), 1, \%states);
	}
}

# Sanity Check
#use DataBrowser; browse(\%states);

# Output Emission Counts
print "GENOMIC_S\n";	emission_table(\%downst_em);
print "\nSTART\n";		emission_table(\%start_em);
print "\nCODING\n";		emission_table(\%exon_em);
print "\nDONOR\n";		emission_table(\%donor_em);
print "\nINTRON\n";		emission_table(\%intron_em);
print "\nACCEPTOR\n";	emission_table(\%acceptor_em);
print "\nSTOP\n";		emission_table(\%stop_em);
print "\nGENOMIC_E\n";	emission_table(\%upst_em);




#-------------------#
#	Generate HMM	#
#-------------------#
my %model = (
	'Genomic_S'		=>	1,
	'Start0'		=>	2,
	'Start1'		=>	3,
	'Start2'		=>	4,
	'Coding'		=>	5,
	'Don0'			=>	6,
	'Don1'			=>	7,
	'Intron'		=>	8,
	'Accep0'		=>	9,
	'Accep1'		=>	10,
	'Stop0'			=>	11,
	'Stop1'			=>	12,
	'Stop2'			=>	13,
	'Genomic_E'		=>	14,
);


open(OUT, ">test_gene_pred.hmm") or die "Could not write into OUT\n";
print OUT "\
#STOCHHMM MODEL FILE

<MODEL INFORMATION>
======================================================
NAME:	Genesmith
DESCRIPTION:	14 state model for gene prediction
CREATION_DATE:	$DATE

<TRACK SYMBOL DEFINITIONS>
======================================================
SEQ:	A,C,G,T

AMBIGUOUS SYMBOL DEFINITIONS
======================================================
SEQ:	N[A,C,G,T], R[A,G], Y[C,T]


<STATE DEFINITIONS>
######################################################
STATE:	
	NAME:	INIT
TRANSITION:	STANDARD:	P(X)
	Genomic_S:	1
";

foreach my $st (sort {$model{$a} <=> $model{$b}} keys(%model)) {
	print OUT "##################################################\n";
	print OUT "STATE:\n";
	print OUT "\tNAME:\t$st\n";
	print OUT "\tGFF_DESC:\t$st\n";
	print OUT "\tPATH_LABEL: $states{$st}{path_label}\n";
	print OUT "TRANSITION:	STANDARD:	P(X)\n";
	print OUT "EMISSION:	SEQ:	COUNTS\n";
	print OUT "\tORDER:	$states{$st}{order}	AMBIGUOUS:	AVG\n";
}


print OUT "##################################################\n";
print OUT "//END\n";

close OUT;




#===================================================================#
# SUBROUTINES														#
#===================================================================#

sub get_state_info{
	my ($name, $length, $order, $states) = @_;
	my $label;
	if ($order == 0) {
		$label = 'C' if ($name =~ /Start|Stop/);
		$label = 'I' if ($name =~ /Accep|Don/);
		for (my $i = 0; $i < $length; $i++) {
			my $st_name = $name . "$i";
			$states->{$st_name}->{order} = $order;
			$states->{$st_name}->{path_label} = $label;
		}
	} else {
		$label = 'C' if ($name =~ /Coding/);
		$label = 'I' if ($name =~ /Intron/);
		$label = 'S' if ($name =~ /Genomic_S/);
		$label = 'E' if ($name =~ /Genomic_E/);
		
		$states->{$name}->{order} = $order;
		$states->{$name}->{path_label} = $label;
	}
	return;
}


sub get_emission_counts{
	my ($seq, $length, $order, $states) = @_;
	
	if ($order == 0) {
		for (my $i = 0; $i < $length; $i++) {
			$states->{$i}->{uc substr($seq, $i, 1)}++;
		}
	} else {
		for (my $i = $order; $i < $length; $i++) {
			my $ctx = uc substr($seq, $i - $order, $order);
			my $nt = uc substr($seq, $i, 1);
			next unless $nt =~ /[ACGT]/;
			$states->{$ctx}->{$nt}++;
		}
	}
	return;
}


sub emission_table{
	my ($states) = @_;
	my @alphabet = qw(A C G T);
	
	foreach my $i (sort keys %$states) {
		if ($i =~ /\d/) {print $i, "\t";}
		my $emission = "";
		foreach my $letter (@alphabet) {
			my $count;
			if (defined $states->{$i}->{$letter}) {
				$count = $states->{$i}->{$letter};
			} else {
				$count = 0;
			}
			$emission .= "$count\t";
		}
		print $emission, "\n";
	}
	return;
}