#!/usr/bin/perl
use strict;
use warnings;
use fatal;
use HMMstar;
use Getopt::Std;
our ($opt_h, $opt_v);
getopts('hv');

# Global Variables
my $GENOMIC_TRANS	= "0.00021739";
my $EXON_TRANS		= "0.00034963";
my $INTRON_TRANS	= "0.00044256";
my $TRACEBACK		= "[FUNCTION:	HMMER	TRACK:	SEQ	COMBINE_LABEL:	C	TO_LABEL:	S]";

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

# Format Emissions for HMM
emission_table(\%downst_em, 'Genomic_S', 1, \%states);
emission_table(\%start_em, 'Start', 0, \%states);
emission_table(\%exon_em, 'Coding', 2, \%states);
emission_table(\%donor_em, 'Don', 0, \%states);
emission_table(\%intron_em, 'Intron', 1, \%states);
emission_table(\%acceptor_em, 'Accep', 0, \%states);
emission_table(\%stop_em, 'Stop', 0, \%states);
emission_table(\%upst_em, 'Genomic_E', 1, \%states);



#-------------------------------#
#	Generate Transition Matrix	#
#-------------------------------#
my @model_states = (
	'Genomic_S',
	'Start0',
	'Start1',
	'Start2',
	'Coding',
	'Don0',
	'Don1',
	'Intron',
	'Accep0',
	'Accep1',
	'Stop0',
	'Stop1',
	'Stop2',
	'Genomic_E',
);

for (my $i = 0; $i < scalar(@model_states); $i++) {
	my $order = $states{$model_states[$i]}{order};
	 
	if ($order == 0) {
		$states{$model_states[$i]}{transition}{$model_states[$i+1]} = 1;
	} elsif ($order == 1) {
		if ($model_states[$i] =~ /_S$/) {
			$states{$model_states[$i]}{transition}{$model_states[$i+1]}	= $GENOMIC_TRANS;
			$states{$model_states[$i]}{transition}{$model_states[$i]}	= 1 - $GENOMIC_TRANS;
		} elsif ($model_states[$i] =~ /^I/) {
			$states{$model_states[$i]}{transition}{$model_states[$i+1]}	= $INTRON_TRANS;
			$states{$model_states[$i]}{transition}{$model_states[$i]}	= 1 - $INTRON_TRANS;
		} elsif ($model_states[$i] =~ /_E$/) {
			$states{$model_states[$i]}{transition}{'END'}				= 1;
			$states{$model_states[$i]}{transition}{$model_states[$i]}	= 1;
		}
	} else {
		$states{$model_states[$i]}{transition}{$model_states[$i+1]}	= $EXON_TRANS;
		$states{$model_states[$i]}{transition}{'Stop0'}				= $EXON_TRANS;
		$states{$model_states[$i]}{transition}{$model_states[$i]}	= 1 - (2 * $EXON_TRANS);
	}
}

# Sanity Check
use DataBrowser; browse(\%states);



#-------------------#
#	Generate HMM	#
#-------------------#
my $DATE = `date`; chomp($DATE);

open(OUT, ">gene_pred14.hmm") or die "Could not write into OUT\n";
print OUT "#STOCHHMM MODEL FILE

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

foreach my $s1 (@model_states) {
	print OUT "##################################################\n";
	print OUT "STATE:\n";
	print OUT "\tNAME:\t$s1\n";
	print OUT "\tGFF_DESC:\t$s1\n";
	print OUT "\tPATH_LABEL: $states{$s1}{path_label}\n";
	print OUT "TRANSITION:	STANDARD:	P(X)\n";
	foreach my $s2 (sort keys %{$states{$s1}{transition}}) {
		print OUT "\t$s2:\t$states{$s1}{transition}{$s2}";
		if ($s1 eq 'Don1')						{print OUT "\t$TRACEBACK";}
		if ($s1 eq 'Genomic_E' and $s1 eq $s2)	{print OUT "\t$TRACEBACK";}
		print OUT "\n";
	}
	
	print OUT "EMISSION:	SEQ:	COUNTS\n";
	print OUT "\tORDER:	$states{$s1}{order}	AMBIGUOUS:	AVG\n";
	print OUT "$states{$s1}{emission}\n";
}


print OUT "##################################################\n";
print OUT "//END\n";

close OUT;




#===================================================================#
# SUBROUTINES														#
#===================================================================#

# Stores Order and Path Label info for each State
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

# Generates hash of Emission counts for each nucleotide
sub get_emission_counts{
	my ($seq, $length, $order, $counts) = @_;
	
	if ($order == 0) {
		for (my $i = 0; $i < $length; $i++) {
			$counts->{$i}->{uc substr($seq, $i, 1)}++;
		}
	} else {
		for (my $i = $order; $i < $length; $i++) {
			my $ctx = uc substr($seq, $i - $order, $order);
			my $nt = uc substr($seq, $i, 1);
			next unless $nt =~ /[ACGT]/;
			$counts->{$ctx}->{$nt}++;
		}
	}
	return;
}

# Uses emission counts to generate final State Emission for HMM
sub emission_table{
	my ($em_counts, $name, $order, $states) = @_;
	my @alphabet = qw(A C G T);
	my $final_emission = "";
	
	foreach my $i (sort keys %$em_counts) {
		my $emission = "";
		foreach my $letter (@alphabet) {
			my $count;
			if (defined $em_counts->{$i}->{$letter}) {
				$count = $em_counts->{$i}->{$letter};
			} else {
				$count = 0;
			}
			$emission .= "$count\t";
		}
		chop($emission);
		if ($i =~ /\d/ and $order == 0) {
			my $st_name = $name . "$i";
			$states->{$st_name}->{emission} = $emission;
		} else {
			$final_emission .= "$emission\n";
		}
	}
	chop($final_emission);
	if ($order != 0) {
		$states->{$name}->{emission} = $final_emission;
	}
	return;
}