#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use FAlite;
use DataBrowser;
use PredictionEval;


die "usage: $0 <exp_fasta> <exp_gff> <predicted_gff>\n" unless @ARGV == 3;
my ($FASTA, $EXP, $OBS) = @ARGV;

my $bp_sum    = new PredictionEval();          # Nucleotide Counts (TP, TN, FP, FN)
my $cds_sum   = new PredictionEval();          # CDS Counts
my $gene_sum  = new PredictionEval();          # Core Gene Counts
my $exp       = gff_coords($EXP);              # Expected GFF CDS Coordinates
my $obs       = gff_coords($OBS);              # Predicted GFF CDS Coordinates
my $slengths  = seq_lengths($FASTA);           # Expected lengths of all FASTA sequences
my %evid;                                      # Contains Evidence Counts for each gene


#--------------------------------------#
# Get Nucleotide, Exon and Gene Counts #
#--------------------------------------#
foreach my $id (keys %$exp) {
	my $slen = $slengths->{$id};               # FASTA sequence length
	my $cds_exp = scalar(@{$exp->{$id}});      # Expected quantity of CDS per gene 
	my $cds_obs;                               # Predicted quantity of CDS per gene
	if (defined @{$obs->{$id}}) {
		$cds_obs = scalar(@{$obs->{$id}});
	} else {
		$cds_obs = 0;
	}
	my $lb_exp = "0" x $slen;
	my $lb_obs = "0" x $slen;	
	
	# Use CDS coordinates to re-label Coding regions of Labeled Seq
	for(my $i = 0; $i < $cds_exp; $i++) {
		my $exp_start = $exp->{$id}->[$i]->[0];
		my $exp_end   = $exp->{$id}->[$i]->[-1];
		my $exp_len   = $exp_end - $exp_start;
		substr($lb_exp, $exp_start, $exp_len) = "1" x $exp_len;
	}
	
	if ($cds_obs > 0) {
		for(my $i = 0; $i < $cds_obs; $i++) {
			my $obs_start = $obs->{$id}->[$i]->[0];
			my $obs_end   = $obs->{$id}->[$i]->[-1];
			my $obs_len = $obs_end - $obs_start;
			substr($lb_obs, $obs_start, $obs_len) = "1" x $obs_len;
		}
	}
	my $bp_counts  = bp_eval($lb_exp, $lb_obs, $slen);
	my $cds_counts = cds_eval($id, $lb_exp, $lb_obs, $cds_obs, $exp, $obs);
		
	# Nucleotide Level Analysis
	my $tp     = $bp_counts->{TP};
	my $tn     = $bp_counts->{TN};
	my $fp     = $bp_counts->{FP};
	my $fn     = $bp_counts->{FN};
	if (!defined $tp)  {$tp = 0;}
	if (!defined $tn)  {$tn = 0;}
	if (!defined $fp)  {$fp = 0;}
	if (!defined $fn)  {$fn = 0;}
	$bp_sum->true_pos($tp);
	$bp_sum->true_neg($tn);
	$bp_sum->false_pos($fp);
	$bp_sum->false_neg($fn);
	
	# CDS Level Analysis
	my $match = $cds_counts->{MATCH};
	my $mis   = $cds_counts->{MISMATCH};
	my $none  = $cds_counts->{MISSING};
	if (!defined $match)  {$match = 0;}
	if (!defined $mis)    {$mis   = 0;}
	if (!defined $none)   {$none  = 0;}
	$cds_sum->match($match);
	$cds_sum->mismatch($mis);
	$cds_sum->missing($none);
			
	# Gene Level Analysis
	my $gene_status = gene_eval($cds_counts, $cds_exp);
	$gene_sum->match(1)    if $gene_status eq 'MATCH';
	$gene_sum->mismatch(1) if $gene_status eq 'MISMATCH';
	$gene_sum->missing(1)  if $gene_status eq 'MISSING';
}


#---------------------------------#
# Calculate Evaluation Parameters #
#---------------------------------#
#browse($bp_sum);
my $tp = $bp_sum->{true_pos};    # Total True Pos.
my $tn = $bp_sum->{true_neg};    # Total True Neg.
my $fp = $bp_sum->{false_pos};   # Total False Pos.
my $fn = $bp_sum->{false_neg};   # Total False Neg.

my $tpr = PredictionEval::calc_tpr($tp, $fn);           # Sensitivity
my $spc = PredictionEval::calc_spc($tn, $fp);           # Specificity
my $acc = PredictionEval::calc_acc($tp, $tn, $fp, $fn); # Accuracy
my $mcc = PredictionEval::calc_mcc($tp, $tn, $fp, $fn); # Matthews Correlation Coefficient
print  "MCC\tACC\tTPR\tSPC\n";
printf "%.3f\t%.3f\t%.3f\t%.3f\n", $mcc, $acc, $tpr, $spc;
#printf "MCC: %.3f\nACC: %.3f\nTPR: %.3f\nSPC: %.3f\n", $mcc, $acc, $tpr, $spc;


#===================================================================#
# SUBROUTINES														#
#===================================================================#

# Stores all CDS coordinates for all IDs in a GFF file
sub gff_coords{
	my ($gff) = @_;
	my %coords;
	
	open(IN, "<$gff") or die "Could not read GFF\n";
	while (<IN>) {
		next unless /\S/;
		my @line = split;
		my ($id, $start, $end);
		if (@line == 8) {
			$id    = $line[0];
			$start = $line[3];
			$end   = $line[4];
			
		} else {
			die "Incorrect GFF format\n";
		}
		push (@{$coords{$id}}, [$start, $end]);
	}
	return \%coords;
}

# Stores FASTA sequence lengths for all IDs
sub seq_lengths{
	my ($fh) = @_;
	my %slen;
	
	open (IN, "<$fh") or die "Error reading FASTA\n";
	my $fasta = new FAlite(\*IN);
	while (my $entry = $fasta->nextEntry) {
		my ($id)   = $entry->def =~ /^>(\S+)$/;
		$slen{$id} = length($entry->seq);
	}
	close IN;
	return \%slen;
}

# Nucleotide Level Analysis (TP, TN, FP, FN)
sub bp_eval{
	my ($lb_exp, $lb_obs, $slen) = @_;
	my %counts;
	
	for(my $i=0; $i < $slen; $i++) {
		my $exp = substr($lb_exp, $i, 1);
		my $obs = substr($lb_obs, $i, 1);
		if    ($exp eq "0" and $obs eq "0") {$counts{TN}++;}
		elsif ($exp eq "1" and $obs eq "1") {$counts{TP}++;}
		elsif ($exp eq "0" and $obs eq "1") {$counts{FP}++;}
		else                                {$counts{FN}++;}
	}
	return \%counts;
}

# CDS Level Analysis
sub cds_eval{
	my ($id, $lb_exp, $lb_obs, $cds_obs, $exp, $obs) = @_;
	my $cds_exp = scalar(@{$exp->{$id}});      # Expected quantity of CDS
	my %counts;
	
	if ($cds_obs == 0) {
		$counts{MISSING} += $cds_exp;
	} else {
		for(my $i = 0; $i < $cds_exp; $i++) {
			my $exp_start = $exp->{$id}->[$i]->[0];
			my $exp_end   = $exp->{$id}->[$i]->[-1];
			my $exp_len   = $exp_end - $exp_start;
			my $obs_start = $obs->{$id}->[$i]->[0];
			my $obs_end   = $obs->{$id}->[$i]->[-1];
			my $obs_len;
			
			if (defined $obs_start and defined $obs_end) {
				$obs_len   = $obs_end - $obs_start;
			} else {
				$obs_start = "NA";
				$obs_end   = "NA";
				$obs_len   = "NA";
			}
			my $exp_seq = substr($lb_exp, $exp_start, $exp_len);
			my $obs_seq = "NA";
			$obs_seq = substr($lb_obs, $obs_start, $obs_len) if $obs_len !~ /NA/;
			
			if (($exp_seq   eq $obs_seq)   and 
			    ($exp_start == $obs_start) and 
			    ($exp_len   == $obs_len))          {$counts{MATCH}++;}  
			elsif ($obs_seq =~ /NA/)               {$counts{MISSING}++;} 
			else                                   {$counts{MISMATCH}++;}
		}
		if ($cds_obs > $cds_exp) {
			$counts{MISMATCH} += ($cds_obs - $cds_exp);
		}
	}
	return \%counts;
}

# Gene Level Analysis
sub gene_eval{
	my ($counts, $cds_exp) = @_;
	my $gene_status = "";
	
	foreach my $ct (keys %$counts) {
		if ($ct eq 'MATCH') {
			my $ex_count = $counts->{$ct};
			$gene_status = 'MATCH' if $ex_count == $cds_exp;
		} elsif ($ct eq 'MISSING') {
			my $ex_count = $counts->{$ct};
			$gene_status = 'MISSING' if $ex_count == $cds_exp;
		} else {
			$gene_status = 'MISMATCH';
		}
	}
	return $gene_status;
}