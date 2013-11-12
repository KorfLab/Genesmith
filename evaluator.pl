#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use DataBrowser;
use PredictionEval;


die "usage: $0 <exp_gff> <predicted_gff>\n" unless @ARGV == 2;
my ($EXP, $OBS) = @ARGV;

my $bp_eval   = new PredictionEval();          # Nucleotide Counts (TP, TN, FP, FN)
my $cds_eval  = new PredictionEval();          # CDS Counts
my $gene_eval = new PredictionEval();          # Core Gene Counts
my $exp       = gff_coords($EXP);              # Expected GFF CDS Coordinates
my $obs       = gff_coords($OBS);              # Predicted GFF CDS Coordinates
my %evid;                                      # Contains Evidence Counts for each gene

#--------------------------------------#
# Get Nucleotide, Exon and Gene Counts #
#--------------------------------------#
foreach my $id (keys %$exp) {
	my $cds_exp = scalar(@{$exp->{$id}});      # Expected quantity of CDS per gene 
	my $cds_obs;                               # Predicted quantity of CDS per gene
	if (defined @{$obs->{$id}}) {
		$cds_obs = scalar(@{$obs->{$id}});
	} else {
		$cds_obs = 0;
	}
	
	my $first = $exp->{$id}->[0]->[0];         # Start position of First CDS
	my $last  = $exp->{$id}->[-1]->[0];        # Start position of Last CDS
		
	# Evaluate Labeled Sequence (1 = CDS, 0 = non-CDS)
	if ($first > $last) {
		my $lb_exp = "0" x ($exp->{$id}->[0]->[-1] + 1);  # Labeled Expected Sequence
		my $lb_obs = "NA";                                # Labeled Predicted Sequence
		if (defined $obs->{$id}->[0]->[-1]) {
			$lb_obs = "0" x ($obs->{$id}->[0]->[-1] + 1);
		}
		
		for(my $i = -1; $i >= -$cds_exp; $i--) {
			my $exp_start = $exp->{$id}->[$i]->[0];
			my $exp_end   = $exp->{$id}->[$i]->[-1];
			my $exp_len   = $exp_end - $exp_start;
			my $obs_start = $obs->{$id}->[$i]->[0];
			my $obs_end   = $obs->{$id}->[$i]->[-1];
			my $obs_len;
			if (defined $obs_start) {
				$obs_len   = $obs_end - $obs_start;
			} else {
				$obs_start = "NA";
				$obs_end   = "NA";
				$obs_len   = "NA";
			}
			
			# Use CDS coordinates to re-label Coding regions of Labeled Seq
			substr($lb_exp, $exp_start, $exp_len) = "1" x $exp_len;
			substr($lb_obs, $obs_start, $obs_len) = "1" x $obs_len if $lb_obs !~ /NA/ and $obs_start !~ /NA/;
		} 
		
		# Nucleotide Level Analysis
		my $counts = evid_counts($id, $lb_exp, $lb_obs, $exp, $obs);
		my $tp     = $counts->{bp}{TP};
		my $tn     = $counts->{bp}{TN};
		my $fp     = $counts->{bp}{FP};
		my $fn     = $counts->{bp}{FN};
		if (!defined $tp)  {$tp = 0;}
		if (!defined $tn)  {$tn = 0;}
		if (!defined $fp)  {$fp = 0;}
		if (!defined $fn)  {$fn = 0;}
		$bp_eval->true_pos($tp);
		$bp_eval->true_neg($tn);
		$bp_eval->false_pos($fp);
		$bp_eval->false_neg($fn);
		
		# CDS Level Analysis
		my $match = $counts->{exon}{MATCH};
		my $mis   = $counts->{exon}{MISMATCH};
		my $none  = $counts->{exon}{MISSING};
		if (!defined $match)  {$match = 0;}
		if (!defined $mis)    {$mis   = 0;}
		if (!defined $none)   {$none  = 0;}
		$cds_eval->match($match);
		$cds_eval->mismatch($mis);
		$cds_eval->missing($none);

		# Gene Level Analysis
		my $gene_status = gene_eval($counts->{exon}, $cds_exp);
		$gene_eval->match(1)    if $gene_status eq 'MATCH';
		$gene_eval->mismatch(1) if $gene_status eq 'MISMATCH';
		$gene_eval->missing(1)  if $gene_status eq 'MISSING';
	} else {
		my $lb_exp = "0" x $exp->{$id}->[-1]->[-1]; # Labeled Expected Sequence
		my $lb_obs = "NA";                          # Labeled Predicted Sequence
		if (defined $obs->{$id}->[0]->[-1]) {
			$lb_obs = "0" x $obs->{$id}->[-1]->[-1];
		}
		
		for(my $i = 0; $i < $cds_exp; $i++) {
			my $exp_start = $exp->{$id}->[$i]->[0];
			my $exp_end   = $exp->{$id}->[$i]->[-1];
			my $exp_len   = $exp_end - $exp_start;
			my $obs_start = $obs->{$id}->[$i]->[0];
			my $obs_end   = $obs->{$id}->[$i]->[-1];
			my $obs_len;
			if (defined $obs_start) {
				$obs_len   = $obs_end - $obs_start;
			} else {
				$obs_start = "NA";
				$obs_end   = "NA";
				$obs_len   = "NA";
			}
			
			# Use CDS coordinates to re-label Coding regions of Labeled Seq
			substr($lb_exp, $exp_start, $exp_len) = "1" x $exp_len;
			substr($lb_obs, $obs_start, $obs_len) = "1" x $obs_len if $lb_obs !~ /NA/ and $obs_start !~ /NA/;
		}
		
		my $counts = evid_counts($id, $lb_exp, $lb_obs, $exp, $obs);
		
		# Nucleotide Level Analysis
		my $tp     = $counts->{bp}{TP};
		my $tn     = $counts->{bp}{TN};
		my $fp     = $counts->{bp}{FP};
		my $fn     = $counts->{bp}{FN};
		if (!defined $tp)  {$tp = 0;}
		if (!defined $tn)  {$tn = 0;}
		if (!defined $fp)  {$fp = 0;}
		if (!defined $fn)  {$fn = 0;}
		$bp_eval->true_pos($tp);
		$bp_eval->true_neg($tn);
		$bp_eval->false_pos($fp);
		$bp_eval->false_neg($fn);
		
		# CDS Level Analysis
		my $match = $counts->{exon}{MATCH};
		my $mis   = $counts->{exon}{MISMATCH};
		my $none  = $counts->{exon}{MISSING};
		if (!defined $match)  {$match = 0;}
		if (!defined $mis)    {$mis   = 0;}
		if (!defined $none)   {$none  = 0;}
		$cds_eval->match($match);
		$cds_eval->mismatch($mis);
		$cds_eval->missing($none);
				
		# Gene Level Analysis
		my $gene_status = gene_eval($counts->{exon}, $cds_exp);
		$gene_eval->match(1)    if $gene_status eq 'MATCH';
		$gene_eval->mismatch(1) if $gene_status eq 'MISMATCH';
		$gene_eval->missing(1)  if $gene_status eq 'MISSING';
	}
}


#---------------------------------#
# Calculate Evaluation Parameters #
#---------------------------------#
#browse($cds_eval);
my $tp = $bp_eval->{true_pos};    # Total True Pos.
my $tn = $bp_eval->{true_neg};    # Total True Neg.
my $fp = $bp_eval->{false_pos};   # Total False Pos.
my $fn = $bp_eval->{false_neg};   # Total False Neg.

my $tpr = PredictionEval::calc_tpr($tp, $fn);           # Sensitivity
my $spc = PredictionEval::calc_spc($tn, $fp);           # Specificity
my $mcc = PredictionEval::calc_mcc($tp, $tn, $fp, $fn); # Matthews Correlation Coefficient
printf "MCC: %.3f\nTPR: %.3f\nSPC: %.3f\n", $mcc, $tpr, $spc;


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

# Compare Labeled Sequences and generate count 
# - Nucleotide Level =>(TP, TN, FP, FN)
# - Exon Level       =>(Match, Missing, Mismatch)
sub evid_counts{
	my ($id, $lb_exp, $lb_obs, $exp, $obs) = @_;
	my %evid;
	
	my $cds_exp = scalar(@{$exp->{$id}});      # Expected quantity of CDS per gene 
	my $first = $exp->{$id}->[0]->[0];         # Start position of First CDS
	my $last  = $exp->{$id}->[-1]->[0];        # Start position of Last CDS
	
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
		my $obs_exp = "NA";
		my $obs_seq = "NA";
		$obs_exp = substr($lb_obs, $exp_start, $exp_len) if $lb_obs !~ /NA/;
		$obs_seq = substr($lb_obs, $obs_start, $obs_len) if $lb_obs !~ /NA/;
		
		# Nucleotide Level Analysis
		for(my $n=0; $n < $exp_len; $n++) {
			my $exp_bp = substr($exp_seq, $n, 1);
			my $obs_bp = "NA";
			$obs_bp = substr($obs_exp, $n, 1) if $obs_exp !~ /NA/;
			
			# separate from Exon evidence =>   ex: $evid{bp}{TN}++
			if    ($exp_bp eq "0" and $obs_bp eq "0") {$evid{bp}{TN}++;}
			elsif ($exp_bp eq "1" and $obs_bp eq "1") {$evid{bp}{TP}++;}
			elsif ($exp_bp eq "0" and $obs_bp eq "1") {$evid{bp}{FP}++;}
			else                                      {$evid{bp}{FN}++;}
		}
		
		# Exon Level Analysis
		# - separate from Nuc. evidence =>   ex: $evid{exon}{MATCH}++
		if (($exp_seq eq $obs_seq) and ($exp_seq eq $obs_exp)) {$evid{exon}{MATCH}++;}  
		elsif ($obs_seq =~ /NA/)                               {$evid{exon}{MISSING}++;} 
		else                                                   {$evid{exon}{MISMATCH}++;}
	}
	return \%evid;
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
