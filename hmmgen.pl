#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use STATE;
use HMMstar;
use DataBrowser;
use Getopt::Std;
use vars qw($opt_h $opt_f $opt_s $opt_C $opt_D $opt_i $opt_A $opt_e $opt_t);
getopts('h:f:s:C:D:i:A:e:t:');


# Help <-h> option
my $HELP;
while (my $line = <DATA>) {$HELP .= $line;}
die "$HELP" if $opt_h;

# DEFAULT setting
# <#states>:<Order>
my $UP      = "1:1";
my $START   = "3:0"; 
my $EXON    = "3:2";
my $DON     = "2:1";
my $INTRON  = "3:2";
my $ACCEP   = "2:1";
my $STOP    = "3:2";
my $DOWN    = "1:1";

die "
usage: $0 [options] <GFF> <FASTA>

options:
  -f <order>          upstream state info       Default = $UP
  -s <order>          start codon state info    Default = $START
  -C <#states:order>  coding state info         Default = $EXON
  -D <#states:order>  donor site state info     Default = $DON
  -i <order>          intron body state info    Default = $INTRON
  -A <#states:order>  acceptor site state info  Default = $ACCEP
  -e <order>          stop codon state info     Default = $STOP
  -t <order>          downstream state info     Default = $DOWN
  -h                  help (format and details)
" unless @ARGV == 2;

$UP	    = edit_cmd_opt('UP',     $UP,     $opt_f) if $opt_f;
$START  = edit_cmd_opt('START',  $START,  $opt_s) if $opt_s;
$EXON   = edit_cmd_opt('EXON',   $EXON,   $opt_C) if $opt_C;
$DON    = edit_cmd_opt('DON',    $DON,    $opt_D) if $opt_D;
$INTRON = edit_cmd_opt('INTRON', $INTRON, $opt_i) if $opt_i;
$ACCEP  = edit_cmd_opt('ACCEP',  $ACCEP,  $opt_A) if $opt_A;
$STOP   = edit_cmd_opt('STOP',   $STOP,   $opt_e) if $opt_e;
$DOWN   = edit_cmd_opt('DOWN',   $DOWN,   $opt_t) if $opt_t;

my ($GFF, $FASTA) = @ARGV;


#---------------------#
# Important Variables #
#---------------------#
my ($i_size)  = $INTRON =~ /^(\d+):/;             # Total Intron States
my ($d_size)  = $DON    =~ /^(\d+):/;             # Total Donor States
my ($a_size)  = $ACCEP  =~ /^(\d+):/;             # Total Acceptor States
my ($e_size)  = $EXON   =~ /^(\d+):/;             # Total Exon States

my ($e_order) = $EXON   =~ /:(\d+)$/;
my $pad_quant;                                    # Total Padded CDS States
$pad_quant = $e_order - 2 if $e_order > 2;
$pad_quant = 0            if $e_order <= 2;

my @init_info = ($UP, $START, $EXON, $DON, $INTRON, $ACCEP, $STOP, $DOWN);
my @init_st	  = qw(GU start cds D i A stop GD);
my @labels    = qw(U C C I I I C D);
my @type      = qw(Upstream Exons Introns Downstream);


#-------------------------------#
# Create Array of State Objects #
#-------------------------------#
my @states;		# will contain all state object refs

# Initialize objects with State name, label, #states and order
for (my $i=0; $i < scalar(@init_st); $i++) {
	if ($init_st[$i] =~ /^[DAi]/) {
		my ($quant, $order) = split(":", $INTRON);
		for (my $c=0; $c < $quant; $c++) {
			my $name = $init_st[$i] . "$c";
			if ($name =~ /^[DA]/) {
				my ($s_quant, $s_order) = split(":", $init_info[$i]);
				$name .= "_";
				for (my $s=0; $s < $s_quant; $s++) {
					$name .= "$s";
					my $st_obj = new STATE($name, $labels[$i], $s_order, $s_quant);
					push(@states, $st_obj);
					chop($name);
				}
			} else {
				my $st_obj = new STATE($name, $labels[$i], $order, $quant);
				push(@states, $st_obj);
				chop($name);
			}	
		}
	} else {
		my ($quant, $order) = split(":", $init_info[$i]);
		if ($quant == 1) {
			my $name = $init_st[$i] . "0";
			my $st_obj = new STATE($name, $labels[$i], $order, $quant);
			push(@states, $st_obj);
			# Pad States
			if ($name =~ /cds/ and $order >= 3) {
				my $pad_sts = pad_states($name, $order);
				foreach my $st (sort keys %$pad_sts) {
					my $pad_order = $pad_sts->{$st}->{order};
					my $pad_quant = $pad_sts->{$st}->{st_quant};
					my $st_obj = new STATE($st, $labels[$i], $pad_order, $pad_quant);
					push(@states, $st_obj);
				}
			}
		} else {
			for (my $c=0; $c < $quant; $c++) {
				my $name = $init_st[$i] . "$c";
				my $st_obj = new STATE($name, $labels[$i], $order, $quant);
				push(@states, $st_obj);
				# Pad States
				if ($name =~ /cds/ and $order >= 3) {
					my $pad_sts = pad_states($name, $order);
					foreach my $st (sort keys %$pad_sts) {
						my $pad_order = $pad_sts->{$st}->{order};
						my $pad_quant = $pad_sts->{$st}->{st_quant};
						my $st_obj = new STATE($st, $labels[$i], $pad_order, $pad_quant);
						push(@states, $st_obj);
					}
				}
			}
		}
	}
}
#foreach my $obj (@states) {browse($obj);}


#---------------------------------#
# Extract Gene Sequence Structure #
#---------------------------------#
my $genome = new HMMstar::Genome($FASTA, $GFF);
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
#-------------------------------------#
my $alt_splice = validate(\%gene_struc, \@type);
foreach my $id (keys %$alt_splice) {delete $gene_struc{$id};}
#browse($alt_splice);
#browse(\%gene_struc);


#----------------------------------------------------------------------#
# Calculate Transition Probabilities                                   #
#----------------------------------------------------------------------#
# Label Gene Seqs based on Structure                                   #
# -  label key:                                                        #
#    GU = 0, start = 1, cds = 2, D = 5, i = 6, A = 7, stop = 8, GD = 9 # 
#----------------------------------------------------------------------#
my %trans_freq;

# Label Sequence
foreach my $id (keys %gene_struc) {
	my $seq    = "";                             # Full gene seq
	my $labseq = "";                             # Labeled gene seq
	my $ex_max = $gene_struc{$id}{tot_exons};    # Total exons
	my $in_max = $gene_struc{$id}{tot_introns};  # Total introns
	foreach my $struc (@type) {
		if ($struc =~ /Upstream/) {
			my $lab  = label('up', $gene_struc{$id}{$struc});
			$labseq .= $lab;
			$seq    .= $gene_struc{$id}{$struc};
		} elsif ($struc =~ /Exons/) {
			for (my $i=1; $i < $ex_max + 1; $i++) {
				my $ex_seq = $gene_struc{$id}{$struc}{$i};
				if ($i == 1 and defined $in_max) {
					my $in_seq   = $gene_struc{$id}{Introns}{$i};
					my $ex_lab   = label('first_ex', $ex_seq);
					my $in_lab   = label('in', $in_seq);
					$labseq     .= $ex_lab . $in_lab;
					$seq        .= $ex_seq . $in_seq; 
				} elsif (!defined $in_max) {
					my $ex_lab   = label('only_ex', $ex_seq);
					$labseq     .= $ex_lab;
					$seq        .= $ex_seq;
				} elsif ($i == $ex_max) {
					my $ex_lab   = label('last_ex', $ex_seq);
					$labseq     .= $ex_lab;
					$seq        .= $ex_seq;
				} else {
					my $in_seq   = $gene_struc{$id}{Introns}{$i};
					my $ex_lab   = label('ex', $ex_seq);
					my $in_lab   = label('in', $in_seq);
					$labseq     .= $ex_lab . $in_lab;
					$seq        .= $ex_seq . $in_seq;
				}
			}
		} elsif ($struc =~ /Downstream/) {
			my $lab  = label('down', $gene_struc{$id}{$struc});
			$labseq .= $lab;
			$seq    .= $gene_struc{$id}{$struc};
		}
	}
	$labseq .= 'END';
	my @labels = split(" ", $labseq);
	
	# Store Transition Counts
	for (my $i=0; $i < scalar(@labels)-1; $i++) {
		my $lab     = $labels[$i];
		my $nxt_lab = $labels[$i+1];
		$trans_freq{$lab}{$nxt_lab}++;
	}
}
#browse(\%trans_freq);



#-------------------#
# Transition Matrix #
#-------------------#
for (my $i=0; $i < scalar(@states); $i++) {
	my $st    = $states[$i]->name;
	my $order = $states[$i]->order;
	my $quant = $states[$i]->st_quant;
	
	for (my $t=0; $t < 10; $t++) {
		if ($t == 0 and $st =~ /GU/) {
			my $total  = tot_trans(\%{$trans_freq{$t}});
			my $t_prob = $trans_freq{$t}{$t}/$total;
			$states[$i]->t_matrix($st, $t_prob);
			$states[$i]->t_matrix($states[$i+1]->name, (1 - $t_prob));
		} elsif ($t == 9 and $st =~ /GD/) {
			my $total  = tot_trans(\%{$trans_freq{$t}});
			my $t_prob = $trans_freq{$t}{$t}/$total;
			$states[$i]->t_matrix($st, $t_prob);
			$states[$i]->t_matrix('END', (1 - $t_prob));
		} elsif ($t == 1 and $st =~ /start/) {
			if ($quant == 1) {
				my $total  = tot_trans(\%{$trans_freq{$t}});
				my $t_prob = $trans_freq{$t}{$t}/$total;
				$states[$i]->t_matrix($st, $t_prob);
				$states[$i]->t_matrix($states[$i+1]->name, (1 - $t_prob));				
			} else {
				$states[$i]->t_matrix($states[$i+1]->name, 1);
			}
		} elsif ($t == 8 and $st =~ /stop/) {
			if ($quant == 1) {
				my $total  = tot_trans(\%{$trans_freq{$t}});
				my $t_prob = $trans_freq{$t}{$t}/$total;
				$states[$i]->t_matrix($st, $t_prob);
				$states[$i]->t_matrix($states[$i+1]->name, (1 - $t_prob));				
			} else {
				$states[$i]->t_matrix($states[$i+1]->name, 1);
			}
		} elsif ($st =~ /cds\d+$/) {
			my ($pos) = $st =~ /cds(\d+)/;
			if ($t ==2) {
				if ($pos == $quant - 1) {
					my $trans          = "D" . "$pos" . "_0";
					my $total          = tot_trans(\%{$trans_freq{$t}});
					my $t_prob         = $trans_freq{$t}{$t}/$total;          # CDS->CDS
					my $in_t_prob      = $trans_freq{$t}{5}/$total;           # CDS->DON
					my $s_t_prob       = 1 - (($t_prob + $in_t_prob)/$total); # CDS->STOP
					$states[$i]->t_matrix('cds0', $t_prob);
					$states[$i]->t_matrix($trans, $in_t_prob);
					$states[$i]->t_matrix('stop0', $s_t_prob);
				} else {
					my $trans = "D" . "$pos" . "_0";
					my $total     = tot_trans(\%{$trans_freq{$t}});
					$total       -= $trans_freq{$t}{8};
					my $t_prob    = $trans_freq{$t}{$t}/$total;
					$states[$i]->t_matrix($states[$i+1]->name, $t_prob)          if $order < 3;
					$states[$i]->t_matrix($states[$i+($order-1)]->name, $t_prob) if $order >= 3;
					$states[$i]->t_matrix($trans, (1 - $t_prob));
				}
			}
		} elsif ($st =~ /cds\d+_\d+/) {
			#print $st, "\t", $order, "\n";
			my ($cds_pos) = $st =~ /cds(\d+)_\d+/;
			my ($count)   = $st =~ /cds\d+_(\d+)/;
			if ($order == 2) {
				$cds_pos++   if $cds_pos != 2;
				$cds_pos = 0 if $cds_pos == 2;
				my $trans = "cds" . "$cds_pos";
				$states[$i]->t_matrix($trans, 1);
			} else {
				$states[$i]->t_matrix($states[$i+1]->name, 1);
			}
		} elsif ($st =~ /^[DA]/) {
			my ($pos)    = $st =~ /\w\d+_(\d+)$/;
			my ($in_pos) = $st =~ /\w(\d+)_\d+$/;
			if ($pos == $quant - 1) {
				if ($t == 5 and $st =~ /^D/) {
					my $trans = "i" . "$in_pos";
					$states[$i]->t_matrix($trans, 1);
				}
				if ($t == 7 and $st =~ /^A/) {
					if ($in_pos == $i_size - 1) {
						if ($pad_quant == 0) {
							$states[$i]->t_matrix('cds0', 1);
						} else {
							my $pad_pos = $pad_quant - 1;
							my $trans   = "cds" . "$in_pos" . "_" . "$pad_pos";
							$states[$i]->t_matrix($trans, 1);
						}
					} elsif ($st =~ /^A/) {
						if ($pad_quant == 0) {
							my $pad_pos = $in_pos + 1;
							my $trans   = "cds" . "$pad_pos";
							$states[$i]->t_matrix($trans, 1);
						} else {
							my $pad_pos = $pad_quant - 1;
							my $trans   = "cds" . "$in_pos" . "_" . "$pad_pos";
							$states[$i]->t_matrix($trans, 1);
						}
					}
				}
			} else {
				$states[$i]->t_matrix($states[$i+1]->name, 1);
			}
		} elsif ($st =~ /i/) {
			if ($t == 6) {
				$i_size = $quant;
				my ($pos)         = $st =~ /i(\d+)/;
				my $trans         = "A" . "$pos" . "_0";
				my $total         = tot_trans(\%{$trans_freq{$t}});
				my $t_prob        = $trans_freq{$t}{$t}/$total;
				$states[$i]->t_matrix($st, $t_prob);
				$states[$i]->t_matrix($trans, (1 - $t_prob));
			}
		}
	}
}
#Sanity Check for High Order CDS Transition Matrix
# foreach my $obj (@states) {
# 	my $st    = $obj->name;
# 	my $order = $obj->order;
# 	my $quant = $obj->st_quant;
# 	
# 	print "\n\n", $st, "\t", $order, "\n";
# 	browse($obj->t_matrix);
# }


#--------------------------------#
# Generate State Emission Counts #
#--------------------------------#
foreach my $id (keys %gene_struc) {
	my $ex_max     = $gene_struc{$id}{tot_exons};      # Total exons
	my $in_max     = $gene_struc{$id}{tot_introns};    # Total introns
	my $upstream   = $gene_struc{$id}{Upstream};       # Upstream sequence
	my $downstream = $gene_struc{$id}{Downstream};     # Downstream sequence
	my $first_ex   = $gene_struc{$id}{Exons}{1};       # First exon sequence
	my $last_ex    = $gene_struc{$id}{Exons}{$ex_max}; # Last exon sequence
	
	for (my $s=0; $s < scalar(@states); $s++) {
		my $st = $states[$s]->name;
		my $order = $states[$s]->order;
		my $quant = $states[$s]->st_quant;
		
		# Calc. Number of State Transitions
		my $tot_trans = 0;
		my $trans = $states[$s]->t_matrix;
		foreach my $t (keys %$trans) {$tot_trans++;}
		
		# Upstream
		if ($st =~ /GU/) {
			my $seq = $gene_struc{$id}{Upstream};
			$states[$s]->emission($st, $order, $tot_trans, $seq);
		}
		
		# Exons/Introns
		if ($ex_max == 1) {
			# KOGs with 1 CDS
			my $seq = "";
			if ($st =~ /start/) {
				$seq = substr($upstream, -$order, $order) if $order > 0;
				$seq .= substr($first_ex, 0, 3);
				$states[$s]->emission($st, $order, $tot_trans, $seq);
			} elsif ($st =~ /cds/) {
				if ($order != 0) {
					$seq   .= substr($first_ex, 3 - $order, $order);
					$seq   .= substr($upstream, -($order - 3), $order - 3) if $order > 3; 
				}
				$seq .= substr($first_ex, 3, length($first_ex) - 6);
				$states[$s]->emission($st, $order, $tot_trans, $seq);
			} elsif ($st =~ /stop/) {
				$seq .= substr($first_ex, -($order + 3), $order) if $order > 0;
				$seq .= substr($first_ex, -3, 3);
				$states[$s]->emission($st, $order, $tot_trans, $seq);
			}
		} else {
			# Start Codon
			if ($st =~ /start/) {
				my $seq = "";
				$seq .= substr($upstream, -$order, $order) if $order != 0;
				$seq .= substr($first_ex, 0, 3);
				$states[$s]->emission($st, $order, $tot_trans, $seq);
			}
			
			# Exon Body 
			if ($st =~ /cds\d+$/) {
				for (my $i=1; $i <= $ex_max; $i++) {
					my $seq = "";
					my $offset = $order - $pad_quant;
					my $exon = $gene_struc{$id}{Exons}{$i};
					if ($i == 1) {
						my $start = substr($exon, 0, 3);
						$seq .= substr($upstream, -$order, $order)   if $order > 3;
						$seq .= substr($start, -$order, $order)      if $order > 0 and $order <= 3;
						$seq .= substr($exon, 3, length($exon) - 3);
						$states[$s]->emission($st, $order, $tot_trans, $seq);
					} elsif ($i == $ex_max) {
						$seq .= substr($gene_struc{$id}{Introns}{$i-1}, -$offset, $offset) if $order != 0;
						$seq .= substr($exon, 0, length($exon) - 3);
						$states[$s]->emission($st, $order, $tot_trans, $seq);
					} else {
						$seq .= substr($gene_struc{$id}{Introns}{$i-1}, -$offset, $offset)   if $order != 0;
						$seq .= $exon;
						$states[$s]->emission($st, $order, $tot_trans, $seq);
					}
				}
			}
			if ($st =~ /cds\d+_\d+/) {
				for (my $i=1; $i <= $ex_max; $i++) {
					my $seq = "";
					my $exon = $gene_struc{$id}{Exons}{$i};
					if ($i != 1) {
						$seq .= substr($gene_struc{$id}{Introns}{$i-1}, -$order, $order);
						$seq .= substr($exon, 0, $pad_quant);
						$states[$s]->emission($st, $order, $tot_trans, $seq);
					}
				}
			} 
			
			# Intron
			if ($st =~ /^[DAi]/) {
				for (my $i=1; $i <= $in_max; $i++) {
					my $seq = "";
					my $exon    = $gene_struc{$id}{Exons}{$i};             # Full Exon
					my $intron  = $gene_struc{$id}{Introns}{$i};           # Full Intron
					my $don     = substr($intron, 0, $d_size);             # Donor Site
					my $accep   = substr($intron, -$a_size, $a_size);      # Acceptor Site
					my $in_len  = length($intron) - ($d_size + $a_size);   # Length of Intron Body (without Splice Sites)
					my $in_body = substr($intron, $d_size, $in_len);       # Intron Body
					
					if ($st =~ /^D/) {
						$seq .= substr($exon, -$order, $order) if $order != 0;
						$seq .= $don;
						$states[$s]->emission($st, $order, $tot_trans, $seq);
					} elsif ($st =~ /^i/) {
						if ($order != 0) {
							if ($order > $d_size) {
								my $offset = $order - $d_size;
								$seq .= substr($exon, -$offset, $offset);
								$seq .= $don;
							} else {
								$seq .= substr($don, -$order, $order);
							}
						}
						$seq .= $in_body;
						$states[$s]->emission($st, $order, $tot_trans, $seq);
					} elsif ($st =~ /^A/) {
						$seq .= substr($in_body, -$order, $order) if $order != 0;
						$seq .= $accep;
						$states[$s]->emission($st, $order, $tot_trans, $seq);
					}
				}
			}
			
			# Stop Codon
			if ($st =~ /stop/) {
				my $seq = "";
				$seq   .= substr($last_ex, -($order + 3), $order) if $order != 0;
				$seq   .= substr($last_ex, -3, 3);
				$states[$s]->emission($st, $order, $tot_trans, $seq);
			}
		}
		
		# Downstream
		if ($st =~ /GD/) {
			my $seq = "";
			$seq   .= substr($last_ex, -$order, $order) if $order != 0;
			$seq   .= $downstream;
			$states[$s]->emission($st, $order, $tot_trans, $seq);
		}
	}
}
# Sanity Check
# foreach my $obj (@states) {
# 	my $st    = $obj->name;
# 	my $order = $obj->order;
# 	my $quant = $obj->st_quant;
# 	
# 	print "\n\n", $st, "\t", $order, "\n";
# 	browse($obj->emission);
# }


#-------------------#
# Create Model File #
#-------------------#
my $DATE      = `date`; chomp($DATE);
my $ST_COUNT  = scalar(@states);
my $TRACEBACK = "[FUNCTION:	HMMER	TRACK:	SEQ	COMBINE_LABEL:	C	TO_LABEL:	U]";
my ($FH)      = $FASTA =~ /(\w+.\w+).fa/;
$FH          .= "_$ST_COUNT" . ".hmm";

open(OUT, ">$FH") or die "Could not write into OUT\n";
print OUT "#STOCHHMM MODEL FILE

<MODEL INFORMATION>
======================================================
NAME:	Genesmith
DESCRIPTION:	$ST_COUNT state gene model
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
	GU0:	1
";

foreach my $obj (@states) {
	my $st        = $obj->name;
	my $label     = $obj->label;
	my $order     = $obj->order;
	my $quant     = $obj->st_quant;
	my $trans     = $obj->t_matrix;
	my $em_counts = $obj->emission;
	
	print OUT "##################################################\n";
	print OUT "STATE:\n";
	print OUT "\tNAME:\t$st\n";
	print OUT "\tGFF_DESC:\t$st\n";
	print OUT "\tPATH_LABEL: $label\n";
	print OUT "TRANSITION:	STANDARD:	P(X)\n";
	foreach my $t (keys %$trans) {
		my $t_prob = $trans->{$t};
		print OUT "\t$t:\t$t_prob";
		if ($st =~ /^D/  and $t =~ /^i/) {print OUT "\t$TRACEBACK";}
		if ($st =~ /^GD/ and $st eq $t)  {print OUT "\t$TRACEBACK";}
		print OUT "\n";
	}
	# Get Final Emission
	my %em_rows = get_em_rows($order); 
	my $em_output = em_table($em_counts, \%em_rows, $order);
	
	print OUT "EMISSION:	SEQ:	COUNTS\n";
	print OUT "\tORDER:	$order	AMBIGUOUS:	AVG\n";
	print OUT $em_output;
}
print OUT "##################################################\n";
print OUT "//END\n";

close OUT;



#===================================================================#
# SUBROUTINES														#
#===================================================================#

# Changes default parameter based on command line input
sub edit_cmd_opt{
	my ($type, $DEFAULT, $cmd_opt) = @_;
	
	my $new_opt;
	if ($type =~ /EXON|DON|ACCEP/) {
		die "
		Incorrect Command Option\nFor extra help refer to [options]<-h>
		" unless $cmd_opt =~ /\d+:\d+/;
		$new_opt = $cmd_opt;
	} else {
		die "
		Incorrect Command Option\nFor extra help refer to [options]<-h>
		" unless $cmd_opt =~ /\d+/;
		($new_opt) = $DEFAULT =~ /^(\d+:)/;
		$new_opt  .= $cmd_opt;
	}
	return $new_opt;
}

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

# Get Names and Order of Padded CDS States
sub pad_states{
	my ($st, $order) = @_;
	my %pad_sts;
	my $pad_quant = $order - 2;
	
	for(my $i=0; $i < $pad_quant; $i++) {
		my $name          = $st . "_" . "$i";
		my $pad_order     = $order - (1 + $i);
		$pad_sts{$name}{order}    = $pad_order;
		$pad_sts{$name}{st_quant} = $pad_quant;
	}
	return \%pad_sts;
}

# Generates Labeled Sequence
sub label{
	my ($type, $seq) = @_;
	my $lab = "";
	
	if ($type eq 'first_ex') {
		my $start  = '1 ' x 3;
		my $ex_lab = '2 ' x (length($seq) - 3);
		$lab      .= $start . $ex_lab;
	} elsif ($type eq 'last_ex') {
		my $stop   = '8 ' x 3;
		my $ex_lab = '2 ' x (length($seq) - 3);
		$lab      .= $ex_lab . $stop;
	} elsif ($type eq 'only_ex') {
		my $start  = '1 ' x 3;
		my $stop   = '8 ' x 3;
		my $ex_lab = '2 ' x (length($seq) - 6);
		$lab      .= $start . $ex_lab . $stop;
	} elsif ($type eq 'ex') {
		my $ex_lab = '2 ' x length($seq);
		$lab      .= $ex_lab;
	} elsif ($type eq 'in') {
		my $don    = '5 ' x 2;
		my $accep  = '7 ' x 2;
		my $in_lab = '6 ' x (length($seq) - 4);
		$lab      .= $don . $in_lab . $accep;
	} elsif ($type eq 'up') {
		$lab = '0 ' x length($seq);
	} elsif ($type eq 'down') {
		$lab = '9 ' x length($seq);
	}
	return $lab;
}

# Gets the Total Transitions for a given state
sub tot_trans{
	my ($counts, $type, $quant) = @_;
	my $total;
	
	foreach my $trans (keys %$counts) {
		$total += $counts->{$trans};
	}
	return $total; 
}

# Generates Table of Emissions formatted for Model File
sub em_table{
	my ($em_counts, $em_rows, $order) = @_;
	my @alph = qw(A C G T);
	my $final_em = "";
	
	if ($order == 0) {
		my $em_row = "";
		foreach my $col (@alph) {
			my $count = $em_counts->{$col};
			if (defined $count) {$em_row .= "$count\t";}
			else                {$em_row .= "0\t";}
		}
		chop($em_row);
		$final_em .= "$em_row\n";
	} else {	
		foreach my $row (sort keys %$em_rows) {
			my $em_row = "";
			if (!defined $em_counts->{$row}) {
				$em_row = "0\t0\t0\t0";
			} else {
				foreach my $col (@alph) {
					my $count = $em_counts->{$row}->{$col};
					if (defined $count) {$em_row .= "$count\t";}
					else                {$em_row .= "0\t";}
				}
				chop($em_row);
			}
			$final_em .= "$em_row\n";
		}
	}
	return $final_em;
}

#-------------------------------------------------------#
# Creates Array of all kmer combinations that represent # 
# rows in the Emission tables for each state depending  #
# on order                                              #
#-------------------------------------------------------#
sub get_em_rows{
	my ($order) = @_;
	my @alph = qw(A C G T);
	my %em_rows;
	
	if ($order == 0) {
		$em_rows{$order} = 1;
	}
	
	if ($order == 1) {
		foreach my $c1 (@alph) {
			$em_rows{$c1} = 1;
		}	
	}
	
	if ($order == 2) {
		foreach my $c1 (@alph) {
			foreach my $c2 (@alph) {
				$em_rows{"$c1$c2"} = 1;
			}
		}	
	}
	
	if ($order == 3) {
		foreach my $c1 (@alph) {
			foreach my $c2 (@alph) {
				foreach my $c3 (@alph) {
					$em_rows{"$c1$c2$c3"} = 1;
				}
			}
		}	
	}
	
	if ($order == 4) {
		foreach my $c1 (@alph) {
			foreach my $c2 (@alph) {
				foreach my $c3 (@alph) {
					foreach my $c4 (@alph) {
						$em_rows{"$c1$c2$c3$c4"} = 1;
					}
				}
			}
		}	
	}
	
	if ($order == 5) {
		foreach my $c1 (@alph) {
			foreach my $c2 (@alph) {
				foreach my $c3 (@alph) {
					foreach my $c4 (@alph) {
						foreach my $c5 (@alph) {
							$em_rows{"$c1$c2$c3$c4$c5"} = 1;
						}
					}
				}
			}
		}	
	}
	return %em_rows;
}




__DATA__
=======================================================

#----------------#
# OPTION DETAILS #
#----------------#
-  Command line options customize the Gene Model
   by modifying the DEFAULT settings used to 
   build the HMM

-  DEFAULT
   -- Init Model Components
      - Upstream Seq   [-F]
      - Start Codon    [-s]
      - cds(Exon)      [-c]
      - Donor Site     [-D]
      - Intron body    [-i]
      - Acceptor Site  [-A]
      - Stop Site      [-e]
      - Downstream Seq [-T] 

	
   -- Values for each Component

      <path label>:<#states>:<$order>
 
      - path label  = Letter by which the state is labeled
      - #states     = Total number of states a given component
      - order       = Order assigned to component states

#-------------#
# EXAMPLE RUN #
#-------------#

hmmgen.pl -D 4:2 -s 1 C.elegans.gff C.elegans.fa

   * Changes for Donor Site
        - #states = 4
        - order   = 2
   * Changes for Start Codon
        - order = 1

=======================================================

