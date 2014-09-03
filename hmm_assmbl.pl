#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use STATE;
use HMMstar;
use DataBrowser;
use Getopt::Std;
use vars qw($opt_h $opt_1 $opt_S $opt_t $opt_5 $opt_m $opt_c $opt_d $opt_i $opt_a $opt_s $opt_3 $opt_b $opt_B $opt_D $opt_A $opt_U $opt_E);
getopts('h1St:5:m:c:d:i:a:s:3:b:BD:A:U:E:');

# Change options to use Getopt::Long to allow an input of 0 for branch and polyA options
# DEFAULT settings [options]
my $UP      = 0;     # Order Upstream
my $START   = 0;     # Order Start, includes canonical ATG at the end, length = 3 + order
my $EXON    = 0;     # Order CDS, always 3 states
my $DON     = 0;     # Order Donor, starts with canonical GT, length = 2 + order
my $INTRON  = 0;     # Order Intron, always 3 states
my $ACCEP   = 0;     # Order Acceptor, ends with canonical AG, length = 2 + order
my $STOP    = 0;     # Order Stop, starts with stop codon, length = 3 + order
my $DOWN    = 0;     # Order Downstream
my $BRANCH  = 0;     # Order for the branch motif states
my $L_DON   = 5;     # Quantity of Donor States
my $L_ACCEP = 10;    # Quantity of Acceptor States 
my $L_UP    = 500;   # Length of Upstream region parsed for training
my $L_DOWN  = 500;   # Length of Downstream region parsed for training
my $TRANS   = 1.0;   # Weight multiplied to CDS transitions

die "
usage: $0 [options] <GFF> <FASTA>

universal parameters:
  -5 <order>   upstream state info                        Default = $UP
  -m <order>   start codon state info                     Default = $START
  -c <order>   coding state info                          Default = $EXON
  -d <order>   donor site state info                      Default = $DON
  -i <order>   intron body state info                     Default = $INTRON
  -a <order>   acceptor site state info                   Default = $ACCEP
  -s <order>   stop codon state info                      Default = $STOP
  -3 <order>   downstream state info                      Default = $DOWN
  -b <order>   branch motif state info                    Default = $BRANCH
  -B           Option to include branch states
  -D <length>  Donor Site Length                          Default = $L_DON
  -A <length>  Acceptor Site Length                       Default = $L_ACCEP
  -U <length>  upstream training                          Default = $L_UP
  -E <length>  downtream training                         Default = $L_DOWN
  -t <digit>   Weight for CDS transitions                 Default = $TRANS
  -1           Convert for Standard to Basic HMM with 1 CDS state
  -S           Option to exclude start and stop states from basic model
  -h           help (format and details)
" unless @ARGV == 2;

$UP	     = $opt_5 if $opt_5;
$START   = $opt_m if $opt_m;
$EXON    = $opt_c if $opt_c;
$DON     = $opt_d if $opt_d;
$INTRON  = $opt_i if $opt_i;
$ACCEP   = $opt_a if $opt_a;
$STOP    = $opt_s if $opt_s;
$DOWN    = $opt_3 if $opt_3;
$BRANCH  = $opt_b if $opt_b;
$L_DON   = $opt_D if $opt_D;
$L_ACCEP = $opt_A if $opt_A;
$L_UP    = $opt_U if $opt_U;
$L_DOWN  = $opt_E if $opt_E;
$TRANS   = $opt_t if $opt_t;

my ($GFF, $FASTA) = @ARGV;

# Sanity Check - options
if ($START   !~ /^\d+$/    or
    $DON     !~ /^\d+$/    or 
    $ACCEP   !~ /^\d+$/    or
    $STOP    !~ /^\d+$/    or
    $UP      !~ /^\d+$/    or
    $EXON    !~ /^\d+$/    or
    $INTRON  !~ /^\d+$/    or
    $DOWN    !~ /^\d+$/    or
    $BRANCH  !~ /^\d+$/    or
    $L_DON   !~ /^\d+$/    or
    $L_ACCEP !~ /^\d+$/    or    
    $L_UP    !~ /^\d+$/    or
    $L_DOWN  !~ /^\d+$/    or
    $TRANS   !~ /\d?.?\d+/ or
    $opt_S and !$opt_1       ) {
    die "Invalid input [options]\n";
}

# OPTIONAL states
my $cds_name  = 'cds';
my $cds_quant = 3;                   # Quantity of CDS states for Standard Model
my $i_quant   = 3;                   # Quantity of Intron states for Standard Model
$cds_name     = 'CDS'     if $opt_1 and !$opt_S;
$cds_name     = 'CDSfull' if $opt_1 and $opt_S;
$cds_quant    = 1         if $opt_1;     # 1 CDS for the Basic model [-1]
$i_quant      = 1         if $opt_1;     # 1 Intron state for the Basic model [-1]
my $MOTIF = "TACTAAC";               # Branch Motif discovered in S.cerevisiae only

# PATH of trained State Emissions
my ($result_dir) = $GFF =~ /\/?(\w+\.*\w+\d*)\.gff$/;
my $path = "./HMM/$result_dir\/";
die "Emissions directory PATH does not exist\n" if !-d $path;


# Info for each Group of States
my @st_order  = ($UP,  $START,  $EXON,       $DON,   $INTRON,  $ACCEP,  $STOP,  $DOWN);
my @st_name   = ('GU', 'start', $cds_name,   'don',  'i',      'accep', 'stop', 'GD');
my @st_label  = ('U',  'C',     'C',         'I',    'I',      'I',     'C',    'D');
my @st_quant  = ( 1,    3,       $cds_quant, $L_DON, $i_quant, $L_ACCEP, 3,      1);


#----------------------------------------------------------#
# Create Array of State Objects and Hash of Emission Files #
#----------------------------------------------------------#
my @states;       # array of state objects

for (my $i=0; $i < @st_quant; $i++) {
	my $order = $st_order[$i];
	my $quant = $st_quant[$i];
	for (my $c=0; $c < $quant; $c++) {
		my $name = $st_name[$i] . "$c";
		my $st_length;
		if  ($name =~ /start|stop/) {
			$st_length  = $quant + $order;
			my $st_obj  = new STATE($name, $st_label[$i], $order, $st_length);
			push(@states, $st_obj) if !$opt_S;
		} elsif ($name =~ /don|accep/) {
			# Produce all three sets of Donor and Acceptor States
			$st_length  = $quant + $order;
			for (my $r=0; $r < $i_quant; $r++) {
				my $st_name = $name . "_$r";
				my $st_obj  = new STATE($st_name, $st_label[$i], $order, $st_length);
				push(@states, $st_obj);
			}
		} else {
			$st_length = $quant;
			# Upstream
			if ($name =~ /GU/) {
				my $st_obj = new STATE($name, $st_label[$i], $order, $st_length);
				push(@states, $st_obj);
			}
			
			# CDS States
			if ($name =~ /cds|CDS/) {
				my $st_obj = new STATE($name, $st_label[$i], $order, $st_length);
				push(@states, $st_obj);
			}
			
			# Intron States, Branch States - Conditional
			if ($name =~ /i/ and $opt_B) {
				# Intron Body States before and after Branch Motif
				my @in_branch = qw(branch_i0 branch_i1 branch_i2);
				foreach my $st (@in_branch) {
					my $st_name = $st . "_$c";
					my $obj = new STATE($st_name, 'I', $INTRON, 1);
					push(@states, $obj);
				}
			
				# Branch Motif States
				for (my $m=0; $m < length($MOTIF); $m++) {
					my $st_name = "branch$m\_$c";
					my $st_obj  = new STATE($st_name, 'I', $BRANCH, 1);
					push (@states, $st_obj);
				}
			} elsif ($name =~ /i/ and !$opt_B) {
				my $st_obj = new STATE($name, $st_label[$i], $order, $st_length);
				push(@states, $st_obj);
			}
			
			# Downstream
			if ($name =~ /GD/) {
				my $st_obj = new STATE($name, $st_label[$i], $order, $st_length);
				push(@states, $st_obj);
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

foreach my $contig ($genome->contigs) {
	foreach my $cds ($contig->coding_sequences) {
		my $id        = $cds->name;
		my $ex_count  = 1;
		my $in_count  = 1;
		
		# Upstream
		my $up_start = $cds->start_site->start - $L_UP;
		$up_start = 0 if $up_start < 0;		
		$gene_struc{$id}{'Upstream'} = uc(substr($contig->dna->sequence, $up_start, $cds->start_site->start - $up_start - 1));
		
		# Exons
		foreach my $exon ($cds->exons) {
			$gene_struc{$id}{'Exons'}{$ex_count} = uc($exon->sequence);
			$gene_struc{$id}{'ex_total'}         = $ex_count;
			$ex_count++;
		}
		# Introns
		foreach my $intron ($cds->introns) {
			# some KOGs only have 1 CDS and no Introns
			if (defined $intron->sequence) {
				$gene_struc{$id}{'Introns'}{$in_count} = uc($intron->sequence);
				$gene_struc{$id}{'in_total'}		   = $in_count;
				$in_count++;
			}
		}
		# Downstream 
		my $dn_beg = $cds->stop_site->end + 2;
		my $dn_len = length($contig->dna->sequence) - $dn_beg +1;
		$dn_len = $L_DOWN if $dn_len > $L_DOWN;
		$gene_struc{$id}{'Downstream'} = uc(substr($contig->dna->sequence, $dn_beg, $dn_len));
	}
}


#----------------------------------------------------------------------#
# Calculate Transition Probabilities                                   #
#----------------------------------------------------------------------#
# Label Gene Seqs based on Structure                                   #
# -  label key:                                                        #
#    GU = 0, start = 1, cds = 2, D = 5, i = 6, A = 7, stop = 8, GD = 9 #
# -  Optional States                                                   #
#    * branch = b, branch_i0 = s, branch_i1 = e, branch_i2 = i         #
#----------------------------------------------------------------------#
my %trans_freq;
my @trans_type = qw(Upstream Exons Downstream); 

# Label Sequence
foreach my $id (keys %gene_struc) {
	my $seq    = "";                             # Gene sequence
	my $labseq = "";                             # Labeled Gene sequence
	my $ex_max = $gene_struc{$id}{'ex_total'};   # Total exons
	my $in_max = $gene_struc{$id}{'in_total'};   # Total introns
	
	foreach my $type (@trans_type) {
		if ($type eq 'Upstream') {
			my $upstr_seq = $gene_struc{$id}{$type};
			$labseq      .= label($type, $upstr_seq, $L_DON, $L_ACCEP);
			$seq         .= $upstr_seq;
		}
		if ($type eq 'Exons') {
			for (my $i=1; $i <= $ex_max; $i++) {
				my $ex_seq  = $gene_struc{$id}{$type}{$i};
				my $in_seq  = $gene_struc{$id}{'Introns'}{$i};
				
				if (!defined $in_max) {
					my $ex_lab   = label('only_ex', $ex_seq, $L_DON, $L_ACCEP);
					$labseq     .= $ex_lab;
					$seq        .= $ex_seq;
				} elsif (defined $in_max) {
					if ($i == 1) {
						my $don     = substr($in_seq, 0, $L_DON);                      # Donor Site
						my $accep   = substr($in_seq, -$L_ACCEP, $L_ACCEP);            # Acceptor Site
						my $in_len  = length($in_seq) - ($L_DON + $L_ACCEP);           # Length of Intron Body
						my $in_body = substr($in_seq, $L_DON, $in_len);                # Intron body without Donor and Acceptor Sites
						
						my $ex_lab    = label('first_ex', $ex_seq, $L_DON, $L_ACCEP);
						my $don_lab   = label('Donor', $don, $L_DON, $L_ACCEP);
						my $accep_lab = label('Acceptor', $accep, $L_DON, $L_ACCEP);

						my ($br_start, $br_end) = get_branch_coord($MOTIF, $in_body);
						if (defined $br_start and defined $br_end and $opt_B) {
							my $in_body1_len = $br_start;
							my $motif_len    = ($br_end + 1) - $br_start;
							my $in_body2_len = length($in_body) - ($br_end + 1);
					
							my $pre_motif  = substr($in_body, 0, $in_body1_len);
							my $motif      = substr($in_body, $br_start, $motif_len);
							my $post_motif = substr($in_body, $br_end+1, $in_body2_len);
							
							my $pre_inlab  = label('branch_i0', $pre_motif, $L_DON, $L_ACCEP);
							my $branch_lab = label('branch', $motif, $L_DON, $L_ACCEP);
							my $post_inlab = label('branch_i1', $post_motif, $L_DON, $L_ACCEP);
							$labseq       .= $ex_lab . $don_lab . $pre_inlab . $branch_lab . $post_inlab . $accep_lab;
							$seq          .= $ex_seq . $in_seq;
						} else {
							my $in_lab;
							$in_lab       = label('Introns', $in_body, $L_DON, $L_ACCEP)   if !$opt_B;
							$in_lab       = label('branch_i2', $in_body, $L_DON, $L_ACCEP) if $opt_B;
							$labseq      .= $ex_lab . $don_lab . $in_lab . $accep_lab;
							$seq         .= $ex_seq . $in_seq;
						} 
					} elsif ($i == $ex_max) {
						my $ex_lab   = label('last_ex', $ex_seq, $L_DON, $L_ACCEP);
						$labseq     .= $ex_lab;
						$seq        .= $ex_seq;
					} else {
						my $don     = substr($in_seq, 0, $L_DON);                      # Donor Site
						my $accep   = substr($in_seq, -$L_ACCEP, $L_ACCEP);            # Acceptor Site
						my $in_len  = length($in_seq) - ($L_DON + $L_ACCEP);           # Length of Intron Body
						my $in_body = substr($in_seq, $L_DON, $in_len);                # Intron body without Donor and Acceptor Sites
						my ($br_start, $br_end) = get_branch_coord($MOTIF, $in_body);

						my $ex_lab     = label('Exons', $ex_seq, $L_DON, $L_ACCEP);
						my $don_lab    = label('Donor', $don, $L_DON, $L_ACCEP);
						my $in_lab;
						$in_lab        = label('Introns', $in_body, $L_DON, $L_ACCEP)   if !$opt_B;
						$in_lab        = label('branch_i2', $in_body, $L_DON, $L_ACCEP) if $opt_B;
						my $accep_lab  = label('Acceptor', $accep, $L_DON, $L_ACCEP);
						$labseq       .= $ex_lab . $don_lab . $in_lab . $accep_lab;
						$seq          .= $ex_seq . $in_seq;
					}
				}
			}
		}
		if ($type eq 'Downstream') {
			my $downstr_seq = $gene_struc{$id}{$type};
			$labseq        .= label($type, $downstr_seq, $L_DON, $L_ACCEP);
			$seq           .= $downstr_seq;
		}
	}
	# Store Transition Counts
	$labseq =~ tr/18/22/ if $opt_S;
	my @labels = split('', $labseq);
	push(@labels, 'END');
	for (my $i=1; $i < @labels; $i++) {
		my $lab = $labels[$i-1];
		my $nxt_lab = $labels[$i];
		$trans_freq{$lab}{$nxt_lab}++;
	}
}
# browse(\%trans_freq);


#----------------------------------------------------------------------#
# Transition Matrix                                                    #
#----------------------------------------------------------------------#
# -  label key:                                                        #
#    GU = 0, start = 1, cds = 2, D = 5, i = 6, A = 7, stop = 8, GD = 9 #
# -  Optional States                                                   #
#    * branch = b, branch_i0 = s, branch_i1 = e, branch_i2 = i         #
#----------------------------------------------------------------------#
my @lb_list;
@lb_list = qw(0 1 2 3 4 5 6 7 8 9)       if !$opt_B;
@lb_list = qw(0 1 2 3 4 5 i s b e 7 8 9) if $opt_B;

for (my $i=0; $i < @states; $i++) {
	my $st    = $states[$i]->name;
	my $order = $states[$i]->order;
	my $quant = $states[$i]->st_length;
	
	# Iterate through transition labels
	foreach my $label (@lb_list) {
		if ($label =~ /0/ and $st =~ /GU/) {
			my $total     = tot_trans(\%{$trans_freq{$label}});
			my $t_prob    = $trans_freq{$label}{$label}/$total;
			my $cds_trans = (1 - $t_prob)*$TRANS;
			$states[$i]->t_matrix($st, (1 - $cds_trans));
			$states[$i]->t_matrix($states[$i+1]->name, $cds_trans);
		}
		if ($label =~ /9/ and $st =~ /GD/) {
			my $total  = tot_trans(\%{$trans_freq{$label}});
			my $t_prob = $trans_freq{$label}{$label}/$total;
			$states[$i]->t_matrix($st, $t_prob);
			$states[$i]->t_matrix('END', (1 - $t_prob));
		}
		if ($label =~ /1/ and $st =~ /start/ and !$opt_S) {
			$states[$i]->t_matrix($states[$i+1]->name, 1);
		}
		if ($label =~ /8/ and $st =~ /stop/ and !$opt_S) {
			$states[$i]->t_matrix($states[$i+1]->name, 1);
		}
		# Standard Model with 3 CDS states
		if ($label =~ /2/ and $st =~ /cds\d+/ and $cds_quant == 3) {
			my ($pos) = $st =~ /cds(\d+)/;
			my $cds_trans  = $trans_freq{$label}{$label};
			my $don_trans  = $trans_freq{$label}{'5'};
			my $stop_trans = $trans_freq{$label}{'8'};
			
			if ($pos == 2) {
				my $tot_3trans  = $cds_trans + $don_trans + $stop_trans;
				$states[$i]->t_matrix('cds0',   ($cds_trans/$tot_3trans));
				$states[$i]->t_matrix('don0_2', ($don_trans/$tot_3trans));
				$states[$i]->t_matrix('stop0',  ($stop_trans/$tot_3trans));
			} else {
				my $tot_2trans  = $cds_trans + $don_trans;
				my $don_st = "don0_$pos";
				$pos++;
				my $cds_st = "cds$pos";
				$states[$i]->t_matrix($cds_st, ($cds_trans/$tot_2trans));
				$states[$i]->t_matrix($don_st, ($don_trans/$tot_2trans));
			}
		}
		# Basic Model with 1 CDS state
		if ($label =~ /2/ and $st =~ /^CDS\d+/ and $cds_quant == 1) {
			my ($pos) = $st =~ /CDS(\d+)/;
			my $cds_trans   = $trans_freq{$label}{$label};
			my $don_trans   = $trans_freq{$label}{'5'};
			my $stop_trans  = $trans_freq{$label}{'8'};
			my $total_trans = $cds_trans + $don_trans + $stop_trans;
			$states[$i]->t_matrix('stop0', ($stop_trans/$total_trans));
			$states[$i]->t_matrix($st,      ($cds_trans/$total_trans));
			$states[$i]->t_matrix('don0_0', ($don_trans/$total_trans));
		}
		
		# Basic Model with 1 CDS state and no Start/Stop states
		if ($label =~ /2/ and $st =~ /^CDSfull\d+/ and $cds_quant == 1) {
			my ($pos) = $st =~ /CDSfull(\d+)/;
			my $cds_trans     = $trans_freq{$label}{$label};
			my $don_trans     = $trans_freq{$label}{'5'};
			my $downstr_trans = $trans_freq{$label}{'9'};
			my $total_trans   = $cds_trans + $don_trans + $downstr_trans;
			$states[$i]->t_matrix('GD0', ($downstr_trans/$total_trans));
			$states[$i]->t_matrix($st,      ($cds_trans/$total_trans));
			$states[$i]->t_matrix('don0_0', ($don_trans/$total_trans));
		}

		
		# Process BRANCH state Transitions
		if (!$opt_B) {
			if ($label =~ /5/ and $st =~ /don/) {
				my ($pos) = $st =~ /don(\d+)_\d+/;
				my ($set) = $st =~ /don\d+_(\d+)/;
				if ($pos == $L_DON - 1) {
					my $st_name = "i$set";
					$states[$i]->t_matrix($st_name, 1);
				} else {
					$pos++;
					my $st_name = "don$pos\_$set";
					$states[$i]->t_matrix($st_name, 1); 
				}
			}
			if ($label =~ /6/ and $st =~ /i/) {
				my ($set) = $st =~ /i(\d+)/;
				my $total  = tot_trans(\%{$trans_freq{$label}});
				my $t_prob = $trans_freq{$label}{$label}/$total;
				$states[$i]->t_matrix($st, $t_prob);
				my $accep_name = "accep0_$set";
				$states[$i]->t_matrix($accep_name, (1 - $t_prob));
			}
		} elsif ($opt_B) {
			if ($label =~ /5/ and $st =~ /don/) {
				my ($pos) = $st =~ /don(\d+)_\d+/;
				my ($set) = $st =~ /don\d+_(\d+)/;
				my $br_in_trans = $trans_freq{$label}{'s'};
				my $in_trans    = $trans_freq{$label}{'i'};
				my $total       = $br_in_trans + $in_trans;
				my $st_name;
				
				if ($pos == $L_DON -1) {
					$st_name = "branch_i0_$set"; $states[$i]->t_matrix($st_name, ($br_in_trans/$total));
					$st_name = "branch_i2_$set"; $states[$i]->t_matrix($st_name, ($in_trans/$total));
				} else {
					$pos++;
					$st_name = "don$pos\_$set";
					$states[$i]->t_matrix($st_name, 1);
				}
			}
			# Branch intron/motif States
			if ($label =~ /s/ and $st =~ /branch_i0/) {
				my ($set)  = $st =~ /branch_i0_(\d+)/;
				my $total  = tot_trans(\%{$trans_freq{$label}});
				my $t_prob = $trans_freq{$label}{$label}/$total;
				my $nxt_st = "branch0_$set";
				$states[$i]->t_matrix($st, $t_prob);
				$states[$i]->t_matrix($nxt_st, (1 - $t_prob));
			}
			if ($label =~ /e/ and $st =~ /branch_i1/) {
				my ($set) = $st =~ /branch_i1_(\d+)/;
				my $total  = tot_trans(\%{$trans_freq{$label}});
				my $t_prob = $trans_freq{$label}{$label}/$total;
				my $nxt_st = "accep0_$set";
				$states[$i]->t_matrix($st, $t_prob);
				$states[$i]->t_matrix($nxt_st, (1 - $t_prob));
			}
			if ($label =~ /i/ and $st =~ /branch_i2/) {
				my ($set) = $st =~ /branch_i2_(\d+)/;
				my $total  = tot_trans(\%{$trans_freq{$label}});
				my $t_prob = $trans_freq{$label}{$label}/$total;
				my $nxt_st = "accep0_$set";
				$states[$i]->t_matrix($st, $t_prob);
				$states[$i]->t_matrix($nxt_st, (1 - $t_prob));
			}
			if ($label =~ /b/ and $st =~ /branch\d+_\d+/) {
				my ($pos) = $st =~ /branch(\d+)_\d+/;
				my ($set) = $st =~ /branch\d+_(\d+)/;
				if ($pos == length($MOTIF) - 1) {
					my $nxt_st = "branch_i1_$set";
					$states[$i]->t_matrix($nxt_st, 1);
				} else {
					$pos++;
					my $nxt_st = "branch$pos\_$set";
					$states[$i]->t_matrix($nxt_st, 1);
				}
			}
		}
		if ($label =~ /7/ and $st =~ /accep/) {
			my ($pos) = $st =~ /accep(\d+)_\d+/;
			my ($set) = $st =~ /accep\d+_(\d+)/;
			if ($pos == $L_ACCEP - 1) {
				if ($set == 2) {
					$states[$i]->t_matrix('cds0', 1);
				} else {
					my $st_name;
					if ($cds_quant == 3) {
						$set++;
						$st_name = "cds$set";
						$states[$i]->t_matrix($st_name, 1);
					} else {
						$st_name = "CDS$set"     if $opt_1 and !$opt_S;
						$st_name = "CDSfull$set" if $opt_1 and $opt_S;
						$states[$i]->t_matrix($st_name, 1);
					}
				}
			} else {
				$pos++;
				my $st_name = "accep$pos\_$set";
				$states[$i]->t_matrix($st_name, 1); 
			}
		}
	}
}
#Sanity Check for High Order CDS Transition Matrix
# foreach my $obj (@states) { 
# 	my $st    = $obj->name;
# 	my $order = $obj->order;
# 	my $quant = $obj->st_length;
# 	print "\n\n", $st, "\t", $order, "\n";
# 	browse($obj->t_matrix);
# }


#-------------------#
# Create Model File #
#-------------------#
my $DATE      = `date`; chomp($DATE);
my $ST_COUNT  = scalar(@states);
my $TRACEBACK = "[FUNCTION:	TRANSLATE	TRACK:	SEQ	COMBINE_LABEL:	C	TO_LABEL:	U]";

# my $FH        = $result_dir . "_$ST_COUNT" . ".hmm";
# open(OUT, ">$FH") or die "Could not write into OUT\n";

print "#STOCHHMM MODEL FILE

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
	my $quant     = $obj->st_length;
	my $trans     = $obj->t_matrix;
	my $em_counts = $obj->emission;
	
	print "##################################################\n";
	print "STATE:\n";
	print "\tNAME:\t$st\n";
	print "\tGFF_DESC:\t$st\n";
	print "\tPATH_LABEL: $label\n";
	print "TRANSITION:	STANDARD:	P(X)\n";
	foreach my $t (keys %$trans) {
		my $t_prob = $trans->{$t};
		print "\t$t:\t$t_prob";
		if ($st =~ /^don/  and $t =~ /i|branch_i0|branch_i2/) {print "\t$TRACEBACK";}
		if ($st =~ /^GD/ and $st eq $t)                       {print "\t$TRACEBACK";}
		print "\n";
	}
	# Get Final Emission
	print "EMISSION:	SEQ:	COUNTS\n";
	print "\tORDER:	$order	AMBIGUOUS:	AVG\n";
	my $em_path = $path . $st;
	if    ($st =~ /CDS\d+/     and !$opt_S) {$em_path .= "-1";}
	elsif ($st =~ /CDSfull\d+/ and $opt_S)  {$em_path .= "-1S";}
	else                                    {$em_path .= "-std";}
	$em_path   .= "-$L_UP"             if $st =~ /GU/;
	$em_path   .= "-$L_DOWN"           if $st =~ /GD/;
	$em_path   .= "-$L_DON"            if $st =~ /don/;
	$em_path   .= "-$L_ACCEP"          if $st =~ /accep/;
	$em_path   .= "-$L_DON\_$L_ACCEP"  if $st =~ /i|branch/;
	$em_path   .= "-$label\-$order\.txt";
	
	open (EM, "<$em_path") or die "Error reading $em_path\n";
	while (my $line = <EM>) {
		print $line;
	}
}
print "##################################################\n";
print "//END\n";

# close OUT;



#===================================================================#
# SUBROUTINES   												    #
#===================================================================#

# Get Position Coordinates for Branch Motif
sub get_branch_coord{
	my ($regex, $inseq) = @_;
	return if not $inseq =~ /$regex/;
	return ($-[0], $+[0]);
}

# Gets the Total Transitions for a given state
sub tot_trans{
	my ($counts) = @_;
	my $total;
	
	foreach my $trans (keys %$counts) {
		$total += $counts->{$trans};
	}
	return $total; 
}


# Generates labeled sequence in the context of the basic gene model
sub label {
	my ($type, $seq, $L_DON, $L_ACCEP) = @_;
	my $labseq = "";
	
	if ($type eq 'only_ex') {
		my $start  = '1' x 3;
		my $stop   = '8' x 3;
		my $ex_lab = '2' x (length($seq) - 6);
		$labseq   .= $start . $ex_lab . $stop;
	} elsif ($type eq 'first_ex') {
		my $start  = '1' x 3;
		my $ex_lab = '2' x (length($seq) - 3);
		$labseq   .= $start . $ex_lab;
	} elsif ($type eq 'last_ex') {
		my $stop   = '8' x 3;
		my $ex_lab = '2' x (length($seq) - 3);
		$labseq   .= $ex_lab . $stop;
	}
	
	if ($type eq 'Upstream')   {$labseq = '0' x length($seq);}
	if ($type eq 'Exons')      {$labseq = '2' x length($seq);}
	if ($type eq 'Donor')      {$labseq = '5' x $L_DON;}
	if ($type eq 'Acceptor')   {$labseq = '7' x $L_ACCEP;}
	if ($type eq 'Introns')    {$labseq = '6' x length($seq);}
	if ($type eq 'Downstream') {$labseq = '9' x length($seq);}
	
	### Branch Labels
	if ($type eq 'branch_i0')  {$labseq = 's' x length($seq);}
	if ($type eq 'branch_i1')  {$labseq = 'e' x length($seq);}
	if ($type eq 'branch_i2')  {$labseq = 'i' x length($seq);}
	if ($type eq 'branch')     {$labseq = 'b' x length($seq);}
	
	return $labseq;
}

