#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use STATE;
use HMMstar;
use DataBrowser;
use Getopt::Std;
use vars qw($opt_h $opt_i $opt_b $opt_D $opt_A);
getopts('h:i:b:D:A:');


my $INTRON  = 0;
my $BRANCH  = 0;
my $L_DON   = 2;
my $L_ACCEP = 2;

die "
usage: $0 [options] <GFF> <FASTA>

universal parameters:
  -i <order>   intron body states info  Default = $INTRON
  -b <order>   Branch states info       Default = $BRANCH
  -D <length>  Donor site length        Default = $L_DON
  -A <length> Acceptor site length      Defualt = $L_ACCEP       
  -h           help (format and details)
" unless @ARGV == 2;

$INTRON  = $opt_i if $opt_i;
$BRANCH  = $opt_b if $opt_b;
$L_DON   = $opt_D if $opt_D;
$L_ACCEP = $opt_A if $opt_A;
my ($GFF, $FASTA) = @ARGV;

if ($INTRON  !~ /^\d+$/ or 
    $BRANCH  !~ /^\d+$/ or
    $L_DON   !~ /^\d+$/ or
    $L_ACCEP !~ /^\d+$/   ) {
	die "Invalid input [options]\n";
}

#-------------------------------#
# Create array of State Objects #
#-------------------------------#
my @states;
my $motif = "TACTAAC";

# Intron Body States before and after Branch Motif
my $in_body1    = "branch_i0";
my $in_body2    = "branch_i1";
my $in_nobranch = "branch_i2";
my $in_obj1 = new STATE($in_body1,    'I', $INTRON, 1);
my $in_obj2 = new STATE($in_body2,    'I', $INTRON, 1);
my $in_obj3 = new STATE($in_nobranch, 'I', $INTRON, 1);
push (@states, $in_obj1);
push (@states, $in_obj2);
push (@states, $in_obj3);

# Branch Motif States
for (my $m=0; $m < length($motif); $m++) {
	my $st_name = "branch$m";
	my $st_obj  = new STATE($st_name, 'I', $BRANCH, 1);
	push (@states, $st_obj);
}


#-------------------------------#
# Extract Intron/Exon Structure #
#-------------------------------#
my $genome = new HMMstar::Genome($FASTA, $GFF);
my %gene_struc;

foreach my $contig ($genome->contigs) {
	foreach my $cds ($contig->coding_sequences) {
		my $id        = $cds->name;
		my $ex_count  = 1;
		my $in_count  = 1;
		
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
	}
}

#--------------------------------#
# Generate State Emission Counts #
#--------------------------------#
foreach my $id (keys %gene_struc) {
	my $ex_max = $gene_struc{$id}{'ex_total'};       # Total exons
	my $in_max = $gene_struc{$id}{'in_total'};       # Total introns
	if ($ex_max > 1) {
		for (my $s=0; $s < @states; $s++) {
			my $st    = $states[$s]->name;
			my $order = $states[$s]->order;
			
			# Introns
			for (my $i=1; $i <= $in_max; $i++) {
				my $exon    = $gene_struc{$id}{'Exons'}{$i};
				my $intron  = $gene_struc{$id}{'Introns'}{$i};
				
				my $don     = substr($intron, 0, $L_DON);                      # Donor Site
				my $accep   = substr($intron, -$L_ACCEP, $L_ACCEP);            # Acceptor Site
				my $in_len  = length($intron) - ($L_DON + $L_ACCEP);           # Length of Intron Body
				my $in_body = substr($intron, $L_DON, $in_len);                # Intron body without Donor and Acceptor Sites
				my ($br_start, $br_end) = get_branch_coord($motif, $in_body);
				if (defined $br_start and defined $br_end) {
					my $in_body1_len = $br_start;
					my $motif_len    = ($br_end + 1) - $br_start;
					my $in_body2_len = length($in_body) - ($br_end + 1);
					
					my $pre_motif  = substr($in_body, 0, $in_body1_len);
					my $motif      = substr($in_body, $br_start, $motif_len);
					my $post_motif = substr($in_body, $br_end+1, $in_body2_len);
					my $seq = "";
					if ($st =~ /branch_i0/) {
						my $prev_seq = $exon . $don;
						my $seq     .= substr($prev_seq, -$order, $order);
						$seq        .= $pre_motif;
						$states[$s]->emission($st, $order, $seq);
					} elsif($st =~ /branch_i1/) {
						$seq   .= substr($motif, -$order, $order);
						$seq   .= $post_motif;
						$states[$s]->emission($st, $order, $seq);
					} else {
						$seq .= substr($pre_motif, -$order, $order);
						$seq .= $motif;
						$states[$s]->emission($st, $order, $seq);
					}
				} else {
					if ($st =~ /branch_i2/) {
						my $prev_seq = $exon . $don;
						my $seq     .= substr($prev_seq, -$order, $order);
						$seq        .= $in_body;
						$states[$s]->emission($st, $order, $seq);
					}
				}
			}
		}
	}
}

# Sanity Check
# foreach my $obj (@states) {
# 	print "\n", $obj->name,  "\n> ",
# 	            $obj->label, "\n> ",
# 	            $obj->order, "\n";
# 	browse($obj->emission);
# }


#------------------------------------#
# Create State Emission Output Files #
#------------------------------------#

# Create Output Directory if it does not Exist
my ($result_dir) = $GFF =~ /\/?(\w+\.*\w+\d*)\.gff$/;
my $path = "./HMM/$result_dir\/";
if (!-d $path) {
	`mkdir $path`;
	print "\n$path\tDIRECTORY CREATED\n";
} else {
	print "\n$path\tDIRECTORY EXISTS\n";
}

foreach my $obj (@states) {
	my $st        = $obj->name;
	my $label     = $obj->label;
	my $order     = $obj->order;
	my $em_counts = $obj->emission;
	my %em_rows   = get_em_rows($order);
	my $em_output = em_table($em_counts, \%em_rows, $order);
	
	for (my $i=0; $i < 3; $i++) {
		my $st_name = $st . "_$i";
		$st_name   .= "-" . $label ."-" . $order . ".txt";
		my $fh = $path . $st_name;
		output_st_em($fh, $em_output);
	}
}


#===============================#
#  SUBROUTINES                  #
#===============================#

# Get Position Coordinates for Branch Motif
sub get_branch_coord{
	my ($regex, $inseq) = @_;
	return if not $inseq =~ /$regex/;
	return ($-[0], $+[0]);
}

# Output state emissions in a text file
sub output_st_em {
	my ($fh, $em_output) = @_;
	open (OUT, ">$fh") or die "Error: cannot write into $fh\n";
	print OUT $em_output;
	close OUT;
	return 0;
}

# Generates Table of Emissions (WITH pseudo counts) formatted for Model File
sub em_table_pseudo_ct{
	my ($em_counts, $em_rows, $order) = @_;
	my @alph = qw(A C G T);
	my $final_em = "";
	
	if ($order == 0) {
		my $em_row = "";
		foreach my $col (@alph) {
			my $count = $em_counts->{$col};
			if (defined $count) {$count++; $em_row .= "$count\t";}
			else                {$em_row .= "1\t";}
		}
		chop($em_row);
		$final_em .= "$em_row\n";
	} else {	
		foreach my $row (sort keys %$em_rows) {
			my $em_row = "";
			if (!defined $em_counts->{$row}) {
				$em_row = "1\t1\t1\t1";
			} else {
				foreach my $col (@alph) {
					my $count = $em_counts->{$row}->{$col};
					if (defined $count) {$count++; $em_row .= "$count\t";}
					else                {$em_row .= "1\t";}
				}
				chop($em_row);
			}
			$final_em .= "$em_row\n";
		}
	}
	return $final_em;
}

# Generates Table of Emissions (NO pseudo counts) formatted for Model File
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
	
	if ($order == 6) {
		foreach my $c1 (@alph) {
			foreach my $c2 (@alph) {
				foreach my $c3 (@alph) {
					foreach my $c4 (@alph) {
						foreach my $c5 (@alph) {
							foreach my $c6 (@alph) {
								$em_rows{"$c1$c2$c3$c4$c5$c6"} = 1;
							}
						}
					}
				}
			}
		}	
	}
	
	if ($order == 7) {
		foreach my $c1 (@alph) {
			foreach my $c2 (@alph) {
				foreach my $c3 (@alph) {
					foreach my $c4 (@alph) {
						foreach my $c5 (@alph) {
							foreach my $c6 (@alph) {
								foreach my $c7 (@alph) {
									$em_rows{"$c1$c2$c3$c4$c5$c6$c7"} = 1;
								}
							}
						}
					}
				}
			}
		}	
	}
	
	if ($order == 8) {
		foreach my $c1 (@alph) {
			foreach my $c2 (@alph) {
				foreach my $c3 (@alph) {
					foreach my $c4 (@alph) {
						foreach my $c5 (@alph) {
							foreach my $c6 (@alph) {
								foreach my $c7 (@alph) {
									foreach my $c8 (@alph) {
										$em_rows{"$c1$c2$c3$c4$c5$c6$c7$c8"} = 1;
									}
								}
							}
						}
					}
				}
			}
		}	
	}
	return %em_rows;
}

