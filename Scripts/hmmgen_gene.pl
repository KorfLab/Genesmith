#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use STATE;
use HMMstar;
use DataBrowser;
use Getopt::Std;
use vars qw($opt_h $opt_5 $opt_m $opt_c $opt_d $opt_i $opt_a $opt_s $opt_3 $opt_D $opt_A $opt_U $opt_E);
getopts('h5:m:c:d:i:a:s:3:D:A:U:E:');


# DEFAULT settings [options]
my $UP      = 0;     # Order Upstream
my $START   = 0;     # Order Start, includes canonical ATG at the end, length = 3 + order
my $EXON    = 0;     # Order CDS, always 3 states
my $DON     = 0;     # Order Donor, starts with canonical GT, length = 2 + order
my $INTRON  = 0;     # Order Intron, always 3 states
my $ACCEP   = 0;     # Order Acceptor, ends with canonical AG, length = 2 + order
my $STOP    = 0;     # Order Stop, starts with stop codon, length = 3 + order
my $DOWN    = 0;     # Order Downstream
my $L_DON   = 5;     # Quantity of Donor States
my $L_ACCEP = 10;     # Quantity of Acceptor States 
my $L_UP    = 500;   # Length of Upstream region parsed for training
my $L_DOWN  = 500;   # Length of Downstream region parsed for training


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
  -D <length>  Donor Site Length                          Default = $L_DON
  -A <length>  Acceptor Site Length                       Default = $L_ACCEP
  -U <length>  upstream training                          Default = $L_UP
  -E <length>  downtream training                         Default = $L_DOWN
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
$L_DON   = $opt_D if $opt_D;
$L_ACCEP = $opt_A if $opt_A;
$L_UP    = $opt_U if $opt_U;
$L_DOWN  = $opt_E if $opt_E;

my ($GFF, $FASTA) = @ARGV;

# Sanity Check - options
if ($START   !~ /^\d+$/ or
    $DON     !~ /^\d+$/ or 
    $ACCEP   !~ /^\d+$/ or
    $STOP    !~ /^\d+$/ or
    $UP      !~ /^\d+$/ or
    $EXON    !~ /^\d+$/ or
    $INTRON  !~ /^\d+$/ or
    $DOWN    !~ /^\d+$/ or
    $L_DON   !~ /^\d+$/ or
    $L_ACCEP !~ /^\d+$/ or    
    $L_UP    !~ /^\d+$/ or
    $L_DOWN  !~ /^\d+$/   ) {
    die "Invalid input [options]\n";
}

# Info for each Group of States
my @st_order  = ($UP,  $START,  $EXON, $EXON, $EXON,     $DON,  $INTRON, $ACCEP,  $STOP,  $DOWN);
my @st_name   = ('GU', 'start', 'cds', 'CDS', 'CDSfull', 'don', 'i',     'accep', 'stop', 'GD');
my @st_label  = ('U',  'C',     'C',   'C',   'C',       'I',   'I',     'I',     'C',    'D');
my @st_quant  = ( 1,    3,       3,     1,     1,        $L_DON, 1,      $L_ACCEP, 3,      1);


#---------------------------------------------------------#
# Create Array of State Objects (or a multi-layered hash) #
#---------------------------------------------------------#
my @states;   # array of state objects

for (my $i=0; $i < @st_quant; $i++) {
	my $order = $st_order[$i];
	my $quant = $st_quant[$i];
	for (my $c=0; $c < $quant; $c++) {
		my $name = $st_name[$i] . "$c";
		my $st_length;
		if  ($name =~ /start|stop|don|accep/) {$st_length = $quant + $order;}
		else                                  {$st_length = $quant;}
		my $st_obj = new STATE($name, $st_label[$i], $order, $st_length);
		push(@states, $st_obj);
	}
}


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
		
		# Exons & CDS
		foreach my $exon ($cds->exons) {
			$gene_struc{$id}{'CDS'} .= uc($exon->sequence);
			push @{$gene_struc{$id}{'Exons'}}, uc($exon->sequence);
		}
		
		# Introns
		foreach my $intron ($cds->introns) {
			push @{$gene_struc{$id}{'Introns'}}, uc($intron->sequence);
		}
		
		# Downstream
		my $dn_beg = $cds->stop_site->end + 2;
		my $dn_len = length($contig->dna->sequence) - $dn_beg +1;
		$dn_len = $L_DOWN if $dn_len > $L_DOWN;
		$gene_struc{$id}{'Downstream'} = uc(substr($contig->dna->sequence, $dn_beg, $dn_len));
	}
}

# Sanity Check
#browse(\%gene_struc);

#--------------------------------#
# Generate State Emission Counts #
#--------------------------------#
my %warnings;   # contains IDs of genes with insufficient lengths

foreach my $id (keys %gene_struc) {
	my $CDS = substr($gene_struc{$id}{'CDS'}, 0, length($gene_struc{$id}{'CDS'}) -3);
	my $upstream   = $gene_struc{$id}{'Upstream'};       # Upstream sequence
	my $downstream = $gene_struc{$id}{'Downstream'};     # Downstream sequence
	
	for (my $s=0; $s < @states; $s++) {
		my $st        = $states[$s]->name;
		my $order     = $states[$s]->order;
		my $st_length = $states[$s]->st_length;
		
		# Upstream
		if ($st =~ /^GU/) {
			$states[$s]->emission($st, $order, $upstream);
		}
		
		# Start
		if ($st =~ /^start\d+$/) {
			my $seq = substr($upstream, -$order, $order) . substr($CDS, 0, 3);
			$states[$s]->emission($st, $order, $seq);
		}
		
		# 3 CDS states [Standard Model]
		if ($st =~ /^cds/) {
			my $cds_str = "";
			$cds_str    = substr($upstream, -($order -3), $order -3) if $order > 3;
			$cds_str   .= $CDS; 
			$states[$s]->emission($st, $order, $cds_str);
		}
		
		# 1 CDS state [Basic Model]
		if ($st =~ /^CDS\d+/) {
			my $seq = $gene_struc{$id}{'CDS'};
			$states[$s]->emission($st, $order, substr($seq, 3, length($seq) -6));
		}
		
		# 1 CDS [Basic Model with no Start/Stop states]
		if ($st =~ /^CDSfull\d+/) {
			my $seq = $gene_struc{$id}{'CDS'};
			$states[$s]->emission($st, $order, $seq);
		}
		
		# Donor
		if ($st =~ /^don/) {
			my $count = 0;
			foreach my $intron (@{$gene_struc{$id}{'Introns'}}) {
				my $exon = $gene_struc{$id}{'Exons'}[$count];
				my $seq  = substr($exon, -$order, $order) . substr($intron, 0, $L_DON);
				$states[$s]->emission($st, $order, $seq);
				$count++;
			}
		}
		
		# Acceptor
		if ($st =~ /^accep/) {
			foreach my $intron (@{$gene_struc{$id}{'Introns'}}) {
				my $seq = substr($intron, -$st_length);
				$states[$s]->emission($st, $order, $seq);
			}
		}
		
		# Intron body
		if ($st =~ /^i/) {
			my $count = 0;
			foreach my $intron (@{$gene_struc{$id}{'Introns'}}) {
				my $exon = $gene_struc{$id}{'Exons'}[$count];
				my $previous_seq = $exon . substr($intron, 0, $L_DON);
				my $seq = substr($previous_seq, -$order, $order);  
				$seq   .= substr($intron, $L_DON, length($intron) - ($L_DON + $L_ACCEP));
				$states[$s]->emission($st, $order, $seq);
				$count++;
			}
		}
		
		# Stop
		if ($st =~ /^stop/) {
			my $seq = substr($gene_struc{$id}{'CDS'}, -$st_length);
			$states[$s]->emission($st, $order, $seq);
		}
		
		# Downstream
		if ($st =~ /^GD/) {
			my $seq = substr($gene_struc{$id}{'CDS'}, -$order, $order) . $downstream;
			$states[$s]->emission($st, $order, $seq);
		}
	}
}

# Sanity Check
# browse(\%warnings);  # lists of Genes Filtered out due to insufficient length
# foreach my $obj (@states) {
# 	print ">", $obj->name, "\n";
# 	browse($obj->emission);
# }

#-------------------------------------------------------------------------------------------#
# Create State Emission Output Files                                                        #
# * Output formats:                                                                         #
#   - Standard                                                                              #
#         [state name]-[Model Type]-[State Label]-[State Order].txt                         #
#   - If state is affected by modified state duration or quantity during training           #
#         [state name]-[Model Type]-[State length/quantity]-[State Label]-[State Order].txt #
#            * State length/quantity params:                                                #
#                - Introns    = "$L_DON_$L_ACCEP"                                           #
#                - Donors     = "$L_DON"                                                    #
#                - Acceptors  = "$L_ACCEP"                                                  #
#                - Upstream   = "$L_UP"                                                     #
#                - Downstream = "$L_DOWN"                                                   #
#-------------------------------------------------------------------------------------------#
my $I_QUANT = 3;   # Number of copies generated for each Intron state

#Create Output Directory
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
	
	# Copy 3 sets of Exon & Intron states
	if ($st =~ /don|i|accep/) {
		for (my $i=0; $i < $I_QUANT; $i++) {
			if ($st =~ /i/) {
				my ($st_name) = $st =~ /(\w+)\d+$/;
				$st_name     .= "$i";
				$st_name     .= "-std-$L_DON\_$L_ACCEP\-$label\-$order\.txt";
				my $fh = $path . $st_name;
				output_st_em($fh, $em_output);
			} else {
				my $st_name = $st . "_$i";
				$st_name   .= "-std-$L_DON"       if $st =~ /don/;
				$st_name   .= "-std-$L_ACCEP"     if $st =~ /accep/;
				$st_name   .= "-$label\-$order\.txt";
				my $fh = $path . $st_name;
				output_st_em($fh, $em_output);
			}
		}
	} else {
		my $st_name;
		if    ($st =~ /^CDS\d+/)     {$st_name = "$st\-1-";}
		elsif ($st =~ /^CDSfull\d+/) {$st_name = "$st\-1S-";}
		else                         {$st_name = "$st\-std-";}
		$st_name    .= "$L_UP\-"             if $st =~ /GU/;
		$st_name    .= "$L_DOWN\-"           if $st =~ /GD/;
		$st_name    .= "$label\-$order\.txt";
		my $fh = $path . $st_name;
		output_st_em($fh, $em_output);
	}
}


#===================================================================#
# SUBROUTINES														#
#===================================================================#

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


