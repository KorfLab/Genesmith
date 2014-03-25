#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use FAlite;
use DataBrowser;
use Getopt::Std;
use vars qw($opt_h $opt_t $opt_o);
getopts('h:t:o:');


#==========================#
# Genesmith Wrapper Script #
#==========================#
my $OPTS  = "none";
my $TRANS = "gencode.txt";

die "
usage: run_genesmith.pl [options] <GFF> <FASTA> <PROTEIN>

options:
  -t <file>    translation table                      Default = $TRANS
  -o <string>  genesmith comma separated cmd options  Default = $OPTS
  -h           help (usage details)
" unless @ARGV == 3;
$TRANS = $opt_t            if $opt_t;
$OPTS  = edit_opts($opt_o) if $opt_o;
my ($GFF, $FASTA, $PROTEIN) = @ARGV;
my $start_run = time();

#----------------#
# Pre-processing #
#----------------#
# Get Species two letter abbreviation
my ($TAXA) = $GFF =~ /(\w\.\w)\w+\.gff/;
$TAXA =~ s/\.//;

print "\nGenesmith Performance Evaluation ($TAXA\)\n";
print "-------------------------------------\n\n";

# Format Translation Table
# print ">>> Formatting Tranlation Table\n\n";
# `format_trans_tbl.pl $TRANS`;

# Create Hash of Profile FHs
print ">>> Create Hash of Profile FHs\n\n";
my %profiles;
foreach my $fh (glob("~/scratch/KOGs_Profiles/*.hmm")) {
	my ($kog_id)      = $fh =~ /(\w+\d+).hmm/;
	my $id_ln         = `grep "$kog_id" $FASTA`;
	next              if $id_ln !~ /$kog_id/;
	my ($fa_id)       = $id_ln =~ />(\S+)\n/;
	$profiles{$fa_id} = $fh;
}

# Create Hash of KOG Protein Sequences
print ">>> Create Hash of KOG Protein Sequences\n\n";
my %proteins;
open(IN, "<$PROTEIN") or die "Error reading $PROTEIN\n";
my $prot_fasta = new FAlite(\*IN);
while (my $entry = $prot_fasta->nextEntry) {
	my ($fa_id)       = $entry->def =~ /^>(\S+)$/;
	my $prot_seq      = $entry->seq;
	$proteins{$fa_id} = $prot_seq;
}

#----------------------------#
# Get Test and Training Sets #
#----------------------------#
my $sets = 1;  # Stores Total number of Sets
print ">>> Creating TEST and TRAINING sets\n";
my $opt_set_size = "-s all";
`test_train_sets.pl $opt_set_size $GFF $FASTA`;

my %files;
my @gff_fhs;
my @fa_fhs;
foreach my $fh (glob("./$TAXA\_*")) {
	if    ($fh =~ /.fa$/)  {push(@fa_fhs, $fh);}
	elsif ($fh =~ /.gff$/) {push(@gff_fhs, $fh);}
}

foreach my $gff (@gff_fhs) {
	foreach my $fa (@fa_fhs) {
		my ($fh) = $gff =~ /.(\/\w+).gff/;
		if ($fa =~ /$fh\.fa/) {
			my ($set) = $fa =~ /\w+(\d+).fa$/;
			$sets = $set + 1 if $set >= $sets;
			$files{$gff} = $fa;
		}
	}
}
print "\t$sets sets\n\n";

#-------------#
# Create HMMs #
#-------------#
print ">>> Generating HMMs\n";
for (my $i=0; $i < $sets; $i++) {
	foreach my $fh (keys %files) {
		if ($fh =~ /train$i\.gff/) {
			my $cmd = "hmmgen.pl ";
			$cmd .= $OPTS if $OPTS ne "none";
			$cmd .= " $fh $files{$fh}";
			`$cmd`;
# 			print "\t$cmd\n";
		}
	}
}

#-------------------#
# Running Genesmith #
#-------------------#
print "\n>>> Running Genesmith\n";
foreach my $fh (glob("$TAXA\_*.hmm")) {
	my ($set, $st_quant) = $fh =~ /$TAXA\_\w+(\d+)_(\d+).hmm/;
# 	print "\tSET: $set\tHMM: $st_quant\ states\n";
	foreach my $gff (keys %files) {
		if ($gff =~ /test/) {
			my ($gffset) = $gff =~ /\w+(\d+).gff/;
			if ($gffset eq $set) {
				my $temp_out = "$TAXA\_test$set\_$st_quant\_pred.gff";
				open (OUT, ">$temp_out") or die "Error writing into OUT\n";
				close OUT;
				
				# Create temp FASTA file with one ID/Seq
				open (IN, "<$files{$gff}") or die "Error reading FASTA\n";
				my $fasta = new FAlite(\*IN);
				while (my $entry  = $fasta->nextEntry) {
					my ($fa_id)   = $entry->def =~ /^>(\S+)$/;
					my $pro_hmm = $profiles{$fa_id};
					my $aa_seq  = $proteins{$fa_id};
					open (ONE, ">one_id.fa") or die "Error writing into OUT\n";
					print ONE $entry;
					close ONE;
					
					my $cmd = "genesmith $fh ./one_id.fa $pro_hmm $aa_seq > one_pred.txt";
					`$cmd`;
					`cat one_pred.txt >> $temp_out`;
				}
				close IN;
			}
		}
	}
}

#---------------------------#
# Evaluate Gene Predictions #
#---------------------------#
print "\n>>> Prediction Evaluation\n\n";
my $exp_gff;
my $pred_gff;
my $exp_fasta;
my $exp_set;
my $ST_QUANT;      # Total number of states

my $all_exp_fasta = "$TAXA\_all_exp.fa";
my $all_exp_gff   = "$TAXA\_all_exp.gff";
my $all_pred_gff  = "$TAXA\_all_pred.gff";
`touch $all_exp_fasta $all_exp_gff $all_pred_gff`;

foreach my $fh (glob("*.gff")) {
	if ($fh =~ /$TAXA\_test\d+/) {
		if ($fh =~ /\w+\d+.gff/) {
			($exp_set) = $fh =~ /\w+(\d+).gff/;
			$exp_gff   = "./$fh";
			$exp_fasta = $files{$exp_gff};
		}
		if ($fh =~ /pred.gff/) {
			$pred_gff = "./$fh";
			my ($pred_set, $st_quant) = $pred_gff =~ /\w+(\d+)_(\d+)_pred.gff/;
			$ST_QUANT = $st_quant;
			if ($pred_set eq $exp_set) {
# 				my $ln_ct = `grep -c ">" $exp_fasta`;
# 				print $exp_fasta, "\t", $ln_ct, "\n";
				`cat $exp_fasta >> $all_exp_fasta`;
				`cat $exp_gff >> $all_exp_gff`;
				`cat $pred_gff >> $all_pred_gff`;
				`rm $exp_fasta $exp_gff $pred_gff`;
			}
		}
	}
}

# my $ln_ct = `grep -c ">" $all_exp_fasta`;
# print $all_exp_fasta, "\t", $ln_ct, "\n";

my $cmd = "evaluator.pl $all_exp_fasta $all_exp_gff $all_pred_gff";
my $cmd_output = `$cmd`;
chomp($cmd_output);
my @eval_stats  = split("\n", $cmd_output);

my $nuc_counts = $eval_stats[0];
my $nuc_stats  = $eval_stats[1];
my $cds_counts = $eval_stats[2];
my $kog_counts = $eval_stats[3];

print $TAXA,        "\t",
      $nuc_stats,   "\t",
      $cds_counts,  "\t",
      $kog_counts,  "\t",
      $ST_QUANT,    "\t\n";


#--------------------#
# Remove Extra Files #
#--------------------#
print "\n\n>>> Removing Extra Files\n";
foreach my $fh (glob("$TAXA\_*"))       {`rm $fh`;}
foreach my $fh (glob("one_*"))          {`rm $fh`;}

my $end_run = time();
my $run_time = $end_run - $start_run;
my $minutes  = int($run_time / 60);
my $seconds  = $run_time % 60;
print "\n>>> COMPLETE!\tTime: $minutes min  $seconds sec\n";


#=============#
# SUBROUTINES #
#=============#

# Generates state order option command implemented while creating HMMs
sub edit_opts{
	my ($cmd_opts) = @_;
	my $cmd = "";
	my @cmd_args = split(",", $cmd_opts);
	my @opts = qw(-f -s -C -D -i -A -e -t);
	die "OPT -o requires 8 comma separated elements\n" unless scalar(@cmd_args) == 8;
	
	for (my $i=0; $i < @cmd_args; $i++) {
		my $feat = $opts[$i];
		if ($feat =~ /-C|-D|-A/) {
			if ($cmd_args[$i] =~ /:/) {
				$cmd .= "$feat $cmd_args[$i] " if $feat =~ /-C/;
				$cmd .= "$feat $cmd_args[$i] " if $feat =~ /-D/;
				$cmd .= "$feat $cmd_args[$i] " if $feat =~ /-A/;
			} else {
				$cmd .= "$feat 3:$cmd_args[$i] " if $feat =~ /-C/;   #CDS states    = 3
				$cmd .= "$feat 2:$cmd_args[$i] " if $feat =~ /-D/;   #Donor states  = 2
				$cmd .= "$feat 2:$cmd_args[$i] " if $feat =~ /-A/;   #Accept states = 2
			}
		} else {
			$cmd .= "$feat $cmd_args[$i] ";
		}
	}
	return $cmd;
}
