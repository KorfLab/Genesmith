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
usage: run_genesmith.pl [options] <GFF> <FASTA>

options:
  -t <file>    translation table        Default = $TRANS
  -o <string>  genesmith order options  Default = $OPTS
  -h           help (usage details)
" unless @ARGV == 2;
$TRANS = $opt_t            if $opt_t;
$OPTS  = edit_opts($opt_o) if $opt_o;
my ($GFF, $FASTA) = @ARGV;
my $start_run = time();


#----------------#
# Pre-processing #
#----------------#
# Get Species two letter abbreviation
my $gff_info = `head -n 1 $GFF`;
my ($TAXA) = $gff_info =~ /^(\w\w)\w+/;
print "\nGenesmith Performance Evaluation ($TAXA\)\n";
print "-------------------------------------\n\n";

# Format Translation Table
print ">>> Formatting Tranlation Table\n\n";
`format_trans_tbl.pl $TRANS`;

# Create Hash of Profile FHs
print ">>> Create Hash of Profile FHs\n\n";
my %profiles;
foreach my $fh (glob("./hmmer3_profiles/*.hmm")) {
	my ($kog_id)      = $fh =~ /(\w+\d+).hmm/;
	my $id_ln         = `grep "$kog_id" $FASTA`;
	next              if $id_ln !~ /$kog_id/;
	my ($fa_id)       = $id_ln =~ />(\S+)\n/;
	$profiles{$fa_id} = $fh;
}


#----------------------------#
# Get Test and Training Sets #
#----------------------------#
my $sets = 1;
print ">>> Creating TEST and TRAINING sets\n";
`test_train_sets.pl $GFF $FASTA`;

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
print "\t4 sets\n\n";

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
			print "\t$cmd\n";
		}
	}
}

#-------------------#
# Running Genesmith #
#-------------------#
print "\n>>> Running Genesmith\n";
foreach my $fh (glob("$TAXA\_*.hmm")) {
	my ($set) = $fh =~ /$TAXA\_\w+(\d+)_\d+.hmm/;
	foreach my $gff (keys %files) {
		if ($gff =~ /test/) {
			my ($gffset) = $gff =~ /\w+(\d+).gff/;
			if ($gffset eq $set) {
				my $temp_out = "$TAXA\_test$set\_pred.gff";
				open (OUT, ">$temp_out") or die "Error writing into OUT\n";
				close OUT;
				
				# Create temp FASTA file with one ID/Seq
				open (IN, "<$files{$gff}") or die "Error reading FASTA\n";
				my $fasta = new FAlite(\*IN);
				while (my $entry  = $fasta->nextEntry) {
					my ($fa_id)   = $entry->def =~ /^>(\S+)$/;
					my ($pro_hmm) = $profiles{$fa_id};
					open (ONE, ">one_id.fa") or die "Error writing into OUT\n";
					print ONE $entry;
					close ONE;
					
					my $cmd = "genesmith $fh ./one_id.fa $pro_hmm > one_pred.txt";
					`$cmd`;
					`cat one_pred.txt >> $temp_out`;
					print $set, "\t", $cmd, "\n";
				}
				close IN;
			}
		}
	}
}

#---------------------------#
# Evaluate Gene Predictions #
#---------------------------#
print "\n>>> Prediction Evaluation\n";
my $exp_gff;
my $pred_gff;
my $exp_fasta;
my $exp_set;
foreach my $fh (glob("*.gff")) {
	if ($fh =~ /$TAXA\_test\d+/) {
		if ($fh =~ /\w+\d+.gff/) {
			($exp_set) = $fh =~ /\w+(\d+).gff/;
			$exp_gff   = "./$fh";
			$exp_fasta = $files{$exp_gff};
		}
		if ($fh =~ /pred.gff/) {
			$pred_gff = "./$fh";
			my ($pred_set) = $pred_gff =~ /\w+(\d+)_pred.gff/;
			if ($pred_set eq $exp_set) {
				my $cmd = "evaluator.pl $exp_fasta $exp_gff $pred_gff";
				my $results = `$cmd`;
				print "--|SET: $pred_set\|-----------\n";
				print $results, "\n";
			}
		}
	}
}

#--------------------#
# Remove Extra Files #
#--------------------#
print "\n>>> Removing Extra Files\n";
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
	die "OPT -o requires 8 number string\n" unless length($cmd_opts) == 8;
	my $cmd = "";
	my @orders = split("", $cmd_opts);
	my @opts = qw(-f -s -C -D -i -A -e -t);
	
	for (my $i=0; $i < @orders; $i++) {
		my $feat = $opts[$i];
		if ($feat =~ /-C|-D|-A/) {
			$cmd .= "$feat 3:$orders[$i] " if $feat =~ /-C/;   #CDS states    = 3
			$cmd .= "$feat 2:$orders[$i] " if $feat =~ /-D/;   #Donor states  = 2
			$cmd .= "$feat 2:$orders[$i] " if $feat =~ /-A/;   #Accept states = 2
		} else {
			$cmd .= "$feat $orders[$i] ";
		}
	}
	return $cmd;
}
