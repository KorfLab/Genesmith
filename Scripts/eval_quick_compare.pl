#!/usr/bin/perl
# Performance Evaluator 
# Ian Korf 2014
use strict;
use warnings 'FATAL' => 'all';
use DataBrowser;
	
die "usage: $0 <gff1> <gff2>\n" unless @ARGV == 2;

# file1 = Expected GFF, file2 = Predicted GFF
my ($file1, $file2) = @ARGV;

my $gff1 = read_gff($file1);
my $gff2 = read_gff($file2);

# collect stats
my %gene   = init_counts();
my %nt     = init_counts();
my %exon   = init_counts();
my %intron = init_counts();
my %acc    = init_counts();
my %don    = init_counts();
my %start  = init_counts();
my %stop   = init_counts();


foreach my $sid (keys %$gff1) {
	# If no prediction was generated 0 counts are added to the sum (across all categories)
	if (!defined $gff2->{$sid}) {
		next;
	}
	
	# exon
	my (@exon1, @exon2);
	foreach my $exon (@{$gff1->{$sid}{exon}}) {push @exon1, stringify($exon)}
	foreach my $exon (@{$gff2->{$sid}{exon}}) {push @exon2, stringify($exon)}
	my $estat = compare(\@exon1, \@exon2);
	$exon{match} += $estat->{match};
	$exon{uni1}  += $estat->{uni1};
	$exon{uni2}  += $estat->{uni2};
	
	# gene
	if ($estat->{uni1} == 0 and $estat->{uni2} == 0) {$gene{match}++}
	else {$gene{uni1}++; $gene{uni2}++}

	# intron
	my (@intron1, @intron2);
	foreach my $intron (@{$gff1->{$sid}{intron}}) {push @intron1, stringify($intron)}
	foreach my $intron (@{$gff2->{$sid}{intron}}) {push @intron2, stringify($intron)}
	my $istat = compare(\@intron1, \@intron2);
	$intron{match} += $istat->{match};
	$intron{uni1}  += $istat->{uni1};
	$intron{uni2}  += $istat->{uni2};
	
	# don
	my (@don1, @don2);
	foreach my $don (@{$gff1->{$sid}{don}}) {push @don1, $don}
	foreach my $don (@{$gff2->{$sid}{don}}) {push @don2, $don}
	my $dstat = compare(\@don1, \@don2);
	$don{match} += $dstat->{match};
	$don{uni1}  += $dstat->{uni1};
	$don{uni2}  += $dstat->{uni2};

	# acc
	my (@acc1, @acc2);
	foreach my $acc (@{$gff1->{$sid}{acc}}) {push @acc1, $acc}
	foreach my $acc (@{$gff2->{$sid}{acc}}) {push @acc2, $acc}
	my $astat = compare(\@acc1, \@acc2);
	$acc{match} += $astat->{match};
	$acc{uni1}  += $astat->{uni1};
	$acc{uni2}  += $astat->{uni2};
		
	# start
	if ($gff1->{$sid}{start} == $gff2->{$sid}{start}) {$start{match}++}
	else {$start{uni1}++; $start{uni2}++}
	
	# stop
	if ($gff1->{$sid}{stop} == $gff2->{$sid}{stop}) {$stop{match}++}
	else {$stop{uni1}++; $stop{uni2}++}
	
	# nucleotide (a rough estimate)
	my $nt1 = nt_array($gff1->{$sid}{exon});
	my $nt2 = nt_array($gff2->{$sid}{exon});
	my $limit = @$nt1 > @$nt2 ? @$nt1 : @$nt2;
	my ($nmatch, $nuni1, $nuni2) = (0, 0, 0);
	for (my $i = 0; $i < $limit; $i++) {
		$nt1->[$i] = 0 if not defined $nt1->[$i];
		$nt2->[$i] = 0 if not defined $nt2->[$i];
		if ($nt1->[$i] == $nt2->[$i]) {$nmatch++} 
		else {$nuni1++; $nuni2++}
	}
	$nt{match} += $nmatch;
	$nt{uni1}  += $nuni1;
	$nt{uni2}  += $nuni2;
}

# report stats
report("Nucleotide", %nt);
report("Gene", %gene);
report("Exon", %exon);
report("Intron", %intron);
report("Start", %start);
report("Donor", %don);
report("Acceptor", %acc);
report("Stop", %stop);


#################
## Subroutines ##
#################
sub init_counts {
	return (match => 0, uni1 => 0, uni2 => 2);
}

sub nt_array {
	my ($exons) = @_;
	my @nt;
	foreach my $exon (@$exons) {
		for (my $i = $exon->{beg}; $i <= $exon->{end}; $i++) {
			$nt[$i-1] = 1;
		}
	}
	return \@nt;
}

sub report {
	my ($title, %data) = @_;
	printf "%-8s\tmatch:%d uni1:%d uni2:%d SN:%.3f SP:%.3f\n",
		$title, $data{match}, $data{uni1}, $data{uni2},
		$data{match} / ($data{match} + $data{uni1}),
		$data{match} / ($data{match} + $data{uni2});
}

sub stringify {
	my ($thing) = @_;
	return $thing->{beg} . $thing->{str} . $thing->{end};
}

sub compare {
	my ($a1, $a2) = @_;
	
	my %comp;
	foreach my $string (@$a1) {$comp{$string} += 1}
	foreach my $string (@$a2) {$comp{$string} += 2}
	
	my ($match, $uni1, $uni2) = (0, 0, 0);
	foreach my $item (keys %comp) {
		if    ($comp{$item} == 1) {$uni1++}
		elsif ($comp{$item} == 2) {$uni2++}
		else                      {$match++}
	}
	
	return {
		match => $match,
		uni1  => $uni1,
		uni2  => $uni2,
	};	
}

sub read_gff {
	my ($file) = @_;
	
	# read CDS features
	my %gff;
	open(my $fh, $file) or die;
	while (<$fh>) {
		chomp;
		my @f = split;
		next unless $f[2] eq 'CDS';
		push @{$gff{$f[0]}{exon}}, {
			beg => $f[3],
			end => $f[4],
			str => $f[6],
		};
	}
	
	# organize into gene structure features
	foreach my $sid (keys %gff) {
		# sort exons
		my @exon = sort {$a->{beg} <=> $b->{beg}} @{$gff{$sid}{exon}};
			
		# create introns
		my @intron;
		my @istring;
		for (my $i = 1; $i < @exon; $i++) {
			push @intron, {
				beg => $exon[$i-1]{end} +1,
				end => $exon[$i  ]{beg} -1,
				str => $exon[$i-1]{str},
			};
		}
			
		# create acceptors and donors
		my (@acc, @don);
		foreach my $intron (@intron) {
			if ($intron->{str} eq '+') {
				push @don, $intron->{beg};
				push @acc, $intron->{end};
			} else {
				push @don, $intron->{end};
				push @acc, $intron->{beg};
			}
		}
			
		# create start & stop
		my ($start, $stop);
		if ($exon[0]{str} eq '+') {
			$start = $exon[0 ]{beg};
			$stop  = $exon[-1]{end};
		} else {
			$start = $exon[-1]{end};
			$stop  = $exon[0 ]{beg};
		}
			
		# keep the structure
		$gff{$sid} = {
			exon => \@exon,
			intron => \@intron,
			acc => \@acc,
			don => \@don,
			start => $start,
			stop => $stop,
		}
	}
	
	return \%gff;
}

