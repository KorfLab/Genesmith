package STATE;
use strict;
use warnings 'FATAL' => 'all';

# State Object
sub new{
	my ($class, $name, $label, $order, $st_length) = @_;
	my $t_matrix = {};
	my $emission = {};
	
	my $self = {
		name     => $name,
		label    => $label,
		order    => $order,
		st_length => $st_length,
		t_matrix => $t_matrix,
		emission => $emission,
	};
	bless $self, $class;
	return $self;
}

# Redefine or Get Object Parameter
sub name{
	my ($self, $name) = @_;
	return $self->{name} unless defined $name;
	$self->{name} = $name;
}

sub label{
	my ($self, $label) = @_;
	return $self->{label} unless defined $label;
	$self->{label} = $label;
}

sub order{
	my ($self, $order) = @_;
	return $self->{order} unless defined $order;
	$self->{order} = $order;
}

# Now this is substr length for parsing rather than quantity of states
sub st_length{
	my ($self, $st_length) = @_;
	return $self->{st_length} unless defined $st_length;
	$self->{st_length} = $st_length;
}

sub t_matrix{
	my ($self, $st, $t_prob) = @_;
	return $self->{t_matrix} unless defined($st and $t_prob);
	$self->{t_matrix}->{$st} = $t_prob;
}

sub emission{
	my ($self, $st, $order, $seq) = @_;
	return $self->{emission} unless defined($st and $order and $seq);
	
	if ($st =~ /start|stop|don|accep|branch/) {
		my $pos;
		($pos) = $st =~ /start(\d+)/  if $st =~ /start\d+/;
		($pos) = $st =~ /stop(\d+)/   if $st =~ /stop\d+/;
		($pos) = $st =~ /don(\d+)/    if $st =~ /don\d+/;
		($pos) = $st =~ /accep(\d+)/  if $st =~ /accep\d+/;
		($pos) = $st =~ /branch(\d+)/ if $st =~ /branch\d+/;
		my $ctx = uc substr($seq, $pos, $order);
		my $nt = uc substr($seq, ($order + $pos), 1);
		if ($ctx =~ /\S+/) {
			$self->{emission}->{$ctx}->{$nt}++;
		} else {
			$self->{emission}->{$nt}++;
		}
	} elsif ($st =~ /cds\d/) {
		my ($phase) = $st =~ /cds(\d)/;
		my $v = 3;
		$v    = $order if $order > 3;
		for (my $i = $phase+$v; $i < length($seq); $i+=3) {
			if ($i >= $order) {
				my $ctx = uc substr($seq, $i - $order, $order);
				my $nt = uc substr($seq, $i, 1);
				next unless $nt =~ /[ACGT]/;
				if ($ctx =~ /\S+/) {
					$self->{emission}->{$ctx}->{$nt}++;
				} else {
					$self->{emission}->{$nt}++;
				}
			}
		}
	} else {
		for (my $i = $order; $i < length($seq); $i++) {
			my $ctx = uc substr($seq, $i - $order, $order);
			my $nt = uc substr($seq, $i, 1);
			next unless $nt =~ /[ACGT]/;
			if ($ctx =~ /\S+/) {
				$self->{emission}->{$ctx}->{$nt}++;
			} else {
				$self->{emission}->{$nt}++;
			}
		}
	}
}

1;