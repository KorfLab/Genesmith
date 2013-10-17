package STATE;
use strict;
use warnings 'FATAL' => 'all';

# State Object
sub new{
	my ($class, $name, $label, $order, $st_quant) = @_;
	my $t_matrix = {};
	my $emission = {};
	
	my $self = {
		name     => $name,
		label    => $label,
		order    => $order,
		st_quant => $st_quant,
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

sub st_quant{
	my ($self, $quant) = @_;
	return $self->{st_quant} unless defined $quant;
	$self->{st_quant} = $quant;
}

sub t_matrix{
	my ($self, $st, $t_prob) = @_;
	return $self->{t_matrix} unless defined($st and $t_prob);
	$self->{t_matrix}->{$st} = $t_prob;
}

sub emission{
	my ($self, $st, $order, $tot_trans, $seq) = @_;
	return $self->{emission} unless defined($st and $order and $tot_trans and $seq);
	
	if ($order == 0) {
		if ($tot_trans == 1) {
			my ($pos) = $st =~ /(\d+)$/;
			my $nt    = uc substr($seq, $pos, 1);
			$self->{emission}->{$nt}++;
		} else {
			for (my $i=0; $i < length($seq); $i++) {
				my $nt = uc substr($seq, $i, 1);
				$self->{emission}->{$nt}++;
			}
		}
	} else {
		if ($tot_trans == 1) {
			my $pos;
			($pos) = $st =~ /\w+(\d+)$/  if $st =~ /\w+\d+/;
			($pos) = $st =~ /\d+_(\d+)$/ if $st =~ /\w+\d+_\d+/;
			
			my $ctx = uc substr($seq, $pos, $order);
			my $nt = uc substr($seq, ($order + $pos), 1);
			$self->{emission}->{$ctx}->{$nt}++;
		} else {
			for (my $i = $order; $i < length($seq); $i++) {
				my $ctx = uc substr($seq, $i - $order, $order);
				my $nt = uc substr($seq, $i, 1);
				next unless $nt =~ /[ACGT]/;
				$self->{emission}->{$ctx}->{$nt}++;
			}
		}
	}
}

1;