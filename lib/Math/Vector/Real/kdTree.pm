package Math::Vector::Real::kdTree;

our $VERSION = '0.01';

use 5.010;
use strict;
use warnings;
use Carp;

use Math::Vector::Real;
use Sort::Key::Top qw(nkeypartref ntop);

our $max_per_pole = 12;
our $recommended_per_pole = 6;

sub new {
    my $class = shift;
    my @v = map Math::Vector::Real::clone($_), @_;
    my @ix = (0..$#v);
    my $tree = _build(\@v, \@ix);
    my $self = { vs => \@v,
                 tree => $tree };
    bless $self, $class;
}

sub _build {
    my ($v, $ix) = @_;
    if (@$ix > $recommended_per_pole) {
        my ($b, $t) = Math::Vector::Real->box(@$v[@$ix]);
        my $axis = ($t - $b)->max_component_index;
        my $bstart = @$ix >> 1;
        my ($l, $r) = nkeypartref { $v->[$_][$axis] } $bstart => @$ix;
        my $lc = ntop -1 => map $v->[$_][$axis], @$l;
        my $rc = ntop  1 => map $v->[$_][$axis], @$r;
        my $median = 0.5 * ($lc + $rc);
        [$axis, _build($v, $l), _build($v, $r), $median, $b->[$axis], $t->[$axis]];
    }
    else {
        [undef, @$ix];
    }
}

sub at {
    my ($self, $ix) = @_;
    $self->{vs}[$ix]
}

sub path {
    my ($self, $vix) = @_;
    use Data::Dumper;
    # print Dumper $self->{tree};
    my $p = _path($self->{tree}, $vix);
    # print "path length: ", scalar(@$p), "\n";
    my $l = 1;
    $l = (($l << 1) | $_) for @$p;
    $l
}

sub _path {
    my ($t, $vix) = @_;
    if (defined $t->[0]) {
        for (0, 1) {
            my $p = _path($t->[1+$_], $vix);
            return [$_, @$p] if $p;
        }
        return undef;
    }
    #print "is $vix in @$t[1..$#$t]?\n";
    ((grep $_ == $vix, @$t[1..$#$t]) ? [] : ());
}

sub all {
    my ($self) = @_;
    sort { $a <=> $b } _all($self->{tree});
}

sub _all {
    my $t = shift;
    if (defined $t->[0]) {
        return (_all($t->[1]), _all($t->[2]))
    }
    return @$t[1..$#$t];
}

sub find {
    my ($self, $v) = @_;
    _find($self->{vs}, $self->{tree}, $v)
}

sub _find {
    my ($vs, $t, $v) = @_;
    while (1) {
        if (defined $t->[0]) {
            my ($axis, $l, $r, $median, $min, $max) = @$t;
            my $c = $v->[$axis];
            return if ($min > $c or $c > $max);
            if ($c < $median) {
                $t = $l;
            }
            else {
                if ($c == $median) {
                    my $ix = _find($vs, $l, $v);
                    return $ix if defined $ix;
                }
                $t = $r;
            }
        }
        else {
            for (@$t[1..$#$t]) {
                return $_ if $v == $vs->[$_];
            }
            return undef;
        }
    }
}

sub find_nearest_neighbor {
    my ($self, $v, $but) = @_;
    my $vs = $self->{vs};
    return unless @$vs;
    $but //= 1;
    _find_nearest_neighbor($self->{vs}, $self->{tree}, $v, 0, $vs->[0]->dist2($v), $but);
}

sub find_nearest_neighbor_internal {
    my ($self, $vix) = @_;
    my $vs = $self->{vs};
    $vix >= @$vs and croak "index out of range";
    return unless @$vs > 1;
    my $start = ($vix ? 0 : 1);
    my $v = $vs->[$vix];
    ( wantarray
      ? _find_nearest_neighbor($self->{vs}, $self->{tree}, $v, $start, $vs->[$start]->dist2($v), $vix)
      : (_find_nearest_neighbor($self->{vs}, $self->{tree}, $v, $start, $vs->[$start]->dist2($v), $vix))[0] )
}

sub _find_nearest_neighbor {
    my ($vs, $t, $v, $ix, $d2, $but) = @_;
    while (1) {
        if (defined $t->[0]) {
            my ($axis, $l, $r, $median) = @$t;
            my $c = $v->[$axis];
            my $cm = $c - $median;
            (my ($first), $t) = (($cm <= 0) ? ($l, $r) : ($r, $l));
            ($ix, $d2) = _find_nearest_neighbor($vs, $first, $v, $ix, $d2, $but);
            return ($ix, $d2) if $d2 <= $cm * $cm;
        }
        else {
            for (@$t[1..$#$t]) {
                my $p = $vs->[$_];
                my $d21 = $p->dist2($v);
                if ($d21 < $d2 and $_ != $but) {
                    $d2 = $d21;
                    $ix = $_;
                }
            }
            return ($ix, $d2)
        }
    }
}

sub find_nearest_neighbor_all_internal {
    my $self = shift;
    my $vs = $self->{vs};
    return unless @$vs > 1;
    my @best = ((0) x @$vs );
    my $first = $vs->[0];
    my @d2 = map $first->dist2($_), @$vs;
    $best[0] = 1;
    $d2[0] = $d2[1];
    _find_nearest_neighbor_all_internal($vs, $self->{tree}, \@best, \@d2);
    return @best;
}

sub _find_nearest_neighbor_all_internal {
    my ($vs, $t, $best, $d2) = @_;
    if (defined $t->[0]) {
        my ($axis, $l, $r, $median) = @$t;
        my @r;
        for my $side (0, 1) {
            my @poles = _find_nearest_neighbor_all_internal($vs, $t->[1 + $side], $best, $d2);
            push @r, @poles;
            my $other = $t->[2-$side];
            for my $pole (@poles) {
                for my $ix (@$pole[1..$#$pole]) {
                    my $v = $vs->[$ix];
                    my $md = $v->[$axis] - $median;
                    if ($d2->[$ix] > $md * $md) {
                        ($best->[$ix], $d2->[$ix]) =
                            _find_nearest_neighbor($vs, $other, $v, $best->[$ix], $d2->[$ix], -1);
                    }
                }
            }
        }
        return @r;
    }
    else {
        for my $i (2..$#$t) {
            my $ix = $t->[$i];
            my $iv = $vs->[$ix];
            for my $jx (@$t[1..$i-1]) {
                my $d21 = $iv->dist2($vs->[$jx]);
                if ($d21 < $d2->[$ix]) {
                    $d2->[$ix] = $d21;
                    $best->[$ix] = $jx;
                }
                if ($d21 < $d2->[$jx]) {
                    $d2->[$jx] = $d21;
                    $best->[$jx] = $ix;
                }
            }
        }
        return $t;
    }
}

1;
__END__

=head1 NAME

Math::Vector::Real::kdTree - kd-Tree implementation on top of Math::Vector::Real

=head1 SYNOPSIS

  use Math::Vector::Real::kdTree;



=head1 DESCRIPTION

This module implements a kd-Tree data structure in Perl.

=head2 EXPORT

None by default.

=head1 SEE ALSO

L<http://en.wikipedia.org/wiki/K-d_tree>

L<Math::Vector::Real>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011 by Salvador Fandiño E<lt>sfandino@yahoo.comE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.3 or,
at your option, any later version of Perl 5 you may have available.


=cut
