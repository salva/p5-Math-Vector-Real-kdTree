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
        [$axis, _build($v, $l), _build($v, $r), $median, $b->[$axis], $t->[$axis], scalar(@$l), scalar(@$r)];
    }
    else {
        [undef, @$ix];
    }
}

sub size { scalar @{shift->{vs}} }

sub at {
    my ($self, $ix) = @_;
    Math::Vector::Real::clone($self->{vs}[$ix]);
}

sub insert {
    my ($self, $v) = @_;
    my $vs = $self->{vs};
    $v = Math::Vector::Real::clone($v);
    push @$vs, $v;
    _insert($vs, $self->{tree}, $#$vs)
}

sub _insert {
    my ($vs, $t, $ix) = @_;
    if (defined $t->[0]) {
        my ($axis, $l, $r, $median, $min, $max, $nl, $nr) = @$t;
        my $c = $vs->[$ix][$axis];
        my $pole;
        if ($c < $median) {
            if (2 * $nr + $max_per_pole >= $nl) {
                $t->[6]++;
                $t->[4] = $c if $c < $min;
                return _insert($vs, $l, $ix);
            }
        }
        else {
            if (2 * $nl + $max_per_pole >= $nr) {
                $t->[7]++;
                $t->[5] = $c if $c > $max;
                return _insert($vs, $r, $ix);
            }
        }
        my @store;
        $#store = $nl + $nr;
        @store = ($ix);
        _push_all($t, \@store);
        $_[1] = _build($vs, \@store);
    }
    elsif ($#$t< $max_per_pole) {
        push @$t, $ix;
    }
    else {
        $_[1] = _build($vs, [$ix, @$t[1..$#$t]])
    }
}

sub move {
    my ($self, $ix, $v) = @_;
    my $vs = $self->{vs};
    ($ix >= 0 and $ix < @$vs) or croak "index out of range";
    my $t = $self->{tree};
    _delete($vs, $t, $ix);
    $vs->[$ix] = Math::Vector::Real::clone($v);
    _insert($vs, $t, $ix);
}


sub _delete {
    my ($vs, $t, $ix) = @_;
    if (defined $t->[0]) {
        my ($axis, $l, $r, $median) = @_;
        my $c = $vs->[$ix][$axis];
        if ($c <= $median and _delete($vs, $l, $ix)) {
            $t->[6]--;
            return 1;
        }
        elsif ($c >= $median and _delete($vs, $r, $ix)) {
            $t->[7]--;
            return 1;
        }
    }
    else {
        my $l = scalar @$t;
        @$t = grep { not defined $_ or $_ =! $ix } @$t;
        return @$t < $l;
    }
}



sub _push_all {
    my ($t, $store) = @_;
    if (defined $t->[0]) {
        _push_all($t->[1], $store);
        _push_all($t->[2], $store);
    }
    else {
        push @$store, @$t[1..$#$t]
    }
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
    my ($self, $v, $d, $but) = @_;
    my $vs = $self->{vs};
    $but //= -1;
    my $start = 0;
    if ($but >= 0) {
        $but >= @$vs and croak "index out of range";
        return if @$vs < 2;
        $start = 1 unless $but;
    }
    else {
        return if @$vs < 1;
    }
    my $d2 = (defined $d ? $d * $d : $vs->[$start]->dist2($v));
    my ($rix, $rd2) = _find_nearest_neighbor($self->{vs}, $self->{tree}, $v, $start, $d2, $but);
    $rix // return;
    wantarray ? ($rix, sqrt($rd2)) : $rix;
}

sub find_nearest_neighbor_internal {
    my ($self, $vix, $d) = @_;
    $vix >= 0 or croak "index out of range";
    $self->find_nearest_neighbor($self->{vs}[$vix], $d, $vix);
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

sub find_in_ball {
    my ($self, $z, $d, $but) = @_;
    $but //= -1;
    _find_in_ball($self->{vs}, $self->{tree}, $z, $d * $d, $but);
}

sub _find_in_ball {
    my ($vs, $t, $z, $d2, $but) = @_;
    if (defined $t->[0]) {
        my ($axis, $l, $r, $median) = @$t;
        my $c = $z->[$axis];
        my $dc = $c - $median;
        my ($f, $s) = (($dc < 0) ? ($l, $r) : ($r, $l));
        if ($dc * $dc <= $d2) {
            if (wantarray) {
                return (_find_in_ball($vs, $f, $z, $d2, $but),
                        _find_in_ball($vs, $s, $z, $d2, $but))
            }
            else {
                return (_find_in_ball($vs, $f, $z, $d2, $but) +
                        _find_in_ball($vs, $s, $z, $d2, $but));
            }
        }
        else {
            return _find_in_ball($vs, $f, $z, $d2, $but);
        }
    }
    else {
        grep { $_ != $but and $vs->[$_]->dist2($z) <= $d2 } @$t[1..$#$t]
    }
}

1;
__END__

=head1 NAME

Math::Vector::Real::kdTree - kd-Tree implementation on top of Math::Vector::Real

=head1 SYNOPSIS

  use Math::Vector::Real::kdTree;

  use Math::Vector::Real;
  use Math::Vector::Real::Random;

  my @v = map Math::Vector::Real->random_normal(4), 1..1000;

  my $tree = Math::Vector::Real::kdTree->new(@v);

  my $ix = $tree->find_nearest_neighbor(V(0, 0, 0, 0));

  say "nearest neighbor is $ix, $v[$ix]";

=head1 DESCRIPTION

This module implements a kd-Tree data structure in Perl and some
related algorithms.

The following methods are provided:

=over 4

=item $t = Math::Vector::Real::kdTree->new(@points)

Creates a new kdTree containing the gived points.

=item $t->insert($p)

Inserts the given point into the kdTree.

=item $s = $t->size

Returns the number of points inside the tree.

=item $p = $t->at($ix)

Returns the point at the given index inside the tree.

=item $t->move($ix, $p)

Moves the point at index C<$ix> to the new given position readjusting
the tree structure accordingly.

=item ($ix, $d) = $t->find_nearest_neighbor($p, $max_d, $but_ix)

Find the nearest neighbor for the given point C<$p> and returns its
index and the distance between the two points (in scalar context the
index is returned).

If C<$max_d> is defined, the search is limited to the points within that distance

If C<$but_ix> is defined, the point with the given index is not considered.

=item @ix = $t->find_nearest_neighbor_all_internal

Returns the index of the nearest neighbor for every point inside the tree.

It is equivalent to (though, internally, it uses a better algorithm):

  @ix = map {
            scalar $t->nearest_neighbor($t->at($_), undef, $_)
        } 0..($t->size - 1);

=item @ix = $t->find_in_ball($z, $d, $but)

=item $n = $t->find_in_ball($z, $d, $but)

Finds the points inside the tree contained in the hypersphere with
center C<$z> and radius C<$d>.

In scalar context returns the number of points found. In list context
returns the indexes of the points.

if the extra argument C<$but> provided. The point with that index is
ignored.

=back

=head1 SEE ALSO

L<http://en.wikipedia.org/wiki/K-d_tree>

L<Math::Vector::Real>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011 by Salvador Fandiño E<lt>sfandino@yahoo.comE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.3 or,
at your option, any later version of Perl 5 you may have available.


=cut
