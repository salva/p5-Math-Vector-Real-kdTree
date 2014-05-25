package Math::Vector::Real::kdTree;

our $VERSION = '0.09';

use 5.010;
use strict;
use warnings;
use Carp;

use Math::Vector::Real;
use Sort::Key::Top qw(nkeypartref nhead ntail);

our $max_per_pole = 12;
our $recommended_per_pole = 6;

use constant _n    => 0; # elements on subtree
use constant _c0   => 1; # corner 0
use constant _c1   => 2; # corner 1
use constant _sum  => 3; # centroid * n
use constant _s0   => 4; # subtree 0
use constant _s1   => 5; # subtree 1
use constant _axis => 6; # cut axis
use constant _cut  => 7; # cut point (mediam)

# on leaf nodes:
use constant _ixs   => 4
use constant _leaf_size => _ixs + 1;

sub new {
    my $class = shift;
    my @v = map Math::Vector::Real::clone($_), @_;
    my @ix = (0..$#v);
    my $tree = _build(\@v, \@ix);
    my $self = { vs   => \@v,
                 tree => $tree };
    bless $self, $class;
}

sub clone {
    my $self = shift;
    require Storable;
    my $clone = { vs   => [@{$self->{vs}}],
                  tree => Storable::dclone($self->{tree}) };
    $clone->{hidden} = { %{$self->{hidden}} } if $self->{hidden};
    bless $clone, ref $self;
}

sub _build {
    my ($v, $ix) = @_;
    if (@$ix > $recommended_per_pole) {
        my ($b, $t) = Math::Vector::Real->box(@$v[@$ix]);
        my $axis = ($t - $b)->max_component_index;
        my $bstart = @$ix >> 1;
        my ($p0, $p1) = nkeypartref { $v->[$_][$axis] } $bstart => @$ix;
        my $s0 = _build($v, $p0);
        my $s1 = _build($v, $p1);
        my ($c0, $c1) = Math::Vector::Real->box(@{$s0}[_c0, _c1], @{$s1}[_c0, _c1]);
        my $cut = 0.5 * ($s0->[_c1][$axis] + $s1->[_c0][$axis]);
        # [n sum s0 s1 axis cut]
        [scalar(@$ix), $c0, $c1, $s0->[_sum] + $s1->[_sum], $s0, $s1, $axis, $cut];
    }
    else {
        # [n, sum, ixs]
        my @vs = @{$v}[@$ixs];
        my ($c0, $c1) = Math::Vector::Real->box(@vs);
        [scalar(@$ix), $c0, $c1, Math::Vector::Real->sum(@vs), $ix];
    }
}

sub size { scalar @{shift->{vs}} }

sub at {
    my ($self, $ix) = @_;
    Math::Vector::Real::clone($self->{vs}[$ix]);
}

sub insert {
    my $self = shift;
    return undef unless @_;
    my $vs = $self->{vs};
    my $ix = @$vs;
    for (@_) {
        my $v = Math::Vector::Real::clone($_);
        push @$vs, $v;
        _insert($vs, $self->{tree}, $#$vs)
    }
    $ix;
}

# _insert does not return anything but modifies its $t argument in
# place. This is really ugly but done to improve performance.

sub _insert {
    my ($vs, $t, $ix) = @_;
    my $n = $t->[_n]++;
    my $v = $vs->[$ix];
    @{$t}[_c0, _c1] = Math::Vector::Real->box(@{$t}[_c0, _c1], $v)

    if (defined (my $axis = $t->[_axis])) {
        my $cut = $t->[_cut];
        my $c = $v->[$axis];

        my $n0 = $t->[_s0][_n];
        my $n1 = $t->[_s1][_n];

        if ($c <= $cut) {
            if (2 * $n1 + $max_per_pole >= $n0) {
                _insert($vs, $t->[_s0], $ix);
                return;
            }
        }
        else {
            if (2 * $n0 + $max_per_pole >= $n1) {
                _insert($vs, $t->[_s1], $ix);
                return;
            }
        }

        # tree needs rebalancing
        my @store;
        $#store = $n; # preallocate space
        @store = ($ix);
        _push_all($t, \@store);
        $_[1] = _build($vs, \@store);
    }
    else {
        my $ixs = $t->[_ixs];
        push @$ixs, $ix;
        if ($n > $max_per_pole) {
            $_[1] = _build($vs, $ixs);
        }
    }
}

sub move {
    my ($self, $ix, $v) = @_;
    my $vs = $self->{vs};
    ($ix >= 0 and $ix < @$vs) or croak "index out of range";
    _delete($vs, $self->{tree}, $ix);
    $vs->[$ix] = Math::Vector::Real::clone($v);
    _insert($vs, $self->{tree}, $ix);
}

sub _delete {
    my ($vs, $t, $ix) = @_;
    if (defined $t->[_axis]) {
        my ($axis, $s0, $s1, $median) = @$t;
        #print "axis: $axis, ix: $ix\n";
        my $c = $vs->[$ix][$axis];
        if ($c <= $median and _delete($vs, $s0, $ix)) {
            if ($#$s0 == 0) {
                # when one subnode is empty, promote the other up:
                @$t = @$s1;
            }
            else {
                $t->[_n]--;
            }
            return 1;
        }
        elsif ($c >= $median and _delete($vs, $s1, $ix)) {
            if ($#$s1 == 0) {
                @$t = @$s0;
            }
            else {
                $t->[_n]--;
            }
            return 1;
        }
        return 0;
    }
    else {
        my $s0 = scalar @$t;
        @$t = grep { not (defined($_) and ($_ == $ix)) } @$t;
        return @$t < $s0;
    }
}

sub hide {
    my ($self, $ix) = @_;
    my $vs = $self->{vs};
    ($ix >= 0 and $ix < @$vs) or croak "index out of range";
    _delete($vs, $self->{tree}, $ix);
    ($self->{hidden} //= {})->{$ix} = 1;
}

sub _push_all {
    my ($t, $store) = @_;
    if (defined $t->[_axis]) {
        _push_all($t->[_s0], $store);
        _push_all($t->[_s1], $store);
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
    if (defined $t->[_axis]) {
        for (0, 1) {
            my $p = _path($t->[_s0 + $_], $vix);
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
        if (defined $t->[_axis]) {
            my ($axis, $s0, $s1, $median) = @$t;
            my $c = $v->[$axis];
            if ($c < $median) {
                $t = $s0;
            }
            else {
                if ($c == $median) {
                    my $ix = _find($vs, $s0, $v);
                    return $ix if defined $ix;
                }
                $t = $s1;
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
    my ($self, $v, $d, @but) = @_;
    my $vs = $self->{vs};
    my $but;
    if (@but) {
        if (@but == 1 and ref $but[0] eq 'HASH') {
            $but = $but[0];
        }
        else {
            my %but = map { $_ => 1 } @but;
            $but = \%but;
        }
    }

    my ($start, $d2);
    if (defined $d) {
        $d2 = $d * $d;
    }
    else {
        my $hidden = $self->{hidden};
        for ($start = 0;
             $start < @$vs or return;
             $start++) {
            last unless ( ( $hidden and $hidden->{$start} ) or
                          ( $but    and $but->{$start}    ) );
        }
        $d2 = $vs->[$start]->dist2($v);
    }
    my ($rix, $rd2) = _find_nearest_neighbor($vs, $self->{tree}, $v, $start, $d2, $but);
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
        if (defined $t->[_axis]) {
            my ($axis, $s0, $s1, $median) = @$t;
            my $c = $v->[$axis];
            my $cm = $c - $median;
            (my ($first), $t) = (($cm <= 0) ? ($s0, $s1) : ($s1, $s0));
            ($ix, $d2) = _find_nearest_neighbor($vs, $first, $v, $ix, $d2, $but);
            return ($ix, $d2) if $d2 <= $cm * $cm;
        }
        else {
            for (@$t[1..$#$t]) {
                next if $but and $but->{$_};
                my $p = $vs->[$_];
                my $d21 = $p->dist2($v);
                if ($d21 < $d2) {
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
    if ($self->{hidden}) {
        $best[$_] = undef for keys %{$self->{hidden}};
    }
    return @best;
}

sub _find_nearest_neighbor_all_internal {
    my ($vs, $t, $best, $d2) = @_;
    if (defined $t->[_axis]) {
        my ($axis, $s0, $s1, $median) = @$t;
        my @r;
        for my $side (0, 1) {
            my @poles = _find_nearest_neighbor_all_internal($vs, $t->[_s0 + $side], $best, $d2);
            push @r, @poles;
            my $other = $t->[_s1 - $side];
            for my $pole (@poles) {
                for my $ix (@$pole[1..$#$pole]) {
                    my $v = $vs->[$ix];
                    my $md = $v->[$axis] - $median;
                    if ($d2->[$ix] > $md * $md) {
                        ($best->[$ix], $d2->[$ix]) =
                            _find_nearest_neighbor($vs, $other, $v, $best->[$ix], $d2->[$ix]);
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
    _find_in_ball($self->{vs}, $self->{tree}, $z, $d * $d, $but);
}

sub _find_in_ball {
    my ($vs, $t, $z, $d2, $but) = @_;
    if (defined $t->[_axis]) {
        my ($axis, $s0, $s1, $median) = @$t;
        my $c = $z->[$axis];
        my $dc = $c - $median;
        my ($f, $s) = (($dc < 0) ? ($s0, $s1) : ($s1, $s0));
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
    elsif ($but) {
        return grep { !$but->{$_} and $vs->[$_]->dist2($z) <= $d2 } @$t[1..$#$t]
    }
    else {
        return grep { $vs->[$_]->dist2($z) <= $d2 } @$t[1..$#$t]
    }
}

sub ordered_by_proximity {
    my $self = shift;
    my @r;
    $#r = $#{$self->{vs}}; $#r = -1; # preallocate
    _ordered_by_proximity($self->{tree}, \@r);
    return @r;
}

sub _ordered_by_proximity {
    my $t = shift;
    my $r = shift;
    if (defined $t->[_axis]) {
        _ordered_by_proximity($t->[_s0], $r);
        _ordered_by_proximity($t->[_s1], $r);
    }
    else {
        push @$r, @{$t}[1..$#$t];
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

Creates a new kd-Tree containing the given points.

=item $t2 = $t->clone

Creates a duplicate of the tree. The two trees will share internal
read only data so this method is more efficient in terms of memory
usage than others performing a deep copy.

=item my $ix = $t->insert($p0, $p1, ...)

Inserts the given points into the kd-Tree.

Returns the index assigned to the first point inserted.

=item $s = $t->size

Returns the number of points inside the tree.

=item $p = $t->at($ix)

Returns the point at the given index inside the tree.

=item $t->move($ix, $p)

Moves the point at index C<$ix> to the new given position readjusting
the tree structure accordingly.

=item ($ix, $d) = $t->find_nearest_neighbor($p, $max_d, @but_ix)

=item ($ix, $d) = $t->find_nearest_neighbor($p, $max_d, \%but_ix)

Find the nearest neighbor for the given point C<$p> and returns its
index and the distance between the two points (in scalar context the
index is returned).

If C<$max_d> is defined, the search is limited to the points within that distance

Optionally, a list of point indexes to be excluded from the search can be
passed or, alternatively, a reference to a hash containing the indexes
of the points to be excluded.

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

If the extra argument C<$but> is provided. The point with that index
is ignored.

=item @ix = $t->ordered_by_proximity

Returns the indexes of the points in an ordered where is likely that
the indexes of near vectors are also in near positions in the list.

=back

=head1 SEE ALSO

L<http://en.wikipedia.org/wiki/K-d_tree>

L<Math::Vector::Real>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011-2013 by Salvador FandiE<ntilde>o E<lt>sfandino@yahoo.comE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.3 or,
at your option, any later version of Perl 5 you may have available.


=cut
