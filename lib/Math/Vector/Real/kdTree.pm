package Math::Vector::Real::kdTree;

our $VERSION = '0.01';

use 5.010;
use strict;
use warnings;
use Carp;

use Math::Vector::Real;
use Sort::Key::Top qw(nkeypartref);

my $max_per_pole = 12;
my $recommended_per_pole = 6;


sub new {
    my $class = shift;
    my @v = map Math::Vector::Real::clone($_), @_;
    my @ix = [0..$#v];
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
        my $bstart = @$ix << 1;
        my ($l, $r) = nkeypartref { $v->[$_][$axis] } $bstart => @$ix;
        my $median = 0.5 * ($v->[$l->[-1]][$axis] + $v->[$r->[0]][$axis]);
        [$axis, _build($v, $l), _build($v, $r), $mediam, $b->[$axis], $t->[$axis]];
    }
    else {
        [undef, @$ix];
    }
}

sub at {
    my ($self, $ix) = @_;
    $self->{vs}[$ix]
}

sub find {
    my ($self, $v) = @_;
    _find($self->{vs}, $self->{tree}, $v)
}

sub _find {
    my ($vs, $t, $v) = @_;
    while (1) {
        if (defined $t->[0]) {
            my ($axis, $l, $r, $mediam, $min, $max) = @$t;
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
            for (@$t[1..$#t]) {
                return $_ if $v == $vs->[$_];
            }
            return undef;
        }
    }
}

sub find_nearest_neighbor {
    my ($self, $v) = @_;
    my $vs = $self->{vs};
    return unless @$vs;
    _find_nearest_neighbor($self->{vs}, $self->{tree}, $v, 0, $vs->[0]->dist2($v), -1);
}

sub find_nearest_neighbor_internal {
    my ($self, $vix) = @_;
    my $vs = $self->{vs};
    $vix >= @$vs and croak "index out of range";
    return unless @$vs > 1;
    my $start = ($vix ? 0 : 1);
    _find_nearest_neighbor($self->{vs}, $self->{tree}, $v, $start, $vs->[$start]->dist2($v), $vix);
}

sub _find_nearest_neighbor {
    my ($vs, $t, $v, $ix, $d2, $but) = @_;
    while (1) {
        if (defined $t->[0]) {
            my ($axis, $median, $l, $r) = @$t;
            my $c = $v->[$axis];
            my $cm = $c - $median;
            (my ($first), $t) = (($cm <= 0) ? ($l, $r) : ($r, $l));
            ($ix, $d2) = _find_nearest_neighbor($vs, $first, $v, $ix, $d2, $but);
            return ($ix, $d2) if $d2 <= $cm * $cm;
        }
        else {
            for (@$t[1..$#t]) {
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



1;
__END__

=head1 NAME

Math::Vector::Real::kdTree - Perl extension for blah blah blah

=head1 SYNOPSIS

  use Math::Vector::Real::kdTree;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for Math::Vector::Real::kdTree, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Salvador Fandiño, E<lt>salva@E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011 by Salvador Fandiño

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.3 or,
at your option, any later version of Perl 5 you may have available.


=cut
