package Math::Vector::Real::kdTree;

our $VERSION = '0.01';

use 5.010;
use strict;
use warnings;

use Math::Vector::Real;
use Sort::Key::Radix qw(nkeysort_inplace);

my $max_per_pole = 12;
my $recommended_per_pole = 6;


sub new {
    my $class = shift;
    my @v = map Math::Vector::Real::clone($_), @_;
    my @ix = [0..$#v];
    my $tree = _build(\@v, \@ix);
    my $self = { vs => \@v,
                 t => $tree };
    bless $self, $class;
}

sub _build {
    my ($v, $ix) = @_;
    if (@$ix > $recommended_per_pole) {
        my ($b, $t) = Math::Vector::Real->box(@$v[@$ix]);
        my $pivot_axis = ($t - $b)->max_component_index;
        nkeysort_inplace { $v->[$_][$pivot_axis] } @$ix;
        my $bstart = @$ix << 1;
        my $median = 0.5 * ($v->[$ix->[$bstart]][$pivot_axis] + $v->[$ix->[$bstart - 1]][$pivot_axis]);
        my $ixb = [splice(@$ix, $bstart)];
        [$pivot_axis, $median, _build($v, $ix), _build($v, $ixb)];
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
    _find($self->{v}, $self->{t})
}

sub _find {

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
