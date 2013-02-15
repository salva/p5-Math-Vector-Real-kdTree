#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 6617;

use_ok('Math::Vector::Real::kdTree');

use Math::Vector::Real;

sub neighbors_bruteforce {
    my ($bottom, $top) = Math::Vector::Real->box(@_);
    my $box = $top - $bottom;
    my $v = [map $_ - $bottom, @_];
    my $ixs = [0..$#_];
    my $dist2 = [($box->abs2 * 10 + 1) x @_];
    my $neighbors = [(undef) x @_];
    _neighbors_bruteforce($v, $ixs, $dist2, $neighbors, $box, 0);
    return @$neighbors;
}

sub _neighbors_bruteforce {
    my ($v, $ixs, $dist2, $neighbors) = @_;
    my $ixix = 0;
    for my $i (@$ixs) {
        $ixix++;
        my $v0 = $v->[$i];
        for my $j (@$ixs[$ixix..$#$ixs]) {
            my $d2 = $v0->dist2($v->[$j]);
            if ($dist2->[$i] > $d2) {
                $dist2->[$i] = $d2;
                $neighbors->[$i] = $j;
            }
            if ($dist2->[$j] > $d2) {
                $dist2->[$j] = $d2;
                $neighbors->[$j] = $i;
            }
        }
    }
}

sub test_neighbors {
    my ($o, $n1, $n2, $msg) = @_;
    my (@d1, @d2);
    for my $ix (0..$#$o) {
        my $eo   = $o->[$ix];
        my $ixn1 = $n1->[$ix];
        my $ixn2 = $n2->[$ix];
        my $en1  = $o->[$ixn1];
        my $en2  = $o->[$ixn2];
        push @d1, $eo->dist2($en1);
        push @d2, $eo->dist2($en2);
    }
    is "@d1", "@d2", $msg;
}

my %gen = ( num => sub { rand },
            int => sub { int rand(10) } );

diag "srand: " . srand;
for my $g (keys %gen) {
    for my $d (1, 2, 3, 10) {
        for my $n (2, 10, 50, 250, 500) {
            my $id = "gen: $g, d: $d, n: $n";
            my @o = map V(map $gen{$g}->(), 1..$d), 1..$n;
            my @nbf = neighbors_bruteforce(@o);

            my $t = Math::Vector::Real::kdTree->new(@o);
            my @n = $t->find_nearest_neighbor_all_internal;
            is ($#n, $#o, "find_nearest_neighbor_all_internal - build - $id");
            test_neighbors(\@o, \@n, \@nbf, "find_nearest_neighbor_all_internal - build - $id");

            $t = Math::Vector::Real::kdTree->new;
            for my $ix (0..$#o) {
                $t->insert($o[$ix]);
                my @obp = $t->ordered_by_proximity;
                is ($ix, $#obp, "ordered_by_proxymity - count - $id, ix: $ix");
            }
            @n = $t->find_nearest_neighbor_all_internal;
            test_neighbors(\@o, \@n, \@nbf, "find_nearest_neighbor_all_internal - insert - $id");
        }
    }
}

