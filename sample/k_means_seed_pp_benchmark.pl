#!/usr/bin/perl

use strict;
use warnings;

use Benchmark qw(cmpthese);

use Math::Vector::Real::kdTree;
use Math::Vector::Real::Random;

my $n = shift @ARGV || 10000;

$| = 1;

print "generating random data...\n";
my @v = map Math::Vector::Real->random_normal(2), 1..$n;
my $t = Math::Vector::Real::kdTree->new(@v);

print "benchmarking...\n";

sub k_means_seed_pp_brute_force {
    my ($vs, $n) = @_;
    my $kmv = $vs->[rand @$vs];
    my @kmv = $kmv;
    my @d2 = map $kmv->dist2($_), @$vs;
    while (@kmv < $n) {
        my $w = 0;
        my @w;
        for my $ix (0..$#$vs) {
            my $d2 = $kmv->dist2($vs->[$ix]);
            my $old_d2 = $d2[$ix];
            if ($d2 < $old_d2) {
                $d2[$ix] = $d2;
                $w += $d2;
            }
            else {
                $w += $old_d2;
            }
            $w[$ix] = $w;
        }

        my $dice = rand($w);
        my $i = 0;
        my $j = @w;
        while ($i < $j) {
            my $pivot = (($i + $j) >> 1);
            if ($w[$pivot] < $dice) {
                $i = $pivot + 1;
            }
            else {
                $j = $pivot;
            }
        }
        $kmv = $vs->[$i];
        push @kmv, $kmv;
    }
}

my @sizes = (2, 8, 32);

my %subs = ( bf => sub { k_means_seed_pp_brute_force(\@v, $_) for @sizes });
for my $i (1, 3, 5, 7, 9, 10) {
    my $p = 0.1 * $i;
    $subs{"p$p"} = sub { $t->k_means_seed_pp($_, $p) for @sizes };
}

cmpthese(-1, \%subs);
