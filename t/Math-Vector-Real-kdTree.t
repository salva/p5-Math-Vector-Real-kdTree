#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 14449;

use_ok('Math::Vector::Real::kdTree');

use Math::Vector::Real;
use Math::Vector::Real::Test qw(eq_vector);

sub nearest_vectors_bruteforce {
    my ($bottom, $top) = Math::Vector::Real->box(@_);
    my $box = $top - $bottom;
    my $v = [map $_ - $bottom, @_];
    my $ixs = [0..$#_];
    my $dist2 = [($box->abs2 * 10 + 1) x @_];
    my $neighbors = [(undef) x @_];
    _nearest_vectors_bruteforce($v, $ixs, $dist2, $neighbors, $box, 0);
    return @$neighbors;
}

sub _nearest_vectors_bruteforce {
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

sub farthest_vectors_bruteforce {
    my @best_ix;
    my @best_d2 = ((0) x @_);
    for my $i (1..$#_) {
        my $v = $_[$i];
        for my $j (0..$i - 1) {
            my $d2 = Math::Vector::Real::dist2($v, $_[$j]);
            if ($d2 > $best_d2[$i]) {
                $best_d2[$i] = $d2;
                $best_ix[$i] = $j;
            }
            if ($d2 > $best_d2[$j]) {
                $best_d2[$j] = $d2;
                $best_ix[$j] = $i;
            }
        }
    }
    return @best_ix;
}

sub test_neighbors {
    unshift @_, $_[0];
    goto &test_neighbors_indirect;
}

sub test_neighbors_indirect {
    my ($o1, $o2, $n1, $n2, $msg) = @_;
    my (@d1, @d2);
    for my $ix (0..$#$o1) {
        my $eo   = $o1->[$ix];
        my $ixn1 = $n1->[$ix];
        defined $ixn1 or do {
            fail("expected index for element $ix is undefined");
            goto break_me;
        };
        my $ixn2 = $n2->[$ix];
        defined $ixn2 or do {
            fail("template index for element $ix is undefined");
            goto break_me;
        };
        $ixn1 < @$o2 or do {
            fail("expected index $ixn1 out of range");
            goto break_me;
        };
        $ixn2 < @$o2 or do {
            fail("template index $ixn1 out of range");
            goto break_me;
        };
        my $en1 = $o2->[$ixn1];
        my $en2  = $o2->[$ixn2];
        push @d1, $eo->dist2($en1);
        push @d2, $eo->dist2($en2);
    }
    is "@d1", "@d2", $msg and return 1;

 break_me:
    diag "break me!";
    0;
}

my %gen = ( num => sub { rand },
            int => sub { int rand(10) } );

#srand 318275924;
diag "srand: " . srand;
for my $g (keys %gen) {
    for my $d (1, 2, 3, 10) {
        for my $n (2, 10, 50, 250, 500) {
            my $id = "gen: $g, d: $d, n: $n";
            my @o = map V(map $gen{$g}->(), 1..$d), 1..$n;
            my @nbf = nearest_vectors_bruteforce(@o);

            my $t = Math::Vector::Real::kdTree->new(@o);

            my @n = map scalar($t->find_nearest_vector_internal($_)), 0..$#o;
            is ($#n, $#o, "count find_nearest_vector_internal - build - $id");
            test_neighbors(\@o, \@n, \@nbf, "find_nearest_vector_internal - build - $id");
            is_deeply([map $t->at($_), 0..$#o], \@o , "at - build - after find_nearest_vector_internal - $id");

            @n = $t->find_nearest_vector_all_internal;
            is ($#n, $#o, "count find_nearest_vector_all_internal - build - $id");
            test_neighbors(\@o, \@n, \@nbf, "find_nearest_vector_all_internal - build - $id");
            is_deeply([map $t->at($_), 0..$#o], \@o , "at - build - after find_nearest_vector_all_internal - $id");

            $t = Math::Vector::Real::kdTree->new;
            for my $ix (0..$#o) {
                $t->insert($o[$ix]);
                my @obp = $t->ordered_by_proximity;
                is ($ix, $#obp, "ordered_by_proxymity - count - $id, ix: $ix");
            }
            is_deeply([map $t->at($_), 0..$#o], \@o , "at - insert - after insert - $id");

            @n = map scalar($t->find_nearest_vector_internal($_)), 0..$#o;
            test_neighbors(\@o, \@n, \@nbf, "find_nearest_vector_internal - insert - $id");
            is_deeply([map $t->at($_), 0..$#o], \@o , "at - insert - after find_nearest_vector_internal - $id");

            @n = $t->find_nearest_vector_all_internal;
            test_neighbors(\@o, \@n, \@nbf, "find_nearest_vector_all_internal - insert - $id");
            is_deeply([map $t->at($_), 0..$#o], \@o , "at - insert - after find_nearest_vector_all_internal - $id");

            my @fbf = farthest_vectors_bruteforce(@o);
            @n = map scalar($t->find_farthest_vector_internal($_)), 0..$#o;
            test_neighbors(\@o, \@n, \@fbf, "find_farthest_vector_internal - insert - $id");
            is_deeply([map $t->at($_), 0..$#o], \@o , "at - insert - after find_farthest_vector_internal - $id");

            my $k;
            for ($k = 1; $k < @n; $k *= 2) {
                my @kms = $t->k_means_start($k);
                is (scalar(@kms), $k, "k_means_start generates $k results - $id");
                my @km = $t->k_means_loop(@kms);
                is (scalar(@km), $k, "k_means_loop generates $k results - $id")
                    or do {
                        diag "break me 2";
                    };
                my @kma = $t->k_means_assign(@km);
                my $t1 = Math::Vector::Real::kdTree->new(@km);
                my @n = map scalar($t1->find_nearest_vector($_)), @o;
                test_neighbors_indirect(\@o, \@km, \@kma, \@n, "k-means assign - k: $k, $id");

                my @sum = map V((0) x $d), 1..$k;
                my @count = ((0) x $k);

                for my $ix (0..$#kma) {
                    my $cluster = $kma[$ix];
                    $count[$cluster]++;
                    $sum[$cluster] += $o[$ix];
                }
                for my $cluster (0..$#sum) {
                    if ($count[$cluster]) {
                        $sum[$cluster] /= $count[$cluster];
                    }
                    else {
                        $sum[$cluster] = $km[$cluster];
                    }

                    eq_vector($sum[$cluster], $km[$cluster], "cluster centroid - $cluster - k: $k, $id");
                }
            }

        }
    }
}
