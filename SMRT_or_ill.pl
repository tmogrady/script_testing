#!/usr/bin/perl

#checks to see if SMRT 5' ends are usually upstream, downstream, or exact matches to CAGE 5' ends

#TO'G 11/11/15

use warnings;
use strict;

my ($SMRT_ends, $CAGE_ends) = @ARGV;

open(INF, "<$SMRT_ends");
open(OUT, ">$SMRT_ends.downstream.bed");
open(OUT2, ">$SMRT_ends.upstream.bed");
open(OUT3, ">$SMRT_ends.equal.bed");

my $chrStart_SMRT;

while (my $line = <INF>) {
    next if ($line =~ /^track/);
    chomp($line);
    my @cols = split("\t", $line);
    $chrStart_SMRT = $cols[1];
    open(INF2, "<$CAGE_ends");
    while (my $line2 = <INF2>) {
        chomp($line2);
        my @cols2 = split("\t", $line2);
        if ((abs ($cols2[1] - $chrStart_SMRT)) <= 2) {
            if ($cols[5] eq "+") {
                if ($cols2[1] < $chrStart_SMRT) {
                    print OUT $line, "\n";
                }
                elsif ($cols2[1] > $chrStart_SMRT) {
                    print OUT2 $line, "\n";
                }
                else {
                    print OUT3 $line, "\n";
                }
            }
            elsif ($cols[5] eq "-") {
                if ($cols2[1] < $chrStart_SMRT) {
                    print OUT2 $line, "\n";
                }
                elsif ($cols2[1] > $chrStart_SMRT) {
                    print OUT $line, "\n";
                }
                else {
                    print OUT3 $line, "\n";
                }
            }
        }
    }
    close(INF2);
}
close(INF);
close(OUT);
close(OUT2);
close(OUT3);