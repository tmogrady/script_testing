#!/usr/bin/perl

#TO'G 11/25/15

use warnings;
use strict;

my ($file) = @ARGV;

open(INF, "<$file") or die "couldn't open input file";
open(OUT, ">$file.MHV68_gembo_loxP.bed") or die "couldn't open output file";

my $chr = "MHV68_gembo_loxP";
my $chrStart;
my $chrEnd;
my $thickStart;
my $thickEnd;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    if ($line =~ /^track/) {
        print OUT "track name=MHV68_gembo_loxP_converted_annotation\n";
    }
    else {
        if ($cols[1] <= 103928) {
            $chrStart = $cols[1] - 22;
            $thickStart = $cols[6] - 22;
            if ($cols[2] <= 103928) {
                $chrEnd = $cols[2] - 22;
                $thickEnd = $cols[7] - 22;
            }
            elsif ($cols[2] > 103928) {
                $chrEnd = $cols[2] + 1168;
                $thickEnd = $cols[7] + 1168;
            }
        }
        else {
            $chrStart = $cols[1] + 1168;
            $chrEnd = $cols[2] + 1168;
            $thickStart = $cols[6] + 1168;
            $thickEnd = $cols[7] + 1168;
        }
        print OUT "$chr\t$chrStart\t$chrEnd\t$cols[3]\t$cols[4]\t$cols[5]\t$thickStart\t$thickEnd\t$cols[8]\t$cols[9]\t$cols[10]\t$cols[11]\n";
    }
}

close(INF);
close(OUT);