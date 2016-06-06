#!/usr/bin/perl

use warnings;
use strict;

my ($file, $file2) = @ARGV; #$file is the bed file with conservation information, $file2 is the .txt file of rMATS output (with duplicates removed by rMATS_to_bed_introns.pl)

open(INF, "<$file");

my @conserv;

while (my $line = <INF>) {
    chomp($line);
    push @conserv, $line;
}

close(INF);

open(INF, "<$file2");
open(OUT, ">$file2.motif_conservation.txt");

my $up_one = 0;
my $up_two = 0;
my $up_five = 0;
my $down_one = 0;
my $down_two = 0;
my $down_five = 0;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    next if $cols[0] eq "ID";
    foreach my $conserv (@conserv) {
        my @cols2 = split("\t", $conserv);
        my ($name) = $cols2[3] =~ /_"(.+)"/;
        my ($exonStart) = $cols2[3] =~ /\:(\d+)-/;
        my ($exonEnd) = $cols2[3] =~ /\-(\d+)\:/;
        my ($gene_name) = $cols[1] =~ /"(.+)"/;
        if ((($name eq $gene_name) and ($exonStart == $cols[5])) and ($exonEnd == $cols[6])) {
            if ($cols2[3] =~ /^up/) {
                $up_one = $up_one + 1;
                if ($cols2[3] =~ /hg19.mm10/) {
                    $up_two = $up_two + 1;
                    if ($cols2[3] =~ /hg19.mm10.rn5.canFam3.galGal4/) {
                        $up_five = $up_five + 1;
                    }
                }
            }
            elsif ($cols2[3] =~ /^down/) {
                $down_one = $down_one + 1;
                if ($cols2[3] =~ /hg19.mm10/) {
                    $down_two = $down_two +1;
                    if ($cols2[3] =~ /hg19.mm10.rn5.canFam3.galGal4/) {
                        $down_five = $down_five + 1;
                    }
                }
            }
        }
        else {
            next;
        }
    }
    print OUT "$line\t$up_one\t$up_two\t$up_five\t$down_one\t$down_two\t$down_five\n";
    $up_one = 0;
    $up_two = 0;
    $up_five = 0;
    $down_one = 0;
    $down_two = 0;
    $down_five = 0;
}

close(INF);
close(OUT);