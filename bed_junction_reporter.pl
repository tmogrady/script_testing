#!/usr/bin/perl

#given a bed file, reports the coordinates of all splice junctions

#TO'G 11/30/15

use warnings;
use strict;

my ($file) = @ARGV;

open(INF, "<$file");
open(OUT, ">${file}_junctions.txt");

my @intron_start;
my @intron_end;
my %ann_intron_coord_pair;
my $start;
my $end;

while (my $line = <INF>) {
    chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
    my @cols = split("\t", $line);
    my $intron_number = $cols[9] - 1;
    next if ($intron_number == 0);
    my @block_sizes = split(",", $cols[10]);
    my @block_starts = split(",", $cols[11]);
    for (my $i = 0; $i < $intron_number; $i = $i + 1) { #for the transcript currently in the "while" loop, creates an array of intron start sites relative to the genome
        $start = $cols[1] + $block_sizes[$i] + $block_starts[$i];
        push(@intron_start, $start);
    }
    for (my $i2 = 1; $i2 < $cols[9]; $i2 = $i2 + 1) { #for the transcript currently in the "while" loop, creates an array of intron end sites relative to the genome
        $end = $cols[1] + $block_starts[$i2];
        push(@intron_end, $end);
    }
    for (my $i3 = 0; $i3 < $intron_number; $i3 = $i3 + 1) { #for the transcript currently in the "while" loop, matches up intron start and end sites to create a hash of complete intron coordinates relative to the genome
#        my $intron_coords = "$cols[0]:$intron_start[$i3]:$intron_end[$i3]:$cols[5]";
        print OUT "$cols[0]\t$intron_start[$i3]\t$intron_end[$i3]\t$cols[5]\n";
#        if (exists $ann_intron_coord_pair{$intron_coords}) {
#            $ann_intron_coord_pair{$intron_coords} = $ann_intron_coord_pair{$intron_coords} + 1; #if the intron is already in the hash (from another transcript), increase the count
#        }
#        else {
#            $ann_intron_coord_pair{$intron_coords} = 1; #if the intron is not already in the hash, adds it with a value of 1
#        }
    }
    @intron_start = ();
    @intron_end = (); #intron starts and ends have been assigned to the %ann_intron_pair hash; empty them for the next transcript
}
close(INF);
close(OUT);