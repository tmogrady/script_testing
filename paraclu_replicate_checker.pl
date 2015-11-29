#!/usr/bin/perl

#Input: weighted average bed files from paraclu_pipe.pl (for 4 replicates). Returns files consensus files of clusters.

#USAGE:
# perl PATH/paraclu_filter.pl numerical_max_distance /PATH/inputfiles

#TO'G 10/28/15

use warnings;
use strict;

my ($dist, $file1, $file2, $file3, $file4) = @ARGV;
my %peaks;

open(INF, "<$file1");

print "Processing first file\n";

my $coord_key;
my $values;

while(my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    my @subcols = split("\:", $cols[3]);
    if ($cols[5] eq "+") {
        $coord_key = "$cols[1]\:$cols[5]"; #creates a key for the hash with chrStart and strand
        $values = "$cols[4]\:$subcols[3]\:1"; #creates a value for the hash with tag count, relative density and replicate count
        $peaks{$coord_key} = $values;
    }
    if ($cols[5] eq "-") {
        $coord_key = "$cols[2]\:$cols[5]"; #creates a key for the hash with chrEnd and strand
        $values = "$cols[4]\:$subcols[3]\:1"; #creates a value for the hash with tag count, relative density and replicate count
        $peaks{$coord_key} = $values;
    }
}

#foreach (sort keys %peaks) {
#    print "$_: $peaks{$_}\n";
#}

close(INF);

open(INF, "<$file2");

print "Processing second file\n";

my $file_start;

while(my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    my @subcols = split("\:", $cols[3]);
    my $found_flag = 0;
    if ($cols[5] eq "+") {
        $file_start = $cols[1];
    }
    if ($cols[5] eq "-") {
        $file_start = $cols[2];
    }
    my $lower_limit = $file_start - $dist;
    my $upper_limit = $file_start + $dist;
    foreach my $hash_start (keys %peaks) {
        my ($coord, $strand) = split("\:", $hash_start);
        if (($coord >= $lower_limit) and ($coord <= $upper_limit)) { #checks to see if the cluster is already in the hash (from an earlier file)
            my @hashcols = split("\:", $peaks{$hash_start}); #extracts information about the cluster from the hash
            my $depth = $cols[4] + $hashcols[0]; #adds tag depths from the current file and the hash
            my $reps = $hashcols[2]+1; #adds one to the replicate count
            my $density = (($subcols[3] + ($hashcols[1]*$hashcols[2]))/($hashcols[2]+1));
            my $new_start = (($coord*$hashcols[0]) + ($file_start*$cols[4]))/($hashcols[0] + $cols[4]);

            $coord_key = "$new_start\:$cols[5]";
            $values = "$depth\:$density\:$reps";
            delete $peaks{$hash_start};
            $peaks{$coord_key} = $values;
            $found_flag = 1;
            last;
        }
    }
    if ($found_flag == 0) {
        $coord_key = "$file_start\:$cols[5]";
        $values = "$cols[4]\:$subcols[3]\:1";
        $peaks{$coord_key} = $values; #adds the coordinate to the hash as a key, with the value assigned above as a value
    }
    
}

open(INF, "<$file3");

print "Processing third file\n";

while(my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    my @subcols = split("\:", $cols[3]);
    my $found_flag = 0;
    if ($cols[5] eq "+") {
        $file_start = $cols[1];
    }
    if ($cols[5] eq "-") {
        $file_start = $cols[2];
    }
    my $lower_limit = $file_start - $dist;
    my $upper_limit = $file_start + $dist;
    foreach my $hash_start (keys %peaks) {
        my ($coord, $strand) = split("\:", $hash_start);
        if (($coord >= $lower_limit) and ($coord <= $upper_limit)) { #checks to see if the cluster is already in the hash (from an earlier file)
            my @hashcols = split("\:", $peaks{$hash_start}); #extracts information about the cluster from the hash
            my $depth = $cols[4] + $hashcols[0]; #adds tag depths from the current file and the hash
            my $reps = $hashcols[2]+1; #adds one to the replicate count
            my $density = (($subcols[3] + ($hashcols[1]*$hashcols[2]))/($hashcols[2]+1));
            my $new_start = (($coord*$hashcols[0]) + ($file_start*$cols[4]))/($hashcols[0] + $cols[4]);
            
            $coord_key = "$new_start\:$cols[5]";
            $values = "$depth\:$density\:$reps";
            delete $peaks{$hash_start};
            $peaks{$coord_key} = $values;
            $found_flag = 1;
            last;
        }
    }
    if ($found_flag == 0) {
        $coord_key = "$file_start\:$cols[5]";
        $values = "$cols[4]\:$subcols[3]\:1";
        $peaks{$coord_key} = $values; #adds the coordinate to the hash as a key, with the value assigned above as a value
    }
    
}

open(INF, "<$file4");

print "Processing fourth file\n";

while(my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    my @subcols = split("\:", $cols[3]);
    my $found_flag = 0;
    if ($cols[5] eq "+") {
        $file_start = $cols[1];
    }
    if ($cols[5] eq "-") {
        $file_start = $cols[2];
    }
    my $lower_limit = $file_start - $dist;
    my $upper_limit = $file_start + $dist;
    foreach my $hash_start (keys %peaks) {
        my ($coord, $strand) = split("\:", $hash_start);
        if (($coord >= $lower_limit) and ($coord <= $upper_limit)) { #checks to see if the cluster is already in the hash (from an earlier file)
            my @hashcols = split("\:", $peaks{$hash_start}); #extracts information about the cluster from the hash
            my $depth = $cols[4] + $hashcols[0]; #adds tag depths from the current file and the hash
            my $reps = $hashcols[2]+1; #adds one to the replicate count
            my $density = (($subcols[3] + ($hashcols[1]*$hashcols[2]))/($hashcols[2]+1));
            my $new_start = (($coord*$hashcols[0]) + ($file_start*$cols[4]))/($hashcols[0] + $cols[4]);
            
            $coord_key = "$new_start\:$cols[5]";
            $values = "$depth\:$density\:$reps";
            delete $peaks{$hash_start};
            $peaks{$coord_key} = $values;
            $found_flag = 1;
            last;
        }
    }
    if ($found_flag == 0) {
        $coord_key = "$file_start\:$cols[5]";
        $values = "$cols[4]\:$subcols[3]\:1";
        $peaks{$coord_key} = $values; #adds the coordinate to the hash as a key, with the value assigned above as a value
    }
    
}

open(OUT, ">any_replicate.temp");
open(OUT2, ">two_replicates.temp");
open(OUT3, ">three_replicates.temp");
open(OUT4, ">four_replicates.temp");
open(OUT5, ">weird_replicates.temp");

my $chrStart;
my $chrEnd;

foreach my $hash_start (sort keys %peaks) {
    my ($coord, $strand) = split("\:", $hash_start);
    if ($strand eq "+") {
        $chrStart = sprintf("%1.0f", $coord);
        $chrEnd = $chrStart + 1;
    }
    if ($strand eq "-") {
        $chrEnd = sprintf("%1.0f", $coord);
        $chrStart = $chrEnd - 1;
    }
    my @hashcols = split("\:", $peaks{$hash_start});
    my $out_line = "chrEBV(Akata_107955to171322_1to107954)\t$chrStart\t$chrEnd\t$peaks{$hash_start}\t$hashcols[0]\t$strand\n";
    print OUT $out_line;
    if ($hashcols[2] > 1) {
        print OUT2 $out_line;
    }
    if ($hashcols[2] > 2) {
        print OUT3 $out_line;
    }
    if ($hashcols[2] > 3) {
        print OUT4 $out_line;
    }
    if ($hashcols[2] > 4) {
        print OUT5 $out_line;
    }
}

system("sort -k 2n any_replicate.temp > any_replicate.bed");
system("sort -k 2n two_replicates.temp > two_replicates.bed");
system("sort -k 2n three_replicates.temp > three_replicates.bed");
system("sort -k 2n four_replicates.temp > four_replicates.bed");
system("sort -k 2n weird_replicates.temp > weird_replicates.bed");

system("rm any_replicate.temp");
system("rm two_replicates.temp");
system("rm three_replicates.temp");
system("rm four_replicates.temp");
system("rm weird_replicates.temp");

close(OUT);
close(OUT2);
close(OUT3);
close(OUT4);
close(OUT5);