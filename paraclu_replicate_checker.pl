#!/usr/bin/perl

#Input: weighted average bed files from paraclu_filter_weighted_average.pl (for 4 replicates). Returns files consensus files of clusters.

#USAGE:
# perl PATH/paraclu_filter.pl numerical_max_distance /PATH/inputfiles

#TO'G 10/28/15

use warnings;
use strict;

my ($dist, $file1, $file2, $file3, $file4) = @ARGV;
my %peaks;


open(INF, "<$file1");

print "Processing first file\n";

while(my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    my @subcols = split("\:", $cols[3]);
    my $value = "$cols[4]\:$subcols[3]\:1"; #creates a variable including the tag depth, relative density for the cluster and a replicate count
    $peaks{$cols[1]} = $value; #adds the coordinate to the hash as a key, with the value assigned above as a value
}

#foreach (sort keys %peaks) {
#    print "$_: $peaks{$_}\n";
#}

close(INF);

open(INF, "<$file2");

print "Processing second file\n";

while(my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    my @subcols = split("\:", $cols[3]);
    my $lower_limit = $cols[1] - $dist;
    my $upper_limit = $cols[1] + $dist;
    my $found_flag = 0;
    foreach my $coord (keys %peaks) {
        if (($coord >= $lower_limit) and ($coord <= $upper_limit)) { #checks to see if the cluster is already in the hash (from an earlier file)
            my @hashcols = split("\:", $peaks{$coord}); #extracts information about the cluster from the hash
            my $depth = $cols[4] + $hashcols[0]; #adds tag depths from the current file and the hash
            my $reps = $hashcols[2]+1; #adds one to the replicate count
            my $density = (($subcols[3] + ($hashcols[1]*$hashcols[2]))/($hashcols[2]+1));
            my $new_coord = (($coord*$hashcols[0]) + ($cols[1]*$cols[4]))/($hashcols[0] + $cols[4]);
            my $value = "$depth\:$density\:$reps";
            delete $peaks{$coord};
            $peaks{$new_coord} = $value;
            $found_flag = 1;
            last;
        }
    }
    if ($found_flag == 0) {
        my $value = "$cols[4]\:$subcols[3]\:1"; #creates a variable including the tag depth, relative density for the cluster and a replicate count
        $peaks{$cols[1]} = $value; #adds the coordinate to the hash as a key, with the value assigned above as a value
    }
    
}

open(OUT, ">test_comparison.bed");

foreach (sort keys %peaks) {
    print OUT "$_: $peaks{$_}\n";
}

close(OUT);