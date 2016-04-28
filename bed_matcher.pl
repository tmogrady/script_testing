#!/usr/bin/perl

#Returns bed features form two files that match on chr, chrStart, chrEnd and strand

#USAGE:
# perl PATH/bed_matcher.pl /PATH/inputfiles

use warnings;
use strict;

my ($file1, $file2) = @ARGV;
my %features;

open(INF, "<$file1");

print "Processing first file\n";

my $coord_key;
my $values;

while(my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    $coord_key = "$cols[0]\:$cols[1]\:$cols[2]\:$cols[5]"; #creates a key for the hash with chr, chrStart, chrEnd and strand
    $features{$coord_key} = 1;
}

close(INF);

open(INF, "<$file2");
open(OUT, ">overlap.temp");

print "Comparing to second file\n";

while(my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    foreach my $feature (keys %features) {
        my ($hash_chr, $hash_chrStart, $hash_chrEnd, $hash_strand) = split("\:", $feature);
        if ((($cols[0] eq $hash_chr) and ($cols[1] == $hash_chrStart)) and (($cols[2] == $hash_chrEnd) and ($cols[5] eq $hash_strand))) {
            print OUT $line, "\n";
            last;
        }
        else {
            next;
        }
    }
}

close(INF);
close(OUT);

system("sort -k1,1 -k2,3n overlap.temp > overlap.bed");