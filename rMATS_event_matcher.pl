#!/usr/bin/perl

#Returns events that match on chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES and downstreamEE in two rMATS SE output files

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
    $coord_key = "$cols[3]\:$cols[4]\:$cols[5]\:$cols[6]\:$cols[7]\:$cols[8]\:$cols[9]\:$cols[10]\:$cols[11]"; #creates a key for the hash with the fields to match
    $features{$coord_key} = $cols[22]; #adds the key to the hash with the DPSI
}

close(INF);

open(INF, "<$file2");
open(OUT, ">overlap.temp");

print "Comparing to second file\n";

while(my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    foreach my $feature (keys %features) {
        my ($hash_chr, $hash_strand, $hash_start, $hash_end, $hash_upstart, $hash_upend, $hash_downstart, $hash_downend) = split("\:", $feature);
        if ((($cols[3] eq $hash_chr) and ($cols[4] eq $hash_strand)) and
            (($cols[5] == $hash_start) and ($cols[6] == $hash_end)) and
            (($cols[7] == $hash_upstart) and ($cols[8] == $hash_upend)) and
            (($cols[9] == $hash_downstart) and ($cols[10] == $hash_downend))) {
            print OUT "$cols[2]\t$hash_chr\t$hash_strand\t$hash_start\t$hash_end\t$hash_upstart\t$hash_upend\t$hash_downstart\t$hash_downend\t$features{$feature}\t $cols[22]\n";
            last;
        }
        else {
            next;
        }
    }
}

close(INF);
close(OUT);

system("sort -k1,1 -k3,4n overlap.temp > event_overlap.bed");
system("rm overlap.temp");