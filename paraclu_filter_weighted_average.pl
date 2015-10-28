#!/usr/bin/perl

#Filters Paraclu output according to user specifications and reformats as a bed file. Also calculates a weighted average for each cluster to report a single 5' end location in a separate bed file.

#USAGE:
# perl PATH/paraclu_filter.pl /PATH/Paraclu_output_file.bed /PATH/Paraclu_input_file.txt

#TO'G 10/27/15

use warnings;
use strict;

die "USAGE: 'perl <PATH/paraclu_filter.pl> </PATH/paraclu_output_file> </PATH/paraclu_input_file>'" unless @ARGV == 2;

my ($CAGE_file, $paraclu_prep_file) = @ARGV;

print "Enter minimum tags per cluster (e.g. 30): ";
my $min_tags = <STDIN>;
chomp $min_tags;

print "Enter minimum density fold change (e.g. 2): ";
my $min_dens = <STDIN>;
chomp $min_dens;

print "Enter minimum cluster length (e.g. 2): ";
my $min_length = <STDIN>;
chomp $min_length;

print "Enter maximum cluster length (e.g. 200): ";
my $max_length = <STDIN>;
chomp $max_length;

my $length;
my $dens;
my $prev_start = 0;
my $prev_end = 0;
my $prev_chr = "start";
my $chrStart;
my $chrEnd;


open(INF, "<$CAGE_file") or die "couldn't open file";
open(OUT, ">$CAGE_file.filtered.$min_tags.$min_dens.$min_length.$max_length.bed") or die "couldn't open file";

while (my $line = <INF>) {
    chomp($line);
    next if ($line =~ /^#/); #skips the header line
    my @cols = split("\t", $line);
    next if ($cols[5] < $min_tags); #skips clusters with too few tags
    $length = $cols[3] - $cols[2] + 1; #calculates the length of the cluster
    if (($length >= $min_length) and ($length <= $max_length)) {
        $dens = $cols[7] / $cols[6]; #calculates the cluster density fold change
        if ($dens >= $min_dens) {
            next if (($cols[0] eq $prev_chr) and ($cols[2] >= $prev_start) and ($cols[2] <= $prev_end)); #removes subclusters
            $chrStart = $cols[2] - 1; #converts from 1-based (from SAM) to 0-based (for BED)
            $chrEnd = $cols[3] - 1; #converts from 1-based (from SAM) to 0-based (for BED)
            if ($dens < 100) { #makes the ouput more readable
                printf OUT "%s\t%d\t%d\t%d%s%.1f\t%d\t%s\n", $cols[0], $chrStart, $chrEnd, $cols[5], ":", $dens, $cols[5], $cols[1];   #limits the density output to 1 decimal place, but doesn't change huge numbers to exponents
            }
            else {
                printf OUT "%s\t%d\t%d\t%d%s%.1e\t%d\t%s\n", $cols[0], $chrStart, $chrEnd, $cols[5], ":", $dens, $cols[5], $cols[1]; #changes large numbers to exponents
            }
            $prev_start = $cols[2]; #resets variables needed to eliminate subclusters
            $prev_end = $cols[3];
            $prev_chr = $cols[0];
        }
    }
}
close(INF);
close(OUT);

#getting weighted averages of Paraclu peaks:

my $chrStart_CAGE;
my $chrEnd_CAGE;
my $strand_CAGE;
my $CAGE_weighted_sum = 0;
#my $tag_depth = 0;
my $CAGE_weighted_average;

open(INF, "<$CAGE_file.filtered.$min_tags.$min_dens.$min_length.$max_length.bed") or die "couldn't open file";
open(OUT, ">$CAGE_file.filtered.$min_tags.$min_dens.$min_length.$max_length.weighted_average.bed") or die "couldn't open file"; #later, reduce output to just one file containing both the weighted average and the cluster extent

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    $chrStart_CAGE = $cols[1];
    $chrEnd_CAGE = $cols[2];
    $strand_CAGE = $cols[5];
    #$tag_depth = 0;
    open(INF2, "<$paraclu_prep_file") or die "couldn't open file";
    while (my $line2 = <INF2>) {
        chomp($line2);
        my @cols2 = split("\t", $line2);
        if ((($cols2[2]-1) >= $chrStart_CAGE) and (($cols2[2]-1) <= $chrEnd_CAGE) and ($cols2[1] eq $strand_CAGE)) {
            $CAGE_weighted_sum = $CAGE_weighted_sum + ($cols2[2]*$cols2[3]);
        }
    }
    $CAGE_weighted_average = ($CAGE_weighted_sum/$cols[4]) - 1;
    printf OUT "%s\t%1.0f\t%1.0f\t%s%s%s%s%s\t%s\t%s\n", $cols[0], $CAGE_weighted_average, $CAGE_weighted_average, $chrStart_CAGE, ":", $chrEnd_CAGE, ":", $cols[3], $cols[4], $cols[5];
    $CAGE_weighted_sum = 0;
    close(INF2);
}

close(INF);
close(OUT);