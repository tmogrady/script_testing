#!/usr/bin/perl

#Filters Paraclu output according to user specifications and reformats as a bed file

#USAGE:
# perl PATH/paraclu_filter.pl /PATH/inputfile(s)

#TO'G 9/9/15

use warnings;
use strict;

die "USAGE: 'perl <PATH/paraclu_filter.pl> </PATH/paraclu_output_file>'" unless @ARGV == 1;

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

foreach my $file(@ARGV) {
    open(INF, "<$file") or die "couldn't open file";
    open(OUT, ">$file.filtered.$min_tags.$min_dens.$min_length.$max_length.bed") or die "couldn't open file";
    
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
                $chrEnd = $cols[3]; #chrEnd stays 1-based (BED is 0-based half-open)
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
}