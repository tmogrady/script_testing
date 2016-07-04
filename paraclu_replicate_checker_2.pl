#!/usr/bin/perl

#Input: weighted average bed files from paraclu_pipe.pl (for 4 replicates). Returns files consensus files of clusters.

#USAGE:
# perl PATH/paraclu_filter.pl numerical_max_distance /PATH/inputfiles

#TO'G 1/7/2016

use warnings;
use strict;

my ($dist, $file1, $file2) = @ARGV;
my %peaks;

open(INF, "<$file1");
open(OUT, ">replicates.temp");

while(my $line = <INF>) {
    chomp($line);
    print OUT "$line\t1\n";
}

close(INF);

open(INF, "<$file2");
while(my $line = <INF>) {
    chomp($line);
    print OUT "$line\t2\n";
}

close(INF);
close(OUT);

system("sort -k1,1 -k2n replicates.temp > replicates_sorted.temp");
system("rm replicates.temp");

open(INF, "<replicates_sorted.temp");
open(OUT, ">any_replicate.txt");
open(OUT2, ">both_replicates.txt");

my $prev_chrom_plus = 0;
my $prev_coord_plus = 0;
my $peak_sum_plus = 1;
my $coordinate_sum_plus = 1;
my $tag_sum_plus = 1;
my $dens_sum_plus= 1;
my $average_coordinate_plus = 1;
my $average_dens_plus = 1;
my $first_plus = 1;
my $prev_sample_plus = 0;
my $multiple_plus = 0;
my $prev_chrom_minus = 0;
my $prev_coord_minus = 0;
my $peak_sum_minus = 1;
my $coordinate_sum_minus = 1;
my $tag_sum_minus = 1;
my $dens_sum_minus= 1;
my $average_coordinate_minus = 1;
my $average_dens_minus = 1;
my $first_minus = 1;
my $prev_sample_minus = 0;
my $multiple_minus = 0;

while(my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    my @subcols = split("\:", $cols[3]);
    if ($cols[5] eq "+") { #plus strand
        if (($cols[0] eq $prev_chrom_plus) and ($cols[1] <= ($prev_coord_plus + $dist))) { #if the coordinate is within the specified distance of the previous one
            $peak_sum_plus = $peak_sum_plus + 1; #counts number of peaks within the window to use in average calculation
            $coordinate_sum_plus = $coordinate_sum_plus + $cols[1]; #adds the coordinate to the sum used to calculate the average position
            $tag_sum_plus = $tag_sum_plus + $cols[4]; #sums the number of CAGE tags contributing to the site
            $dens_sum_plus = $dens_sum_plus + $subcols[3]; #sums the densities of the peaks within to window to eventually get average
            $prev_coord_plus = $cols[1];
            if ($cols[6] != $prev_sample_plus) { #if this peak is from a different sample than the previous one, flips the switch to note that
                $multiple_plus = 1;
            }
            $prev_sample_plus = $cols[6];
        }
        else {
            if ($first_plus == 1) { #"first" flag avoids wonkiness if the first coordinate is far from coordinate 1 (don't need to print out a feature yet)
                $prev_chrom_plus = $cols[0];
                $prev_coord_plus = $cols[1];
                $peak_sum_plus = 1;
                $coordinate_sum_plus = $cols[1];
                $tag_sum_plus = $cols[4];
                $dens_sum_plus = $subcols[3];
                $prev_sample_plus = $cols[6];
                $first_plus = 0;
            }
            else {
                $average_coordinate_plus = sprintf("%1.0f", ($coordinate_sum_plus/$peak_sum_plus)); #calculates coordinate average
                $average_dens_plus = ($dens_sum_plus/$peak_sum_plus); #calculates average density
                if ($multiple_plus == 1) {
                    print OUT $prev_chrom_plus, "\t", $average_coordinate_plus, "\t", $average_coordinate_plus+1, "\t", $average_dens_plus, ":", $tag_sum_plus, ":2\t", $tag_sum_plus, "\t+\n";
                    print OUT2 $prev_chrom_plus, "\t", $average_coordinate_plus, "\t", $average_coordinate_plus+1, "\t", $average_dens_plus, ":", $tag_sum_plus, ":2\t", $tag_sum_plus, "\t+\n";
                }
                else {
                    print OUT $prev_chrom_plus, "\t", $average_coordinate_plus, "\t", $average_coordinate_plus+1, "\t", $average_dens_plus, ":", $tag_sum_plus, ":1\t", $tag_sum_plus, "\t+\n";
                }
                $prev_chrom_plus = $cols[0];
                $prev_coord_plus = $cols[1];
                $peak_sum_plus = 1;
                $coordinate_sum_plus = $cols[1];
                $tag_sum_plus = $cols[4];
                $dens_sum_plus = $subcols[3];
                $prev_sample_plus = $cols[6];
                $multiple_plus = 0;
            }
        }
    }
    elsif ($cols[5] eq "-"){
        if (($cols[0] eq $prev_chrom_minus) and ($cols[2] <= ($prev_coord_minus + $dist))) {
            $peak_sum_minus = $peak_sum_minus + 1; #counts number of peaks within the window to use in average calculation
            $coordinate_sum_minus = $coordinate_sum_minus + $cols[1]; #adds the coordinate to the sum used to calculate the average position
            $tag_sum_minus = $tag_sum_minus + $cols[4]; #sums the number of CAGE tags contributing to the site
            $dens_sum_minus = $dens_sum_minus + $subcols[3]; #sums the densities of the peaks within to window to eventually get average
            $prev_coord_minus = $cols[2];
            if ($cols[6] != $prev_sample_minus) { #if this peak is from a different sample than the previous one, flips the switch to note that
                $multiple_minus = 1;
            }
            $prev_sample_minus = $cols[6];
        }
        else {
            if ($first_minus == 1) { #"first" flag avoids wonkiness if the first coordinate is far from coordinate 1 (don't need to print out a feature yet)
                $prev_chrom_minus = $cols[0];
                $prev_coord_minus = $cols[2];
                $peak_sum_minus = 1;
                $coordinate_sum_minus = $cols[2];
                $tag_sum_minus = $cols[4];
                $dens_sum_minus = $subcols[3];
                $prev_sample_minus = $cols[6];
                $first_minus = 0;
            }
            else {
                $average_coordinate_minus = sprintf("%1.0f", ($coordinate_sum_minus/$peak_sum_minus)); #calculates coordinate average
                $average_dens_minus = ($dens_sum_minus/$peak_sum_minus); #calculates average density
                if ($multiple_minus == 1) {
                    print OUT $prev_chrom_minus, "\t", $average_coordinate_minus-1, "\t", $average_coordinate_minus, "\t", $average_dens_minus, ":", $tag_sum_minus, ":2\t", $tag_sum_minus, "\t-\n";
                    print OUT2 $prev_chrom_minus, "\t", $average_coordinate_minus-1, "\t", $average_coordinate_minus, "\t", $average_dens_minus, ":", $tag_sum_minus, ":2\t", $tag_sum_minus, "\t-\n";
                }
                else {
                    print OUT $prev_chrom_minus, "\t", $average_coordinate_minus-1, "\t", $average_coordinate_minus, "\t", $average_dens_minus, ":", $tag_sum_minus, ":1\t", $tag_sum_minus, "\t-\n";
                }
                $prev_chrom_minus = $cols[0];
                $prev_coord_minus = $cols[2];
                $peak_sum_minus = 1;
                $coordinate_sum_minus = $cols[2];
                $tag_sum_minus = $cols[4];
                $dens_sum_minus = $subcols[3];
                $prev_sample_minus = $cols[6];
                $multiple_minus = 0;
            }
        }
    }
}

#need to deal with last features too:

if ($multiple_plus == 1) {
    print OUT $prev_chrom_plus, "\t", $average_coordinate_plus, "\t", $average_coordinate_plus+1, "\t", $average_dens_plus, ":", $tag_sum_plus, ":2\t", $tag_sum_plus, "\t+\n";
    print OUT2 $prev_chrom_plus, "\t", $average_coordinate_plus, "\t", $average_coordinate_plus+1, "\t", $average_dens_plus, ":", $tag_sum_plus, ":2\t", $tag_sum_plus, "\t+\n";
}
else {
    print OUT $prev_chrom_plus, "\t", $average_coordinate_plus, "\t", $average_coordinate_plus+1, "\t", $average_dens_plus, ":", $tag_sum_plus, ":1\t", $tag_sum_plus, "\t+\n";
}

if ($multiple_minus == 1) {
    print OUT $prev_chrom_minus, "\t", $average_coordinate_minus-1, "\t", $average_coordinate_minus, "\t", $average_dens_minus, ":", $tag_sum_minus, ":2\t", $tag_sum_minus, "\t-\n";
    print OUT2 $prev_chrom_minus, "\t", $average_coordinate_minus-1, "\t", $average_coordinate_minus, "\t", $average_dens_minus, ":", $tag_sum_minus, ":2\t", $tag_sum_minus, "\t-\n";
}
else {
    print OUT $prev_chrom_minus, "\t", $average_coordinate_minus-1, "\t", $average_coordinate_minus, "\t", $average_dens_minus, ":", $tag_sum_minus, ":1\t", $tag_sum_minus, "\t-\n";
}

close(INF);
close(OUT);
close(OUT2);

system("sort -k1,1 -k2n any_replicate.txt > any_replicate.bed");
system("sort -k1,1 -k2n both_replicates.txt > both_replicates.bed");

system("rm replicates_sorted.temp");
system("rm any_replicate.txt");
system("rm both_replicates.txt");