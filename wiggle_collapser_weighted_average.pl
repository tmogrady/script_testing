# !/usr/bin/perl
# wiggle_collapser.pl: collapses wiggle counts within a certain number of bp of each other into bed features. Each bed feature is 1 bp long and is average of the coordinates within the specified distance, weighted by the wiggle counts at each coordinate.
#input wiggle file must be sorted by coordinate. This version can handle wiggle files with stranded counts.
# Usage: perl <path to perl script> number_of_bases_allowance <path to wiggle file>
# TO'G, 5/18/15

use warnings;
use strict;

die "Wrong number of arguments. Usage: perl <path to perl script>  <min_distance_between_wig_features> <path_to_wiggle_file>\n" unless @ARGV == 2;

my ($distance_between_peaks, $file) = @ARGV;

open(INF, "<$file") or die "couldn't open file";	 	
open(OUT, ">$file.distance_between_$distance_between_peaks.temp.bed") or die "couldn't open file";

die "Not enough arguments. Usage: perl <path to perl script> number_of_bases_allowance <path to wiggle file>\n" unless @ARGV == 2;

my $prev_coord_plus = 1;
my $prev_coord_minus = 1;
my $count_sum_plus = 0;
my $count_sum_minus = 0;
my $weighted_coordinate_sum_plus = 0;
my $weighted_coordinate_sum_minus = 0;
my $weighted_average_plus;
my $weighted_average_minus;
my $first_plus = 1;
my $first_minus = 1;
my @coords_plus;
my @coords_minus;
my $chrStart_plus;
my $chrEnd_plus;
my $chrStart_minus;
my $chrEnd_minus;

while (my $line = <INF>) {
	chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
	my @cols = split("\t", $line);
	if ($cols[3] > 0) { #if this coordinate has a positive count...
		if ($cols[1] < $prev_coord_plus + ($distance_between_peaks+1)) { #if the coordinate is within the specified number of bp of the previous coordinate
			$count_sum_plus = $count_sum_plus + $cols[3]; #adds to the sums to eventually calculate the weighted average
			$weighted_coordinate_sum_plus = $weighted_coordinate_sum_plus + ($cols[1]*$cols[3]);
            push (@coords_plus, $cols[1]);
			$prev_coord_plus = $cols[1]; #sets the current coordinate as the "previous coordinate" before moving on
		}
		else { #if the present coordinate is not within the specified number of bp of the previous coordinate, need to print out a feature
			if ($first_plus == 1) { #"first" flag avoids wonkiness if the first coordinate is far from coordinate 1 (don't need to print out a feature yet)
				$count_sum_plus = $cols[3];
				$weighted_coordinate_sum_plus = $cols[1]*$cols[3];
				$prev_coord_plus = $cols[1];
                push (@coords_plus, $cols[1]);
				$first_plus = 0;
			}
			else {
				$weighted_average_plus = $weighted_coordinate_sum_plus/$count_sum_plus; #calculates weighted average
                $chrStart_plus = $coords_plus[0];
                $chrEnd_plus = pop(@coords_plus);
				printf OUT "%s\t%1.0f\t%1.0f\t%d%s%d%s%d\t%d\t%s\n", "chrEBV(Akata_107955to171322_1to107954)", $weighted_average_plus, $weighted_average_plus, $chrStart_plus, ":", $chrEnd_plus, ":", $count_sum_plus, $count_sum_plus, "+"; #prints out weighted average for plus strand features. Use printf to round the weighted average.
                @coords_plus = ($cols[1]);
				$count_sum_plus = $cols[3]; #sets "previous coordinate", count and sum of counts for the current coordinate
				$weighted_coordinate_sum_plus = $cols[1]*$cols[3];
				$prev_coord_plus = $cols[1];
			}
		}
	}
	elsif ($cols[3] < 0) { #if this coordinate has a negative count...
		if ($cols[1] < $prev_coord_minus + ($distance_between_peaks+1)) { #if the coordinate is within the specified number of bp of the previous coordinate
			$count_sum_minus = $count_sum_minus + $cols[3]; #adds to the sums to eventually calculate the weighted average
			$weighted_coordinate_sum_minus = $weighted_coordinate_sum_minus + ($cols[1]*$cols[3]);
            push (@coords_minus, $cols[1]);
			$prev_coord_minus = $cols[1]; #sets the current coordinate as the "previous coordinate" before moving on
		}
		else { #if the present coordinate is not within the specified number of bp of the previous coordinate, need to print out a feature
			if ($first_minus == 1) { #"first" flag avoids wonkiness if the first coordinate is far from coordinate 1 (don't need to print out a feature yet)
				$count_sum_minus = $cols[3];
				$weighted_coordinate_sum_minus = $cols[1]*$cols[3];
				$prev_coord_minus = $cols[1];
                push (@coords_minus, $cols[1]);
				$first_minus = 0;
			}
			else {
				$weighted_average_minus = $weighted_coordinate_sum_minus/$count_sum_minus; #calculates weighted average
                $chrStart_minus = $coords_minus[0];
                $chrEnd_minus = pop(@coords_minus);
				printf OUT "%s\t%1.0f\t%1.0f\t%d%s%d%s%d\t%d\t%s\n", "chrEBV(Akata_107955to171322_1to107954)", $weighted_average_minus, $weighted_average_minus, $chrStart_minus, ":", $chrEnd_minus, ":", $count_sum_minus, $count_sum_minus, "-"; #prints out weighted average for plus strand features. Use printf to round the weighted average.
                @coords_minus = ($cols[1]);
				$count_sum_minus = $cols[3]; #sets "previous coordinate", count and sum of counts for the current coordinate
				$weighted_coordinate_sum_minus = $cols[1]*$cols[3];
				$prev_coord_minus = $cols[1];
			}
		}
	}
}

if ($count_sum_plus > 0) {#calculates and prints out weighted average for the last feature (plus strand)
	$weighted_average_plus = $weighted_coordinate_sum_plus/$count_sum_plus;
    $chrStart_plus = $coords_plus[0];
    $chrEnd_plus = pop(@coords_plus);
	printf OUT "%s\t%1.0f\t%1.0f\t%d%s%d%s%d\t%d\t%s\n", "chrEBV(Akata_107955to171322_1to107954)", $weighted_average_plus, $weighted_average_plus, $chrStart_plus, ":", $chrEnd_plus, ":", $count_sum_plus,  $count_sum_plus, "+";
}

if ($count_sum_minus < 0) {#calculates and prints out weighted average for the last feature (minus strand)
	$weighted_average_minus = $weighted_coordinate_sum_minus/$count_sum_minus;
    $chrStart_minus = $coords_minus[0];
    $chrEnd_minus = pop(@coords_minus);
	printf OUT "%s\t%1.0f\t%1.0f\t%d%s%d%s%d\t%d\t%s\n", "chrEBV(Akata_107955to171322_1to107954)", $weighted_average_minus, $weighted_average_minus, $chrStart_minus, ":", $chrEnd_minus, ":", $count_sum_minus,  $count_sum_minus, "-";
}

close(INF);
close(OUT);

system("sort -k 1,1 -k 2,2n \Q$file.distance_between_$distance_between_peaks.temp.bed\E > \Q$file.distance_between_$distance_between_peaks.bed\E");
system("rm \Q$file.distance_between_$distance_between_peaks.temp.bed\E");