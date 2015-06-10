#!/usr/bin/perl

#Counts the number of non-clipped reads with 3' ends at each genomic position and estimates consensus locations of clusters of 3 ends. Output is a wiggle file of all 3' ends and a bed file of the weighted centers of end clusters. 

#Written to work with PacBio SMRT reads from full-length consensus files, where the read name is in the format putative_isoform_id/number_of_reads/length.

#USAGE:
# perl <PATH/read_end_finder.pl> </PATH/input_sam_file> <maximum_distance_between_ends(optional)>

#TO'G 6/10/2015

use warnings;
use strict;

my $file = shift;
my $distance_between_peaks = shift // 8;
	
system("sort -k 3,3 -k 4,4n \Q$file\E > \Q$file\E.sorted.temp");
system("awk '\$2==0' \Q$file\E.sorted.temp > \Q$file\E.sorted.plus.sam.temp");
system("awk '\$2==16' \Q$file\E.sorted.temp > \Q$file\E.sorted.minus.sam.temp");
system("rm \Q$file\E.sorted.temp");

#processing of PLUS sam file
open(INF, "<$file.sorted.plus.sam.temp") or die "couldn't open file";
open(OUT1, ">$file.sorted.plus.sam.soft_clipped_reads.sam.temp") or die "couldn't open file";
open(OUT2, ">$file.sorted.plus.sam.non_clipped_reads.sam.temp") or die "couldn't open file";
open(OUT3, ">$file.sorted.plus.sam.read_ends.wig.temp") or die "couldn't open file";

my @dist;
my $sum;
my %plus_ends;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    if ($cols[5] =~ m/\d+S$/) {	#removes reads soft-clipped at the 3' end
        print OUT1 $line, "\n";
    }
    else {
        print OUT2 $line, "\n";
        while ($cols[5] =~ /(\d+)[DMNX=]/g) { #these lines use the CIGAR string to determine the downstream coordinate
                push (@dist, $1);
        }
        $sum += $_ for @dist;
        my $end_coord = $cols[3] + $sum - 1;
        my $chr_end_coord = "$cols[2]\:$end_coord"; #combines the chromosome and 3' end coordinate into a key to use for the hash
        $sum = 0;
        @dist = ();
        my @split_id = split("\/", $cols[0]); #extracts the read depth for this putative isoform from its id			
        
        if (exists $plus_ends{$chr_end_coord}) { #if the key is already in the hash, increases the value (count) by 1
            $plus_ends{$chr_end_coord} = $plus_ends{$chr_end_coord} + $split_id[1];	
        }
        else {
            $plus_ends{$chr_end_coord} = $split_id[1]; #if the key is not already in the hash, adds it with a value (count) of the read depth
        }
    }
}

foreach my $chr_end_coord (sort keys %plus_ends) { #prints out a(n inadequately) sorted wiggle file
    my @split_keys = split("\:", $chr_end_coord); 
    print OUT3 "$split_keys[0]\t$split_keys[1]\t$split_keys[1]\t$plus_ends{$chr_end_coord}\n";
}	
close(INF);
close(OUT1);
close(OUT2);
close(OUT3);

system("sort -k 1,1 -k 2,2n \Q$file\E.sorted.plus.sam.read_ends.wig.temp > \Q$file\E.sorted.plus.sam.read_ends.wig"); #properly sorts the wiggle file
system("rm \Q$file\E.sorted.plus.sam.read_ends.wig.temp");
system("rm \Q$file\E.sorted.plus.sam.soft_clipped_reads.sam.temp");
system("rm \Q$file\E.sorted.plus.sam.non_clipped_reads.sam.temp");
system("rm \Q$file\E.sorted.plus.sam.temp");

#processing of MINUS sam file
open(INF, "<$file.sorted.minus.sam.temp") or die "couldn't open file";
open(OUT1, ">$file.sorted.minus.sam.soft_clipped_reads.sam.temp") or die "couldn't open file";
open(OUT2, ">$file.sorted.minus.sam.non_clipped_reads.sam.temp") or die "couldn't open file";
open(OUT3, ">$file.sorted.minus.sam.read_ends.wig") or die "couldn't open file";

my $previous_coordinate=1;
my $count=0;
my $previous_chr = "start";
while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    if ($cols[5] =~ m/^\d+S/) { #removes reads clipped at the 3' end
        print OUT1 $line, "\n";
    }
    else {
        print OUT2 $line, "\n";
        #uses Erik's counting method that should take less memory than the hash method, but only works for coordinates that can be obtained directly from the sorted sam file without having to calculate from the CIGAR string
        my @split_id = split("\/", $cols[0]); #extracts the read depth for this putative isoform from its id
        if ($cols[2] eq $previous_chr && $cols[3] == $previous_coordinate) { 
            $count = $count + $split_id[1]; #increases the count by the read depth for the putative isoform		
        }
        else {
            if ($previous_chr eq "start") { #doesn't print out the placeholder first line.
                
                $previous_chr = $cols[2];	#sets the previous chromosome, previous coordinate and count values			
                $previous_coordinate = $cols[3];				
                $count = $split_id[1];
            }
            else {
                print OUT3 $previous_chr, "\t", $previous_coordinate, "\t", $previous_coordinate, "\t-", $count, "\n"; #prints to output file
                $previous_chr = $cols[2];				
                $previous_coordinate = $cols[3];				
                $count = $split_id[1];
            }
        }
    }
}	
print OUT3 $previous_chr, "\t", $previous_coordinate, "\t", $previous_coordinate, "\t-", $count, "\n"; #prints the last start coordinates to output file
close(INF);
close(OUT1);
close(OUT2);
close(OUT3);

system("cat \Q$file\E.sorted.plus.sam.read_ends.wig \Q$file\E.sorted.minus.sam.read_ends.wig | sort -k2,3n > \Q$file\E.all_read_ends.wig.temp");
system("rm \Q$file\E.sorted.minus.sam.soft_clipped_reads.sam.temp");
system("rm \Q$file\E.sorted.minus.sam.non_clipped_reads.sam.temp");
system("rm \Q$file\E.sorted.minus.sam.temp");

open(INF, "<$file.all_read_ends.wig.temp") or die "couldn't open file";
open(OUT, ">$file.distance_between_$distance_between_peaks.temp.bed") or die "couldn't open file";

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
	printf OUT "%s\t%1.0f\t%1.0f\t%d%s%d%s%d\t%d\t%s\n", "chrEBV(Akata_107955to171322_1to107954)", $weighted_average_plus, $weighted_average_plus, $chrStart_plus, ":", $chrEnd_plus, ":", $count_sum_plus, $count_sum_plus, "+";
}

if ($count_sum_minus < 0) {#calculates and prints out weighted average for the last feature (minus strand)
	$weighted_average_minus = $weighted_coordinate_sum_minus/$count_sum_minus;
    $chrStart_minus = $coords_minus[0];
    $chrEnd_minus = pop(@coords_minus);
	printf OUT "%s\t%1.0f\t%1.0f\t%d%s%d%s%d\t%d\t%s\n", "chrEBV(Akata_107955to171322_1to107954)", $weighted_average_minus, $weighted_average_minus, $chrStart_minus, ":", $chrEnd_minus, ":", $count_sum_minus, $count_sum_minus, "-";
}

close(INF);
close(OUT);

system("sort -k 1,1 -k 2,2n \Q$file.distance_between_$distance_between_peaks.temp.bed\E > \Q$file.distance_between_$distance_between_peaks.bed\E");
system("rm \Q$file.distance_between_$distance_between_peaks.temp.bed\E");
system("rm \Q$file\E.all_read_ends.wig.temp");
