# !/usr/bin/perl
# This script takes sam files of TruSeq stranded RNA-seq reads containing putative polyA tails (output from polyA_read_extractor.pl) and converts them into bed files of putative polyadenylation sites.
# Usage: perl <path to perl script> number_of_bases_allowance <path to sam file>
# EKF and TO'G, 12/12/14

use warnings;
use strict;
use List::Util qw(min max);

my ($distance_between_reads, $polyA_sam_file) = @ARGV;

open(INF, "<$polyA_sam_file") or die "couldn't open file";	 	
open(OUT, ">$polyA_sam_file.distance_between_$distance_between_reads.temp.bed") or die "couldn't open file";	

my $chrom_plus;
my $end_coordinate_plus = 0;
my $first_coordinate_plus;
my $second_coordinate_plus;	
my $cigar_calc;
my @plus_ends;
my @read_dist;
my $cigar_sum;
my $chrom_minus;
my $end_coordinate_minus = 0;	
my @minus_ends;

while (my $line = <INF>) {
	my @split_line = split("\t", $line);
#reads on the positive strand:
	if ($split_line[1] == 81 || $split_line[1] == 83 || $split_line[1] == 89 || $split_line[1] == 16) { 
		if ($chrom_plus) { #checks to see if $chrom_plus has been defined (i.e. if there is a previous plus read)
			while ($split_line[5] =~ /(\d+)[MNX=]/g) {
				push (@read_dist, $1);
			}
			$cigar_sum += $_ for @read_dist;
			$cigar_calc = $split_line[3] + $cigar_sum - 1;
			$cigar_sum = 0;
			@read_dist = ();	
			if (($split_line[2] eq $chrom_plus) and (abs($cigar_calc - $end_coordinate_plus) < ($distance_between_reads))){ #checks to see if the current plus strand read matches the previous plus strand read, within the allowance
				push (@plus_ends, $cigar_calc);
				next;
			}
			else {
				$first_coordinate_plus = min(@plus_ends);
				$second_coordinate_plus = max(@plus_ends);
				print OUT $chrom_plus, "\t", $first_coordinate_plus, "\t", $second_coordinate_plus, "\t", scalar @plus_ends, "\t", scalar @plus_ends, "\t\+\n";
				$chrom_plus = $split_line[2];
				@plus_ends = ($cigar_calc);
				$end_coordinate_plus = $cigar_calc;
				next;
			}
		}
		else { #if $chrom_plus has not been defined (i.e. there is no previous plus strand read)
			while ($split_line[5] =~ /(\d+)[MNX=]/g) {
				push (@read_dist, $1);
			}
			$cigar_sum += $_ for @read_dist;
			$cigar_calc = $split_line[3] + $cigar_sum - 1;
			$cigar_sum = 0;
			@read_dist = ();	
			$chrom_plus = $split_line[2];
			$end_coordinate_plus = $cigar_calc;
			@plus_ends = ($cigar_calc);
			next;
		}
	}

#reads on the negative strand
	elsif ($split_line[1] == 73 || $split_line[1] == 97 || $split_line[1] == 99 || $split_line[1] == 0) {
		if ($chrom_minus) {
			if (($split_line[2] eq $chrom_minus) and (($split_line[3] - $end_coordinate_minus) < ($distance_between_reads))){
				push (@minus_ends, $split_line[3]);
				next;
			}
			else {
				my $first_coordinate_minus = min(@minus_ends);
				my $second_coordinate_minus = max(@minus_ends);
				print OUT $chrom_minus, "\t", $first_coordinate_minus, "\t", $second_coordinate_minus, "\t", scalar @minus_ends, "\t", scalar @minus_ends, "\t-\n";
				$chrom_minus = $split_line[2];
				@minus_ends = ($split_line[3]);
				$end_coordinate_minus = $split_line[3];
				next;
			}
		}
		else {
			$chrom_minus = $split_line[2];
			$end_coordinate_minus = $split_line[3];
			@minus_ends = ($split_line[3]);
			next;
		}
	} 
}

#below to output the last feature on each strand
$first_coordinate_plus = min(@plus_ends);
$second_coordinate_plus = max(@plus_ends);
print OUT $chrom_plus, "\t", $first_coordinate_plus, "\t", $second_coordinate_plus, "\t", scalar @plus_ends, "\t", scalar @plus_ends, "\t\+\n";
my $first_coordinate_minus = min(@minus_ends);
my $second_coordinate_minus = max(@minus_ends);
print OUT $chrom_minus, "\t", $first_coordinate_minus, "\t", $second_coordinate_minus, "\t", scalar @minus_ends, "\t", scalar @minus_ends, "\t-\n";

close(INF);
close(OUT);

system("sort -k 1,1 -k 2,2n \Q$polyA_sam_file.distance_between_$distance_between_reads.temp.bed\E > \Q$polyA_sam_file.distance_between_$distance_between_reads.bed\E");
system("rm \Q$polyA_sam_file.distance_between_$distance_between_reads.temp.bed\E");