# !/usr/bin/perl
# This script takes sam files of TruSeq stranded RNA-seq reads containing putative polyA tails (output from polyA_read_extractor.pl) and converts them into wiggle files of putative polyadenylation sites.
# Usage: perl <path to perl script>

#This works fine for just the EBV genome but has problems when applying to the cellular genome
# TO'G, 12/12/14

use warnings;
use strict;

my ($polyA_sam_file) = @ARGV;

open(INF, "<$polyA_sam_file") or die "couldn't open file";	 	
open(OUT, ">$polyA_sam_file.polyA_sites.temp") or die "couldn't open file";

#create a file with the corresponding to the polyA ends of the reads, and sort it by those coordinates

my $cigar_sum;
my $cigar_calc;
my @plus_ends;
my @read_dist;

while (my $line = <INF>) {
	my @split_line = split("\t", $line);
	if ($split_line[1] == 73 || $split_line[1] == 97 || $split_line[1] == 99 || $split_line[1] == 0) {
		print OUT "$split_line[2]\t$split_line[3]\t0\n";
	}
	elsif ($split_line[1] == 81 || $split_line[1] == 83 || $split_line[1] == 89 || $split_line[1] == 16) {
		while ($split_line[5] =~ /(\d+)[MNX=]/g) {
				push (@read_dist, $1);
		}
		$cigar_sum += $_ for @read_dist;
		$cigar_calc = $split_line[3] + $cigar_sum - 1;
		$cigar_sum = 0;
		@read_dist = ();
		print OUT "$split_line[2]\t$cigar_calc\t1\n";
	}
}
close(INF);
close(OUT);

system("sort -k 1,1 -k 2,2n \Q$polyA_sam_file.polyA_sites.temp\E > \Q$polyA_sam_file.polyA_sites.temp\E.sorted");

#create a wiggle file from the sorted coordinates file

open(INF2, "<$polyA_sam_file.polyA_sites.temp") or die "couldn't open file";	 	
open(OUT2, ">$polyA_sam_file.polyA_sites.temp.wig") or die "couldn't open file";

my $chrom_minus;
my $previous_coordinate_m=0;
my $count_m=0;
my $coordinate_m;
my $chrom_plus;
my $previous_coordinate_p=0;
my $count_p=0;
my $coordinate_p;

while (my $line = <INF2>) {
	
	my @split_line = split("\t", $line);
	
	#reads on the plus strand:	
	if ($split_line[2] == 1) {
		if ($chrom_plus) { #if $chrom_plus has been defined (i.e. there is a previous plus strand read)
			if (($split_line[0] eq $chrom_plus) and ($split_line[1] == $previous_coordinate_p)) {
				$count_p++;
			}
			else {
				$coordinate_p = $previous_coordinate_p - 1;
				print OUT2 "$chrom_plus\t$coordinate_p\t$coordinate_p\t$count_p\n";
				$previous_coordinate_p = $split_line[1];
				$count_p = 1;
			}
		}
		else { #if $chrom_plus has not been defined (i.e. there is no previous plus strand read)
			$chrom_plus = $split_line[0];
			$previous_coordinate_p = $split_line[1];
			$count_p = 1;
		}
	}
	
	#reads on the minus strand:
	elsif ($split_line[2] == 0) {
		if ($chrom_minus) {
			if (($split_line[0] eq $chrom_minus) and ($split_line[1] == $previous_coordinate_m)) {
				$count_m++;
			}
			else {
				$coordinate_m = $previous_coordinate_m - 1;
				print OUT2 "$chrom_minus\t$coordinate_m\t$coordinate_m\t-$count_m\n";
				$chrom_minus = $split_line[0];
				$previous_coordinate_m = $split_line[1];
				$count_m = 1;
			}
		}
		else {
			$chrom_minus = $split_line[0];
			$previous_coordinate_m = $split_line[1];
			$count_m = 1;
		}
	}
}

$coordinate_p = $previous_coordinate_p - 1;
print OUT2 "$chrom_plus\t$coordinate_p\t$coordinate_p\t$count_p\n";
$coordinate_m = $previous_coordinate_m - 1;
print OUT2 "$chrom_minus\t$coordinate_m\t$coordinate_m\t-$count_m\n";

close(INF2);
close(OUT2);

system("sort -k 1,1 -k 2,2n \Q$polyA_sam_file.polyA_sites.temp.wig\E > \Q$polyA_sam_file.polyA_sites.wig\E");
system("rm \Q$polyA_sam_file.polyA_sites.temp.wig\E");
system("rm \Q$polyA_sam_file.polyA_sites.temp\E.sorted");
system("rm \Q$polyA_sam_file.polyA_sites.temp\E");