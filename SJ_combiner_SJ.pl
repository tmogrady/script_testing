#!/usr/bin/perl

#This script takes multiple STAR sj.out files and combines them, removing duplicates and reporting the total uniquely mapped reads for each splice junction in an SJ.tab file.

#USAGE:
# perl PATH/SJ_combiner.pl /PATH/inputfile(s)

# Flemington Lab TO'G 6/22/15

use warnings;
use strict;

my @combined;

foreach my $file(@ARGV) { #concatenates all the files to make one big file
	open(INF, "<$file") or die "couldn't open file";
	
	while (my $line = <INF>) {
		my @line_cols = split("\t", $line);
		#next if $line_cols[0] ne "chrEBV(Akata_107955to171322_1to107954)";
		push (@combined, $line);
	}
	
	close(INF);
}

open(OUT, ">combined_files.txt");

print OUT @combined;

close(OUT);

open(INF2, "<combined_files.txt") or die "couldn't open file";
open(OUT2, ">combined_files.temp");

my %junctions; #creates a hash that will be used to store junctions and their read depths
my $chr_start_end_strand;

while (my $line = <INF2>) {
	chomp($line);
		
	my @cols = split("\t", $line);
	
	next if $cols[3] == 0; #skips the junction if its strand is undefined
	#next if $cols[0] ne "chrEBV(Akata_107955to171322_1to107954)";
	
	if ($cols[3] == 1) { #if the junction is on the plus strand:
		
		my $chr_start_end_strand = "$cols[0]\:$cols[1]\:$cols[2]\:$cols[3]"; #creates a hash key using the chromosome, coordinates and strand
		
		if (exists $junctions{$chr_start_end_strand}) { #if the key is already in the hash, increases the value (count) by the read depth for that junction
			$junctions{$chr_start_end_strand} = $junctions{$chr_start_end_strand} + $cols[6];	
		}
					
		else {
			$junctions{$chr_start_end_strand} = $cols[6]; #if the key is not already in the hash, adds it with a value (count) of the read depth for that junction
				
		}
	}
		if ($cols[3] == 2) { #if the junction is on the minus strand:
		
		my $chr_start_end_strand = "$cols[0]\:$cols[1]\:$cols[2]\:$cols[3]"; #creates a hash key using the chromosome, coordinates and strand
		
		if (exists $junctions{$chr_start_end_strand}) { #if the key is already in the hash, increases the value (count) by the read depth for that junction
			$junctions{$chr_start_end_strand} = $junctions{$chr_start_end_strand} + $cols[6];	
		}
					
		else {
			$junctions{$chr_start_end_strand} = $cols[6]; #if the key is not already in the hash, adds it with a value (count) of the read depth for that junction
				
		}
	}
		
}

foreach my $key (sort keys %junctions) { #prints out a(n inadequately) sorted wiggle file
		my @split_keys = split("\:", $key); 
        print OUT2 "$split_keys[0]\t$split_keys[1]\t$split_keys[2]\t$split_keys[3]\t\.\t0\t$junctions{$key}\t\.\t\.\n";
}		

close(INF2);
close(OUT2);

system("sort -k1,1 -k2,3n combined_files.temp > combined_files_SJ.out.tab");
system("rm combined_files.temp");
system("rm combined_files.txt");
