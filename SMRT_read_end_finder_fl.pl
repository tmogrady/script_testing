#!/usr/bin/perl

#This script will count the number of non-clipped reads with 3' ends at each genomic position.

#Written to work with PacBio SMRT reads from full-length consensus files, where the read name is in the format putative_isoform_id/number_of_reads/length

#USAGE:
# perl PATH/read_end_finder.pl /PATH/inputfile(s)

#TO'G 3/10/2015

use warnings;
use strict;

	
foreach my $file(@ARGV) {
	
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
	
	system("rm \Q$file\E.sorted.minus.sam.soft_clipped_reads.sam.temp");
	system("rm \Q$file\E.sorted.minus.sam.non_clipped_reads.sam.temp");
	system("rm \Q$file\E.sorted.minus.sam.temp");
}