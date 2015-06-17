# !/usr/bin/perl
# This program takes SAM files generated from Pacbio alignments and makes a corresponding .bed annotation for each read. The only issue that it doesn't deal with is where the thick lines should start. This must be added manually.


# USAGE: 
# perl /PATH/pacbio_sam_to_bed_converter.pl /PATH/file_directory/*
	
	
# EKF and T'OG (Flemington lab) 05/12/2015


use warnings;
use strict;

foreach my $file(@ARGV) {
	open(INF, "<$file") or die "couldn't open input file";
	open(OUT, ">$file.bed") or die "couldn't open output introns individual lines file";

	
	while (my $line = <INF>) {
		$line =~ s/\r//g;
		chomp($line);
		
		next if ($line =~ m/\@/); # this is to get rid of header lines in sam files. 
		
		my @split_line = split("\t", $line);
		
		my $strand;
		my $chr = $split_line[2];
		my $chr_start = $split_line[3] - 1;
		my $chr_end = 0;
		my $feature_name = $split_line[0];
		my $score = 1000;
		my $color = "133,0,33";
		
		
		# identify strand
		if ($split_line[1] == 0) {
			$strand = "+";
		}
		elsif ($split_line[1] == 16) {
			$strand = "-";
		}
		
		# This is an interesting way to split. If you put the split marker in parenthesis, the marker will be included in the array. Importantly, however, what is in between each of these array elements (which will be an empty value), is also put into the array. So this has to be dealt with...see bottom portion of this section.
		my @split_CIGAR_temp = split(/(\d+\D)/, $split_line[5]);
		
		my @split_CIGAR;
		
		foreach my $temporary(@split_CIGAR_temp) {
			if ($temporary =~ m/\d+\D/) {
				push(@split_CIGAR, $temporary);
			}
		}
		
		
		
		my $exon_sum = 0;
		
		my @exon_lengths = ();
		my @block_starts = (0);
		
		my $count = 0;
		
		foreach my $split_CIGAR(@split_CIGAR) {
			
			$count++;
			
			# if the first array element is soft clipping, then we need to ignore it.
			if (($count == 1) && (my ($five_prime_clipped_bases) = $split_CIGAR =~ m/(\d+)S$/)) {
				
			}
			
			# if the last array element is soft clipping, then we need to ignore it.
			elsif (($count > 1) && (my ($three_prime_clipped_bases) = $split_CIGAR =~ m/(\d+)S$/)) {
				
			}
				
			
			# if the array element is an intron (indicated by a number followed by an "N"), then we want to 1. stop the exon value summation and put the last value into the exon_lengths array, 2. grab the value of the intron length and add it to the last exon_length plus the last block_start array values to give the new_block_start value (note that you can specify the last element of an array by indicating the array[-1]) and push the new_block_start value into the block_starts array. Lastly, we need to re-set the exon_sum back to 0 so we can start adding up for next exon.
			elsif ($split_CIGAR =~ m/N$/) {
				
				push(@exon_lengths, $exon_sum);
				
				my ($intron_length) = $split_CIGAR =~ m/(\d+)/;
				
				my $new_block_start = $exon_lengths[-1] + $block_starts[-1] + $intron_length;
				
				push(@block_starts, $new_block_start);

				$exon_sum = 0;
			}
			
			# First we grab the numerical value from the element. If the value is associated with an insertion, we don't want to add it to the exon_sum. If it is not an insertion, then we want to add it to the exon_sum and keep counting the length until we hit an intron (previous elsif statement).
			else {
				
				my ($value) = $split_CIGAR =~ m/(\d+)/;
				
				if ($split_CIGAR =~ m/I$/) {
					
					$exon_sum = $exon_sum - 0;
				}
				
				else {
					$exon_sum = $exon_sum + $value;
				}
				

			}
			
			
			
		}
		
		# Once we're done looping through an array, there will be no intron to stop the loop. So, we just push the last exon_sum into the exon_lengths array, we determine the end of the entire gene structure (i.e. chr_end), and we count the number of exons.
		push(@exon_lengths, $exon_sum);
		
		$chr_end = $chr_start + $block_starts[-1] + $exon_lengths[-1];
		
		my $exon_number = @exon_lengths;
		
		# Here we just print out the accumulated data into a bed formatted line.
		print OUT $chr, "\t", $chr_start, "\t", $chr_end, "\t", $feature_name, "\t", $score, "\t", $strand, "\t", $chr_start, "\t", $chr_end, "\t", $color, "\t", $exon_number, "\t", join("\,", @exon_lengths), "\t", join("\,", @block_starts), "\n";
		
	}
	close(INF);
	close(OUT);
}

