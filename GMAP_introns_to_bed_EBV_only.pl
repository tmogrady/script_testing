#!/usr/bin/perl
#GMAP_introns_to_bed.pl by TO'G. Takes junctions files from GMAP (generated with the -f introns argument) and converts them to intron bed files.
use warnings;
use strict;

foreach my $file (@ARGV) {
	open(INF, "<$file" ) or die "couldn't open file";
	open(OUT, ">$file.temp");

	my $line;
	while( $line = <INF> ) {
		chomp($line);
			
		my ($id) = $line =~ /\>(.+)\.i/;
		my ($chr) = $line =~ /\s(.+):/;
		my ($score) = $line =~ /\>.+\/(\d+)\//;
		my ($donor, $acceptor) = $line =~ /:(\d+)\.\.(\d+)/;
		next if $chr ne "chrEBV(Akata_107955to171322_1to107954)";
		if ($acceptor > $donor) {
			print OUT $chr, "\t", $donor, "\t", $acceptor - 1, "\t", $id, "\t", $score, "\t+\n";
		}
		else {
			print OUT $chr, "\t", $acceptor, "\t", $donor - 1, "\t", $id, "\t", $score, "\t-\n";
		}
	}
	close(OUT);
	close(INF);
		
	system("sort -k2,3n \Q$file\E.temp > \Q$file\E.sorted.temp"); #sorts so that duplicate introns will be next to each other.
		
	open(INF2, "<$file.sorted.temp" ) or die "couldn't reopen file";
	open(OUT2, ">$file.EBV_only.bed.temp");
		
	my $plus_previous_chr = "start";
	my $plus_count = 0;
	my $plus_previous_start = 0;
	my $plus_previous_end = 0;
	my $minus_previous_chr = "start";
	my $minus_count = 0;
	my $minus_previous_start = 0;
	my $minus_previous_end = 0;

	while (my $line = <INF2>) {
		chomp($line);
		my @cols = split("\t", $line);
		if ($cols[5] eq "+") { #plus and minus need to be treated separately in case of introns with the same starts and ends annotated on opposite strands
			if ($cols[0] eq $plus_previous_chr && $cols[1] == $plus_previous_start && $cols[2] == $plus_previous_end) { #checks to see if the intron matches the previous intron
				$plus_count = $plus_count + $cols[4];
				$plus_previous_chr = $cols[0];
				$plus_previous_start = $cols[1];
				$plus_previous_end = $cols[2];
			}
			else {
				if ($plus_previous_chr eq "start") { #prevents the initial placeholder value from printing out as a line, and sets the values of the first intron
					$plus_previous_chr = $cols[0];
					$plus_previous_start = $cols[1];
					$plus_previous_end = $cols[2];
					$plus_count = $cols[4];
				}
				else {
					print OUT2 "$plus_previous_chr\t$plus_previous_start\t$plus_previous_end\t$plus_count\t$plus_count\t+\n";
					$plus_count = $cols[4];
					$plus_previous_chr = $cols[0];
					$plus_previous_start = $cols[1];
					$plus_previous_end = $cols[2];
				}
			}
		}
		if ($cols[5] eq "-") {
			if ($cols[0] eq $minus_previous_chr && $cols[1] == $minus_previous_start && $cols[2] == $minus_previous_end) {
				$minus_count = $minus_count + $cols[4];
				$minus_previous_chr = $cols[0];
				$minus_previous_start = $cols[1];
				$minus_previous_end = $cols[2];
			}
			else {
				if ($minus_previous_chr eq "start") {
					$minus_previous_chr = $cols[0];
					$minus_previous_start = $cols[1];
					$minus_previous_end = $cols[2];
					$minus_count = $cols[4];
				}
				else {
					print OUT2 "$minus_previous_chr\t$minus_previous_start\t$minus_previous_end\t$minus_count\t$minus_count\t-\n"; #prints out in bed format
					$minus_count = $cols[4];
					$minus_previous_chr = $cols[0];
					$minus_previous_start = $cols[1];
					$minus_previous_end = $cols[2];
				}
			}
		}
	}
		
	print OUT2 "$plus_previous_chr\t$plus_previous_start\t$plus_previous_end\t$plus_count\t$plus_count\t+\n"; #adds the last plus strand feature
	print OUT2 "$minus_previous_chr\t$minus_previous_start\t$minus_previous_end\t$minus_count\t$minus_count\t-\n"; #adds the last plus strand feature	
	close(OUT2);
	close(INF2);
		
	system("sort -k2,3n \Q$file\E.EBV_only.bed.temp > \Q$file\E.EBV_only.bed");
	system("rm \Q$file\E.temp");
	system("rm \Q$file\E.sorted.temp");
	system("rm \Q$file\E.EBV_only.bed.temp");
}