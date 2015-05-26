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
			if ($acceptor > $donor) {
				print OUT $chr, "\t", $donor, "\t", $acceptor - 1, "\t", $id, "\t", $score, "\t+\n";
			}
			else {
				print OUT $chr, "\t", $acceptor, "\t", $donor - 1, "\t", $id, "\t", $score, "\t-\n";
			}
		}
		close(OUT);
		close(INF);
		
		system("sort -k2,3n \Q$file\E.temp > \Q$file\E.sorted.temp"); #sorts so that duplicate introns will be next to each other
		
		open(INF2, "<$file.sorted.temp" ) or die "couldn't reopen file";
		open(OUT2, ">$file.bed");
		
		my $count = 0;
		my $previous_chr = "start";
		my $previous_donor = 0;
		my $previous_acceptor = 0;
		my $previous_strand = "+";

		while (my $line = <INF2>) {
			chomp($line);
			my @cols = split("\t", $line);
			if ($cols[0] eq $previous_chr && $cols[1] == $previous_donor && $cols[2] == $previous_acceptor && $cols[5] eq $previous_strand) { #checks to see if the intron matches the previous intron
				$count = $count + $cols[4];
				$previous_chr = $cols[0];
				$previous_donor = $cols[1];
				$previous_acceptor = $cols[2];
				$previous_strand = $cols[5];
			}
			else {
				if ($previous_chr eq "start") { #prevents the initial placeholder value from printing out as a line, and sets the values of the first intron
					$previous_chr = $cols[0];
					$previous_donor = $cols[1];
					$previous_acceptor = $cols[2];
					$previous_strand = $cols[5];
					$count = $cols[4];
				}
				else {				
					print OUT2 "$previous_chr\t$previous_donor\t$previous_acceptor\t$count\t$count\t$previous_strand\n"; #prints out in bed format
					$count = $cols[4];
					$previous_chr = $cols[0];
					$previous_donor = $cols[1];
					$previous_acceptor = $cols[2];
					$previous_strand = $cols[5];
				}
			}
		}
		
		print OUT2 "$previous_chr\t$previous_donor\t$previous_acceptor\t$count\t$count\t$previous_strand\n"; #adds the last feature
				
		close(OUT2);
		close(INF2);
		
	system("rm \Q$file\E.temp");
	system("rm \Q$file\E.sorted.temp");
	}