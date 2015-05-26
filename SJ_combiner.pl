#!/usr/bin/perl

#This script takes multiple STAR sj.out files and combines them, removing duplicates and reporting the total uniquely mapped reads for each splice junction

#USAGE:
# perl PATH/SJ_combiner.pl /PATH/inputfile(s)

# Flemington Lab TO'G 5/26/15

use warnings;
use strict;

my @combined;

foreach my $file(@ARGV) {
	open(INF, "<$file") or die "couldn't open file";
	
	while (my $line = <INF>) {
		push (@combined, $line);
	}
	
	close(INF);
}

open(OUT, ">combined.files.txt");

print OUT @combined;

close(OUT);