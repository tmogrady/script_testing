#!/usr/bin/perl
#bed_length_calculator.pl by TO'G. Determines the length of bed features
use warnings;
use strict;

my ($file) = @ARGV;

open(INF, "<$file") or die "couldn't open file";
open(OUT, ">$file.lengths.txt");

my @blockSizes;
my $length;

while (my $line = <INF>) {
	chomp($line);
	next if $line =~ /^track/;
	my @cols = split("\t", $line);
	if ($cols[10] =~ /^[^,]+$/) {
		$length = $cols[10];
		print OUT "$cols[0]\t$cols[3]\t$length\n";
		$length = 0;
	}
	else {
		while ($cols[10] =~ /(\d+)/g) { 
            push (@blockSizes, $1);
        }
        $length += $_ for @blockSizes;
        print OUT "$cols[0]\t$cols[3]\t$length\n";
        @blockSizes = ();
        $length = 0;
	}
}

close(INF);