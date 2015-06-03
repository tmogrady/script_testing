#!/usr/bin/perl
#bed_matcher.pl by TO'G. Takes two bed files and returns the features that match on chromosome, start coordinate, end coordinate and strand.
#3/5/2015
use warnings;
use strict;

my ($file1, $file2) = (@ARGV);
open(INF, "<$file1" ) or die "couldn't open file";
open(INF2, "<$file2" ) or die "couldn't open file";
open(OUT, ">$file1.overlap.bed");
open(OUT2, ">$file2.nooverlap.bed");
open(OUT3, ">$file1.nooverlap.bed");

my $line;
my %features;
	
while( $line = <INF> ) {
	chomp($line);
	my @cols = split("\t", $line);
	my $key_combo = "$cols[0]$cols[1]$cols[2]$cols[5]"; #for each line in file 1, creates a key for the hash combining chromosome, start coordinate, end coordinate and strand
	$features{$key_combo} = $cols[4]; #enters a count value for the key into the hash
}

close(INF);

my $line2;
my %features2;

while( $line2 = <INF2>) {
	chomp($line2);
	my @cols = split("\t", $line2);
	my $key_combo2 = "$cols[0]$cols[1]$cols[2]$cols[5]"; #for each line in file 2, creates a variable/key combining chromosome, start coordinate, end coordinate and strand
	$features2{$key_combo2} = $cols[4]; #enters a count value for the key into the hash
	if (exists $features{$key_combo2}) { #checks to see if the key exists in the hash from file 1: if so, prints it out in the overlap file
		print OUT $cols[0], "\t", $cols[1], "\t", $cols[2], "\t$cols[4]Il_$features{$key_combo2}Pb\t", $cols[4] + $features{$key_combo2}, "\t", $cols[5], "\n";
	}
	else {
		print OUT2 $line2, "\n"; #if the key from file 2 does not exist in the hash for file 1, prints the feature out in the file2.nooverlap file
	}
}

open(INF3, "<$file1" ) or die "couldn't reopen file";

my $line3;

while( $line3 = <INF3>) {
	chomp($line3);
	my @cols = split("\t", $line3);
	my $key_combo3 = "$cols[0]$cols[1]$cols[2]$cols[5]"; #for each line in file 1, creates a key for the hash combining chromosome, start coordinate, end coordinate and strand
	if (exists $features2{$key_combo3}) { #checks to see if the key exists in the hash from file 2: if so, skips (because it is already in file1.overlap)
		next;
	}
	else {
		print OUT3 $line3, "\n"; #if the key from file 1 does not exist in the hash for file 2, prints the feature out in the file1.nooverlap file
	}
}
	
#foreach my $key_combo (keys %features) {
#	print "$key_combo $features{$key_combo}\n";
#}

close(INF3);
close(OUT3);		
close(OUT2);
close(OUT);
close(INF2);