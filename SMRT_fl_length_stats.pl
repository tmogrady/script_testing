#!/usr/bin/perl
#SMRT_fl_length_stats.pl by TO'G. Given a GMAP-aligned SAM file of SMRT fl.consensus reads, returns information about the distribution of read lengths
use warnings;
use strict;

my ($file) = @ARGV;

open(INF, "<$file") or die "couldn't open file";

my $iso_length_sum = 0;
my $iso_number = 0;
my $read_length_sum = 0;
my $read_number = 0;

while (my $line = <INF>) {
	chomp($line);
	next if $line =~ /^@/;
	my @cols = split("\t", $line);
	my @split_id = split("\/", $cols[0]);
	$iso_length_sum = $iso_length_sum + $split_id[2];
	$iso_number = $iso_number + 1;
	$read_length_sum = $read_length_sum + ($split_id[1]*$split_id[2]);
	$read_number = $read_number + $split_id[1];
}
my $iso_mean = $iso_length_sum/$iso_number;
my $read_mean = $read_length_sum/$read_number;

print "Mean of isoform lengths: $iso_mean\nMean of read lengths: $read_mean\n";	
	
close(INF);