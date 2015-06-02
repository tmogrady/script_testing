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
my $un_iso_length_sum = 0;
my $un_iso_number = 0;
my $un_read_length_sum = 0;
my $un_read_number = 0;
my $mapped_iso_length_sum = 0;
my $mapped_iso_number = 0;
my $mapped_read_length_sum = 0;
my $mapped_read_number = 0;

while (my $line = <INF>) {
	chomp($line);
	next if $line =~ /^@/;
	my @cols = split("\t", $line);
	my @split_id = split("\/", $cols[0]);
	$iso_length_sum = $iso_length_sum + $split_id[2];
	$iso_number = $iso_number + 1;
	$read_length_sum = $read_length_sum + ($split_id[1]*$split_id[2]);
	$read_number = $read_number + $split_id[1];
	if ($cols[1] == 4) {
		$un_iso_length_sum = $un_iso_length_sum + $split_id[2];
		$un_iso_number = $un_iso_number + 1;
		$un_read_length_sum = $un_read_length_sum + ($split_id[1]*$split_id[2]);
		$un_read_number = $un_read_number + $split_id[1];
	}
	else {
		$mapped_iso_length_sum = $mapped_iso_length_sum + $split_id[2];
		$mapped_iso_number = $mapped_iso_number + 1;
		$mapped_read_length_sum = $mapped_read_length_sum + ($split_id[1]*$split_id[2]);
		$mapped_read_number = $mapped_read_number + $split_id[1];		
	}
}

my $iso_mean;
my $read_mean;
my $un_iso_mean;
my $un_read_mean;
my $mapped_iso_mean;
my $mapped_read_mean;

if ($iso_number == 0) {
	$iso_mean = "NA";
}
else {
	$iso_mean = $iso_length_sum/$iso_number;
}
if ($read_number == 0) {
	$read_mean = "NA";
}
else {
	$read_mean = $read_length_sum/$read_number;
}
if ($un_iso_number == 0) {
	$un_iso_mean = "NA";
}
else {
	$un_iso_mean = $un_iso_length_sum/$un_iso_number;
}
if ($un_read_number == 0) {
	$un_read_mean = "NA";
}
else {
	$un_read_mean = $un_read_length_sum/$un_read_number;
}
if ($mapped_iso_number == 0) {
	$mapped_iso_mean = "NA";
}
else {
	$mapped_iso_mean = $mapped_iso_length_sum/$mapped_iso_number;
}
if ($mapped_read_number == 0) {
	$mapped_read_mean = "NA";
}
else {
	$mapped_read_mean = $mapped_read_length_sum/$mapped_read_number;
}

print "Mean of total isoform lengths: $iso_mean\nMean of total read lengths: $read_mean\nMean of mapped isoform lengths: $mapped_iso_mean\nMean of mapped read lengths: $mapped_read_mean\nMean of unmapped isoform lengths: $un_iso_mean\nMean of total read lengths: $un_read_mean\n";	
	
close(INF);