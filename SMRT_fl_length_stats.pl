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
my $EBV_iso_length_sum = 0;
my $EBV_iso_number = 0;
my $EBV_read_length_sum = 0;
my $EBV_read_number = 0;
my $cell_iso_length_sum = 0;
my $cell_iso_number = 0;
my $cell_read_length_sum = 0;
my $cell_read_number = 0;

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
		if ($cols[2] eq "chrEBV(Akata_107955to171322_1to107954)") {
			$EBV_iso_length_sum = $EBV_iso_length_sum + $split_id[2];
			$EBV_iso_number = $EBV_iso_number + 1;
			$EBV_read_length_sum = $EBV_read_length_sum + ($split_id[1]*$split_id[2]);
			$EBV_read_number = $EBV_read_number + $split_id[1];
		}
		else {
			$cell_iso_length_sum = $cell_iso_length_sum + $split_id[2];
			$cell_iso_number = $cell_iso_number + 1;
			$cell_read_length_sum = $cell_read_length_sum + ($split_id[1]*$split_id[2]);
			$cell_read_number = $cell_read_number + $split_id[1];
		}
	}
}

my $iso_mean;
my $read_mean;
my $un_iso_mean;
my $un_read_mean;
my $mapped_iso_mean;
my $mapped_read_mean;
my $EBV_iso_mean;
my $EBV_read_mean;
my $cell_iso_mean;
my $cell_read_mean;

#-------------------
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
#-------------------
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
#--------------------
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
#--------------------
if ($EBV_iso_number == 0) {
	$EBV_iso_mean = "NA";
}
else {
	$EBV_iso_mean = $EBV_iso_length_sum/$EBV_iso_number;
}
if ($EBV_read_number == 0) {
	$EBV_read_mean = "NA";
}
else {
	$EBV_read_mean = $EBV_read_length_sum/$EBV_read_number;
}
#--------------------
if ($cell_iso_number == 0) {
	$cell_iso_mean = "NA";
}
else {
	$cell_iso_mean = $cell_iso_length_sum/$cell_iso_number;
}
if ($cell_read_number == 0) {
	$cell_read_mean = "NA";
}
else {
	$cell_read_mean = $cell_read_length_sum/$cell_read_number;
}

print "Total number of isoforms\t$iso_number\tMean of total isoform lengths\t$iso_mean\nTotal number of reads\t$read_number\tMean of total read lengths\t$read_mean\nNumber of mapped isoforms\t$mapped_iso_number\tMean of mapped isoform lengths\t$mapped_iso_mean\nNumber of mapped reads\t$mapped_read_number\tMean of mapped read lengths\t$mapped_read_mean\nNumber of unmapped isoforms\t$un_iso_number\tMean of unmapped isoform lengths\t$un_iso_mean\nNumber of unmapped reads\t$un_read_number\tMean of unmapped read lengths\t$un_read_mean\nNumber of EBV isoforms\t$EBV_iso_number\tMean of EBV isoform lengths\t$EBV_iso_mean\nNumber of EBV reads\t$EBV_read_number\tMean of EBV read lengths\t$EBV_read_mean\nNumber of cellular isoforms\t$cell_iso_number\tMean of cellular isoform lengths: $cell_iso_mean\nNumber of cellular reads\t$cell_read_number\tMean of cellular read lengths: $cell_read_mean\n";	
	
close(INF);