#!/usr/bin/perl

#To get information about full-length coverage of transcripts by SMRT Iso-seq data

#TO'G 7/31/15

use warnings;
use strict;

my ($SMRT_file, $ann_file, $ill_file) = @ARGV;

print "Enter cushion distance for 5' ends: ";
my $start_cushion = <STDIN>;
chomp $start_cushion;

print "Enter cushion distance for 3' ends: ";
my $end_cushion = <STDIN>;
chomp $end_cushion;

my $datestring = localtime();
print "-------------\nstart at $datestring\n-------------\n";

#First, open a SMRT alignment SAM file. Extract the 5' and 3' ends and create a hash of chromosome-strand-start-end sets. Collapse duplicates and keep track of depth.

open(INF, "<$SMRT_file");
open(OUT, ">$SMRT_file.hash.txt");

my @dist;
my $sum;
my %SMRT_coord_hash;

while (my $line = <INF>) {
    $line =~ s/\r//g; #from Erik's sam-to-bed converter. Are there weird carriage returns in here that chomp doesn't get?
    chomp($line);
    next if ($line =~ m/^@/); # ignores SAM header lines
    my @cols = split("\t", $line);
    #my @split_id = split("\/", $cols[0]);
    my $strand;
    my $chr = $cols[2];
    my $chrStart = $cols[3] - 1;
    
    
    while ($cols[5] =~ /(\d+)[DMNX=]/g) { #these lines use the CIGAR string to determine the downstream coordinate
        push (@dist, $1);
    }
    $sum += $_ for @dist;
    my $chrEnd = $cols[3] + $sum - 1;
    
    # identify strand
    if ($cols[1] == 0) {
        $strand = "+";
    }
    elsif ($cols[1] == 16) {
        $strand = "-";
    }
    else {
        $strand = "NA";
    }
    
    my $SMRT_coords = "$chr:$strand:$chrStart:$chrEnd"; #combines chromosome, strand, start and end for the hash
    $sum = 0;
    @dist = ();
    
    if (exists $SMRT_coord_hash{$SMRT_coords}) { #if the key is already in the hash, moves on
        next;
    }
    else {
        $SMRT_coord_hash{$SMRT_coords} = $cols[0]; #if not, adds the key to the hash with the read name as its value
    }
}

foreach my $key (keys %SMRT_coord_hash) {
    print OUT "$key $SMRT_coord_hash{$key}\n";
}

print "SMRT file hash built\n-------------\n";

close(INF);
close(OUT);

#Then, open the bed annotation file. Use a loop to check if the ends are within x bp of ends in the SMRT array. Probably let the user set x so we can check it a few ways.

open(INF, "<$ann_file");
open(OUT, ">$ann_file.no_coverage.bed");
open(OUT2, ">$ann_file.end_coverage.bed");
open(OUT3, ">$ann_file.full_coverage.bed");

my $chr_report = "start";
my @size;
my $length;
my $end_found_flag=0;
my $both_found_flag=0;
my $SMRT_read;
my $plus_boundary;
my $minus_boundary;
my $ill_found_flag=0;

while (my $line = <INF>) {
    chomp($line);
    next if ($line =~ m/^track/);
    my @cols = split("\t", $line);
    my $chr = $cols[0];
    my $chrStart = $cols[1];
    my $chrEnd = $cols[2];
    my $strand = $cols[5];
    my $upper_end_limit;
    my $lower_end_limit;
    my $upper_start_limit;
    my $lower_start_limit;
    
    if ($chr ne $chr_report) { #reports progress (Script can be slow)
        print "processing $chr\n-------------\n";
        $chr_report = $chr;
    }
    
    $end_found_flag=0;
    $both_found_flag=0;
    foreach my $SMRT_coords (keys %SMRT_coord_hash) {
        my @SMRT_cols = split(":", $SMRT_coords);
        if (($chr eq $SMRT_cols[0]) and ($strand eq $SMRT_cols[1])) {
            if ($strand eq "+") {
                $lower_end_limit = $chrEnd - $end_cushion;
                $upper_end_limit = $chrEnd + $end_cushion;
                if (($SMRT_cols[3] > $lower_end_limit) and ($SMRT_cols[3] < $upper_end_limit)){
                    $end_found_flag=1;
                    $lower_start_limit = $chrStart - $start_cushion;
                    $upper_start_limit = $chrStart + $start_cushion;
                    if (($SMRT_cols[2] > $lower_start_limit) and ($SMRT_cols[2] < $upper_start_limit)){
                        $both_found_flag=1;
                        $SMRT_read = $SMRT_coord_hash{$SMRT_coords};
                        last;
                    }
                }
            }
            if ($strand eq "-") {
                $lower_end_limit = $chrStart - $end_cushion;
                $upper_end_limit = $chrStart + $end_cushion;
                if (($SMRT_cols[2] > $lower_end_limit) and ($SMRT_cols[2] < $upper_end_limit)){
                    $end_found_flag=1;
                    $lower_start_limit = $chrEnd - $start_cushion;
                    $upper_start_limit = $chrEnd + $start_cushion;
                    if (($SMRT_cols[3] > $lower_start_limit) and ($SMRT_cols[3] < $upper_start_limit)){
                        $both_found_flag=1;
                        $SMRT_read = $SMRT_coord_hash{$SMRT_coords};
                        last;
                    }
                }
            }
        }
    }
    
    while ($cols[10] =~ /(\d+),/g) {
        push (@size, $1);
    }
    
    $length += $_ for @size;
      
    my $out_line = "$chr\t$chrStart\t$chrEnd\t$cols[3]\t$length\t$strand\t$cols[6]\t$cols[7]\t$cols[8]\t$cols[9]\t$cols[10]\t$cols[11]\n";
 
    if ($both_found_flag == 1) {
        print OUT3 $out_line;
    }
    if (($end_found_flag == 1) and ($both_found_flag == 0)) {
        #at this point, check for Illumina support w/in 100 bp of the annotated 5' site.
        open(INF2, "<$ill_file");
        while (my $line2 = <INF2>) { #or should this be a foreach?
            chomp($line2);
            my @cols2 = split("\t", $line2);
            $ill_found_flag = 0;
            if (($chr eq $cols2[0]) and ($strand eq $cols2[5])) {
                if ($cols2[5] eq "+") {
                    $plus_boundary = ($chrStart + 100);
                    if (($cols2[1] >= $chrStart) and ($cols2[1] <= $plus_boundary)) {
                        print OUT2 $out_line;
                        $ill_found_flag = 1;
                        last;
                    }
                }
                if ($cols2[5] eq "-") {
                    $minus_boundary = ($chrEnd - 100);
                    if (($cols2[2] >= $minus_boundary) and ($cols2[2] <= $chrEnd)) {
                        print OUT2 $out_line;
                        $ill_found_flag = 1;
                        last;
                    }
                }
            }
        }
        if ($ill_found_flag == 0) {
            print OUT $out_line;
        }
       close(INF2);
    }
    if ($end_found_flag == 0) {
        print OUT $out_line;
    }
    
    $end_found_flag = 0;
    $both_found_flag = 0;
    @size = ();
    $length = 0;
}

close(INF);
close(OUT);
close(OUT2);
close(OUT3);


open(INF, "<$ann_file.full_coverage.bed");

my $full_length_sum = 0;
my $full_transcript_count = 0;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    $full_length_sum = $full_length_sum + $cols[4];
    $full_transcript_count = $full_transcript_count + 1;
}

my $full_length = $full_length_sum / $full_transcript_count;

print "-------------\nAverage length of transcript with full-length coverage: $full_length\n";

close(INF);

open(INF, "<$ann_file.end_coverage.bed");

my $part_length_sum = 0;
my $part_transcript_count = 0;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    $part_length_sum = $part_length_sum + $cols[4];
    $part_transcript_count = $part_transcript_count + 1;
}

my $part_length = $part_length_sum / $part_transcript_count;

my $datestring_end = localtime();

print "Average length of transcript with 3' end coverage only: $part_length\n";
print "-------------\nend at $datestring_end\n-------------\n";

close(INF);