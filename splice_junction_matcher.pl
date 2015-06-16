#!/usr/bin/perl

#Accepts a junctions files from GMAP/SMRT (generated with the -f introns argument) and an SJ.out.tab files from STAR/Illumina. Returns 3 bed files: one of SMRT splice junctions, one of Illumina splice junctions and one of junctions detected by both methods.

#SMRT fl read names must be formatted as putative_isoform_id/number_of_reads/length.

#USAGE:
# perl <PATH/splice_junction_matcher.pl> </PATH/SMRT_introns_file> </PATH/Illumina_SJ.out.tab_file>

#TO'G 6/11/2015

use warnings;
use strict;

my ($SMRT_jfile, $ill_jfile, $ann_file) = @ARGV;

print "Enter name of viral chromosome [e.g. chrEBV(Akata_107955to171322_1to107954)]: ";
my $viral_chr = <STDIN>;
chomp $viral_chr;

print "Enter minimum SMRT read depth to report a splice junction (e.g. 1): ";
my $min_SMRTj = <STDIN>;
chomp $min_SMRTj;

print "Enter minimum Illumina read depth to report a splice junction (e.g. 1): ";
my $min_illj = <STDIN>;
chomp $min_illj;

print "------------------------------------------------\n";

#####----------GMAP/SMRT FILE CONVERSION-------------######
open(INF, "<$SMRT_jfile");
open(OUT, ">$SMRT_jfile.temp");

print "Processing SMRT splice junctions...\n";

while(my $line = <INF> ) {
    chomp($line);
    my ($id) = $line =~ /\>(.+)\.i/;
    my ($chr) = $line =~ /\s(.+):/;
    my ($score) = $line =~ /\>.+\/(\d+)\//;
    my ($donor, $acceptor) = $line =~ /:(\d+)\.\.(\d+)/;
    next if $chr ne $viral_chr;
    if ($acceptor > $donor) {
        print OUT $chr, "\t", $donor, "\t", $acceptor - 1, "\t", $id, "\t", $score, "\t+\n";
    }
    else {
        print OUT $chr, "\t", $acceptor, "\t", $donor - 1, "\t", $id, "\t", $score, "\t-\n";
    }
}
close(OUT);
close(INF);

system("sort -k2,3n \Q$SMRT_jfile\E.temp > \Q$SMRT_jfile\E.sorted.temp"); #sorts so that duplicate introns will be next to each other.

open(INF, "<$SMRT_jfile.sorted.temp" ) or die "couldn't reopen file";
open(OUT, ">$SMRT_jfile.bed.temp");

my $plus_previous_chr = "start";
my $plus_count = 0;
my $plus_previous_start = 0;
my $plus_previous_end = 0;
my $minus_previous_chr = "start";
my $minus_count = 0;
my $minus_previous_start = 0;
my $minus_previous_end = 0;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    if ($cols[5] eq "+") { #plus and minus need to be treated separately in case of introns with the same starts and ends annotated on opposite strands
        if ($cols[0] eq $plus_previous_chr && $cols[1] == $plus_previous_start && $cols[2] == $plus_previous_end) { #checks to see if the intron matches the previous intron
            $plus_count = $plus_count + $cols[4];
        }
        else {
            if ($plus_previous_chr eq "start") { #prevents the initial placeholder value from printing out as a line, and sets the values of the first intron
                $plus_previous_chr = $cols[0];
                $plus_previous_start = $cols[1];
                $plus_previous_end = $cols[2];
                $plus_count = $cols[4];
            }
            else {
                print OUT "$plus_previous_chr\t$plus_previous_start\t$plus_previous_end\t$plus_count\t$plus_count\t+\n";
                $plus_previous_chr = $cols[0];
                $plus_previous_start = $cols[1];
                $plus_previous_end = $cols[2];
                $plus_count = $cols[4];
            }
        }
    }
    if ($cols[5] eq "-") {
        if ($cols[0] eq $minus_previous_chr && $cols[1] == $minus_previous_start && $cols[2] == $minus_previous_end) {
            $minus_count = $minus_count + $cols[4];
        }
        else {
            if ($minus_previous_chr eq "start") {
                $minus_previous_chr = $cols[0];
                $minus_previous_start = $cols[1];
                $minus_previous_end = $cols[2];
                $minus_count = $cols[4];
            }
            else {
                print OUT "$minus_previous_chr\t$minus_previous_start\t$minus_previous_end\t$minus_count\t$minus_count\t-\n"; #prints out in bed format
                $minus_count = $cols[4];
                $minus_previous_chr = $cols[0];
                $minus_previous_start = $cols[1];
                $minus_previous_end = $cols[2];
            }
        }
    }
}

print OUT "$plus_previous_chr\t$plus_previous_start\t$plus_previous_end\t$plus_count\t$plus_count\t+\n"; #adds the last plus strand feature
print OUT "$minus_previous_chr\t$minus_previous_start\t$minus_previous_end\t$minus_count\t$minus_count\t-\n"; #adds the last plus strand feature	
close(OUT);
close(INF);

system("sort -k2,3n \Q$SMRT_jfile\E.bed.temp > \Q$SMRT_jfile\E.\Q$viral_chr\E.bed.noheader");
system("rm \Q$SMRT_jfile\E.temp");
system("rm \Q$SMRT_jfile.sorted\E.temp");
system("rm \Q$SMRT_jfile\E.bed.temp");


#add header to bed file
open(INF, "<$SMRT_jfile.$viral_chr.bed.noheader") or die "couldn't open file";
open(OUT, ">$SMRT_jfile.$viral_chr.bed") or die "couldn't open file";

print OUT "track type=bed name=\"$SMRT_jfile.$viral_chr.bed\" description=\"SMRT introns from splice_junction_matcher.pl\"\n";
while (my $line = <INF>) {
    print OUT $line;
}
close(OUT);
close(INF);

system("rm \Q$SMRT_jfile\E.\Q$viral_chr\E.bed.noheader");

#####----------STAR/ILLUMINA FILE CONVERSION-------------######

open(INF, "<$ill_jfile" ) or die "couldn't open file";
open(OUT, ">$ill_jfile.$viral_chr.bed");

print "Processing Illumina splice junctions...\n";
print OUT "track type=bed name=\"$ill_jfile.$viral_chr.bed\" description=\"Illumina STAR introns from splice_junction_matcher.pl\"\n";

while(my $line = <INF> ) {
    chomp($line);
    my @cols = split( "\t", $line );
    tr/12/+-/ foreach ($cols[3]);	#change the numeric strand indicators to + or -
    next if $cols[0] ne $viral_chr; #skip lines that aren't viral
    my $chrStart = $cols[1]-1; #fixes the start coordinate, which is off by 1
    print OUT "$cols[0]\t$chrStart\t$cols[2]\t$cols[4]\t$cols[6]\t$cols[3]\n";
}

close(OUT);
close(INF);

#####----------GMAP/ILLUMINA COMPARISON-------------######

open(INF, "<$ill_jfile.$viral_chr.bed" ) or die "couldn't open file";

print "Checking for matching splice junctions...\n";

my %ill_junctions;

while(my $line = <INF> ) {
	chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
	my @cols = split("\t", $line);
    next if ($cols[4] < $min_illj);
	my $ill_key_combo = "$cols[0]$cols[1]$cols[2]$cols[5]"; #for each line in the Illumina file, creates a key for the hash combining chromosome, start coordinate, end coordinate and strand
	$ill_junctions{$ill_key_combo} = $cols[4]; #enters a count value for the key into the hash
}

close(INF);

open(INF, "<$SMRT_jfile.$viral_chr.bed" ) or die "couldn't open file";
open(OUT, ">$SMRT_jfile.$viral_chr.illumina_support.bed.temp");

while(my $line = <INF>) {
	chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
	my @cols = split("\t", $line);
    next if ($cols[4] < $min_SMRTj);
    my $SMRT_key_combo = "$cols[0]$cols[1]$cols[2]$cols[5]"; #for each line in the SMRT file, creates a variable/key combining chromosome, start coordinate, end coordinate and strand
	if (exists $ill_junctions{$SMRT_key_combo}) { #checks to see if the key exists in the Illumina hash: if so, prints it out
        my $junction_depth = $cols[4] + $ill_junctions{$SMRT_key_combo};
		print OUT "$cols[0]\t$cols[1]\t$cols[2]\t$cols[4]SMRT_$ill_junctions{$SMRT_key_combo}Ill\t$junction_depth\t$cols[5]\n";
	}
    else {
        print OUT "$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]SMRT\t$cols[4]\t$cols[5]\n";
    }
}
close(INF);
close(OUT);

#####----------ANNOTATION FILE COMPARISON-------------######

#First extract intron coordinates from the annotation file
open(INF, "<$ann_file");

print "Processing annotation file...\n";

my @intron_start;
my @intron_end;
my @ann_intron_coord_pair = ();
my $start;
my $end;

while (my $line = <INF>) {
    chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
    my @cols = split("\t", $line);
    my $intron_number = $cols[9] - 1;
    next if ($intron_number == 0);
    my @block_sizes = split(",", $cols[10]);
    my @block_starts = split(",", $cols[11]);
    for (my $i = 0; $i < $intron_number; $i = $i + 1) { #for the transcript currently in the "while" loop, creates an array of intron start sites relative to the genome
        $start = $cols[1] + $block_sizes[$i] + $block_starts[$i];
        push(@intron_start, $start);
    }
    for (my $i2 = 1; $i2 < $cols[9]; $i2 = $i2 + 1) { #for the transcript currently in the "while" loop, creates an array of intron end sites relative to the genome
        $end = $cols[1] + $block_starts[$i2];
        push(@intron_end, $end);
    }
    for (my $i3 = 0; $i3 < $intron_number; $i3 = $i3 + 1) { #for the transcript currently in the "while" loop, matches up intron start and end sites to create an array of complete intron coordinates relative to the genome
        my $intron_coords = "$cols[0]:$intron_start[$i3]:$intron_end[$i3]:$cols[5]";
        push (@ann_intron_coord_pair, $intron_coords);
    }
    @intron_start = ();
    @intron_end = (); #intron starts and ends have been assigned to the @intron_coords array; empty them for the next transcript
}

close(INF);

#Compare introns in the altered (with Illumina data) SMRT file to annotated introns

open(INF, "<$SMRT_jfile.$viral_chr.illumina_support.bed.temp");
open(OUT, ">$SMRT_jfile.$viral_chr.validated_introns.bed");

print "Comparing SMRT junctions to annotation file...\n";

print OUT "track type=bed name=\"$SMRT_jfile.$viral_chr.validated_introns.bed\" description=\"Introns detected by SMRT with read depth at least $min_SMRTj supported by Illumina-detected junctions with read depth at least $min_illj and/or annotation. From splice_junction_matcher.pl\"\n";

while (my $line = <INF>) {
    chomp($line);
    my @SMRT_cols = split("\t", $line);
    my $found_flag=0;
    foreach my $ann_intron_coords (@ann_intron_coord_pair) {
        my @ann_cols = split (":", $ann_intron_coords);
        if (($SMRT_cols[0] eq $ann_cols[0]) and ($SMRT_cols[1] == $ann_cols[1]) and ($SMRT_cols[2] == $ann_cols[2]) and ($SMRT_cols[5] eq $ann_cols[3])) {
            print OUT "$SMRT_cols[0]\t$SMRT_cols[1]\t$SMRT_cols[2]\tann_$SMRT_cols[3]\t$SMRT_cols[4]\t$SMRT_cols[5]\n";
            $found_flag = 1;
            last;
        }
    }
    if ($found_flag == 0) {
        if ($SMRT_cols[3] =~ /.+SMRT_.+Ill/) {
            print OUT "$SMRT_cols[0]\t$SMRT_cols[1]\t$SMRT_cols[2]\tnov_$SMRT_cols[3]\t$SMRT_cols[4]\t$SMRT_cols[5]\n";
        }
    }
}
close(OUT);
close(INF);

system ("rm \Q$SMRT_jfile\E.\Q$viral_chr\E.illumina_support.bed.temp");