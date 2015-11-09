#!/usr/bin/perl

#Accepts a SAM file using SMRT fl data, a SAM file using CAGE data, and a bed file of annotated polyadenylated transcripts. Counts the number of non-clipped SMRT reads with 5' starts at each genomic position and estimates consensus locations of clusters of 5' starts. Uses Paraclu to identify clusters of 5' starts in the CAGE data. Output includes bedgraph files of SMRT 5' starts, a bed file of the weighted centers of SMRT start clusters, a bed file of Paraclu-identified CAGE 5' start clusters, and a bed file of SMRT 5' starts supported by the CAGE data, with their annotation status noted.

#For now, accepts a Paraclu output file rather than a CAGE SAM file.

#SMRT fl read names must be formatted as putative_isoform_id/number_of_reads/length.

#to check: Are soft-clipped and non-softclipped read files generated unnecessarily?
#       Add headers to all output files

#USAGE:
# perl <PATH/start_finder_sam_to_bed.pl> </PATH/SMRT_sam_file> </PATH/CAGE_file> </PATH/Annotation_bed_file>

#TO'G 9/11/2015

use warnings;
use strict;

#die "USAGE: 'perl <PATH/start_finder_sam_to_bed.pl> </PATH/SMRT_sam_file> </PATH/CAGE_file> </PATH/Annotation_bed_file>'" unless @ARGV == 3;

my ($SMRT_file, $CAGE_file, $paraclu_prep_file, $ann_file) = @ARGV;

print "Enter name of viral chromosome [e.g. chrEBV(Akata_107955to171322_1to107954)]: ";
my $viral_chr = <STDIN>;
chomp $viral_chr;

my $distance_between_SMRT_peaks;
my $min_tags;
my $min_dens;
my $min_length;
my $max_length;
my $dist_SMRT_CAGE;
my $min_SMRT;
my $ann_dist;

print "Use default parameters [y/n]? ";
my $answer = <STDIN>;
chomp $answer;

if ($answer eq "y") {
    $distance_between_SMRT_peaks = 8;
    $min_tags = 15;
    $min_dens = 2;
    $min_length = 1;
    $max_length = 20;
    $dist_SMRT_CAGE = 2;
    $min_SMRT = 1;
    $ann_dist = 20;
}
else {
    print "Enter desired window for collapsing SMRT 5' starts (e.g. 8): ";
    my $distance_between_SMRT_peaks = <STDIN>;
    chomp $distance_between_SMRT_peaks;

    print "Enter minimum tags per CAGE cluster: ";
    my $min_tags = <STDIN>;
    chomp $min_tags;

    print "Enter minimum relative density for CAGE clusters: ";
    my $min_dens = <STDIN>;
    chomp $min_dens;

    print "Enter minimum CAGE cluster length: ";
    my $min_length = <STDIN>;
    chomp $min_length;

    print "Enter maximum CAGE cluster length: ";
    my $max_length = <STDIN>;
    chomp $max_length;

    print "Enter desired maximum allowable distance between SMRT and CAGE 5' starts (e.g. 8): ";
    my $dist_SMRT_CAGE = <STDIN>;
    chomp $dist_SMRT_CAGE;

    print "Enter minimum number of SMRT reads to report a 5' start (e.g. 1): ";
    my $min_SMRT = <STDIN>;
    chomp $min_SMRT;

    print "Enter maximum distance in bp from an annotated start to be called as 'annotated' (e.g. 20): ";
    my $ann_dist = <STDIN>;
    chomp $ann_dist;
}

print "------------------------------------------------\n";

#####----------SMRT FILE PROCESSING-------------######
system("sort -k 3,3 -k 4,4n \Q$SMRT_file\E > \Q$SMRT_file\E.sorted.temp");
system("awk '\$2==0' \Q$SMRT_file\E.sorted.temp > \Q$SMRT_file\E.sorted.plus.sam.temp");
system("awk '\$2==16' \Q$SMRT_file\E.sorted.temp > \Q$SMRT_file\E.sorted.minus.sam.temp");
system("rm \Q$SMRT_file\E.sorted.temp");

#processing of PLUS SMRT sam file
print "Processing SMRT plus strand reads...\n";

open(INF, "<$SMRT_file.sorted.plus.sam.temp") or die "couldn't open file";
open(OUT, ">$SMRT_file.sorted.plus.sam.read_starts.bedgraph") or die "couldn't open file";

my $previous_coordinate=1;
my $count=0;
my $previous_chr = "start";

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    next if $cols[2] ne $viral_chr; #skips reads not mapped to the virus
    next if ($cols[5] =~ m/^\d+S/); #skips reads clipped at the 5' end
    my @split_id = split("\/", $cols[0]); #extracts the read depth for this putative isoform from its id
    if (($cols[2] eq $previous_chr) and ($cols[3] == $previous_coordinate)) {
        $count = $count + $split_id[1]; #increases the count by the read depth for the putative isoform
    }
    else {
        if ($previous_chr eq "start") { #doesn't print out the placeholder first line.
            $previous_chr = $cols[2];	#sets the previous chromosome, previous coordinate and count values
            $previous_coordinate = $cols[3];
            $count = $split_id[1];
        }
        else {
            print OUT $previous_chr, "\t", $previous_coordinate-1, "\t", $previous_coordinate-1, "\t", $count, "\n"; #prints to output file, converting to 0-based bedgraph coordinates
            $previous_chr = $cols[2];
            $previous_coordinate = $cols[3];
            $count = $split_id[1];
        }
    }
}

print OUT $previous_chr, "\t", $previous_coordinate-1, "\t", $previous_coordinate-1, "\t", $count, "\n"; #prints the last start coordinates to output file
close(INF);
close(OUT);

system("rm \Q$SMRT_file\E.sorted.plus.sam.temp");

#processing of MINUS SMRT sam file
open(INF, "<$SMRT_file.sorted.minus.sam.temp") or die "couldn't open file";
open(OUT, ">$SMRT_file.sorted.minus.sam.read_starts.bedgraph.temp") or die "couldn't open file";

my @dist;
my $sum;
my %minus_starts;

print "Processing SMRT minus strand reads...\n";

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    next if $cols[2] ne $viral_chr; #skips reads not mapped to the virus
    next if ($cols[5] =~ m/\d+S$/); #skips reads soft-clipped at the 5' end
    while ($cols[5] =~ /(\d+)[DMNX=]/g) { #these lines use the CIGAR string to determine the downstream coordinate
        push (@dist, $1);
    }
    $sum += $_ for @dist;
    my $start_coord = $cols[3] + $sum - 2; #subtract 1 to account for start/end inclusion and 1 to convert to 0-based bedgraph
    my $chr_start_coord = "$cols[2]\:$start_coord"; #combines the chromosome and 5' end coordinate into a key to use for the hash
    $sum = 0;
    @dist = ();
    my @split_id = split("\/", $cols[0]); #extracts the read depth for this putative isoform from its id
    if (exists $minus_starts{$chr_start_coord}) { #if the key is already in the hash, increases the value (count) by the read depth for that putative isoform
        $minus_starts{$chr_start_coord} = $minus_starts{$chr_start_coord} + $split_id[1];
    }
    else {
        $minus_starts{$chr_start_coord} = $split_id[1]; #if the key is not already in the hash, adds it with a value (count) of the read depth for that putative isoform
    }
}

foreach my $chr_start_coord (sort keys %minus_starts) { #prints out a(n inadequately) sorted temporary bedgraph file
    my @split_keys = split("\:", $chr_start_coord);
    print OUT "$split_keys[0]\t$split_keys[1]\t$split_keys[1]\t-$minus_starts{$chr_start_coord}\n";
}
close(INF);
close(OUT);

system("sort -k 1,1 -k 2,2n \Q$SMRT_file\E.sorted.minus.sam.read_starts.bedgraph.temp > \Q$SMRT_file\E.sorted.minus.sam.read_starts.bedgraph");

system("cat \Q$SMRT_file\E.sorted.plus.sam.read_starts.bedgraph \Q$SMRT_file\E.sorted.minus.sam.read_starts.bedgraph.temp | sort -k2,3n > \Q$SMRT_file\E.\Q$viral_chr\E.all_read_starts.bedgraph.noheader");

system("rm \Q$SMRT_file\E.sorted.minus.sam.read_starts.bedgraph.temp");
system("rm \Q$SMRT_file\E.sorted.minus.sam.read_starts.bedgraph");
system("rm \Q$SMRT_file\E.sorted.minus.sam.temp");
system("rm \Q$SMRT_file\E.sorted.plus.sam.read_starts.bedgraph");

#add header to bedgraph file
open(INF, "<$SMRT_file.$viral_chr.all_read_starts.bedgraph.noheader") or die "couldn't open file";
open(OUT, ">$SMRT_file.$viral_chr.all_read_starts.bedgraph") or die "couldn't open file";

print OUT "track type=bedgraph name=\"$SMRT_file.$viral_chr.all_read_starts.bedgraph\" description=\"5' starts of SMRT reads from start_finder_sam_to_bed.pl\"\n";
while (my $line = <INF>) {
    print OUT $line;
}
close(OUT);
close(INF);

system("rm \Q$SMRT_file\E.\Q$viral_chr\E.all_read_starts.bedgraph.noheader");

#make a bed file from the SMRT bedgraph file:
open(INF, "<$SMRT_file.$viral_chr.all_read_starts.bedgraph") or die "couldn't open file";
open(OUT, ">$SMRT_file.starts.temp.bed") or die "couldn't open file";

print "Combining SMRT 5' starts within $distance_between_SMRT_peaks of each other and calculating consensus 5' starts...\n";
collapse_bedgraph($distance_between_SMRT_peaks);

close(INF);
close(OUT);

system("sort -k 1,1 -k 2,2n \Q$SMRT_file\E.starts.temp.bed > \Q$SMRT_file\E.starts.bed.noheader");
system("rm \Q$SMRT_file.starts.temp.bed\E");

#add header to bed file
open(INF, "<$SMRT_file.starts.bed.noheader") or die "couldn't open file";
open(OUT, ">$SMRT_file.$viral_chr.starts.bed") or die "couldn't open file";

print OUT "track type=bed name=\"$SMRT_file.$viral_chr.starts.bed\" description=\"consensus 5' starts of SMRT reads within $distance_between_SMRT_peaks bp collapsed to weighted center from start_finder_sam_to_bed.pl\"\n";
while (my $line = <INF>) {
    print OUT $line;
}
close(OUT);
close(INF);

system("rm \Q$SMRT_file\E.starts.bed.noheader");


#####----------PROCESSING PARACLU OUTPUT-------------######

print "Extracting clusters from CAGE file (Containing $min_tags tags, density fold change at least $min_dens, from $min_length to $max_length bp long)\n";

#filtering clusters:

my $length;
my $dens;
my $prev_start = 0;
my $prev_end = 0;
my $chrStart;
my $chrEnd;

open(INF, "<$CAGE_file") or die "couldn't open file"; #for now, $CAGE_file is the output from Paraclu
open(OUT, ">$CAGE_file.peaks.$min_tags.$min_dens.$min_length.$max_length.bed") or die "couldn't open file";

while (my $line = <INF>) {
    chomp($line);
    next if ($line =~ /^#/); #skips the header line
    my @cols = split("\t", $line);
    next if ($cols[5] < $min_tags);
    $length = $cols[3] - $cols[2] + 1;
    if (($length >= $min_length) and ($length <= $max_length)) {
        $dens = $cols[7] / $cols[6];
        if ($dens >= $min_dens) {
            next if (($cols[2] >= $prev_start) and ($cols[2] <= $prev_end));
            $chrStart = $cols[2] - 1; #converts from 1-based (from SAM) to 0-based (for BED)
            $chrEnd = $cols[3] - 1; #converts from 1-based (from SAM) to 0-based (for BED)
            if ($dens < 100) {
                printf OUT "%s\t%d\t%d\t%d%s%.1f\t%d\t%s\n", $cols[0], $chrStart, $chrEnd, $cols[5], ":", $dens, $cols[5], $cols[1];   #limits the density output to 1 decimal place, but doesn't change huge numbers to exponents
            }
            #                elsif (($dens >= 100) and ($dens <= 10000)) {
            #                    printf OUT "%s\t%d\t%d\t%d%s%d\t%d\t%s\n", $cols[0], $chrStart, $chrEnd, $cols[5], ":", $dens, $cols[5], $cols[1];
            #                }
            else {
                printf OUT "%s\t%d\t%d\t%d%s%.1e\t%d\t%s\n", $cols[0], $chrStart, $chrEnd, $cols[5], ":", $dens, $cols[5], $cols[1]; #changes large numbers to exponents
                #print OUT "$cols[0]\t$cols[2]\t$cols[3]\t$cols[5]:$dens\t$cols[5]\t$cols[1]\n";
            }
            $prev_start = $cols[2];
            $prev_end = $cols[3];
        }
    }
}
close(INF);
close(OUT);

#getting weighted averages of Paraclu peaks:

my $chrStart_CAGE;
my $chrEnd_CAGE;
my $strand_CAGE;
my $CAGE_weighted_sum = 0;
#my $tag_depth = 0;
my $CAGE_weighted_average;

open(INF, "<$CAGE_file.peaks.$min_tags.$min_dens.$min_length.$max_length.bed") or die "couldn't open file";
open(OUT, ">$CAGE_file.peaks_weighted_average.bed") or die "couldn't open file"; #later, reduce output to just one file containing both the weighted average and the cluster extent

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    $chrStart_CAGE = $cols[1];
    $chrEnd_CAGE = $cols[2];
    $strand_CAGE = $cols[5];
    #$tag_depth = 0;
    open(INF2, "<$paraclu_prep_file") or die "couldn't open file";
    while (my $line2 = <INF2>) { #is it a problem to be looping through 2 files at once? Should I slurp one into an array?
        chomp($line2);
        my @cols2 = split("\t", $line2);
        if ((($cols2[2]-1) >= $chrStart_CAGE) and (($cols2[2]-1) <= $chrEnd_CAGE) and ($cols2[1] eq $strand_CAGE)) {
            $CAGE_weighted_sum = $CAGE_weighted_sum + ($cols2[2]*$cols2[3]);
            #$tag_depth = $tag_depth + $cols2[3];
        }
    }
    $CAGE_weighted_average = ($CAGE_weighted_sum/$cols[4]) - 1;
    printf OUT "%s\t%1.0f\t%1.0f\t%s%s%s%s%s\t%s\t%s\n", $cols[0], $CAGE_weighted_average, $CAGE_weighted_average, $chrStart_CAGE, ":", $chrEnd_CAGE, ":", $cols[3], $cols[4], $cols[5];
    $CAGE_weighted_sum = 0;
    close(INF2);
}

close(INF);
close(OUT);




#####----------SEEKING CAGE SUPPORT FOR SMRT STARTS-------------######

#open(INF, "<$CAGE_file.peaks.$min_tags.$min_dens.$min_length.$max_length.bed" ) or die "couldn't open file";
open(INF, "<$CAGE_file.peaks_weighted_average.bed" ) or die "couldn't open file"; #change this back when you revert to outputting a single Paraclu bed file

print "Extracting SMRT 5' starts within $dist_SMRT_CAGE bases of CAGE clusters...\n";

my %features_CAGE;

while(my $line = <INF> ) {
	chomp($line);
	my @cols = split("\t", $line);
	my $key_combo_CAGE = "$cols[0]:$cols[1]:$cols[2]:$cols[5]"; #for each line in the CAGE bed file, creates a key for the hash combining chromosome, start coordinate, end coordinate and strand. Doesn't really need $cols[2], does it?
	$features_CAGE{$key_combo_CAGE} = $cols[4]; #enters a count value for the key into the hash
}

close(INF);

open(INF, "<$SMRT_file.$viral_chr.starts.bed" ) or die "couldn't open file";
open(OUT, ">$SMRT_file.$viral_chr.starts.bed.CAGE_support.bed.temp");

while(my $line = <INF>) {
	chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
	my @SMRT_cols = split("\t", $line);
    next if (abs $SMRT_cols[4] < $min_SMRT);
    my $found_flag=0;
    foreach my $key_combo_CAGE (keys %features_CAGE) {
        my @CAGE_cols = split(":", $key_combo_CAGE);
        my $lower_limit = $SMRT_cols[1]-$dist_SMRT_CAGE;
        my $upper_limit = $SMRT_cols[1]+$dist_SMRT_CAGE;
        if (($SMRT_cols[5] eq $CAGE_cols[3]) and ($CAGE_cols[1] >= $lower_limit) and ($CAGE_cols[1] <= $upper_limit)) {
            my $name = "$SMRT_cols[4]SMRT_$features_CAGE{$key_combo_CAGE}CAGE";
            my $count = $features_CAGE{$key_combo_CAGE} + $SMRT_cols[4];
            print OUT "$SMRT_cols[0]\t$SMRT_cols[1]\t$SMRT_cols[2]\t$name\t$count\t$SMRT_cols[5]\t$SMRT_cols[3]\n";
            $found_flag = 1;
            last; #if the SMRT start is supported by more than one CAGE cluster, only one is reported
        }
    }
    if ($found_flag == 0) { #clean this up (just remove?) if you don't want to validate based on SMRT + ann without CAGE
        my @range_cols = split (":", $SMRT_cols[3]);
        #print OUT "$SMRT_cols[0]\t$SMRT_cols[1]\t$SMRT_cols[2]\t$range_cols[2]SMRT\t$range_cols[2]\t$SMRT_cols[5]\t$SMRT_cols[3]\n";
    }
}

close(OUT);
close(INF);


#####----------COMPARING TO ANNOTATED STARTS-------------######
open(INF, "<$ann_file" ) or die "couldn't open file";
#open(OUT, ">$ann_file.starts_only.bed");

print "Processing annotation file...\n";

#extract 5' starts from the annotation file:
#annotation file must be sorted by chrStart then chrEnd!
my @annotated_starts;
my $plus_prev_chr = 0;
my $plus_prev_coord = 0;
my $minus_prev_chr = 0;
my $minus_prev_coord = 0;

while(my $line = <INF>) {
    chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
	my @ann_cols = split("\t", $line);
    if ($ann_cols[5] eq "+") {
        if ($ann_cols[1] != $plus_prev_coord) {
            push (@annotated_starts, "$ann_cols[0]:$ann_cols[1]:$ann_cols[3]:$ann_cols[5]");
            $plus_prev_chr = $ann_cols[0];
            $plus_prev_coord = $ann_cols[1];
        }
    }
    elsif ($ann_cols[5] eq "-"){
        if ($ann_cols[2] != $minus_prev_coord) {
            push (@annotated_starts, "$ann_cols[0]:$ann_cols[2]:$ann_cols[3]:$ann_cols[5]");
            $minus_prev_chr = $ann_cols[0];
            $minus_prev_coord = $ann_cols[2];
        }
    }
}

#foreach (@annotated_starts){
#    print OUT "$_\n";
#}

my $annotated = scalar @annotated_starts;

close(INF);
#close(OUT);

#compare starts in the altered SMRT starts file (that already has info about CAGE starts) with annotated starts

open(INF, "<$SMRT_file.$viral_chr.starts.bed.CAGE_support.bed.temp" ) or die "couldn't open file";
open(OUT, ">$SMRT_file.$viral_chr.validated_starts.bed");

print "Comparing SMRT starts to annotated starts...\n";

print OUT "track type=bedDetail name=\"$SMRT_file.$viral_chr.starts.bed.CAGE_support.bed\" description=\"consensus SMRT 5' starts of collapse value 8 supported by at least $min_SMRT reads within $dist_SMRT_CAGE bp of CAGE clusters or within $ann_dist of annotated starts. From start_finder_sam_to_bed.pl\"\n";

my $annotated_found_by_SMRT = 0;
my $novel_found_by_SMRT_CAGE = 0;
my $SMRT_annotated = 0; #this is different than $annotated_found_by_SMRT because depending on input parameters two SMRT starts may correspond to a single annotated start or vice versa.

while(my $line = <INF>) {
    chomp($line);
    my @SMRT_cols = split("\t", $line);
    my $found_flag=0;
    foreach my $ann_start (@annotated_starts) {
        my @ann_cols = split(":", $ann_start);
        my $lower_limit = $ann_cols[1]-$ann_dist;
        my $upper_limit = $ann_cols[1]+$ann_dist;
        if (($SMRT_cols[5] eq $ann_cols[3]) and ($SMRT_cols[1]>=$lower_limit) and ($SMRT_cols[1]<=$upper_limit)) {
            if ($found_flag == 0) {
                print OUT "$SMRT_cols[0]\t$SMRT_cols[1]\t$SMRT_cols[2]\tann_$SMRT_cols[3]\t$SMRT_cols[4]\t$SMRT_cols[5]\t$SMRT_cols[6]\n";
                $found_flag = 1;
                $annotated_found_by_SMRT++;
                $SMRT_annotated++;
            }
            elsif ($found_flag == 1) {
                $annotated_found_by_SMRT++;
            }
        }
    }
    if ($found_flag == 0) {
        if ($SMRT_cols[3] =~ /.+SMRT_.+CAGE/) {
            print OUT "$SMRT_cols[0]\t$SMRT_cols[1]\t$SMRT_cols[2]\tnov_$SMRT_cols[3]\t$SMRT_cols[4]\t$SMRT_cols[5]\t$SMRT_cols[6]\n";
            $novel_found_by_SMRT_CAGE++;
        }
    }
}


my $total_found = $SMRT_annotated + $novel_found_by_SMRT_CAGE;

print "------------------------------------------------\n";

if ($SMRT_annotated != $annotated_found_by_SMRT) {
    print "$total_found 5' starts found. $novel_found_by_SMRT_CAGE are novel, $SMRT_annotated are annotated.  $annotated_found_by_SMRT out of $annotated total annotated 5' starts are found.\nNote that two annotated starts may be within $ann_dist bp of a single SMRT start or vice versa.\n\n";
}
else {
    print "$total_found 5' starts found. $novel_found_by_SMRT_CAGE are novel, $SMRT_annotated are annotated (out of a total of $annotated annotated 5' starts).\n\n";
}

close(INF);
close(OUT);

system("rm \Q$SMRT_file\E.\Q$viral_chr\E.starts.bed.CAGE_support.bed.temp\E");

#########################
sub collapse_bedgraph {
    my ($distance_between_peaks) = shift;
    my $prev_coord_plus = 1;
    my $prev_coord_minus = 1;
    my $count_sum_plus = 0;
    my $count_sum_minus = 0;
    my $weighted_coordinate_sum_plus = 0;
    my $weighted_coordinate_sum_minus = 0;
    my $weighted_average_plus;
    my $weighted_average_minus;
    my $first_plus = 1;
    my $first_minus = 1;
    my @coords_plus;
    my @coords_minus;
    my $chrStart_plus;
    my $chrEnd_plus;
    my $chrStart_minus;
    my $chrEnd_minus;
    
    while (my $line = <INF>) {
        chomp($line);
        next if ($line =~ /^track/); #skips the track definition line
        my @cols = split("\t", $line);
        if ($cols[3] > 0) { #if this coordinate has a positive count...
            if ($cols[1] < $prev_coord_plus + ($distance_between_peaks+1)) { #if the coordinate is within the specified number of bp of the previous coordinate
                $count_sum_plus = $count_sum_plus + $cols[3]; #adds to the sums to eventually calculate the weighted average
                $weighted_coordinate_sum_plus = $weighted_coordinate_sum_plus + ($cols[1]*$cols[3]);
                push (@coords_plus, $cols[1]);
                $prev_coord_plus = $cols[1]; #sets the current coordinate as the "previous coordinate" before moving on
            }
            else { #if the present coordinate is not within the specified number of bp of the previous coordinate, need to print out a feature
                if ($first_plus == 1) { #"first" flag avoids wonkiness if the first coordinate is far from coordinate 1 (don't need to print out a feature yet)
                    $count_sum_plus = $cols[3];
                    $weighted_coordinate_sum_plus = $cols[1]*$cols[3];
                    $prev_coord_plus = $cols[1];
                    push (@coords_plus, $cols[1]);
                    $first_plus = 0;
                }
                else {
                    $weighted_average_plus = ($weighted_coordinate_sum_plus/$count_sum_plus); #calculates weighted average, subtracts 1 to make it 0-based for bed file
                    $chrStart_plus = $coords_plus[0];
                    $chrEnd_plus = pop(@coords_plus);
                    printf OUT "%s\t%1.0f\t%1.0f\t%d%s%d%s%d\t%d\t%s\n", $viral_chr, $weighted_average_plus, $weighted_average_plus, $chrStart_plus, ":", $chrEnd_plus, ":", $count_sum_plus, $count_sum_plus, "+"; #prints out weighted average for plus strand features. Use printf to round the weighted average.
                    @coords_plus = ($cols[1]);
                    $count_sum_plus = $cols[3]; #sets "previous coordinate", count and sum of counts for the current coordinate
                    $weighted_coordinate_sum_plus = $cols[1]*$cols[3];
                    $prev_coord_plus = $cols[1];
                }
            }
        }
        elsif ($cols[3] < 0) { #if this coordinate has a negative count...
            if ($cols[1] < $prev_coord_minus + ($distance_between_peaks+1)) { #if the coordinate is within the specified number of bp of the previous coordinate
                $count_sum_minus = $count_sum_minus + $cols[3]; #adds to the sums to eventually calculate the weighted average
                $weighted_coordinate_sum_minus = $weighted_coordinate_sum_minus + ($cols[1]*$cols[3]);
                push (@coords_minus, $cols[1]);
                $prev_coord_minus = $cols[1]; #sets the current coordinate as the "previous coordinate" before moving on
            }
            else { #if the present coordinate is not within the specified number of bp of the previous coordinate, need to print out a feature
                if ($first_minus == 1) { #"first" flag avoids wonkiness if the first coordinate is far from coordinate 1 (don't need to print out a feature yet)
                    $count_sum_minus = $cols[3];
                    $weighted_coordinate_sum_minus = $cols[1]*$cols[3];
                    $prev_coord_minus = $cols[1];
                    push (@coords_minus, $cols[1]);
                    $first_minus = 0;
                }
                else {
                    $weighted_average_minus = ($weighted_coordinate_sum_minus/$count_sum_minus); #calculates weighted average
                    $chrStart_minus = $coords_minus[0];
                    $chrEnd_minus = pop(@coords_minus);
                    printf OUT "%s\t%1.0f\t%1.0f\t%d%s%d%s%d\t%d\t%s\n", $viral_chr, $weighted_average_minus, $weighted_average_minus, $chrStart_minus, ":", $chrEnd_minus, ":", $count_sum_minus, $count_sum_minus, "-"; #prints out weighted average for plus strand features. Use printf to round the weighted average.
                    @coords_minus = ($cols[1]);
                    $count_sum_minus = $cols[3]; #sets "previous coordinate", count and sum of counts for the current coordinate
                    $weighted_coordinate_sum_minus = $cols[1]*$cols[3];
                    $prev_coord_minus = $cols[1];
                }
            }
        }
    }
    
    if ($count_sum_plus > 0) {#calculates and prints out weighted average for the last feature (plus strand)
        $weighted_average_plus = ($weighted_coordinate_sum_plus/$count_sum_plus);
        $chrStart_plus = $coords_plus[0];
        $chrEnd_plus = pop(@coords_plus);
        printf OUT "%s\t%1.0f\t%1.0f\t%d%s%d%s%d\t%d\t%s\n", $viral_chr, $weighted_average_plus, $weighted_average_plus, $chrStart_plus, ":", $chrEnd_plus, ":", $count_sum_plus, $count_sum_plus, "+";
    }
    
    if ($count_sum_minus < 0) {#calculates and prints out weighted average for the last feature (minus strand)
        $weighted_average_minus = ($weighted_coordinate_sum_minus/$count_sum_minus);
        $chrStart_minus = $coords_minus[0];
        $chrEnd_minus = pop(@coords_minus);
        printf OUT "%s\t%1.0f\t%1.0f\t%d%s%d%s%d\t%d\t%s\n", $viral_chr, $weighted_average_minus, $weighted_average_minus, $chrStart_minus, ":", $chrEnd_minus, ":", $count_sum_minus, $count_sum_minus, "-";
    }
}