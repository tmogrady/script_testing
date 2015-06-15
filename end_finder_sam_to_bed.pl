#!/usr/bin/perl

#Accepts a SAM file using SMRT fl data and a SAM file using Illumina data. Counts the number of non-clipped reads with 3' ends at each genomic position and estimates consensus locations of clusters of 3 ends. Output includes wiggle files of all 3' ends and bed files of the weighted centers of end clusters. 

#SMRT fl read names must be formatted as putative_isoform_id/number_of_reads/length.

#USAGE:
# perl <PATH/read_end_finder.pl> </PATH/SMRT_sam_file> </PATH/Illumina_sam_file>

#TO'G 6/11/2015

use warnings;
use strict;

my ($SMRT_file, $ill_file) = @ARGV;

print "Enter name of viral chromosome [e.g. chrEBV(Akata_107955to171322_1to107954)]: ";
my $viral_chr = <STDIN>;
chomp $viral_chr;

print "Enter desired window for collapsing SMRT 3' ends (e.g. 8): ";
my $distance_between_SMRT_peaks = <STDIN>;
chomp $distance_between_SMRT_peaks;

print "Enter minimum number of As for Illumina polyA tails (e.g. 5): ";
my $min_As = <STDIN>;
chomp $min_As;

print "Enter minimum number of mismatches for Illumina polyA tails (e.g. 2): ";
my $min_softclip = <STDIN>;
chomp $min_softclip;

print "Enter desired window for collapsing Illumina 3' ends (e.g. 8): ";
my $distance_between_ill_peaks = <STDIN>;
chomp $distance_between_ill_peaks;

print "Enter desired maximum allowable distance between SMRT and Illumina 3' ends (e.g. 8): ";
my $dist_SMRT_ill = <STDIN>;
chomp $dist_SMRT_ill;

print "Enter minimum number of SMRT reads to report a 3' end (e.g. 5): ";
my $min_SMRT = <STDIN>;
chomp $min_SMRT;

print "Enter minimum number of Illumina polyA tails to support a 3' end (e.g. 1): ";
my $min_ill = <STDIN>;
chomp $min_ill;

print "------------------------------------------------\n";

#####----------SMRT FILE PROCESSING-------------######
system("sort -k 3,3 -k 4,4n \Q$SMRT_file\E > \Q$SMRT_file\E.sorted.temp");
system("awk '\$2==0' \Q$SMRT_file\E.sorted.temp > \Q$SMRT_file\E.sorted.plus.sam.temp");
system("awk '\$2==16' \Q$SMRT_file\E.sorted.temp > \Q$SMRT_file\E.sorted.minus.sam.temp");
system("rm \Q$SMRT_file\E.sorted.temp");

#processing of PLUS sam file
open(INF, "<$SMRT_file.sorted.plus.sam.temp") or die "couldn't open file";
open(OUT1, ">$SMRT_file.sorted.plus.sam.soft_clipped_reads.sam.temp") or die "couldn't open file";
open(OUT2, ">$SMRT_file.sorted.plus.sam.non_clipped_reads.sam.temp") or die "couldn't open file";
open(OUT3, ">$SMRT_file.sorted.plus.sam.read_ends.wig.temp") or die "couldn't open file";

my @dist;
my $sum;
my %plus_ends;

print "Processing SMRT plus strand reads...\n";

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    next if $cols[2] ne $viral_chr;
    if ($cols[5] =~ m/\d+S$/) {	#removes reads soft-clipped at the 3' end
        print OUT1 $line, "\n";
    }
    else {
        print OUT2 $line, "\n";
        while ($cols[5] =~ /(\d+)[DMNX=]/g) { #these lines use the CIGAR string to determine the downstream coordinate
                push (@dist, $1);
        }
        $sum += $_ for @dist;
        my $end_coord = $cols[3] + $sum - 1;
        my $chr_end_coord = "$cols[2]\:$end_coord"; #combines the chromosome and 3' end coordinate into a key to use for the hash
        $sum = 0;
        @dist = ();
        my @split_id = split("\/", $cols[0]); #extracts the read depth for this putative isoform from its id			
        
        if (exists $plus_ends{$chr_end_coord}) { #if the key is already in the hash, increases the value (count) by 1
            $plus_ends{$chr_end_coord} = $plus_ends{$chr_end_coord} + $split_id[1];	
        }
        else {
            $plus_ends{$chr_end_coord} = $split_id[1]; #if the key is not already in the hash, adds it with a value (count) of the read depth
        }
    }
}

foreach my $chr_end_coord (sort keys %plus_ends) { #prints out a(n inadequately) sorted wiggle file
    my @split_keys = split("\:", $chr_end_coord); 
    print OUT3 "$split_keys[0]\t$split_keys[1]\t$split_keys[1]\t$plus_ends{$chr_end_coord}\n";
}	
close(INF);
close(OUT1);
close(OUT2);
close(OUT3);

system("rm \Q$SMRT_file\E.sorted.plus.sam.soft_clipped_reads.sam.temp");
system("rm \Q$SMRT_file\E.sorted.plus.sam.non_clipped_reads.sam.temp");
system("rm \Q$SMRT_file\E.sorted.plus.sam.temp");

#processing of MINUS sam file
open(INF, "<$SMRT_file.sorted.minus.sam.temp") or die "couldn't open file";
open(OUT1, ">$SMRT_file.sorted.minus.sam.soft_clipped_reads.sam.temp") or die "couldn't open file";
open(OUT2, ">$SMRT_file.sorted.minus.sam.non_clipped_reads.sam.temp") or die "couldn't open file";
open(OUT3, ">$SMRT_file.sorted.minus.sam.read_ends.wig.temp") or die "couldn't open file";

my $previous_coordinate=1;
my $count=0;
my $previous_chr = "start";
print "Processing SMRT minus strand reads...\n";
while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    next if $cols[2] ne $viral_chr;
    if ($cols[5] =~ m/^\d+S/) { #removes reads clipped at the 3' end
        print OUT1 $line, "\n";
    }
    else {
        print OUT2 $line, "\n";
        #uses Erik's counting method that should take less memory than the hash method, but only works for coordinates that can be obtained directly from the sorted sam file without having to calculate from the CIGAR string
        my @split_id = split("\/", $cols[0]); #extracts the read depth for this putative isoform from its id
        if ($cols[2] eq $previous_chr && $cols[3] == $previous_coordinate) { 
            $count = $count + $split_id[1]; #increases the count by the read depth for the putative isoform		
        }
        else {
            if ($previous_chr eq "start") { #doesn't print out the placeholder first line.
                
                $previous_chr = $cols[2];	#sets the previous chromosome, previous coordinate and count values			
                $previous_coordinate = $cols[3];				
                $count = $split_id[1];
            }
            else {
                print OUT3 $previous_chr, "\t", $previous_coordinate, "\t", $previous_coordinate, "\t-", $count, "\n"; #prints to output file
                $previous_chr = $cols[2];				
                $previous_coordinate = $cols[3];				
                $count = $split_id[1];
            }
        }
    }
}	
print OUT3 $previous_chr, "\t", $previous_coordinate, "\t", $previous_coordinate, "\t-", $count, "\n"; #prints the last start coordinates to output file
close(INF);
close(OUT1);
close(OUT2);
close(OUT3);

system("cat \Q$SMRT_file\E.sorted.plus.sam.read_ends.wig.temp \Q$SMRT_file\E.sorted.minus.sam.read_ends.wig.temp | sort -k2,3n > \Q$SMRT_file\E.\Q$viral_chr\E.all_read_ends.wig.noheader");

system("rm \Q$SMRT_file\E.sorted.plus.sam.read_ends.wig.temp");
system("rm \Q$SMRT_file\E.sorted.minus.sam.read_ends.wig.temp");
system("rm \Q$SMRT_file\E.sorted.minus.sam.soft_clipped_reads.sam.temp");
system("rm \Q$SMRT_file\E.sorted.minus.sam.non_clipped_reads.sam.temp");
system("rm \Q$SMRT_file\E.sorted.minus.sam.temp");

#add header to wiggle file
open(INF, "<$SMRT_file.$viral_chr.all_read_ends.wig.noheader") or die "couldn't open file";
open(OUT, ">$SMRT_file.$viral_chr.all_read_ends.wig") or die "couldn't open file";

print OUT "track type=wiggle name=\"$SMRT_file.$viral_chr.all_read_ends.wig\" description=\"3' ends of SMRT reads from end_finder_sam_to_bed.pl\"\n";
while (my $line = <INF>) {
    print OUT $line;
}
close(OUT);
close(INF);

system("rm \Q$SMRT_file\E.\Q$viral_chr\E.all_read_ends.wig.noheader");

#make a bed file from the SMRT wiggle file:
open(INF, "<$SMRT_file.$viral_chr.all_read_ends.wig") or die "couldn't open file";
open(OUT, ">$SMRT_file.ends.temp.bed") or die "couldn't open file";

print "Combining SMRT 3' ends within $distance_between_SMRT_peaks of each other and calculating consensus 3' ends...\n";
collapse_wiggle($distance_between_SMRT_peaks);

close(INF);
close(OUT);

system("sort -k 1,1 -k 2,2n \Q$SMRT_file\E.ends.temp.bed > \Q$SMRT_file\E.ends.bed.noheader");
system("rm \Q$SMRT_file.ends.temp.bed\E");

#add header to bed file
open(INF, "<$SMRT_file.ends.bed.noheader") or die "couldn't open file";
open(OUT, ">$SMRT_file.$viral_chr.ends.bed") or die "couldn't open file";

print OUT "track type=bed name=\"$SMRT_file.$viral_chr.ends.bed\" description=\"consensus 3' ends of SMRT reads within $distance_between_SMRT_peaks bp collapsed to weighted center from end_finder_sam_to_bed.pl\"\n";
while (my $line = <INF>) {
    print OUT $line;
}
close(OUT);
close(INF);

system("rm \Q$SMRT_file\E.ends.bed.noheader");

#####----------ILLUMINA FILE PROCESSING-------------######

open(INF, "<$ill_file") or die "couldn't open input file";
open(OUT, ">$ill_file.polyA_ends.temp") or die "couldn't open output file";

print "Extracting Illumina reads with at least $min_As As and at least $min_softclip mismatches...\n";

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    next if ($cols[0] eq "\@HD" || $cols[0] eq "\@PG" || $cols[0] eq "\@SQ"); #skips SAM file header lines
    next if $cols[2] ne $viral_chr;
    if ($cols[1] == 81 || $cols[1] == 83 || $cols[1] == 89 || $cols[1] == 16) {  #selects reads with FLAG codes indicating they are first in pair on the plus strand
        if (($cols[5] =~ m/\d+S$/) and ($cols[9] =~ m/A{$min_As}$/)) { # selects reads with softclipping and a run of As at the end
            my ($softclips) = $cols[5] =~ m/(\d+)S$/; #pulls out the number of softclipped bases
            if ($softclips > $min_softclip) { #selects reads with at least the specified number of softclipped bases
                print OUT $line, "\n";
            }
        }
    }
    elsif ($cols[1] == 73 || $cols[1] == 97 || $cols[1] == 99 || $cols[1] == 0) {  #selects reads with FLAG codes indicating they are first in pair on the minus strand
        if (($cols[5] =~ m/^\d+S/) and ($cols[9] =~ m/^T{$min_As}/)) { #selects reads with softclipping and a run of Ts at the beginning
            my ($softclips) = $cols[5] =~ m/^(\d+)S/; #pulls out the number of softclipped bases
            if ($softclips > $min_softclip) { #selects reads with at least the specified number of softclipped bases
                print OUT $line, "\n";
            }
        }
    }
}

close(INF);
close(OUT);

system ("sort -k 3,3 -k 4,4n \Q$ill_file\E.polyA_ends.temp > \Q$ill_file\E.polyA_ends.sam");
system ("rm \Q$ill_file\E.polyA_ends.temp");

open(INF, "<$ill_file.polyA_ends.sam") or die "couldn't open file";
open(OUT, ">$ill_file.polyA_sites.temp") or die "couldn't open file";

print "Processing Illumina reads with polyA tails...\n";
#create a file with the corresponding to the polyA ends of the reads, and sort it by those coordinates

my $cigar_sum;
my $cigar_calc;
my @plus_ends;
my @read_dist;

while (my $line = <INF>) {
	my @cols = split("\t", $line);
	if ($cols[1] == 73 || $cols[1] == 97 || $cols[1] == 99 || $cols[1] == 0) {
		print OUT "$cols[2]\t$cols[3]\t0\n";
	}
	elsif ($cols[1] == 81 || $cols[1] == 83 || $cols[1] == 89 || $cols[1] == 16) {
		while ($cols[5] =~ /(\d+)[DMNX=]/g) {
            push (@read_dist, $1);
		}
		$cigar_sum += $_ for @read_dist;
		$cigar_calc = $cols[3] + $cigar_sum - 1;
		$cigar_sum = 0;
		@read_dist = ();
		print OUT "$cols[2]\t$cigar_calc\t1\n";
	}
}
close(INF);
close(OUT);

system("sort -k 1,1 -k 2,2n \Q$ill_file.polyA_sites.temp\E > \Q$ill_file.polyA_sites.temp\E.sorted");

#create a wiggle file from the sorted coordinates file

open(INF, "<$ill_file.polyA_sites.temp") or die "couldn't open file";
open(OUT, ">$ill_file.polyA_sites.temp.wig") or die "couldn't open file";

my $chrom_minus;
my $previous_coordinate_m=0;
my $count_m=0;
my $coordinate_m;
my $chrom_plus;
my $previous_coordinate_p=0;
my $count_p=0;
my $coordinate_p;

while (my $line = <INF>) {
	
	my @cols = split("\t", $line);
	
	#reads on the plus strand:
	if ($cols[2] == 1) {
		if ($chrom_plus) { #if $chrom_plus has been defined (i.e. there is a previous plus strand read)
			if (($cols[0] eq $chrom_plus) and ($cols[1] == $previous_coordinate_p)) {
				$count_p++;
			}
			else {
				$coordinate_p = $previous_coordinate_p - 1;
				print OUT "$chrom_plus\t$coordinate_p\t$coordinate_p\t$count_p\n";
				$previous_coordinate_p = $cols[1];
				$count_p = 1;
			}
		}
		else { #if $chrom_plus has not been defined (i.e. there is no previous plus strand read)
			$chrom_plus = $cols[0];
			$previous_coordinate_p = $cols[1];
			$count_p = 1;
		}
	}
	
	#reads on the minus strand:
	elsif ($cols[2] == 0) {
		if ($chrom_minus) {
			if (($cols[0] eq $chrom_minus) and ($cols[1] == $previous_coordinate_m)) {
				$count_m++;
			}
			else {
				$coordinate_m = $previous_coordinate_m - 1;
				print OUT "$chrom_minus\t$coordinate_m\t$coordinate_m\t-$count_m\n";
				$chrom_minus = $cols[0];
				$previous_coordinate_m = $cols[1];
				$count_m = 1;
			}
		}
		else {
			$chrom_minus = $cols[0];
			$previous_coordinate_m = $cols[1];
			$count_m = 1;
		}
	}
}

$coordinate_p = $previous_coordinate_p - 1;
print OUT "$chrom_plus\t$coordinate_p\t$coordinate_p\t$count_p\n";
$coordinate_m = $previous_coordinate_m - 1;
print OUT "$chrom_minus\t$coordinate_m\t$coordinate_m\t-$count_m\n";

close(INF);
close(OUT);

system("sort -k 1,1 -k 2,2n \Q$ill_file\E.polyA_sites.temp.wig > \Q$ill_file\E.polyA_sites.wig.noheader");
system("rm \Q$ill_file\E.polyA_sites.temp.wig");
system("rm \Q$ill_file\E.polyA_sites.temp.sorted");
system("rm \Q$ill_file\E.polyA_sites.temp");

#add header to wiggle file
open(INF, "<$ill_file.polyA_sites.wig.noheader") or die "couldn't open file";
open(OUT, ">$ill_file.$viral_chr.polyA_sites.wig") or die "couldn't open file";

print OUT "track type=wiggle name=\"$ill_file.$viral_chr.polyA_sites.wig\" description=\"polyA sites in Illumina reads with at least 5As and at least 2 mismatches from end_finder_sam_to_bed.pl\"\n";
while (my $line = <INF>) {
    print OUT $line;
}
close(OUT);
close(INF);

system("rm \Q$ill_file\E.polyA_sites.wig.noheader");

#make a bed file from the Illumina wiggle file:
open(INF, "<$ill_file.$viral_chr.polyA_sites.wig") or die "couldn't open file";
open(OUT, ">$ill_file.$viral_chr.polyA_sites.temp.bed") or die "couldn't open file";

print "Combining Illumina polyA tails within $distance_between_ill_peaks of each other and calculating consensus 3' ends...\n";
collapse_wiggle($distance_between_ill_peaks);

close(INF);
close(OUT);

system("sort -k 1,1 -k 2,2n \Q$ill_file\E.\Q$viral_chr\E.polyA_sites.temp.bed > \Q$ill_file\E.\Q$viral_chr\E.polyA_sites.bed.noheader");
system("rm \Q$ill_file\E.\Q$viral_chr\E.polyA_sites.temp.bed\E");

#add header to bed file
open(INF, "<$ill_file.$viral_chr.polyA_sites.bed.noheader") or die "couldn't open file";
open(OUT, ">$ill_file.$viral_chr.polyA_sites.bed") or die "couldn't open file";

print OUT "track type=bed name=\"$ill_file.$viral_chr.polyA_sites.bed\" description=\"consensus polyA sites of Illumina reads with tails of 5 As with 2 mismatches within $distance_between_ill_peaks bp collapsed to weighted centers from end_finder_sam_to_bed.pl\"\n";
while (my $line = <INF>) {
    print OUT $line;
}
close(OUT);
close(INF);

system("rm \Q$ill_file\E.\Q$viral_chr\E.polyA_sites.bed.noheader");

#####----------SEEKING ILLUMINA SUPPORT FOR SMRT ENDS-------------######

open(INF, "<$ill_file.$viral_chr.polyA_sites.bed" ) or die "couldn't open file";

print "Extracting SMRT 3' ends within $dist_SMRT_ill bases of Illumina polyA tails...\n";

my %features_ill;

while(my $line = <INF> ) {
	chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
	my @cols = split("\t", $line);
	my $key_combo_ill = "$cols[0]:$cols[1]:$cols[2]:$cols[5]"; #for each line in the Illumina polyA reads bed file, creates a key for the hash combining chromosome, start coordinate, end coordinate and strand
	$features_ill{$key_combo_ill} = $cols[4]; #enters a count value for the key into the hash
}

close(INF);

open(INF, "<$SMRT_file.$viral_chr.ends.bed" ) or die "couldn't open file";
open(OUT, ">$SMRT_file.$viral_chr.ends.bed.illumina_support.bed");

print OUT "track type=bed name=\"$SMRT_file.$viral_chr.ends.bed.illumina_support.bed\" description=\"consensus SMRT 3' ends of collapse value 8 supported within $dist_SMRT_ill bp by Illumina polyA sites of 5As, 2 mismatches, collapse window 8 from end_finder_sam_to_bed.pl\"\n";

while(my $line = <INF>) {
	chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
	my @SMRT_cols = split("\t", $line);
    next if (abs $SMRT_cols[4] < $min_SMRT);
    foreach my $key_combo_ill (keys %features_ill) {
        my @ill_cols = split(":", $key_combo_ill);
        next if (abs $features_ill{$key_combo_ill} < $min_ill);
        my $lower_limit = $SMRT_cols[1]-$dist_SMRT_ill;
        my $upper_limit = $SMRT_cols[1]+$dist_SMRT_ill;
        if (($SMRT_cols[5] eq $ill_cols[3]) and ($ill_cols[1] > $lower_limit) and ($ill_cols[1]<$upper_limit)) {
            my $count = $features_ill{$key_combo_ill} + $SMRT_cols[4];
            print OUT "$SMRT_cols[0]\t$SMRT_cols[1]\t$SMRT_cols[2]\t$SMRT_cols[4]SMRT_$features_ill{$key_combo_ill}Ill\t$count\t$SMRT_cols[5]\n";
            last;
        }
    }
}

close(OUT);
close(INF);

#########################
sub collapse_wiggle {
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
                    $weighted_average_plus = $weighted_coordinate_sum_plus/$count_sum_plus; #calculates weighted average
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
                    $weighted_average_minus = $weighted_coordinate_sum_minus/$count_sum_minus; #calculates weighted average
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
        $weighted_average_plus = $weighted_coordinate_sum_plus/$count_sum_plus;
        $chrStart_plus = $coords_plus[0];
        $chrEnd_plus = pop(@coords_plus);
        printf OUT "%s\t%1.0f\t%1.0f\t%d%s%d%s%d\t%d\t%s\n", $viral_chr, $weighted_average_plus, $weighted_average_plus, $chrStart_plus, ":", $chrEnd_plus, ":", $count_sum_plus, $count_sum_plus, "+";
    }
    
    if ($count_sum_minus < 0) {#calculates and prints out weighted average for the last feature (minus strand)
        $weighted_average_minus = $weighted_coordinate_sum_minus/$count_sum_minus;
        $chrStart_minus = $coords_minus[0];
        $chrEnd_minus = pop(@coords_minus);
        printf OUT "%s\t%1.0f\t%1.0f\t%d%s%d%s%d\t%d\t%s\n", $viral_chr, $weighted_average_minus, $weighted_average_minus, $chrStart_minus, ":", $chrEnd_minus, ":", $count_sum_minus, $count_sum_minus, "-";
    }
}