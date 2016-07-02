#!/usr/bin/perl

#Input: weighted average bed files from paraclu_pipe.pl (for 4 replicates). Returns files consensus files of clusters.

#USAGE:
# perl PATH/paraclu_filter.pl numerical_max_distance /PATH/inputfiles

#TO'G 10/28/15

use warnings;
use strict;

my ($dist, $file1, $file2, $file3, $file4) = @ARGV;
my %peaks;

open(INF, "<$file1");

print "Processing first file\n";

my $coord_key;
my $values;
my $prev_plus_chr;
my $prev_plus_coord;
my $coord_sum;
my $prev_plus_count;
my $prev_plus_dens;
my $weighted_dens_sum;
my $prev_plus_peak_count;
my $plus_collapse_count;

while(my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    my @subcols = split("\:", $cols[3]);
    if ($cols[5] eq "+") {
        
        if (($cols[0] eq $prev_chr) and ($cols[1] <= ($prev_coord + $dist))) { #if the peak is less than the specified distance from the previous peak
            
            $prev_plus_count = $prev_plus_count + $cols[4]; #adds tag depths
            $plus_collapse_count = $plus_collapse_count + 1; #will be used to calculate average density & start site
            $prev_plus_dens = (($subcols[3] + ($hashcols[1]*$hashcols[2]))/($hashcols[2]+1));
            my $new_start = (($coord*$hashcols[0]) + ($file_start*$cols[4]))/($hashcols[0] + $cols[4]);
            
        }
        else {
            $coord_key = "$prev_plus_chr\:$prev_plus_coord\:+"; #creates a key for the hash with chr, chrStart and strand
            $values = "$prev_plus_count\:$prev_plus_dens\:1"; #creates a value for the hash with tag count, relative density and replicate count
            $peaks{$coord_key} = $values;
            
            $prev_plus_chr = $cols[0];
            $prev_plus_coord = $cols[1];
            $prev_plus_count = $cols[4];
            $prev_plus_dens = $subcols[3];
        }
    }
    if ($cols[5] eq "-") {
        $coord_key = "$cols[0]\:$cols[2]\:$cols[5]"; #creates a key for the hash with chr, chrEnd and strand
        $values = "$cols[4]\:$subcols[3]\:1"; #creates a value for the hash with tag count, relative density and replicate count
        $peaks{$coord_key} = $values;
    }
}
$coord_key = "$prev_plus_chr\:$prev_plus_coord\:+"; #creates a key for the hash with chr, chrStart and strand
$values = "$prev_plus_count\:$prev_plus_dens\:1"; #creates a value for the hash with tag count, relative density and replicate count
$peaks{$coord_key} = $values;


#foreach (sort keys %peaks) {
#    print "$_: $peaks{$_}\n";
#}

close(INF);

open(INF, "<$file2");

print "Processing second file\n";

my $file_start;

while(my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    my @subcols = split("\:", $cols[3]);
    my $found_flag = 0;
    if ($cols[5] eq "+") {
        $file_start = $cols[1];
    }
    if ($cols[5] eq "-") {
        $file_start = $cols[2];
    }
    my $lower_limit = $file_start - $dist;
    my $upper_limit = $file_start + $dist;
    foreach my $hash_start (keys %peaks) {
        my ($hash_chr, $coord, $strand) = split("\:", $hash_start);
        if ((($hash_chr eq $cols[0]) and ($strand eq $cols[5])) and (($coord >= $lower_limit) and ($coord <= $upper_limit))) { #checks to see if the cluster is already in the hash (from the first file)
            my @hashcols = split("\:", $peaks{$hash_start}); #extracts information about the cluster from the hash
            my $depth = $cols[4] + $hashcols[0]; #adds tag depths from the current file and the hash
            #my $reps = $hashcols[2]+1; #adds one to the replicate count
            my $density = (($subcols[3] + ($hashcols[1]*$hashcols[2]))/($hashcols[2]+1));
            my $new_start = (($coord*$hashcols[0]) + ($file_start*$cols[4]))/($hashcols[0] + $cols[4]);

            $coord_key = "$cols[0]\:$new_start\:$cols[5]";
            $values = "$depth\:$density\:2";
            delete $peaks{$hash_start};
            $peaks{$coord_key} = $values;
            $found_flag = 1;
            last;
        }
    }
    if ($found_flag == 0) {
        $coord_key = "$cols[0]\:$file_start\:$cols[5]";
        $values = "$cols[4]\:$subcols[3]\:1";
        $peaks{$coord_key} = $values; #adds the coordinate to the hash as a key, with the value assigned above as a value
    }
    
}

close(INF);

open(OUT, ">any_replicate.temp");
open(OUT2, ">both_replicates.temp");
open(OUT5, ">weird_replicates.temp");

my $chrStart;
my $chrEnd;
my $out_line;

foreach my $hash_start (sort keys %peaks) {
    my ($chr, $coord, $strand) = split("\:", $hash_start);
    my @hashcols = split("\:", $peaks{$hash_start});
    #print "$chr\t$coord\t$strand\n";
    if ($strand eq "+") {
        $chrStart = sprintf("%1.0f", $coord);
        $chrEnd = $chrStart + 1;
        $out_line = "$chr\t$chrStart\t$chrEnd\t$peaks{$hash_start}\t$hashcols[0]\t$strand\n";
    }
    if ($strand eq "-") {
        $chrEnd = sprintf("%1.0f", $coord);
        $chrStart = $chrEnd - 1;
        $out_line = "$chr\t$chrStart\t$chrEnd\t-$peaks{$hash_start}\t$hashcols[0]\t$strand\n";
    }
    print OUT $out_line;
    if ($hashcols[2] == 2) {
        print OUT2 $out_line;
    }
    if ($hashcols[2] > 2) {
        print OUT5 $out_line;
    }
}

system("sort -k1,1 -k2n any_replicate.temp > any_replicate.bed");
system("sort -k1,1 -k2n both_replicates.temp > both_replicates.bed");
system("sort -k1,1 -k2n weird_replicates.temp > weird_replicates.bed");

system("rm any_replicate.temp");
system("rm both_replicates.temp");
system("rm weird_replicates.temp");

close(OUT);
close(OUT2);
close(OUT5);