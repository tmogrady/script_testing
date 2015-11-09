#!/usr/bin/perl

#Takes a SAM file of deepCAGE (or other Illumina) data and counts reads with 5' ends at each genomic position to prepare an input file for Paraclu
#then runs paraclu
#then filters paraclu output and converts to a bed file
#Strand assignment for paired-end sequencing assumes TruSeq stranded prep. Reports the 5' end of each READ, not each fragment. Includes only reads flagged as primary alignments.

#USAGE:
# perl PATH/paraclu_pipe.pl /PATH/inputfile(s).sam

#TO'G 10/9/2015

use warnings;
use strict;

die "USAGE: 'perl <PATH/paraclu_filter.pl> </PATH/sam_file.sam>'" unless @ARGV == 1;

my ($file) = @ARGV;

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

print "preparing input file...\n";

system("sort -k 3,3 -k 4,4n \Q$file\E > \Q$file\E.sorted.temp");
system("awk '\$2==0 \|\| \$2==81 \|\| \$2==83 \|\| \$2==89 \|\| \$2==137 \|\| \$2==161 \|\| \$2==163' \Q$file\E.sorted.temp > \Q$file\E.sorted.plus.sam.temp");
system("awk '\$2==16 \|\| \$2==73 \|\| \$2==97 \|\| \$2==99 \|\| \$2==145 \|\| \$2==147 \|\| \$2==153' \Q$file\E.sorted.temp > \Q$file\E.sorted.minus.sam.temp");

#processing of plus sam file

print "processing plus strand...\n";

open(INF, "<$file.sorted.plus.sam.temp") or die "couldn't open file";
open(OUT, ">$file.read_starts.txt") or die "couldn't open file";

my $previous_coordinate=1;
my $count=0;
my $previous_chr = "start";

while (my $line = <INF>) {
    chomp($line);
    next if ($line =~ m/^@/); #skips header lines
    my @cols = split("\t", $line);
        
    if (($cols[2] eq $previous_chr) and ($cols[3] == $previous_coordinate)) {
        $count++; #increases the count by 1
        }
        
    else {
        if ($previous_chr eq "start") { #doesn't print out the placeholder first line.
            $previous_chr = $cols[2];	#sets the previous chromosome, previous coordinate and count values
            $previous_coordinate = $cols[3];
            $count = 1;
        }
            
        else {
            print OUT $previous_chr, "\t+\t", $previous_coordinate, "\t", $count, "\n"; #prints to output file
            $previous_chr = $cols[2];
            $previous_coordinate = $cols[3];
            $count = 1;
        }
    }
}
print OUT $previous_chr, "\t+\t", $previous_coordinate, "\t", $count, "\n"; #prints the last start coordinates to output file
close(INF);
#close(OUT);

system("rm \Q$file\E.sorted.plus.sam.temp");


#processing of MINUS sam file

print "processing minus strand...\n";

open(INF, "<$file.sorted.minus.sam.temp") or die "couldn't open file";
#open(OUT, ">$file.sorted.minus.sam.read_starts.txt.temp") or die "couldn't open file";

my @dist;
my $sum;
my %minus_starts;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    
    while ($cols[5] =~ /(\d+)[DMNX=]/g) { #these lines use the CIGAR string to determine the downstream coordinate
        push (@dist, $1);
    }
    $sum += $_ for @dist;
    my $start_coord = $cols[3] + $sum - 1;
        
    my $chr_start_coord = "$cols[2]\:$start_coord"; #combines the chromosome and 5' end coordinate into a key to use for the hash
        
    $sum = 0;
    @dist = ();
        
    if (exists $minus_starts{$chr_start_coord}) { #if the key is already in the hash, increases the value (count) by the read depth for that putative isoform
        $minus_starts{$chr_start_coord} = $minus_starts{$chr_start_coord} + 1;
    }
        
    else {
        $minus_starts{$chr_start_coord} = 1; #if the key is not already in the hash, adds it with a value (count) of the read depth for that putative isoform
    }
}

foreach my $chr_start_coord (sort keys %minus_starts) { #prints out a(n inadequately) sorted file
    my @split_keys = split("\:", $chr_start_coord);
    print OUT "$split_keys[0]\t-\t$split_keys[1]\t$minus_starts{$chr_start_coord}\n";
}
close(INF);
close(OUT);

#system("sort -k 1,1 -k 3,3n \Q$file\E.sorted.minus.sam.read_starts.txt.temp > \Q$file\E.minus_read_starts.txt");
#system("rm \Q$file\E.sorted.minus.sam.read_starts.txt.temp");
system("rm \Q$file\E.sorted.minus.sam.temp");
    system("rm \Q$file\E.sorted.temp");

open(INF, "<$file.read_starts.txt") or die "couldn't open file";
open(OUT, ">$file.paraclu.txt.temp");

# paraclu.pl: perform parametric clustering of data attached to sequences

# Written by Martin C Frith 2006
# Genome Exploration Research Group, RIKEN GSC and
# Institute for Molecular Bioscience, University of Queensland

# This program reads in a list of numeric values attached to positions
# in sequences. The list should have four tab- (or space-) separated
# columns containing: the sequence name, the strand, the position, and
# the value. (Multiple values for the same sequence/strand/position
# will be summed.) It outputs the clusters as eight tab-separated
# columns: sequence name, strand, start, end, number of values, sum of
# values, min d, max d. See below for the meaning of "d".

# An example line of input:
# chr1    +       17689   3
# Clustering is performed separately for different strands (as if each
# strand were a completely different sequence).  It does not matter
# whether the position uses 0-based or 1-based coordinates: the
# program does not care, and the output will be consistent with the
# input.

# The clusters are defined as follows. A cluster is a maximal scoring
# segment, where the score of any segment is: the sum of the values in
# the segment minus d times the size of the segment. Large values of d
# give smaller, tighter clusters and small values of d give larger,
# looser clusters. The program finds all possible clusters for any
# value of d, and annotates each cluster with the maximum and minimum
# values of d that produce it. The ratio max d / min d provides a
# measure of the cluster's "stability".

# The output will include two types of obvious/trivial/degenerate
# clusters: those that cover single positions, and those that cover
# all of the positions in a sequence.  For many purposes, it would be
# best to ignore these cases.

use strict;
use List::Util qw(min max);

my %data;

warn "reading...\n";

while (<INF>) {
    s/#.*//;  # ignore comments
    next unless /\S/;  # skip blank lines
    
    my ($seq, $strand, $pos, $value) = split;
    my $key = "$seq $strand";
    push @{$data{$key}}, [ $pos, $value ];
}

warn "clustering...\n";

print OUT "# sequence, strand, start, end, sites, sum of values, min d, max d\n";

for my $key (sort keys %data) {  # iterate over sequences / strands
    my ($seq, $strand) = split " ", $key;
    my $sites = $data{$key};
    
    @$sites = sort { $$a[0] <=> $$b[0] } @$sites;  # sort by position
    
    my $clusters = all_clusters($sites);
    
    for my $c (@$clusters) {
        my ($beg, $end, $tot, $sit, $min, $max) = @$c;
        my $beg_pos = $$sites[$beg][0];
        my $end_pos = $$sites[$end][0];
        printf OUT "$seq\t$strand\t$beg_pos\t$end_pos\t$sit\t$tot\t%.3g\t%.3g\n",
        $min, $max;
    }
}

### Generic code to find clusters in a sparse sequence of values: ###

sub all_clusters {
    our $inf = 1e100;  # hopefully much bigger than any value in the input
    our $sites = shift;  # input: reference to array of site locations & values
    our $clusters = [];  # output: reference to array of clusters
    get_clusters(0, $#$sites, -$inf);
    return $clusters;
}

# get clusters of sites between beg and end with density > min_density
sub get_clusters {
    our ($clusters, $inf);
    my ($beg, $end, $min_density) = @_;
    
    my ($prefix, $pmin, $ptot, $psit) = weakest_prefix($beg, $end);
    my ($suffix, $smin, $stot, $ssit) = weakest_suffix($beg, $end);
    $ptot == $stot and $psit == $ssit or die "internal error!";
    my $max_density = min $pmin, $smin;
    
    unless ($max_density == $inf) {
        my $break = $pmin < $smin ? $prefix + 1 : $suffix;
        my $new_min = max $min_density, $max_density;
        get_clusters($beg, $break-1, $new_min);
        get_clusters($break, $end, $new_min);
    }
    
    push @$clusters, [ $beg, $end, $ptot, $psit, $min_density, $max_density ]
	if $max_density > $min_density;
}

# get least dense prefix (and total of values & sites)
sub weakest_prefix {
    our ($sites, $inf);
    my ($beg, $end) = @_;
    
    my $beg_pos = $$sites[$beg][0];
    my $min_density = $inf;
    my $min_prefix = $end;
    my $tot = 0;
    my $sit = 0;
    
    for (my $i = $beg; $i < $end; ++$i) {
        $tot += $$sites[$i][1];
        next if $$sites[$i][0] == $$sites[$i+1][0];  # idiot-proofing
        ++$sit;
        my $dist = $$sites[$i+1][0] - $beg_pos;
        my $density = $tot / $dist;
        if ($density < $min_density) {
            $min_prefix = $i;
            $min_density = $density;
        }
    }
    
    $tot += $$sites[$end][1];
    ++$sit;
    return ($min_prefix, $min_density, $tot, $sit);
}

# get least dense suffix (and total of values & sites)
sub weakest_suffix {
    our ($sites, $inf);
    my ($beg, $end) = @_;
    
    my $end_pos = $$sites[$end][0];
    my $min_density = $inf;
    my $min_suffix = $beg;
    my $tot = 0;
    my $sit = 0;
    
    for (my $i = $end; $i > $beg; --$i) {
        $tot += $$sites[$i][1];
        next if $$sites[$i][0] == $$sites[$i-1][0];  # idiot-proofing
        ++$sit;
        my $dist = $end_pos - $$sites[$i-1][0];
        my $density = $tot / $dist;
        if ($density < $min_density) {
            $min_suffix = $i;
            $min_density = $density;
        }
    }
    
    $tot += $$sites[$beg][1];
    ++$sit;
    return ($min_suffix, $min_density, $tot, $sit);
}

close(INF);
close(OUT);

system("sort -k2,2 -k3,3n -k4,4rn \Q$file\E.paraclu.txt.temp > \Q$file\E.paraclu.txt");
system("rm \Q$file\E.paraclu.txt.temp");

#####----------PROCESSING PARACLU OUTPUT-------------######

print "Extracting clusters from CAGE file (Containing $min_tags tags, density fold change at least $min_dens, from $min_length to $max_length bp long)\n";

#filtering clusters:

my $length;
my $dens;
my $prev_start = 0;
my $prev_end = 0;
my $chrStart;
my $chrEnd;

open(INF, "<$file.paraclu.txt") or die "couldn't open file"; #for now, $CAGE_file is the output from Paraclu
open(OUT, ">$file.peaks.$min_tags.$min_dens.$min_length.$max_length.bed") or die "couldn't open file";

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

open(INF, "<$file.peaks.$min_tags.$min_dens.$min_length.$max_length.bed") or die "couldn't open file";
open(OUT, ">$file.peaks_weighted_average.bed") or die "couldn't open file"; #later, reduce output to just one file containing both the weighted average and the cluster extent

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    $chrStart_CAGE = $cols[1];
    $chrEnd_CAGE = $cols[2];
    $strand_CAGE = $cols[5];
    #$tag_depth = 0;
    open(INF2, "<$file.read_starts.txt") or die "couldn't open file";
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