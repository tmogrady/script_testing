#!/usr/bin/perl

#Converts rMATS ouput files to bed format

#USAGE:
# perl PATH/rMATS_to_bed.pl /PATH/inputfile

use warnings;
use strict;

my ($file) = @ARGV;

open(INF, "<$file");
open(OUT, ">$file.temp");

my $chr;
my $up_chrStart;
my $up_chrEnd;
my $down_chrStart;
my $down_chrEnd;
my %dup_remover;

while(my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    next if ($cols[0] eq "ID");
    ($chr) = $cols[3] =~ /chr(.+)$/; #gets the chromosome name without the "chr", to make output compatible with Bedtools
    if ($cols[4] eq "+") { #have to treat plus and minus strand separately to get "upstream" and "downstream" right
        $up_chrStart = $cols[5] - 200; #gets coordinates for 200 nt upstream
        my $up_head_start = $up_chrStart + 1; #fixes 1- vs 0-based problems
        print OUT "$chr\t$up_chrStart\t$cols[5]\t$cols[1]\|up_$cols[1]\|$cols[5]-$cols[6]\|$chr\|$up_head_start\|$cols[5]\|1\t$cols[22]\t$cols[4]\n";
        $down_chrEnd = $cols[6] + 200;
        my $down_head_start = $cols[6] + 1;
        print OUT "$chr\t$cols[6]\t$down_chrEnd\t$cols[1]\|down_$cols[1]\|$cols[5]-$cols[6]\|$chr\|$down_head_start\|$down_chrEnd\|1\t$cols[22]\t$cols[4]\n";
    }
    elsif ($cols[4] eq "-") {
        $up_chrEnd = $cols[6] + 200;
        print OUT "$chr\t$cols[6]\t$up_chrEnd\t$cols[2]\|up_$cols[1]\|$cols[5]-$cols[6]\|$chr\|$cols[6]\|$up_chrEnd\|-1\t$cols[22]\t$cols[4]\n";
        $down_chrStart = $cols[5] - 200;
        print OUT "$chr\t$down_chrStart\t$cols[5]\t$cols[2]\|down_$cols[1]\|$cols[5]-$cols[6]\|$chr\|$down_chrStart\|$cols[5]\|-1\t$cols[22]\t$cols[4]\n";
    }
    my $coord_key = "$cols[3]:$cols[5]:$cols[6]:$cols[4]";
    
    if (! exists $dup_remover{$coord_key}) {
        $dup_remover{$coord_key} = $line;
    }
}

close(INF);
close(OUT);

system("sort -k1,1 -k2,3n $file.temp > $file.introns.bed");
system("rm $file.temp");

open(OUT, ">$file.unique.txt");

foreach my $key (keys %dup_remover) {
    print OUT "$dup_remover{$key}\n";
}

close(OUT);