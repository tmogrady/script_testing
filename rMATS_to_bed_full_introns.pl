#!usr/bin/perl

#takes information from rMATS SE output files to produce a bed file of flanking intron coordinates

#USAGE:
#perl PATH/rMATS_to_bed_full_introns.ok /PATH/input_file

use warnings;
use strict;

my ($file) = @ARGV;

open(INF, "<$file");
open(OUT, ">$file.temp");

my $up_chrStart;
my $up_chrEnd;
my $down_chrStart;
my $down_chrEnd;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    next if ($cols[0] eq "ID");
    my ($chr) = $cols[3] =~ /chr(.+)$/; #gets the chromosome name without the "chr", to make output compatible with Bedtools
    if ($cols[4] eq "+") {#have to treat plus and minus strand separately to get "upstream" and "downstream" right
        $up_chrStart = $cols[8];
        $up_chrEnd = $cols[5];
        $down_chrStart = $cols[6];
        $down_chrEnd = $cols[9];
        if ($cols[22] > 0){
            print OUT "$chr\t$up_chrStart\t$up_chrEnd\t$cols[1]:up_in\t1\t$cols[4]\n";
            print OUT "$chr\t$down_chrStart\t$down_chrEnd\t$cols[1]:down_in\t1\t$cols[4]\n";
        }
        elsif ($cols[22] < 0) {
            print OUT "$chr\t$up_chrStart\t$up_chrEnd\t$cols[1]:up_ex\t-1\t$cols[4]\n";
            print OUT "$chr\t$down_chrStart\t$down_chrEnd\t$cols[1]:down_ex\t\-1\t$cols[4]\n";
        }
    }
    elsif ($cols[4] eq "-") {
        $up_chrStart = $cols[6];
        $up_chrEnd = $cols[9];
        $down_chrStart = $cols[8];
        $down_chrEnd = $cols[5];
        if ($cols[22] > 0) {
            print OUT "$chr\t$up_chrStart\t$up_chrEnd\t$cols[1]:up_in\t1\t$cols[4]\n";
            print OUT "$chr\t$down_chrStart\t$down_chrEnd\t$cols[1]:down_in\t1\t$cols[4]\n";
        }
        elsif ($cols[22] < 0) {
            print OUT "$chr\t$up_chrStart\t$up_chrEnd\t$cols[1]:up_ex\t-1\t$cols[4]\n";
            print OUT "$chr\t$down_chrStart\t$down_chrEnd\t$cols[1]:down_ex\t-1\t$cols[4]\n";
        }
    }
}

close(INF);
close(OUT);

system("sort -k1,1 -k2,3n $file.temp | uniq > $file.flanking_introns.bed");
system("rm $file.temp");