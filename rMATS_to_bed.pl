#!/usr/bin/perl

#Converts rMATS ouput files to bed format

#USAGE:
# perl PATH/rMATS_to_bed.pl /PATH/inputfile

use warnings;
use strict;

my ($file) = @ARGV;

open(INF, "<$file");
open(OUT, ">$file.temp");

while(my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    next if ($cols[0] eq "ID");
    print OUT "$cols[3]\t$cols[5]\t$cols[6]\t$cols[1]\:$cols[18]\:$cols[19]\:$cols[22]\t$cols[22]\t$cols[4]\n";
}

close(INF);
close(OUT);

system("sort -k1,1 -k2,3n $file.temp > $file.bed");
system("rm $file.temp");