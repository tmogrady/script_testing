#!/usr/bin/perl

#Converts fasta files from Ensembl into files appropriate for the MEME suite. For genes that have multiple features/transcripts present, selects the longest one.

#USAGE:
# perl PATH/ensembl_to_meme.pl /PATH/inputfile

use warnings;
use strict;

my ($file) = @ARGV;

open(INF, "<$file");

my $name;

while (my $line = <INF>) {
    chomp($line);
    if ($line =~ m/^\>/) {
        $name = $line;
    }
    else {
        next if ($line eq "Sequence unavailable");
        my $length = length ($line);
        print $name, "\t", $length, "\n";
    }
}

close(INF);