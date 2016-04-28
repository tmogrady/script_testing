#!/usr/bin/perl

#Converts fasta files from Ensembl into files appropriate for the MEME suite. For genes that have multiple features/transcripts present, selects the longest one.

#USAGE:
# perl PATH/ensembl_to_meme.pl /PATH/inputfile

#TO'G

use warnings;
use strict;

my ($file) = @ARGV;

open(INF, "<$file");

my $name;
my %features;

while (my $line = <INF>) {
    chomp($line);
    if ($line =~ m/^\>/) {
        $name = $line;
    }
    else {
        next if ($line eq "Sequence unavailable");
        my $length = length ($line);
        #print $name, "\t", $length, "\n";
        if (exists $features{$name}) {
            if ($length > (length $features{$name})) {
                $features{$name} = $line;
            }
            else {
                next;
            }
        }
        else {
            $features{$name} = $line;
        }
    }
}

close(INF);

open(OUT, ">$file.for_meme.fasta");

foreach my $feature (sort keys %features) {
    print OUT $feature, "\n", $features{$feature}, "\n";
}

close(OUT);