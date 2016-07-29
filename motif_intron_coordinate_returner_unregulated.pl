# !/usr/bin/perl

use warnings;
use strict;

my ($regex, $file) = @ARGV;

open (INF, "<$file");
open (OUT, ">$file.up.$regex.txt");
open (OUT2, ">$file.down.$regex.txt");

my $name;
my $up;

while (my $line = <INF>) {
    chomp($line);
    if ($line =~ m/^\>/) { #determines if this is a header line
        ($name) = $line =~ /\>(.+)/;
        if ($line =~ /up/){
            $up = 1;
        }
        elsif ($line =~ /down/) {
            $up = 0;
        }
    }
    else {
        my @pos = match_all_positions($regex, $line);
        if ($up == 0) {
            foreach my $motif_start (@pos) {
                print OUT2 "$motif_start\t$name\n";
            }
        }
        elsif ($up == 1){
            my $length = length($line);
            foreach my $motif_start (@pos) {
                my $up_start = $length - $motif_start;
                print OUT "-$up_start\t$name\n";
            }
        }
    }
}

close(INF);
close(OUT);
close(OUT2);

sub match_all_positions {
    my ($regex, $string) = @_;
    my @ret;
    while ($string =~ /$regex/g) {
        push @ret, ( $-[0] );
    }
    return @ret;
}