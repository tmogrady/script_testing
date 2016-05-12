# !/usr/bin/perl

use warnings;
use strict;

my ($regex, $file) = @ARGV;

open (INF, "<$file");
my $flag = 0;
my @UTR_start;

while (my $line = <INF>) {
    chomp($line);
    if ($line =~ m/^\>/) {
        my @cols = split(/\|/, $line);
        if ($cols[4]){
            my $gene = "$cols[1]:$cols[2]";
            my $chr = $cols[3];
            push (@UTR_start, $cols[4]);
            $flag = 0;
            print "$gene\t$chr\t$UTR_start[0]\n";
        }
        else {
            $flag = 1;
        }
    }
    else {
        next if ($flag == 1);
        my @pos = match_all_positions($regex, $line);
        print "@pos\n";
        foreach my $start (@pos) {
            my $coord = $start + $UTR_start[0];
            print "$coord\n";
        }
    }
}

close(INF);

sub match_all_positions {
    my ($regex, $string) = @_;
    my @ret;
    while ($string =~ /$regex/g) {
        push @ret, ( $-[0] );
    }
    return @ret;
}

#sub match_positions {
#    my ($regex, $string) = @_;
#    return if not $string =~ /$regex/;
#    return ($-[0], $+[0]);
#}
