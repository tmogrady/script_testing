# !/usr/bin/perl

use warnings;
use strict;

my ($regex, $file) = @ARGV;

open (INF, "<$file");
my $flag = 0;

while (my $line = <INF>) {
    chomp($line);
    #next if ($line eq "Sequence unavailable");
    if ($line =~ m/^\>/) {
        my @cols = split(/\|/, $line);
        if ($cols[4]){
            my $gene = "$cols[1]:$cols[2]";
            my $chr = $cols[3];
            my $UTR_start = $cols[4];
            $flag = 0;
            print "$gene\t$chr\t$UTR_start\n";
        }
        else {
            $flag = 1;
        }
    }
    else {
        next if ($flag == 1);
        my @pos = match_all_positions($regex, $line);
        print "@pos\n";
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
