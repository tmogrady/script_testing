# !/usr/bin/perl

use warnings;
use strict;

my ($regex, $string) = @ARGV;

my @pos = match_all_positions($regex, $string);
#my @pos2 = match_positions($regex, $string);

print "\nmatch all @pos\n\n";

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
