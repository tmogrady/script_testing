# !/usr/bin/perl

use warnings;
use strict;

my ($file) = @ARGV;

open (INF, "<$file");
open (OUT, ">$file.conservation.bed");

my $chr;
my $chrStart;
my $chrEnd;
my $seq;
my $hg19 = 0;
my $mm10 = 0;
my $rn5 = 0;
my $canFam3 = 0;
my $galGal4 = 0;
my $first = 0;

while (my $line = <INF>) {
    chomp($line);
    if ($line =~ /^s/) {
        my @cols = split /\s+/, $line;
        #print "$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$cols[5]\t$cols[6]\n";
        if ($cols[1] =~ /hg19/) {
            if ($first == 1) {
                print "$chr\t$chrStart\t$chrEnd\t$hg19.$mm10.$rn5.$canFam3.$galGal4\n";
                print OUT "$chr\t$chrStart\t$chrEnd\t$hg19.$mm10.$rn5.$canFam3.$galGal4\n";
                $mm10 = 0;
                $rn5 = 0;
                $canFam3 = 0;
                $galGal4 = 0;
            }
            ($chr) = ($cols[1] =~ /(chr.+)$/);
            $chrStart = $cols[2];
            $chrEnd = $chrStart + $cols[3];
            $seq = $cols[6];
            $hg19 = "hg19";
            $first = 1;
        }
        elsif ($cols[1] =~ /mm10/){
            if ($cols[6] eq $seq) {
                $mm10 = "mm10";
            }
            else {
                next;
            }
        }
        elsif ($cols[1] =~ /rn5/){
            if ($cols[6] eq $seq) {
                $rn5 = "rn5";
            }
            else {
                next;
            }
        }
        elsif ($cols[1] =~ /canFam3/){
            if ($cols[6] eq $seq) {
                $canFam3 = "canFam3";
            }
            else {
                next;
            }
        }
        elsif ($cols[1] =~ /galGal4/){
            if ($cols[6] eq $seq) {
                $galGal4 = "galGal4";
            }
            else {
                next;
            }
        }
        else {
            next;
        }

    }
}

print "$chr\t$chrStart\t$chrEnd\t$hg19.$mm10.$rn5.$canFam3.$galGal4\n";
print OUT "$chr\t$chrStart\t$chrEnd\t$hg19.$mm10.$rn5.$canFam3.$galGal4\n";

close (INF);