# !/usr/bin/perl

#usage: perl <script.pl> <MAF file> <BED file used to generate MAF file>

use warnings;
use strict;

my ($file, $file2) = @ARGV;

open (INF, "<$file");
open (OUT, ">$file.conservation.temp");

my $chr;
my $chrStart;
my $chrEnd;
my $seq;
my $mm10 = 0;
my $hg19 = 0;
my $rn5 = 0;
my $canFam3 = 0;
my $galGal4 = 0;
my $first = 0;

while (my $line = <INF>) {
    chomp($line);
    if ($line =~ /^s/) {
        my @cols = split /\s+/, $line;
        #print "$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$cols[5]\t$cols[6]\n";
        if ($cols[1] =~ /mm10/) {
            if ($first == 1) {
                print OUT "$chr\t$chrStart\t$chrEnd\t$mm10\t$hg19\t$rn5\t$canFam3\t$galGal4\n";
                $hg19 = 0;
                $rn5 = 0;
                $canFam3 = 0;
                $galGal4 = 0;
            }
            ($chr) = ($cols[1] =~ /(chr.+)$/);
            $chrStart = $cols[2];
            $chrEnd = $chrStart + $cols[3];
            $seq = $cols[6];
            $mm10 = "mm10";
            $first = 1;
        }
        elsif ($cols[1] =~ /hg19/){
            if ($cols[6] eq $seq) {
                $hg19 = "hg19";
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

print OUT "$chr\t$chrStart\t$chrEnd\t$mm10\t$hg19\t$rn5\t$canFam3\t$galGal4\n";

close (INF);
close (OUT);

my %coords_names;

open (INF, "<$file2");

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    my $coord = "$cols[0]\t$cols[1]";
    $coords_names{$coord} = $cols[3];
}

close (INF);

#foreach my $key (keys %coords_names) {
#    print "$key\t$coords_names{$key}\n";
#}

open (INF, "<$file.conservation.temp");
open (OUT, ">$file.conservation.bed");

my $prev_chr = 0;
my $prev_chrStart = 0;
my $prev_chrEnd = 0;
my $prev_hg19 = 0;
my $prev_mm10 = 0;
my $prev_rn5 = 0;
my $prev_canFam3 = 0;
my $prev_galGal4 = 0;
my $coord_check;
my $coord_key;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
#    my $coord_check = "$cols[0]\t$cols[1]";
#    if (exists $coords_names{$coord_check}) {
#        print OUT "$cols[0]\t$cols[1]$cols[2]\t$coords_names{$coord_check}.$cols[3].$cols[4].$cols[5].$cols[6].$cols[7]\n";
#    }
    
    
    if (($cols[0] eq $prev_chr) and ($cols[1] == $prev_chrEnd)) {
        #then it needs to be added to the previous one, unless its the start of a second consecutive motif
        $coord_check = "$cols[0]\t$cols[1]";
        if (exists $coords_names{$coord_check}) {
            $coord_key = "$prev_chr\t$prev_chrStart";
            print OUT "$prev_chr\t$prev_chrStart\t$prev_chrEnd\t$coords_names{$coord_key}.$prev_mm10.$prev_hg19.$prev_rn5.$prev_canFam3.$prev_galGal4\n";
            $prev_chr = $cols[0];
            $prev_chrStart = $cols[1];
            $prev_chrEnd = $cols[2];
            $prev_mm10 = $cols[3];
            $prev_hg19 = $cols[4];
            $prev_rn5 = $cols[5];
            $prev_canFam3 = $cols[6];
            $prev_galGal4 = $cols[7];
        }
        else {
            $prev_chrEnd = $cols[2];
            if ($cols[3] eq $prev_mm10) {
                $prev_mm10 = $cols[3];
            }
            else {
                $prev_mm10 = 0;
            }
            if ($cols[4] eq $prev_hg19) {
                $prev_hg19 = $cols[4];
            }
            else {
                $prev_hg19 = 0;
            }
            if ($cols[5] eq $prev_rn5) {
                $prev_rn5 = $cols[5];
            }
            else {
                $prev_rn5 = 0;
            }
            if ($cols[6] eq $prev_canFam3) {
                $prev_canFam3 = $cols[6];
            }
            else {
                $prev_canFam3 = 0;
            }
            if ($cols[7] eq $prev_galGal4) {
                $prev_galGal4 = $cols[7];
            }
            else {
                $prev_galGal4 = 0;
            }
        }
    }
    else {
        if ($prev_chr eq 0) {
            $prev_chr = $cols[0];
            $prev_chrStart = $cols[1];
            $prev_chrEnd = $cols[2];
            $prev_mm10 = $cols[3];
            $prev_hg19 = $cols[4];
            $prev_rn5 = $cols[5];
            $prev_canFam3 = $cols[6];
            $prev_galGal4 = $cols[7];
        }
        else {
            $coord_check = "$prev_chr\t$prev_chrStart";
            print OUT "$prev_chr\t$prev_chrStart\t$prev_chrEnd\t$coords_names{$coord_check}.$prev_mm10.$prev_hg19.$prev_rn5.$prev_canFam3.$prev_galGal4\n";
            $prev_chr = $cols[0];
            $prev_chrStart = $cols[1];
            $prev_chrEnd = $cols[2];
            $prev_mm10 = $cols[3];
            $prev_hg19 = $cols[4];
            $prev_rn5 = $cols[5];
            $prev_canFam3 = $cols[6];
            $prev_galGal4 = $cols[7];
        }
    }
}
$coord_check = "$prev_chr\t$prev_chrStart";
print OUT "$prev_chr\t$prev_chrStart\t$prev_chrEnd\t$coords_names{$coord_check}.$prev_mm10.$prev_hg19.$prev_rn5.$prev_canFam3.$prev_galGal4\n";

close (INF);
close (OUT);

system("rm $file.conservation.temp");