# !/usr/bin/perl

#usage: perl <script.pl> <MAF file> <BED file used to generate MAF file>

use warnings;
use strict;

my ($file, $file2) = @ARGV;

open (INF, "<$file");
open (OUT, ">$file.conservation.temp");
#open (OUT2, ">$file.temp.bed");

my $chr;
my $chrStart;
my $chrEnd;
my $seq;
my $hg38 = 0;
my $mm10 = 0;
my $rn6 = 0;
my $canFam3 = 0;
my $galGal4 = 0;
my $first = 0;

while (my $line = <INF>) {
    chomp($line);
    if ($line =~ /^s/) {
        my @cols = split /\s+/, $line;
        #print "$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$cols[5]\t$cols[6]\n";
        if ($cols[1] =~ /hg38/) {
            if ($first == 1) {
                print OUT "$chr\t$chrStart\t$chrEnd\t$hg38\t$mm10\t$rn6\t$canFam3\t$galGal4\n";
                #print OUT2 "$chr\t$chrStart\t$chrEnd\t$hg38.$mm10.$rn6.$canFam3.$galGal4\n";
                $mm10 = 0;
                $rn6 = 0;
                $canFam3 = 0;
                $galGal4 = 0;
            }
            ($chr) = ($cols[1] =~ /(chr.+)$/);
            $chrStart = $cols[2];
            $chrEnd = $chrStart + $cols[3];
            $seq = $cols[6];
            $hg38 = "hg38";
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
        elsif ($cols[1] =~ /rn6/){
            if ($cols[6] eq $seq) {
                $rn6 = "rn6";
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

print OUT "$chr\t$chrStart\t$chrEnd\t$hg38\t$mm10\t$rn6\t$canFam3\t$galGal4\n";
#print OUT2 "$chr\t$chrStart\t$chrEnd\t$hg38.$mm10.$rn6.$canFam3.$galGal4\n";

close (INF);
close (OUT);
#close (OUT2);

open (INF, "<$file2"); #this section creates a hash of chromosomes and start coordinates in the bed file for use in condensing the blocks in the maf file into single motifs, while still keeping consecutive motifs as separate features

my %coords_names;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    my $coord = "$cols[0]\t$cols[1]";
    $coords_names{$coord} = $cols[3];
}

close (INF);


open (INF, "<$file.conservation.temp");
open (OUT, ">$file.conservation_collapsed.temp");

my $prev_chr = 0;
my $prev_chrStart = 0;
my $prev_chrEnd = 0;
my $prev_hg38 = 0;
my $prev_mm10 = 0;
my $prev_rn6 = 0;
my $prev_canFam3 = 0;
my $prev_galGal4 = 0;
my $coord_check;
my $coord_key;

my %cons;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    if (($cols[0] eq $prev_chr) and ($cols[1] == $prev_chrEnd)) {
        #then it needs to be added to the previous one, unless its the start of a second consecutive motif
        $coord_check = "$cols[0]\t$cols[1]";
        if (exists $coords_names{$coord_check}) { #checks to see if the start coordinate of the feature exists as a motif start coordinate. If so, this is a subsequent occurence of the motif, not a blocksplit part of the previous motif
            print OUT "$prev_chr\t$prev_chrStart\t$prev_chrEnd\t$prev_hg38.$prev_mm10.$prev_rn6.$prev_canFam3.$prev_galGal4\n"; #prints the full motif line to a (possibly unnecessary) file
            $coord_key = "$prev_chr\t$prev_chrStart";
            $cons{$coord_key} = "$prev_hg38.$prev_mm10.$prev_rn6.$prev_canFam3.$prev_galGal4"; #puts the full motif in a hash to check against the input bed file later
            $prev_chr = $cols[0];
            $prev_chrStart = $cols[1];
            $prev_chrEnd = $cols[2];
            $prev_hg38 = $cols[3];
            $prev_mm10 = $cols[4];
            $prev_rn6 = $cols[5];
            $prev_canFam3 = $cols[6];
            $prev_galGal4 = $cols[7];
        }
        else {
            $prev_chrEnd = $cols[2];
            if ($cols[3] eq $prev_hg38) {
                $prev_hg38 = $cols[3];
            }
            else {
                $prev_hg38 = 0;
            }
            if ($cols[4] eq $prev_mm10) {
                $prev_mm10 = $cols[4];
            }
            else {
                $prev_mm10 = 0;
            }
            if ($cols[5] eq $prev_rn6) {
                $prev_rn6 = $cols[5];
            }
            else {
                $prev_rn6 = 0;
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
            $prev_hg38 = $cols[3];
            $prev_mm10 = $cols[4];
            $prev_rn6 = $cols[5];
            $prev_canFam3 = $cols[6];
            $prev_galGal4 = $cols[7];
        }
        else {
            print OUT "$prev_chr\t$prev_chrStart\t$prev_chrEnd\t$prev_hg38.$prev_mm10.$prev_rn6.$prev_canFam3.$prev_galGal4\n";
            $coord_check = "$prev_chr\t$prev_chrStart";
            $cons{$coord_check} = "$prev_hg38.$prev_mm10.$prev_rn6.$prev_canFam3.$prev_galGal4";
            $prev_chr = $cols[0];
            $prev_chrStart = $cols[1];
            $prev_chrEnd = $cols[2];
            $prev_hg38 = $cols[3];
            $prev_mm10 = $cols[4];
            $prev_rn6 = $cols[5];
            $prev_canFam3 = $cols[6];
            $prev_galGal4 = $cols[7];
        }
    }
}
$coord_check = "$prev_chr\t$prev_chrStart";
print OUT "$prev_chr\t$prev_chrStart\t$prev_chrEnd\t$prev_hg38.$prev_mm10.$prev_rn6.$prev_canFam3.$prev_galGal4\n";
$cons{$coord_check} = "$prev_hg38.$prev_mm10.$prev_rn6.$prev_canFam3.$prev_galGal4";

close (INF);
close (OUT);

system("rm $file.conservation.temp");
system("rm $file.conservation_collapsed.temp");

open (INF, "<$file2"); #this section adds the conservation information to the input bed file.
open (OUT, ">$file.conservation.bed");

my $check;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    $check = "$cols[0]\t$cols[1]";
    print OUT "$line:$cons{$check}\n";
}

close(INF);
close(OUT);