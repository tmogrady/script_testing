# !/usr/bin/perl
# Converts an annotation file corresponding to the Akata genome to correspond to the inverted (107955to171322_1to107954) genome

#USAGE: perl /PATH/Akata_inverter.pl /PATH/Akata_file(s).bed

use warnings;
use strict;

foreach my $file(@ARGV) {
    open(INF, "<$file") or die "couldn't open input file";
    open(OUT, ">${file}_inverted") or die "couldn't open output file";
    open(OUT2, ">${file}_split_spanning") or die "couldn't open output file";
    
    my $chrStart;
    my $chrEnd;
    my $thickStart;
    my $thickEnd;
    
    while (my $line = <INF>) {
        chomp($line);
        if ($line =~ /^track/) {
            print OUT "$line\n";
            next;
        }
        my @cols = split("\t", $line);
        if ($cols[1] <= 107954) {
            if ($cols[2] >= 107954) {
                print OUT2 "$line\n";
            }
            else {
                $chrStart = $cols[1] + 63369;
                $chrEnd = $cols[2] + 63369;
                $thickStart = $cols[6] + 63369;
                $thickEnd = $cols[7] + 63369;
                print OUT "chrEBV_Akata_inverted\t$chrStart\t$chrEnd\t$cols[3]\t$cols[4]\t$cols[5]\t$thickStart\t$thickEnd\t$cols[8]\t$cols[9]\t$cols[10]\t$cols[11]\n";
            }
        }
        if ($cols[1] >= 107954) {
            $chrStart = $cols[1] - 107954;
            $chrEnd = $cols[2] - 107954;
            $thickStart = $cols[6] - 107954;
            $thickEnd = $cols[7] - 107954;
            print OUT "chrEBV_Akata_inverted\t$chrStart\t$chrEnd\t$cols[3]\t$cols[4]\t$cols[5]\t$thickStart\t$thickEnd\t$cols[8]\t$cols[9]\t$cols[10]\t$cols[11]\n";
        }
    }
    close(INF);
    close(OUT);
    close(OUT2);
    
    system ("sort -k 2,2n -k 3,3n \Q$file\E_inverted > \Q$file\E_inverted.bed");
    system ("rm \Q$file\E_inverted");
    system ("sort -k 2,2n -k 3,3n \Q$file\E_split_spanning > \Q$file\E_split_spanning.bed");
    system ("rm \Q$file\E_split_spanning");
    
}