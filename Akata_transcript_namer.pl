# !/usr/bin/perl
# Given a bed file of transcripts aligned to the inverted Akata genome, names transcripts according to the BamHI fragment (or other supplied fragment scheme) transcription starts in, and the direction

# USAGE:
# perl /PATH/Akata_transcript_namer.pl /PATH/BamHI_fragment_bed_file /PATH/transcript_bed_file

#TO'G, 12/15/15

use warnings;
use strict;

my ($frag_file, $tran_file) = @ARGV;

system ("awk '\$6==\"+\"' \Q$tran_file\E | sort -k 2,2n -k 3,3n > \Q$tran_file\E.plus.temp");
system ("awk '\$6==\"-\"' \Q$tran_file\E | sort -k 3,3rn -k 2,2rn > \Q$tran_file\E.minus.temp");
system ("cat \Q$tran_file\E.plus.temp \Q$tran_file\E.minus.temp > \Q$tran_file\E.sorted.temp");

open(INF, "<$frag_file") or die "couldn't open fragment file";

my %frag_coords;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    my $key = "$cols[1]\:$cols[2]";
    $frag_coords{$key} = $cols[3]; #creates a hash with the fragment names and coordinates
}
close(INF);

open(INF, "<$tran_file.sorted.temp") or die "couldn't open transcript file";
open(OUT, ">${tran_file}_named.bed") or die "couldn't open output file";

my $name;
my %names;

while (my $line = <INF>) {
    chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
    my @cols = split("\t", $line);
    foreach my $key (keys %frag_coords) {
        my ($frag_start, $frag_end) = split(":", $key);
        if ($cols[5] eq "+") {
            if (($cols[1] >= $frag_start) and ($cols[1] <= $frag_end)) {
                $name = "B$frag_coords{$key}RT";
                if (exists $names{$name}) {
                    $names{$name} = $names{$name} + 1;
                }
                else {
                    $names{$name} = 1;
                }
            }
        }
        if ($cols[5] eq "-") {
            if (($cols[2] > $frag_start) and ($cols[2] <= $frag_end)) {
                $name = "B$frag_coords{$key}LT";
                if (exists $names{$name}) {
                    $names{$name} = $names{$name} + 1;
                }
                else {
                    $names{$name} = 1;
                }
            }
        }
    }
    print OUT $cols[0], "\t", $cols[1], "\t", $cols[2], "\t", "$name$names{$name}_$cols[3]", "\t", $cols[4], "\t", $cols[5], "\t", $cols[6], "\t", $cols[7], "\t", $cols[8], "\t", $cols[9], "\t", $cols[10], "\t", $cols[11], "\n";
}
close(INF);
close(OUT);


system ("rm \Q$tran_file\E.plus.temp");
system ("rm \Q$tran_file\E.minus.temp");
system ("rm \Q$tran_file\E.sorted.temp");