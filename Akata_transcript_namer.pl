# !/usr/bin/perl
# Given a bed file of transcripts aligned to the inverted Akata genome, names transcripts according to the BamHI fragment (or other supplied fragment scheme) transcription starts in, and the direction

# USAGE:
# perl /PATH/Akata_transcript_namer.pl /PATH/BamHI_fragment_bed_file /PATH/annotation_bed_file /PATH/transcript_bed_file

#TO'G, 12/15/15

use warnings;
use strict;

my ($frag_file, $ann_file, $tran_file) = @ARGV;

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

open(INF, "<$ann_file") or die "couldn't open annotation file";

my %names;
my $name;
my $name2;

while (my $line = <INF>) {
    chomp($line);
#    next if ($line =~ /^track/); #skips the track definition line
#    next if ($line =~ /^LMP/);
#    next if ($line =~ /^EBNA/);
#    next if ($line =~ /^RPMS/);
#    
    if ($line =~ /B\D\DF\d/) {
        if ($line =~ /B\D\DF\d\/B\D\DF\d/) {
            my ($prefix1, $number1, $prefix2, $number2) = $line =~ /(B\D\D)F(\d)\/(B\D\D)F(\d)/;
            $name = "${prefix1}T";
            $name2 = "${prefix2}T";
            if (exists $names{$name}) {
                if ($names{$name} > $number1) {
                    next;
                }
                else {
                    $names{$name} = $number1;
                }
            }
            else {
                $names{$name} = $number1;
            }
            if (exists $names{$name2}) {
                if ($names{$name2} > $number2) {
                    next;
                }
                else {
                    $names{$name2} = $number2;
                }
            }
            else {
                $names{$name2} = $number2;
            }
        }
        else {
            my ($prefix, $number) = $line =~ /(B\D\D)F(\d)/;
            $name = "${prefix}T";
            if (exists $names{$name}) {
                if ($names{$name} > $number) {
                    next;
                }
                else {
                    $names{$name} = $number;
                }
            }
            else {
                $names{$name} = $number;
            }
        }
    }
}

#foreach my $hi (keys %names) {
#    print "$hi $names{$hi}\n";
#}

close(INF);

open(INF, "<$tran_file.sorted.temp") or die "couldn't open transcript file";
open(OUT, ">${tran_file}_named_test.bed") or die "couldn't open output file";

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