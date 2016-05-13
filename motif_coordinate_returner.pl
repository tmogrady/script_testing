# !/usr/bin/perl

use warnings;
use strict;

my ($regex, $file) = @ARGV;

open (INF, "<$file");
open (OUT, ">$file.$regex.temp");

my $flag = 0;
my $name;
my $chr;
my $strand;
my $exon;
my @UTR_starts;
my @UTR_ends;
my $chrStart;
my $chrEnd;
my $motif_length = length($regex); #this probably doesn't work for a lot of regexes

while (my $line = <INF>) {
    chomp($line);
    if ($line =~ m/^\>/) { #determines if this is a header line
        my @cols = split(/\|/, $line);
        if ($cols[4]){ #checks to see if there is coordinate information for a UTR
            $name = "$cols[1]:$cols[2]"; #creates an identifier with the transcript ID and the gene name
            $chr = $cols[3];
            $strand = $cols[6];
            if ($cols[4] =~ /;/) { #determines if the UTR is spliced
                @UTR_starts = split(";", $cols[4]); #extracts start positions of each UTR block
                @UTR_ends = split(";", $cols[5]); #extracts end positions of each UTR block
            }
            else { #if the UTR is not spliced, can just put the start and end coordinates in their arrays (which will have length 1)
                push (@UTR_starts, $cols[4]);
                push (@UTR_ends, $cols[5]);
            }
            $flag = 0; #signals that everything is ok
            print "$name\t$chr\t$UTR_starts[0]\n";
        }
        else {
            $flag = 1; #signals that coordinate information for the UTR is missing: need to skip the "sequence" on the next line
        }
    }
    else { #if it isn't a header line, it's a sequence line
        next if ($flag == 1); #if there wasn't coordinate information on the line above, need to skip this line
        my @pos = match_all_positions($regex, $line); #uses subroutine to create an array of start positions, relative to the sequence, of the regex motif
        print "@pos\n";
        foreach my $motif_start (@pos) {
        my $c_exon_size = 0;
            if ($strand == 1) {
                for (my $i = 0; $i < scalar @UTR_starts; $i = $i + 1) { #goes through the set of exons
                    my $temp_coord = $motif_start - $c_exon_size + $UTR_starts[$i]; #calculates a temporary coordinate for the motif, assuming there is not a splice junction
                    print "UTR_start: $UTR_starts[$i]\n";
                    print "temp_coord: $temp_coord\n";
                    if ($temp_coord < $UTR_ends[$i]) { #checks to see if the coordinate is in the exon
                        print "success: $temp_coord\n\n";
                        $chrStart = $temp_coord - 1; #convert to 0-based for bed
                        $chrEnd = $chrStart + $motif_length;
                        print "$chr\t$chrStart\t$chrEnd\t$name\n";
                        print OUT "$chr\t$chrStart\t$chrEnd\t$name\n";
                        last;
                    }
                    else {
                        $c_exon_size = $c_exon_size + $UTR_ends[$i] - $UTR_starts[$i]; #if the coordinate is not in the exon, adds the exon size to the cumulative exon size
                        print "c_exon_size: $c_exon_size\n";
                    }
                }
            }
            elsif ($strand == -1) { #detects genes on the negative strand
                for (my $i = (scalar @UTR_starts)-1; $i >= 0; $i = $i - 1) { #goes through the set of exons, starting with the last one
                    my $temp_coord = $UTR_ends[$i] - $c_exon_size - $motif_start; #calculates a temporary coordinate for the motif, assuming there is not splice junction
                    print "UTR_start: $UTR_ends[$i]\n";
                    print "temp_coord: $temp_coord\n";
                    if ($temp_coord > $UTR_starts[$i]) { #checks to see if the coordinate is in the exon
                        print "success: $temp_coord\n\n";
                        $chrEnd = $temp_coord;
                        $chrStart = $chrEnd - $motif_length;
                        print "$chr\t$chrStart\t$chrEnd\t$name\n";
                        print OUT "$chr\t$chrStart\t$chrEnd\t$name\n";
                        last;
                    }
                    else {
                        $c_exon_size = $c_exon_size + $UTR_ends[$i] - $UTR_starts[$i]; #if the coordinate is not in the exon, adds the exon size to the cumulative exon size
                        print "c_exon_size: $c_exon_size\n";
                    }
                    
                }
            }
            else { #if there is no strand information, skip it.
                next;
            }
        }
        @UTR_starts = ();
        @UTR_ends = ();
    }
}

close(INF);
close(OUT);

system("sort -k1 -k2,3n $file.$regex.temp | uniq > $file.$regex.bed" ); #can get duplicates because of transcripts. Want to change this or not?

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
