# !/usr/bin/perl

use warnings;
use strict;

my ($regex, $file) = @ARGV;

open (INF, "<$file");
my $flag = 0;
my $exon;
#my @UTR;
my @UTR_starts;
my @UTR_ends;

while (my $line = <INF>) {
    chomp($line);
    if ($line =~ m/^\>/) { #determines if this is a header line
        my @cols = split(/\|/, $line);
        if ($cols[4]){ #checks to see if there is coordinate information for a UTR
            my $gene = "$cols[1]:$cols[2]"; #creates an identifier with the transcript ID and the gene name
            my $chr = $cols[3];
            if ($cols[4] =~ /;/) { #determines if the UTR is spliced
                @UTR_starts = split(";", $cols[4]); #extracts start positions of each UTR block
                @UTR_ends = split(";", $cols[5]); #extracts end positions of each UTR block
                #for (my $i = 0; $i < scalar @starts; $i = $i + 1) { #goes through each UTR block, joining the start and end coordinates, and adds the joined coordinates to an array
                    #$exon = "$starts[$i]:$ends[$i]";
                    #push (@UTR, $exon);
                #}
            }
            else { #if the UTR is not spliced, can just put the start and end coordinates in the array (which will have length 1)
                #$exon = "$cols[4]:$cols[5]";
                #push (@UTR, $exon);
                push (@UTR_starts, $cols[4]);
                push (@UTR_ends, $cols[5]);
            }
            $flag = 0; #signals that everything is ok
            print "$gene\t$chr\t$UTR_starts[0]\n";
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
            for (my $i = 0; $i < scalar @UTR_starts; $i = $i + 1) {
                my $temp_coord = $motif_start - $c_exon_size + $UTR_starts[$i];
                print "UTR_start: $UTR_starts[$i]\n";
                print "temp_coord: $temp_coord\n";
                if ($temp_coord < $UTR_ends[$i]) {
                    print "success: $temp_coord\n\n";
                    last;
                }
                else {
                    $c_exon_size = $c_exon_size + $UTR_ends[$i] - $UTR_starts[$i];
                    print "c_exon_size: $c_exon_size\n";
                }
                
            }
            #my $coord = $start + $UTR_start[0];
            #print "$coord\n";
        }
        @UTR_starts = ();
        @UTR_ends = ();
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
